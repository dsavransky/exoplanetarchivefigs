import numpy as np
import matplotlib
import warnings
import re
import matplotlib.pyplot as plt
import pandas
import astropy.units as u
import astropy.constants as const
from EXOSIMS.util.getExoplanetArchive import (
    getExoplanetArchivePSCP,
    cacheExoplanetArchiveQuery,
)
from EXOSIMS.PlanetPhysicalModel.ForecasterMod import ForecasterMod


def p2sma(mu: u.Quantity, T: u.Quantity):
    """Compute semi-major axis from period and gravitational parameter

    Args:
        mu (~astropy.units.Quantity(np.ndarray(float))):
            Gravitational parameters
        T (~astropy.units.Quantity(np.ndarray(float))):
            Orbital periods

    Returns:
        ~astropy.units.Quantity(np.ndarray(float)):
            semi-major axes

    """

    return ((mu * T**2 / (4 * np.pi**2)) ** (1 / 3.0)).to("AU")


def get_data(min_num_discoveries: int = 30, forceNew: bool = False) -> pandas.DataFrame:
    """
    Retrieve planet data from the Exoplanet Archive and combine with spectral data to
    determine which planets have any photometric measurements

    Args:
        min_num_discoveries (int):
            Any methods with fewer than this number of discoveries are all lumped
            together into a single 'other' category.  To show all methods, set to 0.
            Defaults to 30.
        forceNew (bool):
            Run a fresh query even if results exist on disk.


    """
    data = getExoplanetArchivePSCP(forceNew=forceNew)

    # lump detections with few discoveries into a single 'other' category
    methods, methods_counts = np.unique(
        data["discoverymethod"].values, return_counts=True
    )
    if min_num_discoveries > 0 and np.any(methods_counts < min_num_discoveries):
        othermeths = methods[methods_counts < min_num_discoveries]
        data.loc[data["discoverymethod"].isin(othermeths), "discoverymethod"] = "Other"

    # add column for photometric info
    data["has_photometry"] = False
    data["first_photometry_year"] = np.nan

    # blanket assumption that direct imaging discoveries have some photometric info
    imrows = data["discoverymethod"] == "Imaging"
    data.loc[imrows, "has_photometry"] = True
    data.loc[imrows, "first_photometry_year"] = data.loc[imrows, "disc_year"]

    # technically, we'd want to also do this everywhere where ima_flag = 1, but the two
    # current (circa 2023) non-imaging discoveries with this flag appear to be false
    # positives in the archive. It's also not clear how to figure out *when* the imaging
    # detections took place. TODO: come back to this after response from IPAC support.

    # get spectroscopy data
    basestr = "exoplanetArchiveSpectra"
    querystring = r"select+*+from+spectra"
    specdata = cacheExoplanetArchiveQuery(basestr, querystring, forceNew=forceNew)

    pnames = specdata["pl_name"].values
    upnames = np.unique(pnames)
    if len(set(upnames) - set(data["pl_name"])):
        warnings.warn(
            "Some planets in spectra table do not exist in the composite systems table."
        )

    pubs = specdata["bibcode"].values
    pubdates = np.array([re.match(r"(20\d\d)", d).groups(0)[0] for d in pubs]).astype(
        float
    )

    firstpubdate = np.zeros(upnames.shape)
    for j, n in enumerate(upnames):
        firstpubdate[j] = pubdates[pnames == n].min()

    for n, d in zip(upnames, firstpubdate):
        row = data["pl_name"] == n
        data.loc[row, "first_photometry_year"] = d
        data.loc[row, "has_photometry"] = True

    # final consistency check:
    if not np.all(
        data.loc[(data["pl_nespec"] > 0) | (data["pl_ntranspec"] > 0), "has_photometry"]
    ):
        warnings.warn(
            "There are planets with emission or transmission spectra indicated in the "
            "composite systems table that do not appear in the spectra table."
        )

    return data


def mass_sma_plot(
    data: pandas.DataFrame,
    fig: matplotlib.figure.Figure = None,
    earthmassonly: bool = False,
    title: str = None,
    photom: bool = False,
    syms: str = "osD^pv<vh",
    clrs: list = ["blue", "silver", "crimson", "orange", "darkmagenta"],
    xlim: list = [1e-2, 1e3],
    ylim: list = [1e-4, 40],
) -> None:
    """
    Generate mass vs sma planet plot

    Args:
        data (pandas.DataFrame):
            Output of get_data
        fig (matplotlib.figure.Figure, optional):
            Figure to use. If None, a new figure is craeted.  Defaults None.
            If a figure is passed in, it will be cleared.
        earthmassonly (bool):
            If True, plot only Earth masses. Otherwise show both Earth and Jovian
            masses. Defaults False.
        title (str, optional):
            Plot title.  Defaults to None.
        photom (bool):
            If True, shade only those planets with some photometric measurements and
            leave the rest transparent. Defaults False
        syms (str):
            Plot symbols to use.  If fewer symbols are given than there are methods in
            the data set, the sequence is repeated.
        clrs (list(str)):
            List of named colors to use.  If fewer colors are given than there are
            methods in the data set, the sequence is repeated.
        xlim (list(float)):
            [min, max] x axis values in AU.  Defaults to [1e-2, 1e3]
        ylim (list(float)):
            [min, max] y axis values in Jovian masses. Defaults to [1e-4, 40]

    """

    # define some helper methods
    fmod = ForecasterMod()
    MfromR = fmod.calc_mass_from_radius

    # create or clear the figure and axes
    if fig is None:
        fig, ax = plt.subplots(figsize=(10, 5))
    else:
        fig.clf()
        ax = fig.add_subplot()

    # fill figure with axis
    fig.subplots_adjust(bottom=0.15, top=0.95, left=0.125, right=0.9)

    # figure out distinct methods and set zorder by increasing discovery count
    methods, methods_counts = np.unique(
        data["discoverymethod"].values, return_counts=True
    )
    methodorder = np.zeros(len(methods), dtype=int)
    methodorder[np.argsort(methods_counts)[::-1]] = np.arange(len(methods))

    # ensure that you have enough symbols/colors for plotting
    if len(syms) < len(methods):
        syms = "".join([syms] * int(np.ceil(len(methods) / len(syms))))
    if len(clrs) < len(methods):
        clrs = list(np.hstack([clrs] * int(np.ceil(len(methods) / len(clrs)))))

    # define solar system planets
    # data from https://ssd.jpl.nasa.gov/astro_par.html
    planetnames = [
        "Mercury",
        "Venus",
        "Earth",
        "Mars",
        "Jupiter",
        "Saturn",
        "Uranus",
        "Neptune",
        "Pluto",  # yes, its here. deal with it
    ]
    GMs = [
        22031.868551,  # Me
        324858.592000,  # v
        398600.435507,  # E
        42828.375816,  # Ma
        126712764.100000,  # J
        37940584.841800,  # S
        5794556.400000,  # U
        6836527.100580,  # N
        975.500000,  # P
    ]  # km^3/s^2
    Ms = np.array(GMs) / GMs[4]  # Jupiter Masses
    Mse = np.array(GMs) / GMs[2]  # Earth Masses

    # semi-major axes from: https://ssd.jpl.nasa.gov/planets/approx_pos.html
    smas = np.array(
        [
            0.38709927,  # Me
            0.72333566,  # V
            1.00000261,  # E
            1.52371034,  # Ma
            5.20288700,  # J
            9.53667594,  # S
            19.18916464,  # U
            30.06992276,  # N
            39.482,  # P
        ]
    )  # AU

    # pick an offset for each label
    offs = [
        (5, -5),  # Me
        (-45, -15),  # V
        (0, -17),  # E
        (4, -4),  # Ma
        (6, -4),  # J
        (5, -2),  # S
        (-5, -16),  # U
        (-5, 5),  # N
        (-5, 5),  # P
    ]

    solid_alpha = 0.95 if photom else 0.75

    # loop by method and plot
    for m, s, c, o in zip(methods, syms, clrs, methodorder):
        inds = data["discoverymethod"] == m
        mj = data.loc[inds, "pl_bmassj"]
        sma = data.loc[inds, "pl_orbsmax"]
        # determine which points should be plotted as solid based on phot input
        make_solid = data.loc[inds, "has_photometry"]
        if not photom:
            make_solid.loc[~make_solid] = True

        # fill in missing masses from radii
        radj = data.loc[mj.index[mj.isna()]]["pl_radj"]
        mfills = MfromR(radj.values * u.R_jup).to(u.M_jup).value
        mj.loc[mj.isna()] = mfills

        # fill in missing smas from period
        orbper = data.loc[sma.index[sma.isna()]]["pl_orbper"]
        stmass = data.loc[sma.index[sma.isna()]]["st_mass"]

        GMs = const.G * (stmass.values * u.solMass)  # units of solar mass
        T = orbper.values * u.day
        smafill = p2sma(GMs, T)
        sma.loc[sma.isna()] = smafill.value

        # if only plotting earth masses, scale accordingly
        if earthmassonly:
            mj /= Ms[2]

        # plot solid symbols
        ax.scatter(
            sma.loc[make_solid],
            mj.loc[make_solid],
            marker=s,
            s=60,
            zorder=o,
            facecolors=c,
            edgecolors="k",
            alpha=solid_alpha,
            label=m,
        )
        # plot transparent symbols
        ax.scatter(
            sma.loc[~make_solid],
            mj.loc[~make_solid],
            marker=s,
            s=60,
            zorder=-1,
            facecolors=c,
            edgecolors="k",
            alpha=0.1,
            label=None,
        )

    # now plot the solar system planets
    tmp = Mse if earthmassonly else Ms
    ax.scatter(
        smas,
        tmp,
        marker="o",
        s=60,
        facecolors="yellow",
        edgecolors="k",
        alpha=1,
        zorder=methodorder.max(),
    )
    for a, m, n, off in zip(smas, tmp, planetnames, offs):
        ax.annotate(n, (a, m), ha="left", xytext=off, textcoords="offset points")

    # let's log scale and trim
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(xlim)
    ax.set_xlabel("Semi-Major Axis (AU)")

    # label y axis and add another one, as needed
    if earthmassonly:
        ax.set_ylim(np.array(ylim) / Ms[2])
        ax.set_ylabel("(Minimum) Mass (Earth Masses)")
        plt.subplots_adjust(right=0.95)
    else:
        ax.set_ylim(ylim)
        ax.set_ylabel("(Minimum) Mass (M$_J$)")
        ax2 = ax.twinx()
        ax2.set_yscale("log")
        ax2.set_ylim(np.array(ax.get_ylim()) / Ms[2])
        ax2.set_ylabel(r"M$_\oplus$")
        plt.subplots_adjust(right=0.88)

    ax.legend(loc="lower right", scatterpoints=1, fancybox=True, prop={"size": 14})

    if title is not None:
        plt.title(title)
        plt.subplots_adjust(top=0.93)
