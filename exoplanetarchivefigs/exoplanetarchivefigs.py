import numpy as np
import matplotlib
import warnings
import re
import matplotlib.pyplot as plt
import pandas
from EXOSIMS.util.getExoplanetArchive import (
    getExoplanetArchivePSCP,
    cacheExoplanetArchiveQuery,
)


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
    syms: str = None,
):
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

    """

    # create or clear the figure and axes
    if fig is None:
        fig, ax = plt.subplots(figsize=(10, 5))
    else:
        fig.clf()
        ax = fig.add_subplot()

    # fill figure with axis
    fig.subplots_adjust(bottom=0.15, top=0.95, left=0.125, right=0.9)


    methods, methods_counts = np.unique(
        data["discoverymethod"].values, return_counts=True
    )
    methodorder = np.argsort(methods_counts)[::-1]

    syms = "os^pvD<"


    for m, s, c, o in zip(methods, syms, cmap, methodorder):
        inds = data["discoverymethod"] == m
        mj = data.loc[inds, "pl_bmassj"]
        sma = data.loc[inds, "pl_orbsmax"]

        # fill in missing masses from radii
        radj = data.iloc[mj.index[mj.isna()]]["pl_radj"]
        mfills = MfromR(radj.values * u.R_jup).to(u.M_jup).value
        mj.loc[mj.isna()] = mfills

        # fill in missing smas from period
        orbper = data.iloc[sma.index[sma.isna()]]["pl_orbper"]
        stmass = data.iloc[sma.index[sma.isna()]]["st_mass"]

        GMs = const.G * (stmass.values * u.solMass)  # units of solar mass
        T = orbper.values * u.day
        smafill = p2sma(GMs, T)
        sma.loc[sma.isna()] = smafill.value

        if earthmassonly:
            mj /= Ms[2]
        ax.scatter(
            sma,
            mj,
            marker=s,
            s=60,
            zorder=o,
            facecolors=c,
            edgecolors="k",
            alpha=0.75,
            label=m,
        )

    if earthmassonly:
        tmp = Mse
    else:
        tmp = Ms
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
    for a, m, n, ha, off in zip(smas, tmp, planetnames, has, offs):
        ax.annotate(n, (a, m), ha=ha, xytext=off, textcoords="offset points")

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim([1e-2, 1e3])
    ax.set_xlabel("Semi-Major Axis (AU)")

    if earthmassonly:
        ax.set_ylim(np.array([1e-4, 40]) / Ms[2])
        ax.set_ylabel("(Minimum) Mass (Earth Masses)")
    else:
        ax.set_ylim([1e-4, 40])
        ax.set_ylabel("(Minimum) Mass (M$_J$)")
    ax.legend(loc="lower right", scatterpoints=1, fancybox=True, prop={"size": 14})

    if not earthmassonly:
        ax2 = ax.twinx()
        ax2.set_yscale("log")
        ax2.set_ylim(np.array(ax.get_ylim()) / Ms[2])
        ax2.set_ylabel("M$_\oplus$")
        plt.subplots_adjust(right=0.88)
    else:
        plt.subplots_adjust(right=0.95)

    if title:
        plt.title(title)
        plt.subplots_adjust(top=0.93)
