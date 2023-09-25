import numpy as np
from EXOSIMS.util.getExoplanetArchive import getExoplanetArchivePSCP, cacheExoplanetArchiveQuery


def get_data(min_num_discoveries: int = 30, forceNew: bool = False) -> None:
    """
    Retrieve planet data from the Exoplanet Archive and filter

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

        methods, methods_counts = np.unique(
            data["discoverymethod"].values, return_counts=True
        )

    # figure out plotting order (by decreasing counts
    methodorder = np.argsort(methods_counts)[::-1]

    # get spectroscopy data
    basestr = "exoplanetArchiveAtmospheres"
    querystring = r"select+*+from+atmospheres"

    specdata = cacheExoplanetArchiveQuery(basestr, querystring, forceNew=forceNew)



