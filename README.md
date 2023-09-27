# exoplanetarchivefigs
Generate figures based on data from NASA's Exoplanet Archive

The code in this repository is a slightly repackaged version of the code used to make Figure 2.5 in the Astro2020 Decadal Survey final report (Pathways to Discovery in Astronomy and Astrophysics for the 2020s - https://nap.nationalacademies.org/read/26141/chapter/4#40).  

## Installation 
To create your own versions of these figures, first clone or download a copy of this repository.  Then, in the top-level directory of the repository, run:

    pip install --user .
   
This will install the backend routines used in the example Jupyter notebook, along with all required python packages.  If you plan on modifying the backend code, you should install in developer mode:

    pip install --user -e .
    
In order to run the example figures notebook, you will also need jupyter/jupyter-lab and ipympl.  You can install them as:

    pip install --user jupyterlab ipympl

See here for more details on using Jupyter: https://jupyter.org/

## What the code is doing

The method `get_data` downloads and caches (or reads from an existing cache) all of the data from the NASA Exoplanet Archive's composite planetary systems and atmospheric spectra tables.  Details on these can be found at: https://exoplanetarchive.ipac.caltech.edu/docs/help.html

Detection methods with fewer than `min_num_discoveries` (default value = 30) all get lumped into a single 'Other' category.  For direclty imaged planets, it is assumed that some photometric information is available.  Similalry, anything with any transmission or emission spectra is similarly assumed to have been photometrically characterized. This information is stored in the boolean column `has_photometry`. For objects found via imaging, the year of the first photometric characterization is taken to be the discovery year.  For objects with spectral information found by other means,  the year of the first photometric characterization is taken to be the year of the earliest publication associated with the planet in the spectra table.  This information is stored in the column `first_photometry_year`.  See the notebook for examples on filtering by these columns.

## Data Anomalies 

The Exoplanet Archive staff do a tremendous job in curating and updating their data, but no dataset is perfect, and there will inevitably be errors/anomalies in the data set.  When you find these, do the entire community a favor and report them to IPAC by submitting a helpdesk ticket, here: https://exoplanetarchive.ipac.caltech.edu/cgi-bin/Helpdesk/nph-genTicketForm

## Crediting

Feel free to use/modify/redistribute figures based on this code in any way you wish, but please consider crediting this repository, and the NASA Exoplanet Archive (see here: https://exoplanetarchive.ipac.caltech.edu/docs/acknowledge.html). 
