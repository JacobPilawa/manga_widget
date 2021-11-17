# Galaxies Widget, Fall 2021

## Background


I wanted to do some visualization of the MANGA IFU data, and I came across the package [marvin](https://sdss-marvin.readthedocs.io/en/latest/index.html) which helps you do that. The documentation and tutorials are the primary inspiration and guide for this widget. 

There are essentially a bunch of sub-functions which do various analyses on the MANGA data. One set of routines produces plots of the gas and stellar kinematics. The other set of routines finds the spaxels which are labelled "star-forming" from their BPT diagram placement, and uses some stellar population models to infer the stellar mass and stellar surface density of a galaxy. It results in a nice plot of the "spatially resolved stellar surface density-gas phase metallicity relation."

All you need to do is specify the object_name, which in this case refers to the mange "plateifu" number (which is almost universally reported in the literature). 

There are sample plots included in the ``kinematic_plots'' and ``metallicity_plots'' directories in case you don't want to run the functions yourself. They can take a minute to run. 

## Instructions

* Github wouldn't let me upload one of the files you need (though you can still run the kinematic analysis without the file). The file you need is the summary file from the [Firefly stellar population fitting results file](https://data.sdss.org/sas/dr14/manga/spectro/firefly/v1_0_3/manga_firefly-v2_1_2-STELLARPOP.fits) (1.8 GB). Create a directory called "data/" in the same directory as widget.ipynb. Place the .fits Firefly summary file in data. 
* "widget.ipydnb" contains a sample run of the functions with the outputs. The main functionality is:
    * run_kinematics(object_name, save_dir): 
        * object_name: manga plateifu designation; if object_name = None, a random galaxy is queried and displayed which can make for some fun images!
        * save_dir: name of directory to be generated for plot storage. automatically creates directory if it does not exist. 
        * purpose: generates three plots in directory "save_dir": (i) SDSS image of galaxy, (b) 6 kinematic maps of the stellar and gas velocities, dispersions, and differences, and (iii) plots of the Halpha emission lines and Calcium H and K absorption features used to meausre the velocities and dispersions. Note that the third plot is sometimes inaccurate due to inaccurate redshifts in the reported table/non-detections in the particular spectra I choose to display. A more careful treatment of this would be nice!
    * run_metallicity():
        * object_name and save_dir as same as above, including random query
        * purpose: generates five plots: (i) SDSS image of galaxy, (ii) bpt diagram of the spaxels to show which are star forming spaxels, (iii) plot of the star forming spaxel metallicity using calibration from Petini and Pagel+2004, (iv) stellar mass and stellar surface density maps from the Firefly stellar population fitting results, (v) plot of stellar surface density vs gas-phase metallicity for our spaxels and the fit to the galaxy version from Barerra-Ballestros+2016
