# Galaxies Widget, Fall 2021

## Background


I wanted to do some visualization of the MANGA IFU data, and I came across the package [marvin](https://sdss-marvin.readthedocs.io/en/latest/index.html) which helps you do that. The documentation and tutorials are the primary inspiration and guide for this widget. 

The purpose of the widget is to visualize (i) the kinematics of stars and gas in galaxies, and (ii) plot the mass-metallicity relation of star-forming spaxels (based on BPT diagram placement) compared to the mass-metallicity relation of Barerra-Ballestros+2016. There are two functions that do these tasks, called ``run_kinematics'' and ``run_metallicity'' respectively. The widget.ipynb shows two example function calls. All you need to do is specify the object_name, which in this case refers to the mange "plateifu" number (which is almost universally reported in the literature). There are sample plots included in the ``kinematic_plots'' and ``metallicity_plots'' directories in case you don't want to run the functions yourself. They can take a minute to run. 

I think this widget could be useful in a variety of ways. First, I am a big fan of the plots that are made as a result, and I think they could be useful in lecture or presentations (maybe with slight modification). I think they also can help students make connections between morphology and kinematics, as well as morphology and star formation given the spatially resolved spectral maps that are plotted as a result. In creating this widget, I gained a better understanding of BPT diagrams as well, which might transfer to folks using this widget.

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
