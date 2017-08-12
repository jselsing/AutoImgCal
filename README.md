# AutoImgCal
Automatic astrometric - and flux calibration of astronomical images

Rutine to automatically do astrometric calibration and photometry of detected sources. Uses astrometry.net to correct the astrometric solution of the image. Input images need to be larger than ~10 arcmin for this to work. This correction includes image distortions. Queries  Pan-STARRS, SDSS and USNO in that order for coverage for reference photometry against which to do the calibration. This is achieved with gr_cat.py developed by Thomas Krühler which can be consulted for additional documentation. Sextractor is run on the astrometrically calibrated image using the function sextract, heavily inspired by autoastrometry.py by Daniel Perley and available at http://www.dark-cosmology.dk/~dperley/code/code.html. Handling of the entire sextractor interfacing is heavily based on autoastrometry.py. The two lists of images are then matched with a k-d tree algorithm and sextracted magntiudes can be calibrated against the chosen catalog.

## Depedencies

We use a range of astronomical software packages for this, but mostly astropy- and astropy affiliated packages that can be easily installed using either conda or pip. 

python dependencies

  - numpy
  - scipy
  - pandas
  - astropy
  - astropy: photutils
  - astropy: astroscrappy
  
Additionally astrometry-net and sextractor are stand-alone packages that are required, but these can easily be install through homebrew. 
