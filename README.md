# AutoImgCal
Automatic astrometric - and flux calibration of astronomical images

Rutine to automatically do astrometric calibration and photometry of detected sources. Uses Astrometry-net to correct the astrometric solution of the image. Input images need to be larger than ~10 arcmin for this to work. Queries Pan-STARRS, SDSS and USNO in that order for coverage for reference photometry against which to do the calibration. This is achieved with a modified version of gr_cat.py developed by Thomas Krühler which can be consulted for additional documentation. Handling of the sextractor interfacing is inspired by autoastrometry.py, written by Daniel Perley and available at http://www.dark-cosmology.dk/~dperley/code/code.html. The two lists of images are then matched with a k-d tree algorithm and sextracted magntiudes can be calibrated against the chosen catalog.

## Depedencies

We use a range of astronomical software packages for this, but mostly astropy- and astropy affiliated packages that can be easily installed using either conda or pip. 

python dependencies

  - numpy
  - scipy
  - pandas
  - astropy
  - astropy: photutils
  - astropy: astroscrappy
  
Additionally Astrometry-net and Sextractor are stand-alone packages that are required, but these can easily be install through homebrew. Additionally, Astrometry-net needs index files which can be downloaded from http://data.astrometry.net/4200/ - please see the README for Astrometry-net at http://astrometry.net/doc/readme.html. Index files need to go in the Astrometry-net dir, e.g. /usr/local/Cellar/astrometry-net/HEAD-99d4344/data/.
