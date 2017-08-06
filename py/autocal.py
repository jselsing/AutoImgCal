#!/usr/local/anaconda3/envs/py36 python
# -*- coding: utf-8 -*-

# Plotting
import matplotlib; matplotlib.use('TkAgg')
import matplotlib.pyplot as pl
import seaborn; seaborn.set_style('ticks')

# Imports
import numpy as np
import scipy.stats
import astroscrappy
import math
import subprocess
import os
import sys
import logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)
from astropy.io import fits


sexpath = ''  # if "sex" works in any directory, leave blank

defaulttolerance = 0.01  # these defaults should generally not be altered.
defaultpatolerance = 1.4   
defaultminfwhm = 1.5
defaultmaxfwhm = 40

fastmatch = 1
showmatches = 0


def writeparfile():
    params = '''X_IMAGE
    Y_IMAGE
    ALPHA_J2000
    DELTA_J2000
    MAG_AUTO
    MAGERR_AUTO
    ELLIPTICITY
    FWHM_IMAGE
    FLAGS'''
    pf = open('temp.param','w')
    pf.write(params)
    pf.close()


def writeconfigfile(satlevel=55000.):
    configs='''
    #-------------------------------- Catalog ------------------------------------

    CATALOG_NAME     temp_sex.cat   # name of the output catalog
    CATALOG_TYPE     ASCII_HEAD     # NONE,ASCII,ASCII_HEAD, ASCII_SKYCAT,
                                    # ASCII_VOTABLE, FITS_1.0 or FITS_LDAC
    PARAMETERS_NAME  temp.param     # name of the file containing catalog contents

    #------------------------------- Extraction ----------------------------------

    DETECT_TYPE      CCD            # CCD (linear) or PHOTO (with gamma correction)
    DETECT_MINAREA   5              # minimum number of pixels above threshold
    DETECT_THRESH    3              # <sigmas> or <threshold>,<ZP> in mag.arcsec-2
    ANALYSIS_THRESH  3              # <sigmas> or <threshold>,<ZP> in mag.arcsec-2

    FILTER           Y              # apply filter for detection (Y or N)?
    FILTER_NAME      sex.conv       # name of the file containing the filter

    DEBLEND_NTHRESH  16             # Number of deblending sub-thresholds
    DEBLEND_MINCONT  0.02           # Minimum contrast parameter for deblending

    CLEAN            Y              # Clean spurious detections? (Y or N)?
    CLEAN_PARAM      1.0            # Cleaning efficiency

    MASK_TYPE        CORRECT        # type of detection MASKing: can be one of
                                    # NONE, BLANK or CORRECT

    #------------------------------ Photometry -----------------------------------

    PHOT_APERTURES   5              # MAG_APER aperture diameter(s) in pixels
    PHOT_AUTOPARAMS  2.5, 3.5       # MAG_AUTO parameters: <Kron_fact>,<min_radius>
    PHOT_PETROPARAMS 2.0, 3.5       # MAG_PETRO parameters: <Petrosian_fact>,
                                    # <min_radius>



    MAG_ZEROPOINT    0.0            # magnitude zero-point
    MAG_GAMMA        4.0            # gamma of emulsion (for photographic scans)
    GAIN             0.0            # detector gain in e-/ADU
    PIXEL_SCALE      1.0            # size of pixel in arcsec (0=use FITS WCS info)

    #------------------------- Star/Galaxy Separation ----------------------------

    SEEING_FWHM      1.2            # stellar FWHM in arcsec
    STARNNW_NAME     default.nnw    # Neural-Network_Weight table filename

    #------------------------------ Background -----------------------------------

    BACK_SIZE        64             # Background mesh: <size> or <width>,<height>
    BACK_FILTERSIZE  3              # Background filter: <size> or <width>,<height>

    BACKPHOTO_TYPE   GLOBAL         # can be GLOBAL or LOCAL

    #------------------------------ Check Image ----------------------------------

    CHECKIMAGE_TYPE  NONE           # can be NONE, BACKGROUND, BACKGROUND_RMS,
                                    # MINIBACKGROUND, MINIBACK_RMS, -BACKGROUND,
                                    # FILTERED, OBJECTS, -OBJECTS, SEGMENTATION,
                                    # or APERTURES
    CHECKIMAGE_NAME  check.fits     # Filename for the check-image

    #--------------------- Memory (change with caution!) -------------------------

    MEMORY_OBJSTACK  3000           # number of objects in stack
    MEMORY_PIXSTACK  300000         # number of pixels in stack
    MEMORY_BUFSIZE   1024           # number of lines in buffer

    #----------------------------- Miscellaneous ---------------------------------

    VERBOSE_TYPE     QUIET          # can be QUIET, NORMAL or FULL
    WRITE_XML        N              # Write XML file (Y/N)?
    XML_NAME         sex.xml        # Filename for XML output
    '''
    #SATUR_LEVEL      '''+str(satlevel)+'''        # level (in ADUs) at which arises saturation
    pf = open('sex.config','w')
    pf.write(configs)
    pf.close()

    convol='''CONV NORM
    # 3x3 ``all-ground'' convolution mask with FWHM = 2 pixels.
    1 2 1
    2 4 2
    1 2 1
    '''
    if not os.path.exists('sex.conv'): 
        cf = open('sex.conv','w')
        cf.write(convol)
        cf.close()

class Obj:
    ra = 0.0
    dec = 0.0
    mag = 0.0

    ra_rad = 0.0
    dec_rad = 0.0

    def __init__(self, inra, indec, inmag):
        self.ra = inra
        self.dec = indec
        self.ra_rad = inra * math.pi/180
        self.dec_rad = indec * math.pi/180
        self.mag = inmag

    def rotate(self, dpa_deg, ra0, dec0):
        dpa_rad = dpa_deg * math.pi/180
        sindpa = sin(dpa_rad)
        cosdpa = cos(dpa_rad)
        rascale = cos(dec0*math.pi/180)

        #this is only valid for small fields away from the pole.
        x = (self.ra  - ra0 ) * rascale
        y = (self.dec - dec0)

        xrot = cosdpa * x - sindpa * y
        yrot = sindpa * x + cosdpa * y

        self.ra   = (xrot / rascale) + ra0
        self.dec  =  yrot + dec0
        self.ra_rad  = self.ra  * math.pi/180
        self.dec_rad =  self.dec * math.pi/180

class SexObj(Obj):
    x = 0.
    y = 0.
    mag = 0.0
    magerr = 0.0
    ellip = 0.0
    fwhm = 0.0
    flag = 0

    def __init__(self, inline):
        inlinearg = inline.split()

        if len(inlinearg) < 8: return # maybe throw an error?
        self.x = float(inlinearg[0])
        self.y = float(inlinearg[1])
        self.ra = float(inlinearg[2])
        self.dec = float(inlinearg[3])
        self.mag = float(inlinearg[4])
        self.magerr = float(inlinearg[5])
        self.ellip = float(inlinearg[6])
        self.fwhm = float(inlinearg[7])
        if len(inlinearg) >= 9: self.flag = int(inlinearg[8])

        self.ra_rad = self.ra * math.pi/180
        self.dec_rad = self.dec * math.pi/180


def writetextfile(filename, objlist):
    out = open(filename,'w')
    for ob in objlist:
      out.write("%11.7f %11.7f %5.2f\n" % (ob.ra, ob.dec, ob.mag))
    out.close()


def writeregionfile(filename, objlist, color="green",sys=''):
    if sys == '': sys = 'wcs'
    out = open(filename,'w')
    i = -1
    out.write('# Region file format: DS9 version 4.0\nglobal color='+color+' font="helvetica 10 normal" select=1 highlite=1 edit=1 move=1 delete=1 include=1 fixed=0 source\n')
    if sys == 'wcs':
      out.write('fk5\n')
      for ob in objlist:
        i += 1
        out.write("point(%.7f,%.7f) # point=boxcircle text={%i}\n" % (ob.ra, ob.dec, i))
    if sys == 'img':
      out.write('image\n')
      for ob in objlist:
        i += 1
        out.write("point(%.3f,%.3f) # point=boxcircle text={%i}\n" % (ob.x, ob.y, i))
    out.close()


def sextract(sexfilename, nxpix, nypix, border=3, corner=12, minfwhm=1.5, maxfwhm=25, maxellip=0.5, saturation=-1):

    if maxellip == -1: maxellip = 0.5
    if saturation > 0:
       sexsaturation = saturation
    else:
       sexsaturation = 1e10

    try:
       # Sextract the image !
       os.system(sexpath + "sex " + sexfilename + " -c sex.config -SATUR_LEVEL "+str(sexsaturation))
    except (OSError, IOError):
       logger.warn("Sextractor failed to be executed.", exc_info=1)
       sys.exit(1)

    # Read in the sextractor catalog
    try:
       cat = open("temp_sex.cat",'r')
       catlines = cat.readlines()
       cat.close()
    except:
        logger.warn("Cannot load sextractor output file!", exc_info=1)
        sys.exit(1)

    if len(catlines) == 0:
        logger.warn("Sextractor catalog is empty: try a different catalog?", exc_info=1)
        sys.exit(1)

    minx = border
    miny = border
    maxx = nxpix - border    # This should be generalized
    maxy = nypix - border


    l = -1
    nsexinit = 0
    nsexpass = 0
    xlist = []
    ylist = []
    sexlist = []
    fwhmlist = []
    elliplist = []
    flaglist = []
    while l < len(catlines)-1:
        l += 1
        if (len(catlines[l]) <= 1 or catlines[l][0] == '#'):
            continue

        iobj = SexObj(catlines[l]) #process the line into an object
        nsexinit += 1

        #Initial filtering
        if iobj.ellip > maxellip : continue
        if iobj.fwhm < minfwhm: continue
        if iobj.fwhm > maxfwhm: continue
        if iobj.x < minx: continue
        if iobj.y < miny: continue
        if iobj.x > maxx: continue
        if iobj.y > maxy: continue
        if iobj.x + iobj.y < corner: continue
        if iobj.x + (nypix-iobj.y) < corner: continue
        if (nxpix-iobj.x) < corner: continue
        if (nxpix-iobj.x) + (nypix-iobj.y) < corner: continue
        if saturation > 0:
           if iobj.flag > 0: continue  # this will likely overdo it for very deep fields.

        sexlist.append(iobj)
        xlist.append(iobj.x)
        ylist.append(iobj.y)
        fwhmlist.append(iobj.fwhm)
        elliplist.append(iobj.ellip)
        flaglist.append(iobj.flag)
        nsexpass += 1

    print(nsexinit, 'raw sextractor detections')
    print(nsexpass, 'pass initial critiera')

     # Remove detections along bad columns
    threshprob = 0.0001
    ctbadcol = 0
    for i in range(5):
        txp = 1.0
        xthresh = 1
        while txp > threshprob:
          txp *= min((len(sexlist)*1.0/nxpix),0.8) # some strange way of estimating the threshold.
          xthresh += 1                          #what I really want is a general analytic expression for
        removelist = []                         #the 99.99% prob. threshold for value of n for >=n out 
        modex = scipy.stats.mode(xlist)[0]                 #of N total sources to land in the same bin (of NX total bins)
        for j in range(len(sexlist)):
           if (sexlist[j].x > modex-1) and (sexlist[j].x < modex+1):
             removelist.append(j)
        removelist.reverse()
        if len(removelist) > xthresh:
         #print removelist
         for k in removelist:
           del xlist[k]
           del ylist[k]
           del sexlist[k]
           del fwhmlist[k]
           del elliplist[k]
           del flaglist[k]
           ctbadcol += 1

        typ = 1.0
        ythresh = 1
        while typ > threshprob:
          typ *= min((len(sexlist)*1.0/nypix),0.8)
          ythresh += 1
        removelist = []
        modey = scipy.stats.mode(ylist)[0]
        for j in range(len(sexlist)):
           if (sexlist[j].y > modey-1) and (sexlist[j].y < modey+1):
             removelist.append(j)
        removelist.reverse()
        if len(removelist) > ythresh:
         for k in removelist:
           del xlist[k]
           del ylist[k]
           del sexlist[k]
           del fwhmlist[k]
           del elliplist[k]
           del flaglist[k]
           ctbadcol += 1
    if ctbadcol > 0: print(' Removed ', ctbadcol, ' detections along bad columns.')

    # Remove galaxies and cosmic rays
    if len(fwhmlist) > 5:
       # fwhmlist.sort()
       fwhm20 = np.percentile(fwhmlist, 0.2)
       fwhm25 = np.percentile(fwhmlist, 0.25)
       fwhm50 = np.percentile(fwhmlist, 0.50)     #percentile values
       fwhm75 = np.percentile(fwhmlist, 0.75)
       fwhmmode = scipy.stats.mode(fwhmlist)[0]
    else:
       fwhmmode = minfwhm
       fwhm20 = minfwhm
    #hifwhmlist = []
    #for f in fwhmlist:
    #   if f > fwhmmode*1.5: hifwhmlist.append(f)
    #fwhmhimode = mode(hifwhmlist)

 #   ellipmode = mode(elliplist)   # in theory, if this is large we could use a theta cut, too.  (THETA_IMAGE)
 #   ellipstdev = stdev(elliplist)
 #   elliptol = 1.0 #min(0.25, 2.5*ellipstdev)
             # this effectively disables the ellipticity filter.


    # formerly a max, but occasionally a preponderance of long CR's could cause fwhmmode to be bigger than the stars
    refinedminfwhm = np.median([0.75*fwhmmode,0.9*fwhm20,minfwhm]) # if CR's are bigger and more common than stars, this is dangerous...
    print('Refined min FWHM:', refinedminfwhm, 'pix')
    #refinedmaxfwhm = 35


    ngood = 0
    goodsexlist = []
    for sex in sexlist:
       if sex.fwhm > refinedminfwhm: # and sex.ellip < ellipmode +elliptol:
          goodsexlist.append(sex)
          ngood += 1

    # i = 0
    # for sex in goodsexlist:
    #      if i < 1000: print(i, sex.mag, sex.fwhm)
    #      i += 1

    writetextfile('det.init.txt', goodsexlist)
    writeregionfile('det.im.reg', goodsexlist, 'red', 'img')

    print(len(sexlist), 'objects detected in image ('+ str(len(sexlist)-len(goodsexlist)) +' discarded)')


    return goodsexlist


def get_catalog(img_ra, img_dec, img_filt, radius = 1, catalog = "PS"):

  gr_cat_arg = "python gr_cat.py -c %s%s -r %s -s %s -b %s -f temp_cat.dat"%(img_ra, img_dec, radius, catalog, img_filt)

  # Run gr_cat_arg to get catalog around ra and dec
  try:
      os.system(gr_cat_arg)
  except (OSError, IOError):
      logger.warn("gr_cat.py failed to be executed.", exc_info=1)

  # Read in the catalog
  try:
      cat = np.genfromtxt("temp_cat.dat")
  except (OSError, IOError):
      logger.warn("Cannot load catalog file!", exc_info=1)
  # Check for exsistence of targets
  if cat.shape[0] == 0:
      logger.warn("Catalog is empty: try a different catalog?", exc_info=1)
      sys.exit(1)

  return cat



def main():
    """
    Rutine to automatically do astrometric calibration and photometry of detected sources. Uses astrometry.net to correct the astrometric solution of the image. This correction includes image distortions. Queries  Pan-STARRS, SDSS and USNO in that order for coverage for reference photometry against which to do the calibration. This is achieved with gr_cat.py developed by Thomas Krühler which can be consulted for additional documentation. Sextractor is run on the astrometrically calibrated image using the function sextract, written by Daniel Perley and available at http://www.dark-cosmology.dk/~dperley/code/code.html. Handling of the entire sextractor interfacing is heavily based on autoastrometry.py The two lists of images are then matched with a k-d tree algorithm and sextracted magntiudes can be calibrated against the chosen catalog. 
    """

    filename = "../test_data/XSHOO.2016-10-18T00:09:56.640.fits"
    fitsfile = fits.open(filename)
    header = fitsfile[0].header
    img_filt = header["HIERARCH ESO INS FILT1 NAME"][0]

    img_ra, img_dec = header["RA"], header["DEC"]


    cat = get_catalog(img_ra, img_dec, img_filt)

    writeparfile()
    saturation = -1
    if not os.path.exists('sex.config'): writeconfigfile(saturation)


    #Read the header info from the file for sextractor
    try:
        # no longer drawing RA and DEC from here.
        key = 'NAXIS1'
        nxpix = header[key]
        key = 'NAXIS2'
        nypix = header[key]
    except:
        logger.warn("Cannot find necessary WCS header keyword.", exc_info=1)
        sys.exit(1)
    try:
        key = 'CRVAL1'
        cra =  float(header[key])
        key = 'CRVAL2'
        cdec = float(header[key])

        key = 'CRPIX1'
        crpix1 = float(header[key])
        key = 'CRPIX2'
        crpix2 = float(header[key])

        key = 'CD1_1'
        cd11 = float(header[key])
        key = 'CD2_2'
        cd22 = float(header[key])
        key = 'CD1_2'
        cd12 = float(header[key]) # deg / pix
        key = 'CD2_1'
        cd21 = float(header[key])

        equinox = float(header.get('EQUINOX', 2000.))
        if abs(equinox-2000.) > 0.2: print('Warning: EQUINOX is not 2000.0')
    except:
        if pixelscale == -1:
            print('Cannot find necessary WCS header keyword', key)
            logger.warn("Must specify pixel scale (-px VAL) or provide provisional basic WCS info via CD matrix.", exc_info=1)
            #Some images might use CROT parameters, could try to be compatible with this too...?
            sys.exit(1)
    seeing = -1
    if (seeing == -1):
        minfwhm = defaultminfwhm #1.5
        maxfwhm = defaultmaxfwhm #40
    else:
        minfwhm = 0.7 * seeing
        maxfwhm = 2 * seeing
    maxellip = -1

    gain = fitsfile[0].header['HIERARCH ESO DET OUT1 GAIN']
    ron = fitsfile[0].header['HIERARCH ESO DET OUT1 RON']

    # Clean for cosmics
    crmask, clean_arr = astroscrappy.detect_cosmics(fitsfile[0].data, cleantype='medmask', sepmed=True)

    # Replace data array with cleaned image
    fitsfile[0].data = clean_arr

    # Update file
    fitsfile.writeto(filename+"cosmicced.fits", output_verify='fix', clobber=True)

    # Sextract stars to produce image star catalog
    goodsexlist = sextract(filename+"cosmicced.fits", nxpix, nypix, 3, 12, minfwhm=minfwhm, maxfwhm=maxfwhm, maxellip=maxellip, saturation=saturation)



    try:
       os.remove('temp\*')
    except:
       print('Could not remove temp.param for some reason')

if __name__ == '__main__':
    main()