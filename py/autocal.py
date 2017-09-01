#!/usr/local/anaconda3/envs/py36 python
# -*- coding: utf-8 -*-

# Plotting
import matplotlib; matplotlib.use('TkAgg')
import matplotlib.pyplot as pl
# import seaborn; seaborn.set_style('ticks')

# Imports
import numpy as np
import scipy.stats
import astroscrappy
import math
import subprocess
import os
import glob
import sys
import getopt
import logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)
from astropy.io import fits
from astropy import wcs
from astropy.table import Table
from scipy.spatial import cKDTree
import pandas as pd
from upper_limit import limiting_magnitude


sexpath = ''  # if "sex" works in any directory, leave blank

defaulttolerance = 0.01  # these defaults should generally not be altered.
defaultpatolerance = 1.4   
defaultminfwhm = 1.5
defaultmaxfwhm = 40
defaultelliptol = 0.2
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
    FILTER_NAME      sex_temp.conv  # name of the file containing the filter

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
    PIXEL_SCALE      0            # size of pixel in arcsec (0=use FITS WCS info)

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
    pf = open('sex_temp.config','w')
    pf.write(configs)
    pf.close()

    convol='''CONV NORM
    # 3x3 ``all-ground'' convolution mask with FWHM = 2 pixels.
    1 2 1
    2 4 2
    1 2 1
    '''
    if not os.path.exists('sex_temp.conv'): 
        cf = open('sex_temp.conv','w')
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
      out.write("%11.7f %11.7f %5.2f %5.2f %5.2f %5.2f\n" % (ob.ra, ob.dec, ob.mag, ob.magerr, ob.cat_mag, ob.cat_magerr))
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
        out.write("point(%.7f,%.7f) # point=boxcircle text={%.2f +- %0.2f}\n" % (ob.ra, ob.dec, ob.cat_mag, ob.cat_magerr))
    if sys == 'img':
      out.write('image\n')
      for ob in objlist:
        i += 1
        out.write("point(%.3f,%.3f) # point=boxcircle text={%.2f +- %0.2f}\n" % (ob.x, ob.y, ob.cat_mag, ob.cat_magerr))
    out.close()


def sextract(sexfilename, nxpix, nypix, border=3, corner=12, minfwhm=1.5, maxfwhm=25, maxellip=0.5, saturation=-1, zeropoint=0):

    if maxellip == -1: maxellip = 0.5
    if saturation > 0:
       sexsaturation = saturation
    else:
       sexsaturation = 1e10

    try:
       # Sextract the image !
       subprocess.run(['sex', '%s'%sexfilename, '-c', 'sex_temp.config', '-SATUR_LEVEL', '%s'%sexsaturation, '-MAG_ZEROPOINT', '%s'%zeropoint])
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

    # formerly a max, but occasionally a preponderance of long CR's could cause fwhmmode to be bigger than the stars
    refinedminfwhm = np.median([0.75*fwhmmode,0.9*fwhm20,minfwhm]) # if CR's are bigger and more common than stars, this is dangerous...
    print('Refined min FWHM:', refinedminfwhm, 'pix')
    #refinedmaxfwhm = 35


    ngood = 0
    goodsexlist = []
    for sex in sexlist:
       if sex.fwhm > refinedminfwhm and sex.ellip < maxellip:
          goodsexlist.append(sex)
          ngood += 1

    print(len(sexlist), 'objects detected in image ('+ str(len(sexlist)-len(goodsexlist)) +' discarded)')

    return goodsexlist


def get_catalog(img_ra, img_dec, img_filt, radius = 5, catalog = "PS"):

  gr_cat_arg = [sys.executable, 'gr_cat.py', '-c', '%s%s'%(img_ra, img_dec), '-r', '%s'%radius, '-s', '%s'%catalog, '-b', '%s'%img_filt, '-f', 'temp_cat.dat', '-d', 'temp_cat.reg']

  # Run gr_cat_arg to get catalog around ra and dec
  try:
      subprocess.run(gr_cat_arg)
  except (OSError, IOError):
      logger.warn("gr_cat.py failed to be executed.", exc_info=1)

  # Read in the catalog
  try:
      cat = pd.read_csv("temp_cat.dat").values
  except (OSError, IOError):
      logger.warn("Cannot load catalog file!", exc_info=1)
  # Check for exsistence of targets
  if cat.shape[0] == 0:
      logger.warn("Catalog is empty: try a different catalog?", exc_info=1)
      sys.exit(1)
  return cat


def run_astrometry_net(img_name, img_ra, img_dec):
    # Shell command to run astrometry-net
    astrometry_args = ['solve-field', '-g', '-p', '-O', '--fits-image', '%s'%(img_name), '--ra', '%s'%img_ra, '--dec', '%s'%img_dec, '--radius', '%s'%(1/60)]

    # Run astrometry-net on field
    try:
       subprocess.run(astrometry_args)
    except (OSError, IOError):
      logger.warn("astrometry-net failed to be executed.", exc_info=1)

    # Read in the calibrated image
    calib_img_name = img_name.replace("temp", "new")
    try:
        calib_img_name = img_name.replace("temp", "new")
        calib_img = fits.open(calib_img_name)
        img_name = calib_img_name
    except (OSError, IOError):
        # logger.warn("Astrometry solution did not solve! Continuing without astrometric calibration.", exc_info=1)
        logger.warn("Astrometry solution did not solve! Continuing without astrometric calibration.")

    return img_name


def joint_catalog(cat_1, cat_2):
    """
    Small function to match to arrays based on the first two columns, which is assumed to be ra and dec
    """
    # Grow the tree
    tree_data = np.array([cat_1[:, 0], cat_1[:, 1]]).T
    tree = cKDTree(tree_data)

    # Find mapping indices
    idx_map_cat2, idx_map_cat1 = [], []
    tol = 1e-3 # Distance in degrees - This could change depending on the accuracy of the astrometric solution
    for ii, kk in enumerate(cat_2):
        # find the k nearest neighbours
        distance, indice = tree.query(kk[0:2], k=1)
        if distance < tol:
            # Store
            idx_map_cat1.append(indice)
            idx_map_cat2.append(ii)

    cat_1 = cat_1[np.array(idx_map_cat1)]
    cat_2 = cat_2[np.array(idx_map_cat2)]
    # Return joint lists
    return cat_1, cat_2


def autocal(filename = "../test_data/FORS_R_OB_ana.fits", catalog = "SDSS", sigclip = 50, objlim = 75, filter = None, cosmic_rejection = True, astrometry = True):

    """
    Rutine to automatically do astrometric calibration and photometry of detected sources. Uses astrometry.net to correct the astrometric solution of the image. Input images need to be larger than ~10 arcmin for this to work. This correction includes image distortions. Queries  Pan-STARRS, SDSS and USNO in that order for coverage for reference photometry against which to do the calibration. This is achieved with gr_cat.py developed by Thomas Krühler which can be consulted for additional documentation. Sextractor is run on the astrometrically calibrated image using the function sextract, heavily inspired by autoastrometry.py by Daniel Perley and available at http://www.dark-cosmology.dk/~dperley/code/code.html. Handling of the entire sextractor interfacing is heavily based on autoastrometry.py. The two lists of images are then matched with a k-d tree algorithm and sextracted magntiudes can be calibrated against the chosen catalog.
    """

    fitsfile = fits.open(filename)
    header = fitsfile[0].header

    img_ra, img_dec = header["CRVAL1"], header["CRVAL2"]

    # temp_filename = filename
    temp_filename = filename.replace("fits", "")+"temp"

    # Get gain and readnoise
    try:
      gain_key = [x for x in header.keys() if "GAIN" in x][0]
      ron_key = [x for x in header.keys() if "RON" in x or "RDNOISE" in x][0]
      gain = header[gain_key]
      ron = header[ron_key]
    except:
      logger.warn("Gain and RON keys not understood. Setting to default values")
      gain = 2
      ron = 3.3

    if cosmic_rejection:
      # Clean for cosmics
      crmask, clean_arr = astroscrappy.detect_cosmics(fitsfile[0].data, gain=gain, readnoise=ron, sigclip=sigclip, objlim=objlim, cleantype='medmask', sepmed=True, verbose=True)

      # Replace data array with cleaned image
      fitsfile[0].data = clean_arr/gain

    # Save cosmicced file to temp
    fitsfile.writeto(temp_filename, output_verify='fix', clobber=True)

    # Attempt astrometric calibration
    if astrometry:
      temp_filename = run_astrometry_net(temp_filename, img_ra, img_dec)

    # Read in cosmic-ray rejected, possibly astrometrically calibrated image
    fitsfile = fits.open(temp_filename)
    header = fitsfile[0].header

    # Get header keyword for catalog matching
    if filter is None:
      try:
        img_filt = header["HIERARCH ESO INS FILT1 NAME"][0] # image filter name
      except KeyError:
        try:
          img_filt = header["FILTER"][0]
        except KeyError:
          try:
            img_filt = header["NCFLTNM2"][0]
          except KeyError:
            logger.warn("Filter keyword not recognized.", exc_info=1)
            sys.exit(1)

    img_ra, img_dec = header["CRVAL1"], header["CRVAL2"] # ra and dec
    # Ensure sign convention for gr_cat
    if not img_dec < 0:
        img_dec = "+"+str(img_dec)

    w = wcs.WCS(header)
    pixscale = wcs.utils.proj_plane_pixel_scales(w)
    nxpix = header['NAXIS1']
    nypix = header['NAXIS2']

    img_radius = np.sqrt((pixscale[0]*nxpix*60)**2 + (pixscale[1]*nypix*60)**2) # Largest image dimension to use as catalog query radius in arcmin

    # Get the catalog sources
    if img_filt == "I":
      # Get sdss filters for Lupton (2005) tranformations - http://www.sdss3.org/dr8/algorithms/sdssUBVRITransform.php
      cat_i = get_catalog(img_ra, img_dec, "i", catalog=catalog, radius = img_radius)
      cat_z = get_catalog(img_ra, img_dec, "z", catalog=catalog, radius = img_radius)
      cat_i, cat_z = joint_catalog(cat_i, cat_z) # Get joint catalog
      # Do filter transformation
      cat_i[:, 2] = cat_i[:, 2] - 0.3780*(cat_i[:, 2] - cat_z[:, 2]) - 0.3974
      # Account for transformation scatter
      cat_i[:, 3] = np.sqrt(cat_i[:, 3]**2 + 0.0063**2)
      cat = cat_i.copy()
    elif img_filt == "R":
      # Get sdss filters for Lupton (2005) tranformations - http://www.sdss3.org/dr8/algorithms/sdssUBVRITransform.php
      cat_r = get_catalog(img_ra, img_dec, "r", catalog=catalog, radius = img_radius)
      cat_i = get_catalog(img_ra, img_dec, "i", catalog=catalog, radius = img_radius)
      cat_r, cat_i = joint_catalog(cat_r, cat_i) # Get joint catalog
      # Do filter transformation
      cat_r[:, 2] = cat_r[:, 2] - 0.2936*(cat_r[:, 2] - cat_i[:, 2]) - 0.1439
      # Account for transformation scatter
      cat_r[:, 3] = np.sqrt(cat_r[:, 3]**2 + 0.0072**2)
      cat = cat_r.copy()
    else:
      cat = get_catalog(img_ra, img_dec, img_filt, catalog=catalog, radius = img_radius)

    print(cat)
    # Prepare sextractor
    writeparfile()
    saturation = 30000
    writeconfigfile(saturation)

    # Sextract stars to produce image star catalog
    goodsexlist = sextract(temp_filename, nxpix, nypix, border = 3, corner = 12, saturation=saturation)

    # Get sextracted ra, dec list for k-d Tree algoritm
    ra_sex, dec_sex = [], []
    for ii in goodsexlist:
        ra_sex.append(ii.ra)
        dec_sex.append(ii.dec)

    # Grow the tree
    tree_data = np.array([ra_sex, dec_sex]).T
    tree = cKDTree(tree_data)

    # Find mapping indices
    idx_map_sex, idx_map_cat = [], []
    tol = 1e-3 # Distance in degrees - This could change depending on the accuracy of the astrometric solution
    for ii, kk in enumerate(cat):
        # find the k nearest neighbours
        distance, indice = tree.query(kk[0:2], k=1)
        # print(distance, indice)
        if distance < tol:
            # Store
            idx_map_sex.append(indice)
            idx_map_cat.append(ii)

    # Add catalog photometry to sextractor object
    for ii, kk in enumerate(idx_map_sex):
        goodsexlist[kk].cat_mag = cat[idx_map_cat[ii]][2]
        goodsexlist[kk].cat_magerr = cat[idx_map_cat[ii]][3]

    # Bad matches from sextracted star catalog
    idx_bad = [ii for ii in np.arange(len(goodsexlist)) if ii not in idx_map_sex]

    # Remove mismatches
    for ii, kk in enumerate(idx_bad[::-1]):
        goodsexlist.pop(kk)

    # writetextfile('det.init.txt', goodsexlist)
    writeregionfile(temp_filename+'.det.im.reg', goodsexlist, 'red', 'img')

    # Get sextracted magnitudes and equivalent catalog magnitudes
    n_good = len(goodsexlist)
    mag, magerr, cat_mag, cat_magerr = np.zeros(n_good), np.zeros(n_good), np.zeros(n_good), np.zeros(n_good)
    for ii, kk in enumerate(goodsexlist):
        mag[ii] = kk.mag #+ 2.5*np.log10(exptime) # Correct for exposure time
        magerr[ii] = kk.magerr
        cat_mag[ii] = kk.cat_mag
        cat_magerr[ii] = kk.cat_magerr

    # Filter away 5-sigma outliers in the zero point
    zp = cat_mag - mag
    # print(zp)
    zp_l, zp_m, zp_h = np.percentile(zp, [16, 50, 84])

    sig_l = zp_m - zp_l
    sig_h = zp_h - zp_m
    # Filter zp's
    sigma_mask = 3
    mask = (zp > zp_m - sigma_mask * sig_l) & (zp < zp_m + sigma_mask * sig_h)
    zp = zp[mask]
    zp_m, zp_std = np.mean(zp), np.std(zp)
    zp_scatter = np.std(zp)
    print(np.mean(zp), np.std(zp), np.std(zp)/np.sqrt(len(zp)))

    # Fit for zero point
    from scipy.optimize import curve_fit
    from scipy import odr

    def func(p, x):
      b = p
      return x + b

    # Model object
    lin_model = odr.Model(func)

    # Create a RealData object
    data = odr.RealData(mag[mask], cat_mag[mask], sx=magerr[mask], sy=cat_magerr[mask])

    # Set up ODR with the model and data.
    odr = odr.ODR(data, lin_model, beta0=[np.mean(zp)])

    # Run the regression.
    out = odr.run()

    #print fit parameters and 1-sigma estimates
    popt = out.beta
    perr = out.sd_beta
    print('fit parameter 1-sigma error')
    print('———————————-')
    for i in range(len(popt)):
      print(str(popt[i])+' +- '+str(perr[i]))
    zp_m, zp_std = popt[0], perr[0]
    # prepare confidence level curves
    nstd = 5. # to draw 5-sigma intervals
    popt_up = popt + nstd * perr
    popt_dw = popt - nstd * perr

    x_fit = np.linspace(min(mag[mask]), max(mag[mask]), 100)
    fit = func(popt, x_fit)
    fit_up = func(popt_up, x_fit)
    fit_dw= func(popt_dw, x_fit)

    #plot
    pl.errorbar(mag[mask], cat_mag[mask], xerr=magerr[mask], yerr=cat_magerr[mask], fmt = 'k.', label = str(zp_m)+' +- '+str(zp_std))
    pl.plot(x_fit, fit, lw=2, label='best fit curve')
    pl.fill_between(x_fit, fit_up, fit_dw, alpha=.25, label='5-sigma interval')
    pl.legend()
    pl.savefig(filename+".pdf")
    pl.close()
    # Add catalog photometry to sextractor object
    for ii, kk in enumerate(goodsexlist):
        kk.cat_mag = mag[ii] + zp_m
        kk.cat_magerr = np.sqrt(magerr[ii]**2 + zp_std**2)
    writeregionfile(temp_filename+'.cal.im.reg', goodsexlist, 'red', 'img')

    # Get seeing fwhm for catalog object
    fwhm = np.zeros(len(goodsexlist))
    for ii, kk in enumerate(goodsexlist):
        fwhm[ii] = kk.fwhm

    # Filtered mean and std seeing FWHM in pixels
    l_fwhm, m_fwhm, h_fwhm = np.percentile(fwhm, [16, 50, 84])
    sig_l = m_fwhm - l_fwhm
    sig_h = h_fwhm - m_fwhm
    sigma_mask = 3
    mask = (fwhm > m_fwhm - sigma_mask * sig_l) & (fwhm < m_fwhm + sigma_mask * sig_h)
    fwhm = fwhm[mask]
    fwhm, fwhm_std = np.mean(fwhm), np.std(fwhm)

    # Median seeing in arcsec for sextractor
    seeing_fwhm = fwhm*pixscale[0] * 3600 # Seeing in arcsec
    # gain = 1e4
    try:
       # Sextract the image using the derived zero-point and fwhm!
       subprocess.run(['sex', '%s'%temp_filename, '-c', 'sex_temp.config', '-SEEING_FWHM', '%s'%seeing_fwhm, '-SATUR_LEVEL', '%s'%saturation, '-MAG_ZEROPOINT', '%s'%zp_m, '-CATALOG_NAME', 'temp_sex_obj.cat', '-GAIN', '%s'%gain, '-CHECKIMAGE_NAME', '%s_objfree.fits, %s_backrms.fits, %s_aper.fits'%(temp_filename, temp_filename, temp_filename), '-CHECKIMAGE_TYPE', '-OBJECTS, BACKGROUND_RMS, APERTURES', '-DETECT_THRESH', '3', '-BACK_SIZE', '64', '-BACK_FILTERSIZE', '3', '-DEBLEND_NTHRESH', '64', '-DEBLEND_MINCONT', '0.0001'])
    except (OSError, IOError):
       logger.warn("Sextractor failed to be executed.", exc_info=1)
       sys.exit(1)



    # From sextractors background rms image, get variance
    back_rms_image = fits.open("%s_backrms.fits"%temp_filename)
    l_rms, m_rms, h_rms = np.percentile(back_rms_image[0].data, [16, 50, 84])
    sig_l = m_rms - l_rms
    sig_h = h_rms - m_rms
    sigma_mask = 3
    mask = (back_rms_image[0].data > m_rms - sigma_mask * sig_l) & (back_rms_image[0].data < m_rms + sigma_mask * sig_h)
    back_rms_image[0].data = back_rms_image[0].data[mask]
    rms, rms_std = np.mean(back_rms_image[0].data), np.std(back_rms_image[0].data)



    lim_mag = limiting_magnitude(img_rms = rms, img_fwhm = fwhm, img_zp = zp_m, sigma_limit = 5)
    print("Limiting magnitude")
    print(lim_mag)

    fin_img = fits.open('%s_aper.fits'%temp_filename)
    fin_img[0].header["LIMMAG"] = lim_mag[0]
    fin_img.writeto('%s_calibrated.fits'%temp_filename, clobber = True)

    # Read in the sextractor catalog
    try:
       cat = open("temp_sex_obj.cat",'r')
       catlines = cat.readlines()
       cat.close()
    except:
        logger.warn("Cannot load sextractor output file!", exc_info=1)
        sys.exit(1)

    if len(catlines) == 0:
        logger.warn("Sextractor catalog is empty: try a different catalog?", exc_info=1)
        sys.exit(1)

    l = -1
    sexlist = []
    while l < len(catlines)-1:
        l += 1
        if (len(catlines[l]) <= 1 or catlines[l][0] == '#'):
            continue
        iobj = SexObj(catlines[l]) #process the line into an object
        sexlist.append(iobj)

    for ii, kk in enumerate(sexlist):
        if kk.mag <= lim_mag:
          kk.cat_mag = kk.mag
          kk.cat_magerr = np.sqrt(kk.magerr**2 + zp_std**2)
        elif kk.mag > lim_mag:
          kk.cat_mag = lim_mag
          kk.cat_magerr = 9.99
        else:
          sys.exit(1)
    writeregionfile(temp_filename+'.obj.im.reg', sexlist, 'red', 'img')

    try:
        for fl in glob.glob("*temp*"):
            os.remove(fl)
    except:
       print('Could not remove temp files for some reason')


    try:
        for fl in glob.glob("*temp*"):
            os.remove(fl)
    except:
       print('Could not remove temp files for some reason')


def main():



    # gfilelist = glob.glob("/Users/jselsing/Dropbox/SN2017eaw_PHOT/ALFOSC/*g0*.fits")
    # rfilelist = glob.glob("/Users/jselsing/Dropbox/SN2017eaw_PHOT/ALFOSC/*r0*.fits")
    # ifilelist = glob.glob("/Users/jselsing/Dropbox/SN2017eaw_PHOT/ALFOSC/*i0*.fits")
    filelist = glob.glob("/Users/jselsing/Dropbox/SN2017eaw_PHOT/REDUCED_PHOT/*.fits")
    # filelist = gfilelist + rfilelist + ifilelist + zfilelist
    # print(filelist)
    # exit()
    for ii in filelist:
      autocal(filename = ii, catalog = "2MASS", sigclip = 50, objlim = 75, cosmic_rejection = False, astrometry = True)


if __name__ == '__main__':
    main()