#!/usr/local/anaconda3/envs/py36 python
# -*- coding: utf-8 -*-

def limiting_magnitude(img_rms = 100, img_fwhm = 5, img_zp = 30, sigma_limit = 5, profile = "Moffat", return_image = False):
    """
    Small "simulation-like" calculation to get limiting magnitude of image given Background RMS, PSF - FWHM, and zero-point
    """
    from astropy.modeling import models
    from photutils import CircularAperture, aperture_photometry
    import numpy as np
    # Generate image parameters
    tmp_amplitude = 1000

    aperture_radius = img_fwhm * (66/100)
    img_size = int(np.ceil(10 * img_fwhm / 2.) * 2)
    source_pos = [img_size/2, img_size/2]

    # Make synthetic sky with similar noise characteristics as observed image
    coviance = [[img_rms**2, 0], [0, img_rms**2]]
    sky = np.random.multivariate_normal([0, 0], coviance, [img_size, img_size])[:, :, 0]


    # Simulate source as either Gaussian or Moffat
    y, x = np.mgrid[:img_size, :img_size]
    if profile == "Moffat":
        beta = 4.765
        gamma = img_fwhm / (2 * np.sqrt(2**(1/beta) - 1))
        source = models.Moffat2D.evaluate(x, y, tmp_amplitude, source_pos[0], source_pos[1], gamma, beta)
    elif profile == "Gaussian":
        img_sigma = img_fwhm / 2.35
        source = models.Gaussian2D.evaluate(x, y, tmp_amplitude, source_pos[0], source_pos[1], img_sigma, img_sigma, 0)
    else:
        print("Exiting ... Profile needs to be either Moffat or Gaussian")
        sys.exit(1)

    # Define aperture
    apertures = CircularAperture(source_pos, r=aperture_radius)

    # Do aperture photmetry of both source and sky
    source_counts = aperture_photometry(source, apertures)["aperture_sum"].data
    error_counts = aperture_photometry(sky**2, apertures)["aperture_sum"].data

    # Calculate S/N of detection
    SN = source_counts / np.sqrt(error_counts)

    # Rescale source intensity to be detected at N-sigma level
    rescaled_source = source * (sigma_limit / SN)

    # Get magnitude of N-sigma detected source
    source_counts = aperture_photometry(source, apertures)["aperture_sum"].data

    # "Perfect" Aperture correction
    aperture_correction = -2.5*np.log10(np.sum(rescaled_source)/source_counts)

    # Limiting magnitude
    magnitude_limit = -2.5*np.log10(source_counts) + img_zp + aperture_correction

    if return_image:
        return rescaled_source + sky
    elif not return_image:
        return magnitude_limit



def main():

    FORSz = limiting_magnitude(img_rms = 40.74, img_fwhm = 2.74, img_zp = 32.36, sigma_limit=5, profile="Gaussian", return_image=True)

    pl.imshow(FORSz, cmap = "viridis")
    pl.show()


if __name__ == '__main__':
    main()