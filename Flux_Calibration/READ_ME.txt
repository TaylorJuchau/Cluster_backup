This READ_ME is for the Functions.py file and how to use the included functions. See "flux_calibration/Example_script.ipynb" for flushed out examples

First, I will assume that the .py file is imported as follows, simply copy the chunk to the top of your script:


import os
os.chdir('/d/ret1/Taylor/jupyter_notebooks/Research') #TJ change working directory to be the parent directory
from Py_files.Functions import *



#################Major functions:#####################

**************************************************************************************************************************************************************
results = compare_photometry(list_of_IFUs, image_files, filter_files, location, radius)
    Generates a dictionary with the following keys:
    'add_datasets' - list of structured arrays, each entry is its own structured array, with results['add_datasets'][0]['wavelength'] being the wavelength array for the first cube in the list_of_IFUs, and results['add_datasets'][0]['intensity'] being the flux density at each wavelength. There is expected to be some overlap between results['add_datasets'][0]['wavelength'] and results['add_datasets'][1]['wavelength']
    'mult_datasets' - list of structured arrays, same format as add_datasets
    'ifu_files' - list of ifu files used in the dataset, in case you need to reload them when analyzing the results
    'image_files' - list of image files used in the dataset
    'location' - skycoord object for the location of the center of the aperture
    'loc_idx' - index for the location in the list of the 4 locations 
        [[202.5062429, 47.2143358], [202.4335225, 47.1729608], [202.4340450, 47.1732517], [202.4823742, 47.1958589]]
    'radius' - radius of aperture with units of angular size (for example 1.25*u.arcsec)
    'add_correction_values' - list of correction values used for additive correction, zero if the cube was not adjusted, has same length as datasets
    'mult_correction_values' - list of correction values used for multiplicative correction, unity if the cube was not adjusted
    'filter_names' - list of filter names used for the data
    'filter_wavelengths' - mean wavelengths for the filters, in the same order as filter_names
    'add_synthetic_fluxes' - total fluxes calculated for the additive correction methods (used to compare to image fluxes)
    'mult_synthetic_fluxes' - total fluxes calculated for the multiplicative correction methods (used to compare to image fluxes)
    'photo_fluxes' - total fluxes calculated just from summing over the aperture in the image files
    'wavelength' - array of wavelengths (used for plotting)
    'add_intensity' - array of intensities used in the additive method (for plotting)
    'mult_intensity' - array of intensities used in the multiplicative method (for plotting)
**************************************************************************************************************************************************************
plot_results(results, correction = 'mult', show_images = [], color_min_max = [1, 99.5])
    returns None
    Plots the synthetic spectrum and comparison to image data. 
    show_images takes a list of filter_names (max of two) to show both the synthetic and raw images on the plot for clarifying drastic discrepencies.
**************************************************************************************************************************************************************
syntehtic_flux = get_Fnu_transmission(Fnu_array, wl_array, transmission_array, trans_wl_array, warnings = True)
    applies a synthetic filter to a given flux vs wavelength array.
**************************************************************************************************************************************************************
corrected_data, correction_value_used = adjust_spectrum(ifu_file, get_largest_filter_within(ifu_file), image_files, location, radius, adjustment_operation = 'add')
    Used to generate new dataset of flux density as a function of wavelength. Flux densities are added across all pixels in the aperture (radius and center location given as arguments). Flux levels are corrected to match the photometric flux in the image_files. If get_largest_filter_within(original_ifu) argument is None, then the original IFU spectrum is returned without correction.
**************************************************************************************************************************************************************
largest_filter = get_largest_filter_within(IFU_file)
    Used to determine the largest filter name that is entirely within an IFU's wavelength coverage. Returns None if there isnt a single filter that is entirely within the coverage. 
**************************************************************************************************************************************************************
psf = generate_psf(psf_type, fwhm_pix=None, size=51, psf_array=None)
    Used to create a psf to convolve a fits file with. 
    psf_type can be either 'gaussian' or 'array' (the latter is a custom array that you supply yourself)
    fwhm_pix is the pixel number of the FWHM value (note that this is NOT the angular size, but the pixel number)
    size is how large of a convolution array to use (higher number is higher resolution array)
    psf_array is only needed if psf_type is specified as array
**************************************************************************************************************************************************************
output_filepath = convolve_ifu_cube(ifu_fits_path, psf, output_path, ext='SCI')
    Convolve an ifu file with a psf (created using the "generate_psf" function), save output to filepath prescribed

###############Minor functions########################

**************************************************************************************************************************************************************
filter_name = extract_filter_name(file_name)
    Used to extract the upper case filter name from a long file name. This keys on a single upper or lower case "f", followed by some numbers, concluded by a single letter that is an upper or lower case "W", "M", or "N". Returns the first instance of this match.
    For example: 'Data_files/Image_files/ngc5194_nircam_lv3_f560_f115w_i2d.fits' will return "F115W"
    However, the file 'Data_files/Image_files/ngcf5194_f21nircam_lv3_f115w_i2d.fits' will return "F21N"
**************************************************************************************************************************************************************
filter_wavelength, filter_transmission = get_filter_data(filter_name filter_file_list = filter_files)
    Used to generate transmission vs wavelength array for a given filter. Uses the filter files at '/d/crow1/tools/cigale/database_builder/filters/jwst/' by default. If you use another filter file format you will need to redo this function, it relies on the file having a very particular format.
**************************************************************************************************************************************************************





