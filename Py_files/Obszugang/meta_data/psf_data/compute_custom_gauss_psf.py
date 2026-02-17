"""
create all needed data for the HST PSFs
"""
import numpy as np
import os
from astropy.io import fits
# from phangs_data_access import phot_tools, phangs_info

from werkzeugkiste import helper_func, phot_tools
from obszugang import obs_info
from malkasten import plotting_tools

import matplotlib.pyplot as plt
import pickle
from astropy.convolution import Gaussian2DKernel
from astropy.convolution import convolve


# these are the parameters from the original creation of the convolution which is ultimately the
pixel_scale_arcsec = (0.47435208445882804 + 0.47165711155398005) / 2
orig_distance = 0.84
degrade_distance = 10.
target = 'm33'
obs = 'hst'

distance_ratio = degrade_distance/orig_distance
std_orig_img_pix = np.sqrt((distance_ratio*2.2)**2 - 2.2**2)/2.35
print('pixel_scale_arcsec ', pixel_scale_arcsec)
print('std_orig_img_pix ', std_orig_img_pix)


std_psf_pix = std_orig_img_pix / distance_ratio
std_psf_arcsec = std_psf_pix * pixel_scale_arcsec
print('std_psf_pix ', std_psf_pix)
print('std_psf_arcsec ', std_psf_arcsec)

# create kernel at original size and oversampled
kernel_size = 73
oversample_fact = 4
# define the scales
oversampled_pixel_scale = pixel_scale_arcsec / oversample_fact
std_psf_oversampled_pix = std_psf_pix * oversample_fact
std_psf_oversampled_arcsec = std_psf_arcsec * oversample_fact
print('oversampled_pixel_scale ', oversampled_pixel_scale)
print('std_psf_oversampled_pix ', std_psf_oversampled_pix)
print('std_psf_oversampled_arcsec ', std_psf_oversampled_arcsec)


# create the kernael that is our PSF
kernel = Gaussian2DKernel(x_stddev=std_psf_pix, x_size=kernel_size, y_size=kernel_size)
kernel_oversampled = Gaussian2DKernel(x_stddev=std_psf_oversampled_pix, x_size=kernel_size * oversample_fact, y_size=kernel_size * oversample_fact)
psf_array = kernel.array
psf_array_oversampled = kernel_oversampled.array


# encirceled energy values
ee_values = [0.25, 0.30, 0.35, 0.40, 0.45, 0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90, 0.95, 0.99]
ee_str = ['25', '30', '35', '40', '45', '50', '55', '60', '65', '70', '75', '80', '85', '90', '95', '99']

central_x_pos_oversampled = psf_array_oversampled.shape[0]/2
central_y_pos_oversampled = psf_array_oversampled.shape[1]/2
max_rad_oversampled = np.min(psf_array_oversampled.shape) / 2

# get radial profile stats:
rad_profile_stat_dict_oversampled = phot_tools.ProfileTools.get_rad_profile(
    data=psf_array_oversampled, x_pos=central_x_pos_oversampled, y_pos=central_y_pos_oversampled, max_rad=max_rad_oversampled, err=None,
    pix_scale=oversampled_pixel_scale, method='exact')


# plt.plot(rad_profile_stat_dict_oversampled['rad'], rad_profile_stat_dict_oversampled['profile'] / max(rad_profile_stat_dict_oversampled['profile']), color='green')
# plt.show()
# exit()

ee_values_arcsec = phot_tools.ProfileTools.get_src_ee(data=psf_array_oversampled, x_pos=central_x_pos_oversampled, y_pos=central_y_pos_oversampled,
                                            max_rad=max_rad_oversampled, ee_values=ee_values,
                                            pix_scale=oversampled_pixel_scale,
                                            err=None)
ee_values_pix= phot_tools.ProfileTools.get_src_ee(data=psf_array_oversampled, x_pos=central_x_pos_oversampled, y_pos=central_y_pos_oversampled,
                                            max_rad=max_rad_oversampled, ee_values=ee_values,
                                            pix_scale=1/oversample_fact,
                                            err=None)

print(ee_values_arcsec[5], ee_values_arcsec[11], rad_profile_stat_dict_oversampled['gaussian_fwhm'])
psf_dict_custom_gauss = {

    'F275W': {
        # the PSF itself
        'over_sampled_psf': psf_array_oversampled,
        'pixel_scale_psf_over_sampled': oversampled_pixel_scale,
        'n_over_sampled': oversample_fact,
        # parametrization of the radial profile
        'radius_arcsec': rad_profile_stat_dict_oversampled['rad'],
        'psf_profile': rad_profile_stat_dict_oversampled['profile'],
        'gaussian_profile': rad_profile_stat_dict_oversampled['gaussian_profile'],
        'gaussian_fwhm': rad_profile_stat_dict_oversampled['gaussian_fwhm'],
        'gaussian_amp': rad_profile_stat_dict_oversampled['gaussian_amp'],
        'gaussian_mean': rad_profile_stat_dict_oversampled['gaussian_mean'],
        'gaussian_std': rad_profile_stat_dict_oversampled['gaussian_std'],
    },

    'F336W': {
        # the PSF itself
        'over_sampled_psf': psf_array_oversampled,
        'pixel_scale_psf_over_sampled': oversampled_pixel_scale,
        'n_over_sampled': oversample_fact,
        # parametrization of the radial profile
        'radius_arcsec': rad_profile_stat_dict_oversampled['rad'],
        'psf_profile': rad_profile_stat_dict_oversampled['profile'],
        'gaussian_profile': rad_profile_stat_dict_oversampled['gaussian_profile'],
        'gaussian_fwhm': rad_profile_stat_dict_oversampled['gaussian_fwhm'],
        'gaussian_amp': rad_profile_stat_dict_oversampled['gaussian_amp'],
        'gaussian_mean': rad_profile_stat_dict_oversampled['gaussian_mean'],
        'gaussian_std': rad_profile_stat_dict_oversampled['gaussian_std'],
    },

    'F475W': {
        # the PSF itself
        'over_sampled_psf': psf_array_oversampled,
        'pixel_scale_psf_over_sampled': oversampled_pixel_scale,
        'n_over_sampled': oversample_fact,
        # parametrization of the radial profile
        'radius_arcsec': rad_profile_stat_dict_oversampled['rad'],
        'psf_profile': rad_profile_stat_dict_oversampled['profile'],
        'gaussian_profile': rad_profile_stat_dict_oversampled['gaussian_profile'],
        'gaussian_fwhm': rad_profile_stat_dict_oversampled['gaussian_fwhm'],
        'gaussian_amp': rad_profile_stat_dict_oversampled['gaussian_amp'],
        'gaussian_mean': rad_profile_stat_dict_oversampled['gaussian_mean'],
        'gaussian_std': rad_profile_stat_dict_oversampled['gaussian_std'],
    },

    'F814W': {
        # the PSF itself
        'over_sampled_psf': psf_array_oversampled,
        'pixel_scale_psf_over_sampled': oversampled_pixel_scale,
        'n_over_sampled': oversample_fact,
        # parametrization of the radial profile
        'radius_arcsec': rad_profile_stat_dict_oversampled['rad'],
        'psf_profile': rad_profile_stat_dict_oversampled['profile'],
        'gaussian_profile': rad_profile_stat_dict_oversampled['gaussian_profile'],
        'gaussian_fwhm': rad_profile_stat_dict_oversampled['gaussian_fwhm'],
        'gaussian_amp': rad_profile_stat_dict_oversampled['gaussian_amp'],
        'gaussian_mean': rad_profile_stat_dict_oversampled['gaussian_mean'],
        'gaussian_std': rad_profile_stat_dict_oversampled['gaussian_std'],
    }
}

# encircled energy values
for idx_ee, ee in enumerate(ee_str):
    psf_dict_custom_gauss.update({'ee_%s_arcsec' % ee: ee_values_arcsec[idx_ee]})
    psf_dict_custom_gauss.update({'ee_%s_pix' % ee: ee_values_pix[idx_ee]})

# save dictionary
if not os.path.isdir('data_output'):
    os.makedirs('data_output')

with open('data_output/custom_gauss_psf_dict.pickle', 'wb') as file_name:
    pickle.dump(psf_dict_custom_gauss, file_name)


