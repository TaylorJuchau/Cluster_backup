"""
Script to compute profile fitting correction values for HST
"""

import os.path
import pickle

import numpy as np
import matplotlib.pyplot as plt

from werkzeugkiste import helper_func, phot_tools
from obszugang import obs_info
from malkasten import plotting_tools

from astropy.convolution import Gaussian2DKernel
from astropy.convolution import convolve


fact_img_size = 5

degrade_distance = 10.
target = 'm33'
obs = 'hst'

band_suffix = '%s_%s_%s_mpc' % (target, obs, degrade_distance)

band_list = [
    'F275W',
    'F336W',
    'F475W',
    'F814W',
]



custom_gauss_apert_gauss_corr_dict = {}

for band in band_list:
    print('band: ', band)
    psf_dict = phot_tools.PSFTools.load_obs_psf_dict(band=band, instrument='custom_gauss')

    psf_data = psf_dict['over_sampled_psf']
    psf_fwhm_pix = psf_dict['gaussian_fwhm'] / psf_dict['pixel_scale_psf_over_sampled']
    psf_std_pix = psf_dict['gaussian_std'] / psf_dict['pixel_scale_psf_over_sampled']
    print('psf_std_arcsec ', psf_dict['gaussian_std'])
    print('psf_std_pix ', psf_std_pix)

    apert_radii_arcsec_list = np.linspace(psf_dict['gaussian_std'], 4.0, 50)
    apert_radii_pix_list = apert_radii_arcsec_list / psf_dict['pixel_scale_psf_over_sampled']


    # get maximal convolution
    max_convolution_for_psf_img = np.min(psf_data.shape)/5
    # check if this is feasible:
    if max_convolution_for_psf_img > (10 * psf_std_pix):
        max_convolution = 10 * psf_std_pix
    else:
        max_convolution = max_convolution_for_psf_img

    list_convolutions = np.linspace(0, max_convolution, 15)
    list_convolutions[0] = list_convolutions[1]/100
    print('list_convolutions ', list_convolutions)

    measured_std_values = np.zeros(len(list_convolutions))
    measured_flux_frac_in_apert = np.zeros((len(apert_radii_arcsec_list), len(list_convolutions)))


    for conv_idx, conv in enumerate(list_convolutions):
        print('conv_idx, conv ', conv_idx, conv)

        # now get a gaussian function as a primarry source
        x_stddev = conv
        y_stddev = conv
        theta = 0

        data_width = np.rint(np.max((conv, psf_std_pix)) * fact_img_size)
        n_pixels = data_width * 2 + 1
        data_width_pix = data_width
        x_bins = np.linspace(-data_width_pix, data_width_pix, int(n_pixels))
        y_bins = np.linspace(-data_width_pix, data_width_pix, int(n_pixels))
        x_mesh, y_mesh = np.meshgrid(x_bins, y_bins)  # get 2D variables instead of 1D
        x0 = (x_bins[0] + x_bins[-1]) / 2
        y0 = (y_bins[0] + y_bins[-1]) / 2

        moffat_data = helper_func.FuncAndModels.star_cluster_moffat2d(
            x=x_mesh, y=y_mesh, x0=x0, y0=y0, mu_0=1, nu=1.3, fwhm=conv*2)

        # gauss_data = helper_func.FuncAndModels.gauss2d_rot(
        #     x=x_mesh, y=y_mesh, amp=1, x0=x0, y0=y0, sig_x=conv, sig_y=conv, theta=theta)
        data_convolve = convolve(psf_data, moffat_data)
        # norm to 1
        data_convolve = data_convolve/np.sum(data_convolve)

        data_err = np.ones(data_convolve.shape) * np.max(data_convolve) / 100

        x_center_pix = data_convolve.shape[0] / 2
        y_center_pix = data_convolve.shape[1] / 2

        # measure photometry in all radii
        for apert_idx in range(len(apert_radii_pix_list)):
            apert_stats = phot_tools.ApertTools.get_apert_stats(
            data=data_convolve, data_err=data_err, x_pos=x_center_pix, y_pos=y_center_pix,
            aperture_rad_pix=apert_radii_pix_list[apert_idx], mask=None, sig_clip=None, sum_method='exact')

            measured_flux_frac_in_apert[apert_idx, conv_idx] = apert_stats.sum

        profile_dict = phot_tools.ProfileTools.compute_axis_profiles(data=data_convolve,
                                                                     x_pos=x_center_pix,
                                                                     y_pos=y_center_pix, n_slits=12,
                                                                     err=data_err)

        amp_list, mu_list, sig_list, amp_err_list, mu_err_list, sig_err_list = \
            phot_tools.ProfileTools.fit_gauss2rad_profiles(rad_profile_dict=profile_dict, std_pix=psf_std_pix,
                                                upper_sig_fact=10, central_rad_fact=3)
        amp_list = np.array(amp_list)
        mu_list = np.array(mu_list)
        sig_list = np.array(sig_list)
        amp_err_list = np.array(amp_err_list)
        sig_err_list = np.array(sig_err_list)

        # get all the gaussian functions that make sense
        # then need to be central
        mask_mu = np.abs(mu_list) < psf_std_pix * 3
        # they must have a positive amplitude
        mask_amp = (amp_list > 0)
        mean_sigma = np.nanmean(sig_list[mask_mu * mask_amp]) * psf_dict['pixel_scale_psf_over_sampled']

        measured_std_values[conv_idx] = mean_sigma

    custom_gauss_apert_gauss_corr_dict.update(
        {
            band:
                {
                    'apert_radii_arcsec_list': apert_radii_arcsec_list,
                    'measured_std_values': measured_std_values,
                    'measured_flux_frac_in_apert': measured_flux_frac_in_apert
                }})



    # plot_the output
    # interpolate function
    interp_func = helper_func.InterpTools.interp2dgrid(
        x_bins=apert_radii_arcsec_list, y_bins=measured_std_values,func_values=measured_flux_frac_in_apert,
        method='cubic')

    fine_interpolated = helper_func.InterpTools.get2dinterp_fine_grid(
        interp_func=interp_func, x_min=np.min(apert_radii_arcsec_list), x_max=np.max(apert_radii_arcsec_list),
        y_min=np.min(measured_std_values), y_max=np.max(measured_std_values), n_x_bins=100, n_y_bins=100)

    fig = plt.figure(figsize=(30, 12))
    fontsize = 30
    ax_corr_val = fig.add_axes((0.07, 0.07, 0.8, 0.92))
    ax_cbar = fig.add_axes((0.9, 0.3, 0.02, 0.4))
    norm = plotting_tools.ColorBarTools.compute_cbar_norm(vmin_vmax=(np.min(measured_flux_frac_in_apert), np.max(measured_flux_frac_in_apert)))
    ax_corr_val.imshow(fine_interpolated,
               extent=(np.min(apert_radii_arcsec_list), np.max(apert_radii_arcsec_list),
                       np.min(measured_std_values), np.max(measured_std_values)),
               origin='lower',
               norm=norm, cmap='plasma')
    # plot colorbar
    plotting_tools.ColorBarTools.create_cbar(ax_cbar=ax_cbar, cmap='plasma', norm=norm, cbar_label='EE in Aperture', fontsize=fontsize, ticks=None, labelpad=2, tick_width=2,
                    orientation='vertical', top_lable=True, label_color='k',
                    extend='neither')

    # plot also the dots from where we got the interpolation
    x_mesh, y_mesh = np.meshgrid(apert_radii_arcsec_list, measured_std_values)

    ax_corr_val.scatter(x_mesh, y_mesh, c=measured_flux_frac_in_apert.T, norm=norm, cmap='plasma', edgecolors='k')
    ax_corr_val.set_xlabel('Aperture [arcsec]', fontsize=fontsize,)
    ax_corr_val.set_ylabel(r'Measured Gauss. $\sigma$ [arcsec]', fontsize=fontsize,)
    ax_corr_val.tick_params(axis='both', which='both', width=3, length=15, right=False, top=True, direction='in',
                   labelsize=fontsize)
    if not os.path.isdir('plot_output/apert_profile_corr'): os.makedirs('plot_output/apert_profile_corr')
    fig.savefig('plot_output/apert_profile_corr/apert_profile_corr_%s.png' % band)
    plt.clf()
    plt.cla()

    # # lastly check that everything has correct values
    # print('x_val=0.25, y_val=0.1 ', helper_func.InterpTools.get2dinterp_value(interp_func=interp_func, x_val=0.25, y_val=0.1))
    # print('x_val=0.40, y_val=0.1 ', helper_func.InterpTools.get2dinterp_value(interp_func=interp_func, x_val=0.40, y_val=0.1))
    # print('x_val=0.59, y_val=0.1 ', helper_func.InterpTools.get2dinterp_value(interp_func=interp_func, x_val=0.59, y_val=0.1))
    #
    # print('x_val=0.25, y_val=0.2 ', helper_func.InterpTools.get2dinterp_value(interp_func=interp_func, x_val=0.25, y_val=0.25))
    # print('x_val=0.25, y_val=0.3 ', helper_func.InterpTools.get2dinterp_value(interp_func=interp_func, x_val=0.25, y_val=0.35))
    # print('x_val=0.25, y_val=0.4 ', helper_func.InterpTools.get2dinterp_value(interp_func=interp_func, x_val=0.25, y_val=0.41))


# # save dictionary
# if not os.path.isdir('data_output/apert_gauss_corr'):
#     os.makedirs('data_output/apert_gauss_corr')
#
with open('data_output/apert_gauss_corr/custom_gauss_apert_gauss_corr_dict.pickle', 'wb') as file_name:
    pickle.dump(custom_gauss_apert_gauss_corr_dict, file_name)






















#
# original_pixel_size_acs_arcsec = 0.05
# original_pixel_size_uvis_arcsec = 0.039
#
# orig_distance = 0.84
# degrade_distance = 10.
# distance_ratio = degrade_distance/orig_distance
# std_orig_img_pix = np.sqrt((distance_ratio*2.2)**2 - 2.2**2)/2.35
#
# std_reprojected_pix = std_orig_img_pix / distance_ratio
# std_reprojected_acs_arcsec = std_reprojected_pix * original_pixel_size_acs_arcsec
# std_reprojected_uvis_arcsec = std_reprojected_pix * original_pixel_size_uvis_arcsec
#
# print('std_orig_img_pix ', std_orig_img_pix)
# print('std_reprojected_pix ', std_reprojected_pix)
#
# print('std_reprojected_acs_arcsec ', std_reprojected_acs_arcsec)
# print('std_reprojected_uvis_arcsec ', std_reprojected_uvis_arcsec)
#
#
# oversample_fact = 4
# oversampled_pixelscale_acs = original_pixel_size_acs_arcsec / oversample_fact
# oversampled_pixelscale_uvis = original_pixel_size_uvis_arcsec / oversample_fact
# # create the kernael that is our PSF
# kernel = Gaussian2DKernel(x_stddev=std_reprojected_pix * oversample_fact, x_size=250, y_size=250)
#
# psf_array_oversampled = kernel.array
#
#
# # print(kernel.__dict__)
# # plt.imshow(np.log10(psf_array_oversampled))
# # plt.show()
#
#
# # astropy_conv = convolve_fft(hdup.data, kernel, allow_huge=True)
#
# # acs_wfc_band_list = obs_info.acs_wfc_psf_band_list
# # wfc3_uv_band_list = obs_info.wfc3_uv_psf_band_list
# # wfc3_ir_band_list = obs_info.wfc3_ir_psf_band_list
# #
# # acs_wfc_fact_img_size = 5
# # wfc3_uv_fact_img_size = 5
# # wfc3_ir_fact_img_size = 5
# #
#
# acs_wfc_apert_gauss_corr_dict = {}
#
# # psf_dict = phot_tools.PSFTools.load_obs_psf_dict(band=band, instrument='acs')
#
#
#
# # psf_fwhm_pix = psf_dict['gaussian_fwhm'] / psf_dict['pixel_scale_psf_over_sampled']
# # psf_std_pix = psf_dict['gaussian_std'] / psf_dict['pixel_scale_psf_over_sampled']
# # print('psf_std_arcsec ', psf_dict['gaussian_std'])
#
# psf_std_oversampeled_pix = std_reprojected_pix * oversample_fact
# print('psf_std_oversampeled_pix ', psf_std_oversampeled_pix)
#
# apert_radii_acs_arcsec_list = np.linspace(std_reprojected_acs_arcsec, 1.0, 50)
# apert_radii_uvis_arcsec_list = np.linspace(std_reprojected_uvis_arcsec, 1.0, 50)
# apert_radii_acs_pix_list = apert_radii_acs_arcsec_list / oversampled_pixelscale_acs
# apert_radii_uvis_pix_list = apert_radii_uvis_arcsec_list / oversampled_pixelscale_uvis
#
# print('apert_radii_acs_arcsec_list ', apert_radii_acs_arcsec_list)
# print('apert_radii_uvis_arcsec_list ', apert_radii_uvis_arcsec_list)
# print('apert_radii_acs_pix_list ', apert_radii_acs_pix_list)
# print('apert_radii_uvis_pix_list ', apert_radii_uvis_pix_list)
#
#
# # get maximal convolution
# max_convolution_for_psf_img_pix = np.min(psf_array_oversampled.shape)/5
# print('max_convolution_for_psf_img_pix ', max_convolution_for_psf_img_pix)
#
#
# # check if this is feasible:
# if max_convolution_for_psf_img_pix > (10 * psf_std_oversampeled_pix):
#     max_convolution = 10 * psf_std_oversampeled_pix
# else:
#     max_convolution = max_convolution_for_psf_img_pix
#
# list_convolutions = np.linspace(0, max_convolution, 15)
# list_convolutions[0] = list_convolutions[1]/100
# print('list_convolutions ', list_convolutions)
#
# measured_std_values_acs = np.zeros(len(list_convolutions))
# measured_std_values_uvis = np.zeros(len(list_convolutions))
# measured_flux_frac_in_apert_acs = np.zeros((len(apert_radii_acs_arcsec_list), len(list_convolutions)))
# measured_flux_frac_in_apert_uvis = np.zeros((len(apert_radii_uvis_arcsec_list), len(list_convolutions)))
#
#
# for conv_idx, conv in enumerate(list_convolutions):
#     print('conv_idx, conv ', conv_idx, conv)
#
#     # now get a gaussian function as a primarry source
#     x_stddev = conv
#     y_stddev = conv
#     theta = 0
#
#     data_width = np.rint(np.max((conv, std_reprojected_pix)) * oversample_fact)
#     n_pixels = data_width * 2 + 1
#     data_width_pix = data_width
#     x_bins = np.linspace(-data_width_pix, data_width_pix, int(n_pixels))
#     y_bins = np.linspace(-data_width_pix, data_width_pix, int(n_pixels))
#     x_mesh, y_mesh = np.meshgrid(x_bins, y_bins)  # get 2D variables instead of 1D
#     x0 = (x_bins[0] + x_bins[-1]) / 2
#     y0 = (y_bins[0] + y_bins[-1]) / 2
#
#     moffat_data = helper_func.FuncAndModels.star_cluster_moffat2d(
#         x=x_mesh, y=y_mesh, x0=x0, y0=y0, mu_0=1, nu=1.3, fwhm=conv*2)
#
#     # gauss_data = helper_func.FuncAndModels.gauss2d_rot(
#     #     x=x_mesh, y=y_mesh, amp=1, x0=x0, y0=y0, sig_x=conv, sig_y=conv, theta=theta)
#     data_convolve = convolve(psf_array_oversampled, moffat_data)
#     # norm to 1
#     data_convolve = data_convolve/np.sum(data_convolve)
#
#     data_err = np.ones(data_convolve.shape) * np.max(data_convolve) / 100
#
#     x_center_pix = data_convolve.shape[0] / 2
#     y_center_pix = data_convolve.shape[1] / 2
#
#     # measure photometry in all radii
#     for apert_idx in range(len(apert_radii_acs_pix_list)):
#         apert_stats = phot_tools.ApertTools.get_apert_stats(
#         data=data_convolve, data_err=data_err, x_pos=x_center_pix, y_pos=y_center_pix,
#         aperture_rad_pix=apert_radii_acs_pix_list[apert_idx], mask=None, sig_clip=None, sum_method='exact')
#
#         measured_flux_frac_in_apert_acs[apert_idx, conv_idx] = apert_stats.sum
#     for apert_idx in range(len(apert_radii_uvis_pix_list)):
#         apert_stats = phot_tools.ApertTools.get_apert_stats(
#         data=data_convolve, data_err=data_err, x_pos=x_center_pix, y_pos=y_center_pix,
#         aperture_rad_pix=apert_radii_uvis_pix_list[apert_idx], mask=None, sig_clip=None, sum_method='exact')
#
#         measured_flux_frac_in_apert_uvis[apert_idx, conv_idx] = apert_stats.sum
#
#
#     profile_dict = phot_tools.ProfileTools.compute_axis_profiles(data=data_convolve,
#                                                                  x_pos=x_center_pix,
#                                                                  y_pos=y_center_pix, n_slits=12,
#                                                                  err=data_err)
#
#     amp_list, mu_list, sig_list, amp_err_list, mu_err_list, sig_err_list = \
#         phot_tools.ProfileTools.fit_gauss2rad_profiles(rad_profile_dict=profile_dict, std_pix=psf_std_oversampeled_pix,
#                                             upper_sig_fact=10, central_rad_fact=3)
#     amp_list = np.array(amp_list)
#     mu_list = np.array(mu_list)
#     sig_list = np.array(sig_list)
#     amp_err_list = np.array(amp_err_list)
#     sig_err_list = np.array(sig_err_list)
#
#     # get all the gaussian functions that make sense
#     # then need to be central
#     mask_mu = np.abs(mu_list) < psf_std_oversampeled_pix * 3
#     # they must have a positive amplitude
#     mask_amp = (amp_list > 0)
#     mean_sigma_acs = np.nanmean(sig_list[mask_mu * mask_amp]) * oversampled_pixelscale_acs
#     mean_sigma_uvis = np.nanmean(sig_list[mask_mu * mask_amp]) * oversampled_pixelscale_uvis
#
#     measured_std_values_acs[conv_idx] = mean_sigma_acs
#     measured_std_values_uvis[conv_idx] = mean_sigma_uvis
#
#
# hst_custom_psf_apert_gauss_corr_dict = {
#     'acs':
#         {
#             'apert_radii_arcsec_list': apert_radii_acs_arcsec_list,
#             'measured_std_values': measured_std_values_acs,
#             'measured_flux_frac_in_apert': measured_flux_frac_in_apert_acs
#             },
#     'uvis':
#         {
#             'apert_radii_arcsec_list': apert_radii_uvis_arcsec_list,
#             'measured_std_values': measured_std_values_uvis,
#             'measured_flux_frac_in_apert': measured_flux_frac_in_apert_uvis
#             },
# }
#
# # plot_the output
# # interpolate function
# interp_func_acs = helper_func.InterpTools.interp2dgrid(
#     x_bins=apert_radii_acs_arcsec_list, y_bins=measured_std_values_acs, func_values=measured_flux_frac_in_apert_acs,
#     method='cubic')
#
# interp_func_uvis = helper_func.InterpTools.interp2dgrid(
#     x_bins=apert_radii_uvis_arcsec_list, y_bins=measured_std_values_uvis, func_values=measured_flux_frac_in_apert_uvis,
#     method='cubic')
#
# fine_interpolated_acs = helper_func.InterpTools.get2dinterp_fine_grid(
#     interp_func=interp_func_acs, x_min=np.min(apert_radii_acs_arcsec_list), x_max=np.max(apert_radii_acs_arcsec_list),
#     y_min=np.min(measured_std_values_acs), y_max=np.max(measured_std_values_acs), n_x_bins=100, n_y_bins=100)
#
# fine_interpolated_uvis = helper_func.InterpTools.get2dinterp_fine_grid(
#     interp_func=interp_func_uvis, x_min=np.min(apert_radii_uvis_arcsec_list), x_max=np.max(apert_radii_uvis_arcsec_list),
#     y_min=np.min(measured_std_values_uvis), y_max=np.max(measured_std_values_uvis), n_x_bins=100, n_y_bins=100)
#
# fig = plt.figure(figsize=(30, 12))
# fontsize = 30
# ax_corr_val = fig.add_axes((0.07, 0.07, 0.8, 0.92))
# ax_cbar = fig.add_axes((0.9, 0.3, 0.02, 0.4))
# norm = plotting_tools.ColorBarTools.compute_cbar_norm(vmin_vmax=(np.min(measured_flux_frac_in_apert_acs), np.max(measured_flux_frac_in_apert_acs)))
# ax_corr_val.imshow(fine_interpolated_acs,
#            extent=(np.min(apert_radii_acs_arcsec_list), np.max(apert_radii_acs_arcsec_list),
#                    np.min(measured_std_values_acs), np.max(measured_std_values_acs)),
#            origin='lower',
#            norm=norm, cmap='plasma')
# # plot colorbar
# plotting_tools.ColorBarTools.create_cbar(ax_cbar=ax_cbar, cmap='plasma', norm=norm, cbar_label='EE in Aperture', fontsize=fontsize, ticks=None, labelpad=2, tick_width=2,
#                 orientation='vertical', top_lable=True, label_color='k',
#                 extend='neither')
#
# # plot also the dots from where we got the interpolation
# x_mesh, y_mesh = np.meshgrid(apert_radii_acs_arcsec_list, measured_std_values_acs)
#
# ax_corr_val.scatter(x_mesh, y_mesh, c=measured_flux_frac_in_apert_acs.T, norm=norm, cmap='plasma', edgecolors='k')
# ax_corr_val.set_xlabel('Aperture [arcsec]', fontsize=fontsize,)
# ax_corr_val.set_ylabel(r'Measured Gauss. $\sigma$ [arcsec]', fontsize=fontsize,)
# ax_corr_val.tick_params(axis='both', which='both', width=3, length=15, right=False, top=True, direction='in',
#                labelsize=fontsize)
# if not os.path.isdir('plot_output/apert_profile_corr'): os.makedirs('plot_output/apert_profile_corr')
# fig.savefig('plot_output/apert_profile_corr/custom_psf_apert_profile_corr_acs.png')
# plt.clf()
# plt.cla()
#
#
#
# fig = plt.figure(figsize=(30, 12))
# fontsize = 30
# ax_corr_val = fig.add_axes((0.07, 0.07, 0.8, 0.92))
# ax_cbar = fig.add_axes((0.9, 0.3, 0.02, 0.4))
# norm = plotting_tools.ColorBarTools.compute_cbar_norm(vmin_vmax=(np.min(measured_flux_frac_in_apert_uvis), np.max(measured_flux_frac_in_apert_uvis)))
# ax_corr_val.imshow(fine_interpolated_uvis,
#            extent=(np.min(apert_radii_uvis_arcsec_list), np.max(apert_radii_uvis_arcsec_list),
#                    np.min(measured_std_values_uvis), np.max(measured_std_values_uvis)),
#            origin='lower',
#            norm=norm, cmap='plasma')
# # plot colorbar
# plotting_tools.ColorBarTools.create_cbar(ax_cbar=ax_cbar, cmap='plasma', norm=norm, cbar_label='EE in Aperture', fontsize=fontsize, ticks=None, labelpad=2, tick_width=2,
#                 orientation='vertical', top_lable=True, label_color='k',
#                 extend='neither')
#
# # plot also the dots from where we got the interpolation
# x_mesh, y_mesh = np.meshgrid(apert_radii_uvis_arcsec_list, measured_std_values_uvis)
#
# ax_corr_val.scatter(x_mesh, y_mesh, c=measured_flux_frac_in_apert_uvis.T, norm=norm, cmap='plasma', edgecolors='k')
# ax_corr_val.set_xlabel('Aperture [arcsec]', fontsize=fontsize,)
# ax_corr_val.set_ylabel(r'Measured Gauss. $\sigma$ [arcsec]', fontsize=fontsize,)
# ax_corr_val.tick_params(axis='both', which='both', width=3, length=15, right=False, top=True, direction='in',
#                labelsize=fontsize)
# if not os.path.isdir('plot_output/apert_profile_corr'): os.makedirs('plot_output/apert_profile_corr')
# fig.savefig('plot_output/apert_profile_corr/custom_psf_apert_profile_corr_uvis.png')
# plt.clf()
# plt.cla()
#
#
#
# # save dictionary
# if not os.path.isdir('data_output/apert_gauss_corr'):
#     os.makedirs('data_output/apert_gauss_corr')
#
# with open('data_output/apert_gauss_corr/hst_custom_psf_apert_gauss_corr_dict.pickle', 'wb') as file_name:
#     pickle.dump(hst_custom_psf_apert_gauss_corr_dict, file_name)


