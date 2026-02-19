import os
from pathlib import Path
import numpy as np
from astropy.io import fits
import astropy.units as u
from astropy.coordinates import SkyCoord
import matplotlib.pyplot as plt
from sternenfabrik import StarFactory
from obszugang import access_config, obs_info, ObsTools
from werkzeugkiste import helper_func
from malkasten import plotting_tools

################################
#### some useful parameters ####
################################
hst_overview_img_params = {
    'color_r': '#FF4433',
    'color_g': '#0FFF50',
    'color_b': '#1F51FF',
    'min_max_r': (0.3, 99.5),
    'min_max_g': (0.3, 99.5),
    'min_max_b': (0.3, 99.5),
    'gamma_r': 17.5,
    'gamma_g': 17.5,
    'gamma_b': 17.5,
    'gamma_corr_r': 17.5,
    'gamma_corr_g': 17.5,
    'gamma_corr_b': 17.5,
    'combined_gamma': 17.5}



# get data access
# define access parameters
target_name = 'ngc5194'
nircam_data_ver='v0p3p2'
miri_data_ver='v0p3p2'
nirspec_data_ver = 'karin_reduction_v1_oct2024'
miri_mrs_data_ver = 'karin_reduction'
# create an access instance
analysis_access = StarFactory(target_name=target_name,
                              x_target_name=target_name, radio_target_name=target_name,
                              nircam_data_ver=nircam_data_ver, miri_data_ver=miri_data_ver,
                              nirspec_data_ver=nirspec_data_ver,
                              miri_mrs_data_ver=miri_mrs_data_ver)
# load all bands that will beused
print('All available HST bands ', ObsTools.get_hst_obs_band_list(target=target_name))
print('All available NIRcam bands ', ObsTools.get_nircam_obs_band_list(target=target_name, version=nircam_data_ver))
print('All available MIRI bands ', ObsTools.get_miri_obs_band_list(target=target_name, version=miri_data_ver))
band_list_hst_overview = ['F814W', 'F555W', 'F435W']
band_list_zoom_in_hst = ['F275W', 'F435W', 'F658N']
band_list_zoom_in_nircam_1 = ['F115W', 'F164N', 'F200W']
band_list_zoom_in_nircam_2 = ['F115W', 'F187N', 'F200W']
band_list_zoom_in_nircam_3 = ['F300M', 'F335M', 'F360M']
band_list_zoom_in_miri = ['F770W', 'F1000W', 'F1130W']
band_list = list(np.unique(band_list_hst_overview + band_list_zoom_in_hst + band_list_zoom_in_nircam_1 +
                           band_list_zoom_in_nircam_2 + band_list_zoom_in_nircam_3 + band_list_zoom_in_miri))
print('all used bands ', band_list)
# load band
analysis_access.load_obs_bands(band_list=band_list)

# get the obs hulls of the JWST spec observations
nirspec_hull_dict = analysis_access.get_nirspec_obs_coverage_hull_dict()
miri_mrs_hull_dict = analysis_access.get_miri_mrs_obs_coverage_hull_dict()

# get center of hulls
# NE
mean_ra_ne_g140m_f100lp = np.mean(nirspec_hull_dict['NE']['G140M/F100LP'][0]['ra'])
mean_dec_ne_g140m_f100lp = np.mean(nirspec_hull_dict['NE']['G140M/F100LP'][0]['dec'])
# N
mean_ra_n_g140m_f100lp = np.mean(nirspec_hull_dict['N']['G140M/F100LP'][0]['ra'])
mean_dec_n_g140m_f100lp = np.mean(nirspec_hull_dict['N']['G140M/F100LP'][0]['dec'])
# SE
mean_ra_se_g140m_f100lp = np.mean(nirspec_hull_dict['SE']['G140M/F100LP'][0]['ra'])
mean_dec_se_g140m_f100lp = np.mean(nirspec_hull_dict['SE']['G140M/F100LP'][0]['dec'])


# get an HST overview image
img_overview, wcs_overview = analysis_access.get_target_overview_rgb_img(
            red_band='F814W', green_band='F555W', blue_band='F435W',
            overview_img_size=(2000, 2000), **hst_overview_img_params)

# get the zoom-in cutouts
zoom_in_cutout_size = (17, 20)
# get hst zoom_in
img_zoom_in_hst_ne, wcs_zoom_in_hst_ne = analysis_access.get_rgb_zoom_in(
    ra=mean_ra_ne_g140m_f100lp, dec=mean_dec_ne_g140m_f100lp, cutout_size=zoom_in_cutout_size,
    band_red=band_list_zoom_in_hst[2], band_green=band_list_zoom_in_hst[1], band_blue=band_list_zoom_in_hst[0])
img_zoom_in_hst_n, wcs_zoom_in_hst_n = analysis_access.get_rgb_zoom_in(
    ra=mean_ra_n_g140m_f100lp, dec=mean_dec_n_g140m_f100lp, cutout_size=zoom_in_cutout_size,
    band_red=band_list_zoom_in_hst[2], band_green=band_list_zoom_in_hst[1], band_blue=band_list_zoom_in_hst[0])
img_zoom_in_hst_se, wcs_zoom_in_hst_se = analysis_access.get_rgb_zoom_in(
    ra=mean_ra_se_g140m_f100lp, dec=mean_dec_se_g140m_f100lp, cutout_size=zoom_in_cutout_size,
    band_red=band_list_zoom_in_hst[2], band_green=band_list_zoom_in_hst[1], band_blue=band_list_zoom_in_hst[0])
# nircam 1
img_zoom_in_nircam_1_ne, wcs_zoom_in_nircam_1_ne = analysis_access.get_rgb_zoom_in(
    ra=mean_ra_ne_g140m_f100lp, dec=mean_dec_ne_g140m_f100lp, cutout_size=zoom_in_cutout_size,
    band_red=band_list_zoom_in_nircam_1[2], band_green=band_list_zoom_in_nircam_1[1], band_blue=band_list_zoom_in_nircam_1[0])
img_zoom_in_nircam_1_n, wcs_zoom_in_nircam_1_n = analysis_access.get_rgb_zoom_in(
    ra=mean_ra_n_g140m_f100lp, dec=mean_dec_n_g140m_f100lp, cutout_size=zoom_in_cutout_size,
    band_red=band_list_zoom_in_nircam_1[2], band_green=band_list_zoom_in_nircam_1[1], band_blue=band_list_zoom_in_nircam_1[0])
img_zoom_in_nircam_1_se, wcs_zoom_in_nircam_1_se = analysis_access.get_rgb_zoom_in(
    ra=mean_ra_se_g140m_f100lp, dec=mean_dec_se_g140m_f100lp, cutout_size=zoom_in_cutout_size,
    band_red=band_list_zoom_in_nircam_1[2], band_green=band_list_zoom_in_nircam_1[1], band_blue=band_list_zoom_in_nircam_1[0])
# nircam 2
img_zoom_in_nircam_2_ne, wcs_zoom_in_nircam_2_ne = analysis_access.get_rgb_zoom_in(
    ra=mean_ra_ne_g140m_f100lp, dec=mean_dec_ne_g140m_f100lp, cutout_size=zoom_in_cutout_size,
    band_red=band_list_zoom_in_nircam_2[2], band_green=band_list_zoom_in_nircam_2[1], band_blue=band_list_zoom_in_nircam_2[0])
img_zoom_in_nircam_2_n, wcs_zoom_in_nircam_2_n = analysis_access.get_rgb_zoom_in(
    ra=mean_ra_n_g140m_f100lp, dec=mean_dec_n_g140m_f100lp, cutout_size=zoom_in_cutout_size,
    band_red=band_list_zoom_in_nircam_2[2], band_green=band_list_zoom_in_nircam_2[1], band_blue=band_list_zoom_in_nircam_2[0])
img_zoom_in_nircam_2_se, wcs_zoom_in_nircam_2_se = analysis_access.get_rgb_zoom_in(
    ra=mean_ra_se_g140m_f100lp, dec=mean_dec_se_g140m_f100lp, cutout_size=zoom_in_cutout_size,
    band_red=band_list_zoom_in_nircam_2[2], band_green=band_list_zoom_in_nircam_2[1], band_blue=band_list_zoom_in_nircam_2[0])
# nircam 3
img_zoom_in_nircam_3_ne, wcs_zoom_in_nircam_3_ne = analysis_access.get_rgb_zoom_in(
    ra=mean_ra_ne_g140m_f100lp, dec=mean_dec_ne_g140m_f100lp, cutout_size=zoom_in_cutout_size,
    band_red=band_list_zoom_in_nircam_3[2], band_green=band_list_zoom_in_nircam_3[1], band_blue=band_list_zoom_in_nircam_3[0])
img_zoom_in_nircam_3_n, wcs_zoom_in_nircam_3_n = analysis_access.get_rgb_zoom_in(
    ra=mean_ra_n_g140m_f100lp, dec=mean_dec_n_g140m_f100lp, cutout_size=zoom_in_cutout_size,
    band_red=band_list_zoom_in_nircam_3[2], band_green=band_list_zoom_in_nircam_3[1], band_blue=band_list_zoom_in_nircam_3[0])
img_zoom_in_nircam_3_se, wcs_zoom_in_nircam_3_se = analysis_access.get_rgb_zoom_in(
    ra=mean_ra_se_g140m_f100lp, dec=mean_dec_se_g140m_f100lp, cutout_size=zoom_in_cutout_size,
    band_red=band_list_zoom_in_nircam_3[2], band_green=band_list_zoom_in_nircam_3[1], band_blue=band_list_zoom_in_nircam_3[0])
# miri
img_zoom_in_miri_ne, wcs_zoom_in_miri_ne = analysis_access.get_rgb_zoom_in(
    ra=mean_ra_ne_g140m_f100lp, dec=mean_dec_ne_g140m_f100lp, cutout_size=zoom_in_cutout_size,
    band_red=band_list_zoom_in_miri[2], band_green=band_list_zoom_in_miri[1], band_blue=band_list_zoom_in_miri[0])
img_zoom_in_miri_n, wcs_zoom_in_miri_n = analysis_access.get_rgb_zoom_in(
    ra=mean_ra_n_g140m_f100lp, dec=mean_dec_n_g140m_f100lp, cutout_size=zoom_in_cutout_size,
    band_red=band_list_zoom_in_miri[2], band_green=band_list_zoom_in_miri[1], band_blue=band_list_zoom_in_miri[0])
img_zoom_in_miri_se, wcs_zoom_in_miri_se = analysis_access.get_rgb_zoom_in(
    ra=mean_ra_se_g140m_f100lp, dec=mean_dec_se_g140m_f100lp, cutout_size=zoom_in_cutout_size,
    band_red=band_list_zoom_in_miri[2], band_green=band_list_zoom_in_miri[1], band_blue=band_list_zoom_in_miri[0])


# get radio observations
analysis_access.load_radio_data(band='L', map='tt0')
cutout_tt0_ne = helper_func.CoordTools.get_img_cutout(
    img=analysis_access.radio_data['L_tt0_data_img'], wcs=analysis_access.radio_data['L_tt0_wcs_img'],
    coord=SkyCoord(ra=mean_ra_ne_g140m_f100lp*u.deg, dec=mean_dec_ne_g140m_f100lp*u.deg),
    cutout_size=zoom_in_cutout_size)
cutout_tt0_n = helper_func.CoordTools.get_img_cutout(
    img=analysis_access.radio_data['L_tt0_data_img'], wcs=analysis_access.radio_data['L_tt0_wcs_img'],
    coord=SkyCoord(ra=mean_ra_n_g140m_f100lp*u.deg, dec=mean_dec_n_g140m_f100lp*u.deg),
    cutout_size=zoom_in_cutout_size)
cutout_tt0_se = helper_func.CoordTools.get_img_cutout(
    img=analysis_access.radio_data['L_tt0_data_img'], wcs=analysis_access.radio_data['L_tt0_wcs_img'],
    coord=SkyCoord(ra=mean_ra_se_g140m_f100lp*u.deg, dec=mean_dec_se_g140m_f100lp*u.deg),
    cutout_size=zoom_in_cutout_size)

# get X-ray data
# load X-ray data
analysis_access.load_chandra_data(energy='0p5to2')
analysis_access.load_chandra_data(energy='2to7')

# get cutout
cutout_0p5to2_ne = helper_func.CoordTools.get_img_cutout(
    img=analysis_access.x_data['0p5to2_data_img'], wcs=analysis_access.x_data['0p5to2_wcs_img'],
    coord=SkyCoord(ra=mean_ra_ne_g140m_f100lp*u.deg, dec=mean_dec_ne_g140m_f100lp*u.deg), cutout_size=zoom_in_cutout_size)
cutout_2to7_ne = helper_func.CoordTools.get_img_cutout(
    img=analysis_access.x_data['2to7_data_img'], wcs=analysis_access.x_data['2to7_wcs_img'],
    coord=SkyCoord(ra=mean_ra_ne_g140m_f100lp*u.deg, dec=mean_dec_ne_g140m_f100lp*u.deg), cutout_size=zoom_in_cutout_size)
# get x-ray_rgb image
x_ray_rgb_ne = plotting_tools.ImgTools.get_2color_img(data1=cutout_0p5to2_ne.data, data2=cutout_2to7_ne.data)

cutout_0p5to2_n = helper_func.CoordTools.get_img_cutout(
    img=analysis_access.x_data['0p5to2_data_img'], wcs=analysis_access.x_data['0p5to2_wcs_img'],
    coord=SkyCoord(ra=mean_ra_n_g140m_f100lp*u.deg, dec=mean_dec_n_g140m_f100lp*u.deg), cutout_size=zoom_in_cutout_size)
cutout_2to7_n = helper_func.CoordTools.get_img_cutout(
    img=analysis_access.x_data['2to7_data_img'], wcs=analysis_access.x_data['2to7_wcs_img'],
    coord=SkyCoord(ra=mean_ra_n_g140m_f100lp*u.deg, dec=mean_dec_n_g140m_f100lp*u.deg), cutout_size=zoom_in_cutout_size)
# get x-ray_rgb image
x_ray_rgb_n = plotting_tools.ImgTools.get_2color_img(data1=cutout_0p5to2_n.data, data2=cutout_2to7_n.data)

cutout_0p5to2_se = helper_func.CoordTools.get_img_cutout(
    img=analysis_access.x_data['0p5to2_data_img'], wcs=analysis_access.x_data['0p5to2_wcs_img'],
    coord=SkyCoord(ra=mean_ra_se_g140m_f100lp*u.deg, dec=mean_dec_se_g140m_f100lp*u.deg), cutout_size=zoom_in_cutout_size)
cutout_2to7_se = helper_func.CoordTools.get_img_cutout(
    img=analysis_access.x_data['2to7_data_img'], wcs=analysis_access.x_data['2to7_wcs_img'],
    coord=SkyCoord(ra=mean_ra_se_g140m_f100lp*u.deg, dec=mean_dec_se_g140m_f100lp*u.deg), cutout_size=zoom_in_cutout_size)
# get x-ray_rgb image
x_ray_rgb_se = plotting_tools.ImgTools.get_2color_img(data1=cutout_0p5to2_se.data, data2=cutout_2to7_se.data)



fig = plt.figure(figsize=(30, 80))
fontsize_medium = 45
fontsize_large = 60
# add axis
ax_overview = fig.add_axes((0.1, 0.445, 0.8, 0.8), projection=wcs_overview)

ax_zoom_in_hst_ne = fig.add_axes((0.05, 0.545 - 0.02, 0.28, 0.2), projection=wcs_zoom_in_hst_ne)
ax_zoom_in_hst_n = fig.add_axes((0.37, 0.545 - 0.02, 0.28, 0.2), projection=wcs_zoom_in_hst_n)
ax_zoom_in_hst_se = fig.add_axes((0.69, 0.545 - 0.02, 0.28, 0.2), projection=wcs_zoom_in_hst_se)

ax_zoom_in_nircam_1_ne = fig.add_axes((0.05, 0.45 - 0.02, 0.28, 0.2), projection=wcs_zoom_in_nircam_1_ne)
ax_zoom_in_nircam_1_n = fig.add_axes((0.37, 0.45 - 0.02, 0.28, 0.2), projection=wcs_zoom_in_nircam_1_n)
ax_zoom_in_nircam_1_se = fig.add_axes((0.69, 0.45 - 0.02, 0.28, 0.2), projection=wcs_zoom_in_nircam_1_se)

ax_zoom_in_nircam_2_ne = fig.add_axes((0.05, 0.355 - 0.02, 0.28, 0.2), projection=wcs_zoom_in_nircam_2_ne)
ax_zoom_in_nircam_2_n = fig.add_axes((0.37, 0.355 - 0.02, 0.28, 0.2), projection=wcs_zoom_in_nircam_2_n)
ax_zoom_in_nircam_2_se = fig.add_axes((0.69, 0.355 - 0.02, 0.28, 0.2), projection=wcs_zoom_in_nircam_2_se)

ax_zoom_in_nircam_3_ne = fig.add_axes((0.05, 0.26 - 0.02, 0.28, 0.2), projection=wcs_zoom_in_nircam_3_ne)
ax_zoom_in_nircam_3_n = fig.add_axes((0.37, 0.26 - 0.02, 0.28, 0.2), projection=wcs_zoom_in_nircam_3_n)
ax_zoom_in_nircam_3_se = fig.add_axes((0.69, 0.26 - 0.02, 0.28, 0.2), projection=wcs_zoom_in_nircam_3_se)

ax_zoom_in_miri_ne = fig.add_axes((0.05, 0.165 - 0.02, 0.28, 0.2), projection=wcs_zoom_in_miri_ne)
ax_zoom_in_miri_n = fig.add_axes((0.37, 0.165 - 0.02, 0.28, 0.2), projection=wcs_zoom_in_miri_n)
ax_zoom_in_miri_se = fig.add_axes((0.69, 0.165 - 0.02, 0.28, 0.2), projection=wcs_zoom_in_miri_se)

ax_zoom_in_tt0_ne = fig.add_axes((0.05, 0.07 - 0.02, 0.28, 0.2), projection=cutout_tt0_ne.wcs)
ax_zoom_in_tt0_n = fig.add_axes((0.37, 0.07 - 0.02, 0.28, 0.2), projection=cutout_tt0_n.wcs)
ax_zoom_in_tt0_se = fig.add_axes((0.69, 0.07 - 0.02, 0.28, 0.2), projection=cutout_tt0_se.wcs)

ax_zoom_in_x_ray_ne = fig.add_axes((0.05, -0.025 - 0.02, 0.28, 0.2), projection=cutout_0p5to2_ne.wcs)
ax_zoom_in_x_ray_n = fig.add_axes((0.37, -0.025 - 0.02, 0.28, 0.2), projection=cutout_0p5to2_n.wcs)
ax_zoom_in_x_ray_se = fig.add_axes((0.69, -0.025 - 0.02, 0.28, 0.2), projection=cutout_0p5to2_se.wcs)



ax_overview.imshow(img_overview)
# erase axis
plotting_tools.WCSPlottingTools.arr_axis_params(
    ax=ax_overview, ra_tick_label=False, dec_tick_label=False, ra_axis_label=' ', dec_axis_label=' ',)
# add hulls
# nispec
plotting_tools.WCSPlottingTools.plot_obs_hull(
    ax=ax_overview, wcs=wcs_overview, hull_dict=nirspec_hull_dict['NE']['G140M/F100LP'],
    line_width=4, line_style='-', color='tab:red')
plotting_tools.WCSPlottingTools.plot_obs_hull(
    ax=ax_overview, wcs=wcs_overview, hull_dict=nirspec_hull_dict['N']['G140M/F100LP'],
    line_width=4, line_style='-', color='tab:green')
plotting_tools.WCSPlottingTools.plot_obs_hull(
    ax=ax_overview, wcs=wcs_overview, hull_dict=nirspec_hull_dict['SE']['G140M/F100LP'],
    line_width=4, line_style='-', color='tab:orange')
# miri mrs
plotting_tools.WCSPlottingTools.plot_obs_hull(
    ax=ax_overview, wcs=wcs_overview, hull_dict=miri_mrs_hull_dict['NE']['CH3'],
    line_width=4, line_style='-', color='tab:red')
plotting_tools.WCSPlottingTools.plot_obs_hull(
    ax=ax_overview, wcs=wcs_overview, hull_dict=miri_mrs_hull_dict['N']['CH3'],
    line_width=4, line_style='-', color='tab:green')
plotting_tools.WCSPlottingTools.plot_obs_hull(
    ax=ax_overview, wcs=wcs_overview, hull_dict=miri_mrs_hull_dict['SE']['CH3'],
    line_width=4, line_style='-', color='tab:orange')
# add a scale bar for reference
plotting_tools.WCSPlottingTools.plot_img_scale_bar(ax=ax_overview, img_shape=img_overview[:, :, 0].shape,
                                                   wcs=wcs_overview, bar_length=5, length_unit='kpc', target_dist_mpc=None,
                           phangs_target=target_name, bar_color='white', text_color='white', line_width=4, fontsize=fontsize_large,
                           va='bottom', ha='left', x_offset=0.05, y_offset=0.05, text_y_offset_diff=0.01,
                           path_eff=True, path_eff_color='white')
plotting_tools.WCSPlottingTools.plot_img_scale_bar(ax=ax_overview, img_shape=img_overview[:, :, 0].shape,
                                                   wcs=wcs_overview, bar_length=1, length_unit='arcmin', target_dist_mpc=None,
                           phangs_target=target_name, bar_color='white', text_color='white', line_width=4, fontsize=fontsize_large,
                           va='bottom', ha='left', x_offset=0.05, y_offset=0.1, text_y_offset_diff=0.01,
                           path_eff=True, path_eff_color='white')
# add text
plotting_tools.StrTools.display_text_in_corner(ax=ax_overview, text='HST',
                                               fontsize=fontsize_large,
                                               text_color='white', x_frac=0.02, y_frac=0.98,
                                               horizontal_alignment='left', vertical_alignment='top',
                                               path_eff=True, path_err_linewidth=3, path_eff_color='white')
plotting_tools.StrTools.display_text_in_corner(ax=ax_overview, text='I',
                                               fontsize=fontsize_large,
                                               text_color='red', x_frac=0.02, y_frac=0.93,
                                               horizontal_alignment='left', vertical_alignment='top',
                                               path_eff=True, path_err_linewidth=3, path_eff_color='white')
plotting_tools.StrTools.display_text_in_corner(ax=ax_overview, text='V',
                                               fontsize=fontsize_large,
                                               text_color='green', x_frac=0.02, y_frac=0.88,
                                               horizontal_alignment='left', vertical_alignment='top',
                                               path_eff=True, path_err_linewidth=3, path_eff_color='white')
plotting_tools.StrTools.display_text_in_corner(ax=ax_overview, text='B',
                                               fontsize=fontsize_large,
                                               text_color='blue', x_frac=0.02, y_frac=0.83,
                                               horizontal_alignment='left', vertical_alignment='top',
                                               path_eff=True, path_err_linewidth=3, path_eff_color='white')



ax_zoom_in_hst_ne.imshow(img_zoom_in_hst_ne)
ax_zoom_in_hst_ne.set_title('NE', fontsize=fontsize_large)
# erase axis
plotting_tools.WCSPlottingTools.arr_axis_params(
    ax=ax_zoom_in_hst_ne, ra_tick_label=False, dec_tick_label=False, ra_axis_label=' ', dec_axis_label=' ',)
# display hulls
plotting_tools.WCSPlottingTools.plot_obs_hull(
    ax=ax_zoom_in_hst_ne, wcs=wcs_zoom_in_hst_ne, hull_dict=nirspec_hull_dict['NE']['G140M/F100LP'],
    line_width=4, line_style='-', color='tab:red')
plotting_tools.WCSPlottingTools.plot_obs_hull(
    ax=ax_zoom_in_hst_ne, wcs=wcs_zoom_in_hst_ne, hull_dict=miri_mrs_hull_dict['NE']['CH3'],
    line_width=4, line_style='-', color='tab:red')
# add a scale bar for reference
plotting_tools.WCSPlottingTools.plot_img_scale_bar(
    ax=ax_zoom_in_hst_ne, img_shape=img_zoom_in_hst_ne[:, :, 0].shape, wcs=wcs_zoom_in_hst_ne, bar_length=50,
    length_unit='pc', phangs_target=target_name,  x_offset=0.1, y_offset=0.05, fontsize=fontsize_medium)
plotting_tools.WCSPlottingTools.plot_img_scale_bar(
    ax=ax_zoom_in_hst_ne, img_shape=img_zoom_in_hst_ne[:, :, 0].shape, wcs=wcs_zoom_in_hst_ne, bar_length=5,
    length_unit='arcsec', phangs_target=target_name,  x_offset=0.1, y_offset=0.2, fontsize=fontsize_medium,)
# add text
plotting_tools.StrTools.display_text_in_corner(ax=ax_zoom_in_hst_ne, text=r'H$\alpha$',
                                               fontsize=fontsize_medium,
                                               text_color='red', x_frac=0.98, y_frac=0.96,
                                               horizontal_alignment='right', vertical_alignment='top',
                                               path_eff=True, path_err_linewidth=3, path_eff_color='white')
plotting_tools.StrTools.display_text_in_corner(ax=ax_zoom_in_hst_ne, text='B',
                                               fontsize=fontsize_medium,
                                               text_color='green', x_frac=0.98, y_frac=0.88,
                                               horizontal_alignment='right', vertical_alignment='top',
                                               path_eff=True, path_err_linewidth=3, path_eff_color='white')
plotting_tools.StrTools.display_text_in_corner(ax=ax_zoom_in_hst_ne, text='NUV',
                                               fontsize=fontsize_medium,
                                               text_color='blue', x_frac=0.98, y_frac=0.80,
                                               horizontal_alignment='right', vertical_alignment='top',
                                               path_eff=True, path_err_linewidth=3, path_eff_color='white')
# N
ax_zoom_in_hst_n.imshow(img_zoom_in_hst_n)
ax_zoom_in_hst_n.set_title('N', fontsize=fontsize_large)
# erase axis
plotting_tools.WCSPlottingTools.arr_axis_params(
    ax=ax_zoom_in_hst_n, ra_tick_label=False, dec_tick_label=False, ra_axis_label=' ', dec_axis_label=' ',)
# display hulls
plotting_tools.WCSPlottingTools.plot_obs_hull(
    ax=ax_zoom_in_hst_n, wcs=wcs_zoom_in_hst_n, hull_dict=nirspec_hull_dict['N']['G140M/F100LP'],
    line_width=4, line_style='-', color='tab:green')
plotting_tools.WCSPlottingTools.plot_obs_hull(
    ax=ax_zoom_in_hst_n, wcs=wcs_zoom_in_hst_n, hull_dict=miri_mrs_hull_dict['N']['CH3'],
    line_width=4, line_style='-', color='tab:green')
# add a scale bar for reference
plotting_tools.WCSPlottingTools.plot_img_scale_bar(
    ax=ax_zoom_in_hst_n, img_shape=img_zoom_in_hst_n[:, :, 0].shape, wcs=wcs_zoom_in_hst_n, bar_length=50,
    length_unit='pc', phangs_target=target_name,  x_offset=0.1, y_offset=0.05, fontsize=fontsize_medium)
plotting_tools.WCSPlottingTools.plot_img_scale_bar(
    ax=ax_zoom_in_hst_n, img_shape=img_zoom_in_hst_n[:, :, 0].shape, wcs=wcs_zoom_in_hst_n, bar_length=5,
    length_unit='arcsec', phangs_target=target_name,  x_offset=0.1, y_offset=0.2, fontsize=fontsize_medium,)
# SE
ax_zoom_in_hst_se.imshow(img_zoom_in_hst_se)
ax_zoom_in_hst_se.set_title('SE', fontsize=fontsize_large)
# erase axis
plotting_tools.WCSPlottingTools.arr_axis_params(
    ax=ax_zoom_in_hst_se, ra_tick_label=False, dec_tick_label=False, ra_axis_label=' ', dec_axis_label=' ',)
# display hulls
plotting_tools.WCSPlottingTools.plot_obs_hull(
    ax=ax_zoom_in_hst_se, wcs=wcs_zoom_in_hst_se, hull_dict=nirspec_hull_dict['SE']['G140M/F100LP'],
    line_width=4, line_style='-', color='tab:orange')
plotting_tools.WCSPlottingTools.plot_obs_hull(
    ax=ax_zoom_in_hst_se, wcs=wcs_zoom_in_hst_se, hull_dict=miri_mrs_hull_dict['SE']['CH3'],
    line_width=4, line_style='-', color='tab:orange')
# add a scale bar for reference
plotting_tools.WCSPlottingTools.plot_img_scale_bar(
    ax=ax_zoom_in_hst_se, img_shape=img_zoom_in_hst_se[:, :, 0].shape, wcs=wcs_zoom_in_hst_se, bar_length=50,
    length_unit='pc', phangs_target=target_name,  x_offset=0.1, y_offset=0.05, fontsize=fontsize_medium)
plotting_tools.WCSPlottingTools.plot_img_scale_bar(
    ax=ax_zoom_in_hst_se, img_shape=img_zoom_in_hst_se[:, :, 0].shape, wcs=wcs_zoom_in_hst_se, bar_length=5,
    length_unit='arcsec', phangs_target=target_name,  x_offset=0.1, y_offset=0.2, fontsize=fontsize_medium,)




ax_zoom_in_nircam_1_ne.imshow(img_zoom_in_nircam_1_ne)
# erase axis
plotting_tools.WCSPlottingTools.arr_axis_params(
    ax=ax_zoom_in_nircam_1_ne, ra_tick_label=False, dec_tick_label=False, ra_axis_label=' ', dec_axis_label=' ',)
# display hulls
plotting_tools.WCSPlottingTools.plot_obs_hull(
    ax=ax_zoom_in_nircam_1_ne, wcs=wcs_zoom_in_nircam_1_ne, hull_dict=nirspec_hull_dict['NE']['G140M/F100LP'],
    line_width=4, line_style='-', color='tab:red')
plotting_tools.WCSPlottingTools.plot_obs_hull(
    ax=ax_zoom_in_nircam_1_ne, wcs=wcs_zoom_in_nircam_1_ne, hull_dict=miri_mrs_hull_dict['NE']['CH3'],
    line_width=4, line_style='-', color='tab:red')
# add text
plotting_tools.StrTools.display_text_in_corner(
    ax=ax_zoom_in_nircam_1_ne, text='F200W', fontsize=fontsize_medium, text_color='red', x_frac=0.98, y_frac=0.96,
    horizontal_alignment='right')
plotting_tools.StrTools.display_text_in_corner(
    ax=ax_zoom_in_nircam_1_ne, text='F164N', fontsize=fontsize_medium, text_color='green', x_frac=0.98, y_frac=0.88,
    horizontal_alignment='right')
plotting_tools.StrTools.display_text_in_corner(
    ax=ax_zoom_in_nircam_1_ne, text='F115W', fontsize=fontsize_medium, text_color='blue', x_frac=0.98, y_frac=0.80,
    horizontal_alignment='right')
ax_zoom_in_nircam_1_n.imshow(img_zoom_in_nircam_1_n)
# erase axis
plotting_tools.WCSPlottingTools.arr_axis_params(
    ax=ax_zoom_in_nircam_1_n, ra_tick_label=False, dec_tick_label=False, ra_axis_label=' ', dec_axis_label=' ',)
# display hulls
plotting_tools.WCSPlottingTools.plot_obs_hull(
    ax=ax_zoom_in_nircam_1_n, wcs=wcs_zoom_in_nircam_1_n, hull_dict=nirspec_hull_dict['N']['G140M/F100LP'],
    line_width=4, line_style='-', color='tab:green')
plotting_tools.WCSPlottingTools.plot_obs_hull(
    ax=ax_zoom_in_nircam_1_n, wcs=wcs_zoom_in_nircam_1_n, hull_dict=miri_mrs_hull_dict['N']['CH3'],
    line_width=4, line_style='-', color='tab:green')
ax_zoom_in_nircam_1_se.imshow(img_zoom_in_nircam_1_se)
# erase axis
plotting_tools.WCSPlottingTools.arr_axis_params(
    ax=ax_zoom_in_nircam_1_se, ra_tick_label=False, dec_tick_label=False, ra_axis_label=' ', dec_axis_label=' ',)
# display hulls
plotting_tools.WCSPlottingTools.plot_obs_hull(
    ax=ax_zoom_in_nircam_1_se, wcs=wcs_zoom_in_nircam_1_se, hull_dict=nirspec_hull_dict['SE']['G140M/F100LP'],
    line_width=4, line_style='-', color='tab:orange')
plotting_tools.WCSPlottingTools.plot_obs_hull(
    ax=ax_zoom_in_nircam_1_se, wcs=wcs_zoom_in_nircam_1_se, hull_dict=miri_mrs_hull_dict['SE']['CH3'],
    line_width=4, line_style='-', color='tab:orange')



ax_zoom_in_nircam_2_ne.imshow(img_zoom_in_nircam_2_ne)
# erase axis
plotting_tools.WCSPlottingTools.arr_axis_params(
    ax=ax_zoom_in_nircam_2_ne, ra_tick_label=False, dec_tick_label=False, ra_axis_label=' ', dec_axis_label=' ',)
# display hulls
plotting_tools.WCSPlottingTools.plot_obs_hull(
    ax=ax_zoom_in_nircam_2_ne, wcs=wcs_zoom_in_nircam_2_ne, hull_dict=nirspec_hull_dict['NE']['G140M/F100LP'],
    line_width=4, line_style='-', color='tab:red')
plotting_tools.WCSPlottingTools.plot_obs_hull(
    ax=ax_zoom_in_nircam_2_ne, wcs=wcs_zoom_in_nircam_2_ne, hull_dict=miri_mrs_hull_dict['NE']['CH3'],
    line_width=4, line_style='-', color='tab:red')
plotting_tools.StrTools.display_text_in_corner(
    ax=ax_zoom_in_nircam_2_ne, text='F200W', fontsize=fontsize_medium, text_color='red', x_frac=0.98, y_frac=0.96,
    horizontal_alignment='right')
plotting_tools.StrTools.display_text_in_corner(
    ax=ax_zoom_in_nircam_2_ne, text='F187N', fontsize=fontsize_medium, text_color='green', x_frac=0.98, y_frac=0.88,
    horizontal_alignment='right')
plotting_tools.StrTools.display_text_in_corner(
    ax=ax_zoom_in_nircam_2_ne, text='F115W', fontsize=fontsize_medium, text_color='blue', x_frac=0.98, y_frac=0.80,
    horizontal_alignment='right')
ax_zoom_in_nircam_2_n.imshow(img_zoom_in_nircam_2_n)
# erase axis
plotting_tools.WCSPlottingTools.arr_axis_params(
    ax=ax_zoom_in_nircam_2_n, ra_tick_label=False, dec_tick_label=False, ra_axis_label=' ', dec_axis_label=' ',)
# display hulls
plotting_tools.WCSPlottingTools.plot_obs_hull(
    ax=ax_zoom_in_nircam_2_n, wcs=wcs_zoom_in_nircam_2_n, hull_dict=nirspec_hull_dict['N']['G140M/F100LP'],
    line_width=4, line_style='-', color='tab:green')
plotting_tools.WCSPlottingTools.plot_obs_hull(
    ax=ax_zoom_in_nircam_2_n, wcs=wcs_zoom_in_nircam_2_n, hull_dict=miri_mrs_hull_dict['N']['CH3'],
    line_width=4, line_style='-', color='tab:green')
ax_zoom_in_nircam_2_se.imshow(img_zoom_in_nircam_2_se)
# erase axis
plotting_tools.WCSPlottingTools.arr_axis_params(
    ax=ax_zoom_in_nircam_2_se, ra_tick_label=False, dec_tick_label=False, ra_axis_label=' ', dec_axis_label=' ',)
# display hulls
plotting_tools.WCSPlottingTools.plot_obs_hull(
    ax=ax_zoom_in_nircam_2_se, wcs=wcs_zoom_in_nircam_2_se, hull_dict=nirspec_hull_dict['SE']['G140M/F100LP'],
    line_width=4, line_style='-', color='tab:orange')
plotting_tools.WCSPlottingTools.plot_obs_hull(
    ax=ax_zoom_in_nircam_2_se, wcs=wcs_zoom_in_nircam_2_se, hull_dict=miri_mrs_hull_dict['SE']['CH3'],
    line_width=4, line_style='-', color='tab:orange')




ax_zoom_in_nircam_3_ne.imshow(img_zoom_in_nircam_3_ne)
# erase axis
plotting_tools.WCSPlottingTools.arr_axis_params(
    ax=ax_zoom_in_nircam_3_ne, ra_tick_label=False, dec_tick_label=False, ra_axis_label=' ', dec_axis_label=' ',)
# display hulls
plotting_tools.WCSPlottingTools.plot_obs_hull(
    ax=ax_zoom_in_nircam_3_ne, wcs=wcs_zoom_in_nircam_3_ne, hull_dict=nirspec_hull_dict['NE']['G140M/F100LP'],
    line_width=4, line_style='-', color='tab:red')
plotting_tools.WCSPlottingTools.plot_obs_hull(
    ax=ax_zoom_in_nircam_3_ne, wcs=wcs_zoom_in_nircam_3_ne, hull_dict=miri_mrs_hull_dict['NE']['CH3'],
    line_width=4, line_style='-', color='tab:red')
plotting_tools.StrTools.display_text_in_corner(
    ax=ax_zoom_in_nircam_3_ne, text='F360M', fontsize=fontsize_medium, text_color='red', x_frac=0.98, y_frac=0.96,
    horizontal_alignment='right')
plotting_tools.StrTools.display_text_in_corner(
    ax=ax_zoom_in_nircam_3_ne, text='F335M', fontsize=fontsize_medium, text_color='green', x_frac=0.98, y_frac=0.88,
    horizontal_alignment='right')
plotting_tools.StrTools.display_text_in_corner(
    ax=ax_zoom_in_nircam_3_ne, text='F300M', fontsize=fontsize_medium, text_color='blue', x_frac=0.98, y_frac=0.80,
    horizontal_alignment='right')
ax_zoom_in_nircam_3_n.imshow(img_zoom_in_nircam_3_n)
# erase axis
plotting_tools.WCSPlottingTools.arr_axis_params(
    ax=ax_zoom_in_nircam_3_n, ra_tick_label=False, dec_tick_label=False, ra_axis_label=' ', dec_axis_label=' ',)
# display hulls
plotting_tools.WCSPlottingTools.plot_obs_hull(
    ax=ax_zoom_in_nircam_3_n, wcs=wcs_zoom_in_nircam_3_n, hull_dict=nirspec_hull_dict['N']['G140M/F100LP'],
    line_width=4, line_style='-', color='tab:green')
plotting_tools.WCSPlottingTools.plot_obs_hull(
    ax=ax_zoom_in_nircam_3_n, wcs=wcs_zoom_in_nircam_3_n, hull_dict=miri_mrs_hull_dict['N']['CH3'],
    line_width=4, line_style='-', color='tab:green')
ax_zoom_in_nircam_3_se.imshow(img_zoom_in_nircam_3_se)
# erase axis
plotting_tools.WCSPlottingTools.arr_axis_params(
    ax=ax_zoom_in_nircam_3_se, ra_tick_label=False, dec_tick_label=False, ra_axis_label=' ', dec_axis_label=' ',)
# display hulls
plotting_tools.WCSPlottingTools.plot_obs_hull(
    ax=ax_zoom_in_nircam_3_se, wcs=wcs_zoom_in_nircam_3_se, hull_dict=nirspec_hull_dict['SE']['G140M/F100LP'],
    line_width=4, line_style='-', color='tab:orange')
plotting_tools.WCSPlottingTools.plot_obs_hull(
    ax=ax_zoom_in_nircam_3_se, wcs=wcs_zoom_in_nircam_3_se, hull_dict=miri_mrs_hull_dict['SE']['CH3'],
    line_width=4, line_style='-', color='tab:orange')




ax_zoom_in_miri_ne.imshow(img_zoom_in_miri_ne)
# erase axis
plotting_tools.WCSPlottingTools.arr_axis_params(
    ax=ax_zoom_in_miri_ne, ra_tick_label=False, dec_tick_label=False, ra_axis_label=' ', dec_axis_label=' ',)
# display hulls
plotting_tools.WCSPlottingTools.plot_obs_hull(
    ax=ax_zoom_in_miri_ne, wcs=wcs_zoom_in_miri_ne, hull_dict=nirspec_hull_dict['NE']['G140M/F100LP'],
    line_width=4, line_style='-', color='tab:red')
plotting_tools.WCSPlottingTools.plot_obs_hull(
    ax=ax_zoom_in_miri_ne, wcs=wcs_zoom_in_miri_ne, hull_dict=miri_mrs_hull_dict['NE']['CH3'],
    line_width=4, line_style='-', color='tab:red')
plotting_tools.StrTools.display_text_in_corner(
    ax=ax_zoom_in_miri_ne, text='F1130W', fontsize=fontsize_medium, text_color='red', x_frac=0.98, y_frac=0.96,
    horizontal_alignment='right')
plotting_tools.StrTools.display_text_in_corner(
    ax=ax_zoom_in_miri_ne, text='F1000W', fontsize=fontsize_medium, text_color='green', x_frac=0.98, y_frac=0.88,
    horizontal_alignment='right')
plotting_tools.StrTools.display_text_in_corner(
    ax=ax_zoom_in_miri_ne, text='F770W', fontsize=fontsize_medium, text_color='blue', x_frac=0.98, y_frac=0.80,
    horizontal_alignment='right')
ax_zoom_in_miri_n.imshow(img_zoom_in_miri_n)
# erase axis
plotting_tools.WCSPlottingTools.arr_axis_params(
    ax=ax_zoom_in_miri_n, ra_tick_label=False, dec_tick_label=False, ra_axis_label=' ', dec_axis_label=' ',)
# display hulls
plotting_tools.WCSPlottingTools.plot_obs_hull(
    ax=ax_zoom_in_miri_n, wcs=wcs_zoom_in_miri_n, hull_dict=nirspec_hull_dict['N']['G140M/F100LP'],
    line_width=4, line_style='-', color='tab:green')
plotting_tools.WCSPlottingTools.plot_obs_hull(
    ax=ax_zoom_in_miri_n, wcs=wcs_zoom_in_miri_n, hull_dict=miri_mrs_hull_dict['N']['CH3'],
    line_width=4, line_style='-', color='tab:green')
ax_zoom_in_miri_se.imshow(img_zoom_in_miri_se)
# erase axis
plotting_tools.WCSPlottingTools.arr_axis_params(
    ax=ax_zoom_in_miri_se, ra_tick_label=False, dec_tick_label=False, ra_axis_label=' ', dec_axis_label=' ',)
# display hulls
plotting_tools.WCSPlottingTools.plot_obs_hull(
    ax=ax_zoom_in_miri_se, wcs=wcs_zoom_in_miri_se, hull_dict=nirspec_hull_dict['SE']['G140M/F100LP'],
    line_width=4, line_style='-', color='tab:orange')
plotting_tools.WCSPlottingTools.plot_obs_hull(
    ax=ax_zoom_in_miri_se, wcs=wcs_zoom_in_miri_se, hull_dict=miri_mrs_hull_dict['SE']['CH3'],
    line_width=4, line_style='-', color='tab:orange')



# radio data
# get color map norm
norm_tt0 = plotting_tools.ColorBarTools.compute_cbar_norm(
    vmin_vmax=(np.nanmin(cutout_tt0_ne.data), np.nanmax(cutout_tt0_ne.data)), log_scale=True)
ax_zoom_in_tt0_ne.imshow(cutout_tt0_ne.data, cmap='Grays', norm=norm_tt0)
# erase axis
plotting_tools.WCSPlottingTools.arr_axis_params(
    ax=ax_zoom_in_tt0_ne, ra_tick_label=False, dec_tick_label=False, ra_axis_label=' ', dec_axis_label=' ',)
# display hulls
plotting_tools.WCSPlottingTools.plot_obs_hull(
    ax=ax_zoom_in_tt0_ne, wcs=cutout_tt0_ne.wcs, hull_dict=nirspec_hull_dict['NE']['G140M/F100LP'],
    line_width=4, line_style='-', color='tab:red')
plotting_tools.WCSPlottingTools.plot_obs_hull(
    ax=ax_zoom_in_tt0_ne, wcs=cutout_tt0_ne.wcs, hull_dict=miri_mrs_hull_dict['NE']['CH3'],
    line_width=4, line_style='-', color='tab:red')
plotting_tools.StrTools.display_text_in_corner(
    ax=ax_zoom_in_tt0_ne, text='VLA L-band', fontsize=fontsize_medium, text_color='red', x_frac=0.98, y_frac=0.96,
    horizontal_alignment='right')
# get color map norm
norm_tt0 = plotting_tools.ColorBarTools.compute_cbar_norm(
    vmin_vmax=(np.nanmin(cutout_tt0_n.data), np.nanmax(cutout_tt0_n.data)), log_scale=True)
ax_zoom_in_tt0_n.imshow(cutout_tt0_n.data, cmap='Grays', norm=norm_tt0)
# erase axis
plotting_tools.WCSPlottingTools.arr_axis_params(
    ax=ax_zoom_in_tt0_n, ra_tick_label=False, dec_tick_label=False, ra_axis_label=' ', dec_axis_label=' ',)
# display hulls
plotting_tools.WCSPlottingTools.plot_obs_hull(
    ax=ax_zoom_in_tt0_n, wcs=cutout_tt0_n.wcs, hull_dict=nirspec_hull_dict['N']['G140M/F100LP'],
    line_width=4, line_style='-', color='tab:green')
plotting_tools.WCSPlottingTools.plot_obs_hull(
    ax=ax_zoom_in_tt0_n, wcs=cutout_tt0_n.wcs, hull_dict=miri_mrs_hull_dict['N']['CH3'],
    line_width=4, line_style='-', color='tab:green')

# get color map norm
norm_tt0 = plotting_tools.ColorBarTools.compute_cbar_norm(
    vmin_vmax=(np.nanmin(cutout_tt0_se.data), np.nanmax(cutout_tt0_se.data)), log_scale=True)
ax_zoom_in_tt0_se.imshow(cutout_tt0_se.data, cmap='Grays', norm=norm_tt0)
# erase axis
plotting_tools.WCSPlottingTools.arr_axis_params(
    ax=ax_zoom_in_tt0_se, ra_tick_label=False, dec_tick_label=False, ra_axis_label=' ', dec_axis_label=' ',)
# display hulls
plotting_tools.WCSPlottingTools.plot_obs_hull(
    ax=ax_zoom_in_tt0_se, wcs=cutout_tt0_se.wcs, hull_dict=nirspec_hull_dict['SE']['G140M/F100LP'],
    line_width=4, line_style='-', color='tab:orange')
plotting_tools.WCSPlottingTools.plot_obs_hull(
    ax=ax_zoom_in_tt0_se, wcs=cutout_tt0_se.wcs, hull_dict=miri_mrs_hull_dict['SE']['CH3'],
    line_width=4, line_style='-', color='tab:orange')




# x-ray data
ax_zoom_in_x_ray_ne.imshow(x_ray_rgb_ne)
# erase axis
plotting_tools.WCSPlottingTools.arr_axis_params(
    ax=ax_zoom_in_x_ray_ne, ra_tick_label=False, dec_tick_label=False, ra_axis_label=' ', dec_axis_label=' ',)
# display hulls
plotting_tools.WCSPlottingTools.plot_obs_hull(
    ax=ax_zoom_in_x_ray_ne, wcs=cutout_0p5to2_ne.wcs, hull_dict=nirspec_hull_dict['NE']['G140M/F100LP'],
    line_width=4, line_style='-', color='tab:red')
plotting_tools.WCSPlottingTools.plot_obs_hull(
    ax=ax_zoom_in_x_ray_ne, wcs=cutout_0p5to2_ne.wcs, hull_dict=miri_mrs_hull_dict['NE']['CH3'],
    line_width=4, line_style='-', color='tab:red')
plotting_tools.StrTools.display_text_in_corner(
    ax=ax_zoom_in_x_ray_ne, text='0.5-2 eV', fontsize=fontsize_medium, text_color='red', x_frac=0.98, y_frac=0.96,
    horizontal_alignment='right')
plotting_tools.StrTools.display_text_in_corner(
    ax=ax_zoom_in_x_ray_ne, text='2-7 eV', fontsize=fontsize_medium, text_color='blue', x_frac=0.98, y_frac=0.88,
    horizontal_alignment='right')

ax_zoom_in_x_ray_n.imshow(x_ray_rgb_n)
# erase axis
plotting_tools.WCSPlottingTools.arr_axis_params(
    ax=ax_zoom_in_x_ray_n, ra_tick_label=False, dec_tick_label=False, ra_axis_label=' ', dec_axis_label=' ',)
# display hulls
plotting_tools.WCSPlottingTools.plot_obs_hull(
    ax=ax_zoom_in_x_ray_n, wcs=cutout_0p5to2_n.wcs, hull_dict=nirspec_hull_dict['N']['G140M/F100LP'],
    line_width=4, line_style='-', color='tab:green')
plotting_tools.WCSPlottingTools.plot_obs_hull(
    ax=ax_zoom_in_x_ray_n, wcs=cutout_0p5to2_n.wcs, hull_dict=miri_mrs_hull_dict['N']['CH3'],
    line_width=4, line_style='-', color='tab:green')

ax_zoom_in_x_ray_se.imshow(x_ray_rgb_se)
# erase axis
plotting_tools.WCSPlottingTools.arr_axis_params(
    ax=ax_zoom_in_x_ray_se, ra_tick_label=False, dec_tick_label=False, ra_axis_label=' ', dec_axis_label=' ',)
# display hulls
plotting_tools.WCSPlottingTools.plot_obs_hull(
    ax=ax_zoom_in_x_ray_se, wcs=cutout_0p5to2_se.wcs, hull_dict=nirspec_hull_dict['SE']['G140M/F100LP'],
    line_width=4, line_style='-', color='tab:orange')
plotting_tools.WCSPlottingTools.plot_obs_hull(
    ax=ax_zoom_in_x_ray_se, wcs=cutout_0p5to2_se.wcs, hull_dict=miri_mrs_hull_dict['SE']['CH3'],
    line_width=4, line_style='-', color='tab:orange')




fig.savefig('plot_output/overview_jwst_spec.png')


