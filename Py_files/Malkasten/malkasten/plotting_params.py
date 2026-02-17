"""
This script gathers functions to help plotting procedures
"""

# import os
# from pathlib import Path
import astropy.units as u
# from astropy.coordinates import SkyCoord
# from astropy.stats import SigmaClip
# from astropy.visualization.wcsaxes import SphericalCircle, Quadrangle, add_beam
#
# from matplotlib.colors import Normalize, LogNorm
# from matplotlib.colorbar import ColorbarBase
# from matplotlib import patheffects
import matplotlib.pyplot as plt
# from matplotlib.patches import ConnectionPatch, Ellipse
#
# from matplotlib import text as mtext
# import math
#
#
# import decimal
#
# from photutils.segmentation import make_2dgaussian_kernel
# from photutils.background import Background2D, MedianBackground
# from astropy.stats import sigma_clipped_stats
# from astropy.visualization.mpl_normalize import ImageNormalize
# from astropy.visualization import SqrtStretch
#
# from astropy.convolution import convolve
#
# from regions import PixCoord, RectanglePixelRegion
#
# import numpy as np
#
# from malkasten import multicolorfits as mcf
# from werkzeugkiste import helper_func, phys_params, spec_tools
# from obszugang.obs_tools import ObsTools
# from obszugang import gal_access
# from sternenstaub import dust_tools



nuvb_label_dict = {
    1: {'offsets': [0.25, -0.1], 'ha': 'center', 'va': 'bottom', 'label': r'1,2,3 Myr'},
    5: {'offsets': [0.05, 0.1], 'ha': 'right', 'va': 'top', 'label': r'5 Myr'},
    10: {'offsets': [0.1, -0.2], 'ha': 'left', 'va': 'bottom', 'label': r'10 Myr'},
    100: {'offsets': [-0.1, 0.0], 'ha': 'right', 'va': 'center', 'label': r'100 Myr'}
}
ub_label_dict = {
    1: {'offsets': [0.2, -0.1], 'ha': 'center', 'va': 'bottom', 'label': r'1,2,3 Myr'},
    5: {'offsets': [0.05, 0.1], 'ha': 'right', 'va': 'top', 'label': r'5 Myr'},
    10: {'offsets': [0.1, -0.2], 'ha': 'left', 'va': 'bottom', 'label': r'10 Myr'},
    100: {'offsets': [-0.1, 0.0], 'ha': 'right', 'va': 'center', 'label': r'100 Myr'}
}
bv_label_dict = {
    1: {'offsets': [0.2, -0.1], 'ha': 'center', 'va': 'bottom', 'label': r'1,2,3 Myr'},
    5: {'offsets': [0.05, 0.1], 'ha': 'right', 'va': 'top', 'label': r'5 Myr'},
    10: {'offsets': [0.1, -0.1], 'ha': 'left', 'va': 'bottom', 'label': r'10 Myr'},
    100: {'offsets': [-0.1, 0.1], 'ha': 'right', 'va': 'center', 'label': r'100 Myr'}
}

nuvb_annotation_dict = {
    500: {'offset': [-0.5, +0.0], 'label': '500 Myr', 'ha': 'right', 'va': 'center'},
    1000: {'offset': [-0.7, +0.5], 'label': '1 Gyr', 'ha': 'right', 'va': 'center'},
    13750: {'offset': [+0.05, -0.9], 'label': '13.8 Gyr', 'ha': 'left', 'va': 'center'}
}
ub_annotation_dict = {
    500: {'offset': [-0.5, +0.0], 'label': '500 Myr', 'ha': 'right', 'va': 'center'},
    1000: {'offset': [-0.5, +0.5], 'label': '1 Gyr', 'ha': 'right', 'va': 'center'},
    13750: {'offset': [-0.0, -0.7], 'label': '13.8 Gyr', 'ha': 'left', 'va': 'center'}
}
bv_annotation_dict = {
    500: {'offset': [-0.5, +0.3], 'label': '500 Myr', 'ha': 'right', 'va': 'center'},
    1000: {'offset': [-0.5, +0.5], 'label': '1 Gyr', 'ha': 'right', 'va': 'center'},
    13750: {'offset': [-0.0, -0.7], 'label': '13.8 Gyr', 'ha': 'left', 'va': 'center'}
}


# define a color seelction
color_list_dark2 = plt.get_cmap('Dark2')(range(8))
color_list_set2 = plt.get_cmap('Set2')(range(8))
color_list_set3 = plt.get_cmap('Set3')(range(12))
color_list_paired = plt.get_cmap('Paired')(range(12))
color_list_tab10 = plt.get_cmap('tab10')(range(10))
color_list_pastel2 = plt.get_cmap('Pastel2')(range(8))

color_list_rainbow = plt.get_cmap('rainbow')
color_list_hsv = plt.get_cmap('hsv')



ppxf_overview_plot_param_dict = {
    'figsize': (45, 60),

    'fontsize_label': 45,
    'fontsize_title': 55,

    # overview spec
    'overview_spec_left_align': 0.07,
    'overview_spec_bottom_align': 0.8,
    'overview_spec_width': 0.925,
    'overview_spec_height': 0.195,
    #
    'overview_spec_residuals_left_align': 0.07,
    'overview_spec_residuals_bottom_align': 0.73,
    'overview_spec_residuals_width': 0.925,
    'overview_spec_residuals_height': 0.065,
    #
    # light weight frac
    'light_weight_frace_left_align': 0.07,
    'light_weight_frace_bottom_align': 0.575,
    'light_weight_frace_width': 0.60,
    'light_weight_frace_height': 0.125,
    # colorbar
    'cbar_light_weight_frace_left_align': 0.68,
    'cbar_light_weight_frace_bottom_align': 0.595,
    'cbar_light_weight_frace_width': 0.015,
    'cbar_light_weight_frace_height': 0.085,

    'fit_param_text_pos': (0.76, 0.695),

    'light_weight_cmap': 'Reds',

    # units
    'x_unit': u.Angstrom,
    'y_unit': u.erg / (u.Angstrom * u.s * u.cm * u.cm),

    # emission lines
    'line_list_1': ['Hbeta', '[OIII]5007_d'],
    # axis
    'line_axis_1_left_align': 0.07,
    'line_axis_1_bottom_align': 0.415,
    'line_axis_1_width': 0.40,
    'line_axis_1_height': 0.13,
    # residuals
    'line_axis_residuals_1_left_align': 0.07,
    'line_axis_residuals_1_bottom_align': 0.38,
    'line_axis_residuals_1_width': 0.40,
    'line_axis_residuals_1_height': 0.03,

    'line_list_2': ['HeI5876'],
    # axis
    'line_axis_2_left_align': 0.535,
    'line_axis_2_bottom_align': 0.415,
    'line_axis_2_width': 0.2,
    'line_axis_2_height': 0.13,
    # residuals
    'line_axis_residuals_2_left_align': 0.535,
    'line_axis_residuals_2_bottom_align': 0.38,
    'line_axis_residuals_2_width': 0.2,
    'line_axis_residuals_2_height': 0.03,


    'line_list_3': ['[OI]6300_d'],
    # axis
    'line_axis_3_left_align': 0.795,
    'line_axis_3_bottom_align': 0.415,
    'line_axis_3_width': 0.2,
    'line_axis_3_height': 0.13,
    # residuals
    'line_axis_residuals_3_left_align': 0.795,
    'line_axis_residuals_3_bottom_align': 0.38,
    'line_axis_residuals_3_width': 0.2,
    'line_axis_residuals_3_height': 0.03,


    'line_list_4': ['Halpha', '[NII]6583_d'],
    # axis
    'line_axis_4_left_align': 0.07,
    'line_axis_4_bottom_align': 0.24,
    'line_axis_4_width': 0.40,
    'line_axis_4_height': 0.13,
    # residuals
    'line_axis_residuals_4_left_align': 0.07,
    'line_axis_residuals_4_bottom_align': 0.205,
    'line_axis_residuals_4_width': 0.40,
    'line_axis_residuals_4_height': 0.03,


    'line_list_5': ['[SII]6716', '[SII]6731'],
    # axis
    'line_axis_5_left_align': 0.07,
    'line_axis_5_bottom_align': 0.065,
    'line_axis_5_width': 0.40,
    'line_axis_5_height': 0.13,
    # residuals
    'line_axis_residuals_5_left_align': 0.07,
    'line_axis_residuals_5_bottom_align': 0.03,
    'line_axis_residuals_5_width': 0.40,
    'line_axis_residuals_5_height': 0.03,



    'line_data_color': 'black',
    'line_data_line_width': 5,
    'best_fit_color': 'tab:red',
    'best_fit_line_width': 5,
    'line_comp_color_list': ['', 'tab:blue', 'tab:orange'],
    'line_comp_line_width': 5,

    'tick_label_color': 'black',
    'label_color': 'black',
    'line_label_color': 'black',
    'line_label_path_eff_color': 'tab:red',
    'display_x_label': True,



}



