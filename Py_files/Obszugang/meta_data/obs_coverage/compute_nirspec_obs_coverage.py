"""
Script to develop how to check the observational coverage of a PHANGS target
"""

import numpy as np
import os
import pickle
from werkzeugkiste import helper_func
from obszugang import obs_info, spec_access, ObsTools
import matplotlib.pyplot as plt


plot_flag = True


nirspec_data_ver = 'karin_reduction_v1_oct2024'

grating = 'G140M/F100LP'



target_list = list(obs_info.nirspec_obs_target_dict.keys())

print(target_list)
for target in target_list:

    obs_hull_dict = {}

    for region_name in obs_info.nirspec_obs_target_dict[target].keys():
        print(target, region_name)
        obs_hull_dict.update({region_name: {}})

        # get data access
        analysis_access = spec_access.SpecAccess(spec_target_name=target, nirspec_data_ver=nirspec_data_ver)

        print('gratings, available: ', obs_info.nirspec_obs_target_dict[target][region_name][nirspec_data_ver]['gratings'])

        for grating in obs_info.nirspec_obs_target_dict[target][region_name][nirspec_data_ver]['gratings']:
            print('grating ', grating)
            analysis_access.load_nirspec_cube(region_name=region_name, grating=grating)

            spec_cube = analysis_access.nirspec_datacube_data[region_name]['data_cube_%s' % grating]
            spec_cube_img_data = np.nansum(spec_cube, axis=0)
            cube_wcs_2d = analysis_access.nirspec_datacube_data[region_name]['wcs_2d_%s' % grating]
            mask_covered_pixels = spec_cube_img_data > 0

            hull_dict = helper_func.GeometryTools.contour2hull(
                data_array=mask_covered_pixels, level=0, contour_index=0, n_max_rejection_vertice=100)

            # fig = plt.figure(figsize=(20, 20))
            # ax_img = fig.add_axes((0.05, 0.05, 0.945, 0.945), projection=cube_wcs_2d)
            # ax_img.imshow(np.log10(spec_cube_img_data))
            print(grating, ' n of hulls: ', len(hull_dict.keys()))

            hull_coord_dict = {}
            for idx in hull_dict.keys():
                # ax_img.plot(hull_dict[idx]['x_convex_hull'], hull_dict[idx]['y_convex_hull'])
                # transform into coordinates
                coordinates = cube_wcs_2d.pixel_to_world(hull_dict[idx]['x_convex_hull'], hull_dict[idx]['y_convex_hull'])
                ra = coordinates.ra.deg
                dec = coordinates.dec.deg
                hull_coord_dict.update({idx: {'ra': ra, 'dec': dec}})
                # plt.plot(hull_dict[idx]['x_convex_hull'] - 1, hull_dict[idx]['y_convex_hull'] - 1)

            obs_hull_dict[region_name].update({grating: hull_coord_dict})

    with open('data_output/%s_nirspec_obs_hull_dict_%s.pickle' % (target, nirspec_data_ver),
              'wb') as file_name:
        pickle.dump(obs_hull_dict, file_name)
