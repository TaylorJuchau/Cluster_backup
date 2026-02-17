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


miri_mrs_data_version_name = 'karin_reduction'


target_list = list(obs_info.miri_mrs_obs_target_dict.keys())

print(target_list)
for target in target_list:

    obs_hull_dict = {}

    for region_name in obs_info.miri_mrs_obs_target_dict[target].keys():
        print(target, region_name)
        obs_hull_dict.update({region_name: {}})

        # get data access
        analysis_access = spec_access.SpecAccess(spec_target_name=target, miri_mrs_data_ver=miri_mrs_data_version_name)

        print('channels, available: ', obs_info.miri_mrs_obs_target_dict[target][region_name][miri_mrs_data_version_name]['channels'])

        for channel in obs_info.miri_mrs_obs_target_dict[target][region_name][miri_mrs_data_version_name]['channels']:
            print('channel ', channel)
            analysis_access.load_miri_mrs_cube(region_name=region_name, channel=channel)

            spec_cube = analysis_access.miri_mrs_datacube_data[region_name]['data_cube_%s' % channel]
            spec_cube_img_data = np.nansum(spec_cube, axis=0)

            new_spec_cube_img_data = np.zeros((spec_cube_img_data.shape[0] + 2, spec_cube_img_data.shape[1] + 2))
            new_spec_cube_img_data[1:-1, 1:-1] = spec_cube_img_data
            cube_wcs_2d = analysis_access.miri_mrs_datacube_data[region_name]['wcs_2d_%s' % channel]
            mask_covered_pixels = new_spec_cube_img_data > 0

            hull_dict = helper_func.GeometryTools.contour2hull(
                data_array=mask_covered_pixels, level=0, contour_index=0, n_max_rejection_vertice=150)

            # fig = plt.figure(figsize=(20, 20))
            # ax_img = fig.add_axes((0.05, 0.05, 0.945, 0.945), projection=cube_wcs_2d)
            # ax_img.imshow(np.log10(new_spec_cube_img_data))
            print(channel, ' n of hulls: ', len(hull_dict.keys()))

            hull_coord_dict = {}
            for idx in hull_dict.keys():
                # ax_img.plot(hull_dict[idx]['x_convex_hull'], hull_dict[idx]['y_convex_hull'])
                # transform into coordinates
                coordinates = cube_wcs_2d.pixel_to_world(hull_dict[idx]['x_convex_hull'] - 1, hull_dict[idx]['y_convex_hull'] - 1)
                ra = coordinates.ra.deg
                dec = coordinates.dec.deg
                hull_coord_dict.update({idx: {'ra': ra, 'dec': dec}})
                # plt.plot(hull_dict[idx]['x_convex_hull'] - 1, hull_dict[idx]['y_convex_hull'] - 1)

            obs_hull_dict[region_name].update({channel: hull_coord_dict})
            # plt.show()
    with open('data_output/%s_miri_mrs_obs_hull_dict_%s.pickle' % (target, miri_mrs_data_version_name),
              'wb') as file_name:
        pickle.dump(obs_hull_dict, file_name)
