"""
Script to develop how to check the observational coverage of a PHANGS target
"""

import numpy as np
import os
import pickle
from werkzeugkiste import helper_func
from obszugang import obs_info, phot_access, ObsTools

import matplotlib.pyplot as plt

target_list = list(obs_info.full_hst_galaxy_list)


for target in target_list:
    # if os.path.isfile('data_output/%s_hst_obs_hull_dict.pickle' % target):
    #     continue

    # target = 'ngc2835'

    if target == 'ngc1510':
        continue

    # now get hst bands
    phangs_phot = phot_access.PhotAccess(phot_target_name=target)

    # get band list
    band_list = ObsTools.get_hst_obs_band_list(target=target)

    print(target, ' bands, available: ', band_list)

    obs_hull_dict = {}
    for band in band_list:

        exp_file_name = phangs_phot.get_hst_img_file_name(band=band, file_type='wht')
        if phangs_phot.phot_hst_target_name == 'ngc5194': hdu_number = 'WHT'
        else: hdu_number = 0

        data, header, wcs = helper_func.FileTools.load_img(file_name=exp_file_name, hdu_number=hdu_number)
        mask_covered_pixels = data > 0

        # plt.imshow(np.log10(data))

        # fig = plt.figure(figsize=(20, 20))
        # ax_img = fig.add_axes((0.05, 0.05, 0.945, 0.945), projection=wcs)
        # ax_img.imshow(np.log10(data))

        hull_dict = helper_func.GeometryTools.contour2hull(data_array=data, level=0, contour_index=0, n_max_rejection_vertice=1000)

        print(band, ' n of hulls: ', len(hull_dict.keys()))

        # now save the hull points as coordinates
        hull_coord_dict = {}
        for idx in hull_dict.keys():

            # ax_img.plot(hull_dict[idx]['x_convex_hull'], hull_dict[idx]['y_convex_hull'])

            # transform into coordinates
            coordinates = wcs.pixel_to_world(hull_dict[idx]['x_convex_hull'], hull_dict[idx]['y_convex_hull'])
            ra = coordinates.ra.deg
            dec = coordinates.dec.deg
            hull_coord_dict.update({idx: {'ra': ra, 'dec': dec}})
            # plt.plot(hull_dict[idx]['x_convex_hull'] - 1, hull_dict[idx]['y_convex_hull'] - 1)

        obs_hull_dict.update({band: hull_coord_dict})
        # plt.show()
    # save dictionary

    # exit()
    if not os.path.isdir('data_output'):
        os.makedirs('data_output')

    with open('data_output/%s_hst_obs_hull_dict.pickle' % target, 'wb') as file_name:
        pickle.dump(obs_hull_dict, file_name)



