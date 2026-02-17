"""
Script to develop how to check the observational coverage of a PHANGS target
"""

import numpy as np
import os
import pickle
from werkzeugkiste import helper_func
from obszugang import obs_info
from obszugang import spec_access
import matplotlib.pyplot as plt

target_list = obs_info.phangs_muse_galaxy_list


for target in target_list:

    # if os.path.isfile('data_output/%s_muse_obs_hull_dict.pickle' % target):
    #     continue

    # target = 'ngc2835'

    # now get muse bands
    phangs_spec = spec_access.SpecAccess(spec_target_name=target)

    print(target)

    obs_hull_dict = {}
    for res in ['copt', 'native', '150pc', '15asec']:

        # load H-alpha muse DAP
        phangs_spec.load_muse_dap_map(res=res)

        # plt.imshow(phangs_spec.muse_dap_map_data['dap_map_data_copt_fiducial_HA6562_FLUX'])
        # plt.show()

        data = phangs_spec.muse_dap_map_data['dap_map_data_copt_fiducial_HA6562_FLUX']
        wcs = phangs_spec.muse_dap_map_data['dap_map_wcs_copt_fiducial_HA6562_FLUX']

        mask_coverage = np.array(np.invert(np.isnan( data)), dtype=float)

        hull_dict = helper_func.GeometryTools.contour2hull(data_array=mask_coverage, level=0, contour_index=0, n_max_rejection_vertice=100)
        print(res, ' n of hulls: ', len(hull_dict.keys()))


        # fig = plt.figure(figsize=(20, 20))
        # ax_img = fig.add_axes((0.05, 0.05, 0.945, 0.945), projection=wcs)
        # ax_img.imshow(data)



        # now save the hull points as coordinates
        hull_coord_dict = {}
        for idx in hull_dict.keys():

            # ax_img.plot(hull_dict[idx]['x_convex_hull'], hull_dict[idx]['y_convex_hull'])

            # transform into coordinates
            coordinates = wcs.pixel_to_world(hull_dict[idx]['x_convex_hull'], hull_dict[idx]['y_convex_hull'])
            ra = coordinates.ra.deg
            dec = coordinates.dec.deg
            hull_coord_dict.update({idx: {'ra': ra, 'dec': dec}})
        obs_hull_dict.update({res: hull_coord_dict})

        # plt.show()




    # exit()


    # save dictionary
    if not os.path.isdir('data_output'):
        os.makedirs('data_output')

    with open('data_output/%s_muse_obs_hull_dict.pickle' % target, 'wb') as file_name:
        pickle.dump(obs_hull_dict, file_name)



