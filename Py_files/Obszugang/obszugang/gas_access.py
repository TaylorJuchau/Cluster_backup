"""
Construct a data access structure for CO and HI gas observations
"""
from pathlib import Path
import pickle
import numpy as np
import astropy.units as u
from astropy.coordinates import SkyCoord
from werkzeugkiste import helper_func, phys_params
from obszugang import access_config, gal_access


class GasAccess:
    """
    Access class to organize data structure of ALMA observations
    """

    def __init__(self, gas_target_name=None):
        """

        Parameters
        ----------
        gas_target_name : str
            Default None. Target name
        gas_target_name : str
            Default None. Target name used for Hs observation
        """

        # get target specifications
        # check if the target names are compatible

        if gas_target_name is not None:
            gas_target_name = helper_func.FileTools.target_name_no_directions(target=gas_target_name)

        # if (gas_target_name not in phangs_info.phangs_alma_galaxy_list) & (gas_target_name is not None):
        #     # raise AttributeError('The target %s is not in the PHANGS ALMA sample or has not been added to '
        #     #                      'the current package version' % gas_target_name)

        self.gas_target_name = gas_target_name

        # loaded data dictionaries
        self.alma_data = {}

        # get path to observation coverage hulls
        self.path2obs_cover_hull = (Path(__file__).parent.parent.absolute() / 'meta_data' / 'obs_coverage' /
                                    'data_output')

        super().__init__()

    def get_alma_co_mom_map_file_name(self, telescope_config='12m+7m+tp', co_line='co21', line_measure_type='broad',
                                      mom='mom0', res=150):
        """

        Parameters
        ----------
        telescope_config : str
        co_line : str
        line_measure_type : str
        mom : str
        res : int or str
        Returns
        -------
        data_file_path : ``Path``
        """
        if isinstance(res, int):
            res_str = str(res) + 'pc_'
        elif res == 'native':
            res_str = ''
        else:
            raise KeyError('res must be either int with the resolution in pc or native')

        file_path = (Path(access_config.phangs_config_dict['alma_data_path']) /
                     ('delivery_%s' % access_config.phangs_config_dict['alma_data_ver']) / self.gas_target_name)
        file_name = '%s_%s_%s_%s%s_%s.fits' % (self.gas_target_name, telescope_config, co_line, res_str,
                                               line_measure_type, mom)

        return file_path / file_name

    def get_alma_co_cube_file_name(self, telescope_config='12m+7m+tp', co_line='co21', err=False):
        """

        Parameters
        ----------
        telescope_config : str
        co_line : str
        err : bool
        Returns
        -------
        data_file_path : ``Path``
        """

        file_path = (Path(access_config.phangs_config_dict['alma_data_path']) /
                     ('delivery_%s' % access_config.phangs_config_dict['alma_data_ver']) / self.gas_target_name)
        if err:
            file_name = '%s_%s_%s_noise.fits' % (self.gas_target_name, telescope_config, co_line)
        else:
            file_name = '%s_%s_%s.fits' % (self.gas_target_name, telescope_config, co_line)

        return file_path / file_name

    def get_alma_co21_conv_map_file_name(self, alpha_co_method='S20_MUSEGPR'):
        """

        Parameters
        ----------
        alpha_co_method : str
            can be S20_MUSEGPR, S20_scaling, B13_MUSEGPR, B13_scaling,
            see Sun+ (2020, ApJ, 892, 148) or the read me file
        Returns
        -------
        data_file_path : ``Path``
        """

        file_path = (Path(access_config.phangs_config_dict['alma_conv_map_data_path']) /
                     access_config.phangs_config_dict['alma_conv_map_data_ver'])
        file_name = '%s_alphaCO21_%s.fits' % (self.gas_target_name.upper(), alpha_co_method)

        return file_path / file_name

    def load_alma_co_mom_maps(self,  telescope_config='12m+7m+tp', co_line='co21', line_measure_type='broad',
                              mom='mom0', res=150, reload=False):
        """

        Parameters
        ----------
        telescope_config : str
        co_line : str
        line_measure_type : str
        mom : str
        res : int or str
        reload : bool
        """
        # check if data is already loaded
        if (('%s_%s_data_img' % (mom, str(res))) in self.alma_data.keys()) & (not reload):
            return None

        file_name = self.get_alma_co_mom_map_file_name(telescope_config=telescope_config, co_line=co_line,
                                                       line_measure_type=line_measure_type, mom=mom, res=res)

        img_data, img_header, img_wcs = helper_func.FileTools.load_img(file_name=file_name)

        self.alma_data.update({
            '%s_%s_%s_%s_%s_data_img' % (telescope_config, co_line, mom, line_measure_type, str(res)): img_data,
            '%s_%s_%s_%s_%s_header_img' % (telescope_config, co_line, mom, line_measure_type, str(res)): img_header,
            '%s_%s_%s_%s_%s_wcs_img' % (telescope_config, co_line, mom, line_measure_type, str(res)): img_wcs})

    def load_alma_co_data_cubes(self,  telescope_config='12m+7m+tp', co_line='co21', err=False, reload=False):
        """

        Parameters
        ----------
        telescope_config : str
        co_line : str
        reload : bool
        err : bool
        """
        # check if data is already loaded
        if (('%s_%s_data_cube' % (telescope_config, co_line)) in self.alma_data.keys()) & (not reload) & (not err):
            return None
        # in case loading error
        if (('%s_%s_data_cube_err' % (telescope_config, co_line)) in self.alma_data.keys()) & (not reload) & err:
            return None

        file_name = self.get_alma_co_cube_file_name(telescope_config=telescope_config, co_line=co_line, err=err)

        data_cube, x_axis_values, header, wcs3d, wcs2d = helper_func.FileTools.load_alma_cube(file_name=file_name)
        if err:
            self.alma_data.update({'%s_%s_data_cube_err' % (telescope_config, co_line): data_cube,
                                   '%s_%s_x_axis_values_cube_err' % (telescope_config, co_line): x_axis_values,
                                   '%s_%s_header_cube_err' % (telescope_config, co_line): header,
                                   '%s_%s_wcs3d_cube_err' % (telescope_config, co_line): wcs3d,
                                   '%s_%s_wcs2d_cube_err' % (telescope_config, co_line): wcs2d})
        else:
            self.alma_data.update({'%s_%s_data_cube' % (telescope_config, co_line): data_cube,
                                   '%s_%s_x_axis_values_cube' % (telescope_config, co_line): x_axis_values,
                                   '%s_%s_header_cube' % (telescope_config, co_line): header,
                                   '%s_%s_wcs3d_cube' % (telescope_config, co_line): wcs3d,
                                   '%s_%s_wcs2d_cube' % (telescope_config, co_line): wcs2d})

    def get_alma_co_mom_map_beam(self, telescope_config='12m+7m+tp', co_line='co21', line_measure_type='broad',
                                 res=150):
            """

            Parameters
            ----------
            telescope_config : str
            co_line : str
            line_measure_type : str
            mom : str
            res : int or str
            reload : bool
            """
            # make sure the data cube is loaded
            self.load_alma_co_mom_maps(telescope_config=telescope_config, co_line=co_line,
                                       line_measure_type=line_measure_type, mom='mom0', res=res)
            header = self.alma_data[
                '%s_%s_mom0_%s_%s_header_img' % (telescope_config, co_line, line_measure_type, str(res))]

            beam_major = header['BMAJ']
            beam_minor = header['BMIN']
            beam_pa = header['BPA']

            return beam_major, beam_minor, beam_pa

    def get_alma_co_cube_beam(self, telescope_config='12m+7m+tp', co_line='co21'):
            """

            Parameters
            ----------
            telescope_config : str
            co_line : str
            """
            # make sure the data cube is loaded
            self.load_alma_co_data_cubes(telescope_config=telescope_config, co_line=co_line)
            header = self.alma_data['%s_%s_header_cube' % (telescope_config, co_line)]

            beam_major = header['BMAJ']
            beam_minor = header['BMIN']
            beam_pa = header['BPA']

            return beam_major, beam_minor, beam_pa

    def extract_alma_co_spec_circ_app(self, ra, dec, rad_arcsec, telescope_config='12m+7m+tp', co_line='co21'):
        """

        Parameters
        ----------
        ra, dec : float
        rad_arcsec : float
            radius in arcseconds
        wave_range : tuple or None
        res : str

        Return
        ---------
        spec_dict : dict
            Dictionary with all spectra related values, including a ln rebinned spectrum
        """
        # make sure the data cube is loaded
        self.load_alma_co_data_cubes(telescope_config=telescope_config, co_line=co_line)
        # make sure error is loaded
        self.load_alma_co_data_cubes(telescope_config=telescope_config, co_line=co_line, err=True)

        # get the coordinates
        coords_exctract = SkyCoord(ra=ra*u.deg, dec=dec*u.deg)
        obj_coords_alma_pix = (
            self.alma_data['%s_%s_wcs2d_cube' % (telescope_config, co_line)].world_to_pixel(coords_exctract))
        # get the extraction radius in pixel
        selection_radius_pix = helper_func.CoordTools.transform_world2pix_scale(
            length_in_arcsec=rad_arcsec, wcs=self.alma_data['%s_%s_wcs2d_cube' % (telescope_config, co_line)])

        # select spaxels which are inside the selected radius
        x_lin_alma = np.linspace(1, self.alma_data['%s_%s_data_cube' % (telescope_config, co_line)].shape[2],
                                 self.alma_data['%s_%s_data_cube' % (telescope_config, co_line)].shape[2])
        y_lin_alma = np.linspace(1, self.alma_data['%s_%s_data_cube' % (telescope_config, co_line)].shape[1],
                                 self.alma_data['%s_%s_data_cube' % (telescope_config, co_line)].shape[1])
        x_data_alma, y_data_alma = np.meshgrid(x_lin_alma, y_lin_alma)
        mask_spectrum = (np.sqrt((x_data_alma - obj_coords_alma_pix[0]) ** 2 +
                                 (y_data_alma - obj_coords_alma_pix[1]) ** 2) < selection_radius_pix)

        spec_flx = np.sum(self.alma_data['%s_%s_data_cube' % (telescope_config, co_line)][:, mask_spectrum], axis=1)
        spec_flx_err = np.sqrt(np.sum(
            self.alma_data['%s_%s_data_cube_err' % (telescope_config, co_line)][:, mask_spectrum] ** 2, axis=1))

        return {
            # general description of spectrum
            'rad_arcsec': rad_arcsec,
            # spectrum
            'x_axis_values': self.alma_data['%s_%s_x_axis_values_cube' % (telescope_config, co_line)],
            'spec_flx': spec_flx,
            'spec_flx_err': spec_flx_err
        }


    def load_alma_conv_map(self, alpha_co_method='S20_MUSEGPR', reload=False):
        """
        PHANGS-ALMA CO-to-H2 Conversion Factor Maps

        Parameters
        ----------
        alpha_co_method : str
        reload : bool

        """
        if ('alpha_co21_data_img' in self.alma_data.keys()) & (not reload):
            return None
        file_name = self.get_alma_co21_conv_map_file_name(alpha_co_method=alpha_co_method)

        img_data, img_header, img_wcs = helper_func.FileTools.load_img(file_name=file_name)

        self.alma_data.update({'alpha_co21_data_img': img_data,
                               'alpha_co21_header_img': img_header,
                               'alpha_co21_wcs_img': img_wcs})

    def get_alma_h2_map(self, telescope_config='12m+7m+tp', co_line='co21', line_measure_type='broad', res=150,
                        alpha_co_method='S20_MUSEGPR'):
        """

        Parameters
        ----------
        telescope_config : str
        co_line : str
        line_measure_type : str
        res : int or str
        alpha_co_method : str

        """
        self.load_alma_co_mom_maps(telescope_config=telescope_config, co_line=co_line,
                                   line_measure_type=line_measure_type, res=res)
        self.load_alma_conv_map(alpha_co_method=alpha_co_method)

        alma_mom0_img = self.alma_data['%s_%s_mom0_%s_%s_data_img' %
                                       (telescope_config, co_line, line_measure_type, str(res))]
        alma_mom0_wcs = self.alma_data['%s_%s_mom0_%s_%s_wcs_img' %
                                       (telescope_config, co_line, line_measure_type, str(res))]

        # reproject conversion map to alma co map
        conv_map_reprojected = helper_func.CoordTools.reproject_image(
            data=self.alma_data['alpha_co21_data_img'], wcs=self.alma_data['alpha_co21_wcs_img'],
            new_wcs=alma_mom0_wcs, new_shape=alma_mom0_img.shape)

        return alma_mom0_img * conv_map_reprojected, alma_mom0_wcs

    def get_alma_co21_cloud_cat_file_name(self, res=150):
        """
        Parameters
        ----------
        res : int or str
        Returns
        -------
        data_file_path : ``Path``
        """
        # To do: add homogenized maps
        if isinstance(res, int):
            res_folder = 'matched_%spc' % str(res)
            res_str = '%spc_nativenoise' % str(res)
        elif res == 'native':
            res_folder = 'native'
            res_str = 'native'

        else:
            raise KeyError('res must be either int with the resolution in pc or native')

        '/media/benutzer/Extreme Pro/data/phangs_data_products/cloud_catalogs/v4p0_ST1p6/v4p0_gmccats/matched_150pc'
        file_path = (Path(access_config.phangs_config_dict['alma_cloud_cat_data_path']) /
                     ('%s_%s' % (access_config.phangs_config_dict['alma_data_ver'],
                                 access_config.phangs_config_dict['alma_cloud_cat_data_release_ver'])) /
                     ('%s_gmccats' % access_config.phangs_config_dict['alma_data_ver']) / res_folder)
        file_name = '%s_12m+7m+tp_co21_%s_props.fits' % (self.gas_target_name, res_str)

        return file_path / file_name

    def load_alma_cloud_cat(self, res=150, reload=False):
        """
        best description can be found in Rosolowsky+2021 2021MNRAS.502.1218R
        Parameters
        ----------
        res : int or str
        reload : bool

        """
        if ('alpha_cloud_cat_data' in self.alma_data.keys()) & (not reload):
            return None
        file_name = self.get_alma_co21_cloud_cat_file_name(res=res)
        cat_data, cat_header = helper_func.FileTools.load_fits_table(file_name=file_name, hdu_number=1)

        self.alma_data.update({'alpha_cloud_cat_data': cat_data, 'alpha_cloud_cat_header': cat_header})

    def get_cloud_coords(self, res=150):
        """
        Get positions from giant molecular clouds
        Parameters
        ----------
        res : int or str

        """
        # load cloud catalog
        self.load_alma_cloud_cat(res=res)
        return self.alma_data['alpha_cloud_cat_data']['XCTR_DEG'], self.alma_data['alpha_cloud_cat_data']['YCTR_DEG']

    def get_cloud_rad_pc(self, res=150):
        """
        Get radius from giant molecular clouds in pc
        Parameters
        ----------
        res : int or str

        """
        # load cloud catalog
        self.load_alma_cloud_cat(res=res)
        return self.alma_data['alpha_cloud_cat_data']['RAD3D_PC']

    def get_cloud_rad_arcsec(self, res=150):
        """
        Get radius from giant molecular clouds in arcsec
        Parameters
        ----------
        res : int or str

        """
        rad_cloud_pc = self. get_cloud_rad_pc(res=res)
        sample_access = gal_access.PhangsSampleAccess()
        target_dist_mpc = sample_access.get_target_dist(target=self.gas_target_name)
        central_target_pos = helper_func.CoordTools.get_target_central_simbad_coords(target_name=self.gas_target_name,
                                                                                     target_dist_mpc=target_dist_mpc)

        off_set_central_pos = SkyCoord(ra=central_target_pos.ra + 1*u.arcsec, dec=central_target_pos.dec,
                                       distance=target_dist_mpc*u.Mpc)
        value_pc_per_arc_sec = central_target_pos.separation_3d(off_set_central_pos).to(u.pc).value
        return rad_cloud_pc / value_pc_per_arc_sec

    def get_cloud_surf_dens(self, res=150):
        """
        Get surface density from giant molecular clouds
        Parameters
        ----------
        res : int or str

        """
        # load cloud catalog
        self.load_alma_cloud_cat(res=res)
        return self.alma_data['alpha_cloud_cat_data']['SURFDENS']

    def get_alma_obs_coverage_hull_dict(self):
        """
        Function to load the coverage dict of ALMA observations

        Returns
        -------
        coverage_dict : dict
        """
        # return np.load(self.path2obs_cover_hull / ('%s_alma_obs_hull_dict.npy' % self.gas_target_name),
        #                allow_pickle=True).item()
        with open(self.path2obs_cover_hull / ('%s_alma_obs_hull_dict.pickle' % self.gas_target_name), 'rb') as file_name:
            return pickle.load(file_name)

    def check_coords_covered_by_alma(self, ra, dec, res='native', max_dist_dist2hull_arcsec=2):
        """
        Function to check if coordinate points are inside MUSE observation

        Parameters
        ----------
        ra : float or ``np.ndarray``
        dec : float or ``np.ndarray``
        res : str
        max_dist_dist2hull_arcsec : float

        Returns
        -------
        coverage_dict : ``np.ndarray``
        """

        if isinstance(ra, float):
            ra = [ra]
            dec = [dec]

        hull_dict = self.get_alma_obs_coverage_hull_dict()[res]
        coverage_mask = np.zeros(len(ra), dtype=bool)
        hull_data_ra = np.array([])
        hull_data_dec = np.array([])

        for hull_idx in hull_dict.keys():
            ra_hull = hull_dict[hull_idx]['ra']
            dec_hull = hull_dict[hull_idx]['dec']
            hull_data_ra = np.concatenate([hull_data_ra, ra_hull])
            hull_data_dec = np.concatenate([hull_data_dec, dec_hull])
            coverage_mask += helper_func.GeometryTools.check_points_in_polygon(x_point=ra, y_point=dec,
                                                                               x_data_hull=ra_hull,
                                                                               y_data_hull=dec_hull)

        coverage_mask *= helper_func.GeometryTools.flag_close_points2ensemble(
            x_data=ra, y_data=dec, x_data_ensemble=hull_data_ra,y_data_ensemble=hull_data_dec,
            max_dist2ensemble=max_dist_dist2hull_arcsec/3600)

        return coverage_mask
