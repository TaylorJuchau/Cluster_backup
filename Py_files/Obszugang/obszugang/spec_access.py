"""
Data access structure for all kind of spectroscopic data products related to PHANGS
"""

from pathlib import Path
import numpy as np
import pickle
from astropy.wcs import WCS
import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.io import fits
try:
    from TardisPipeline.readData.MUSE_WFM import get_MUSE_polyFWHM
except ImportError:
    print('Spectroscopic analysis of MUSe data relying on the TardisPipeline is not available. '
          'If you want to use this, clone the git repo at https://gitlab.com/francbelf/ifu-pipeline '
          'and add package to your root directory. ')
from obszugang import access_config, obs_info
from werkzeugkiste import helper_func, spec_tools


class SpecAccess:
    """
    Access class to organize all kind of spectroscopic data available
    """

    def __init__(self, spec_target_name=None, nirspec_data_ver=None, miri_mrs_data_ver=None):
        """

        Parameters
        ----------
        spec_target_name : str
            Default None. Target name
        """
        if spec_target_name is not None:
            spec_target_name = helper_func.FileTools.target_name_no_directions(target=spec_target_name)
        self.spec_target_name = spec_target_name

        self.nirspec_data_ver = nirspec_data_ver
        self.miri_mrs_data_ver = miri_mrs_data_ver

        # loaded data dictionaries
        self.muse_dap_map_data = {}
        self.muse_datacube_data = {}
        self.kcwi_datacube_data = {}
        self.nirspec_datacube_data = {}
        self.miri_mrs_datacube_data = {}

        # get path to observation coverage hulls
        self.path2obs_cover_hull = (Path(__file__).parent.parent.absolute() / 'meta_data' / 'obs_coverage' /
                                    'data_output')

        super().__init__()

    def get_muse_data_file_name(self, data_prod='MUSEDAP', res='copt', ssp_model='fiducial'):
        """

        Parameters
        ----------
        data_prod : str
        res : str
        ssp_model : str
            This one is only for the copt DAPMAP data products

        Return
        ---------
        fiel_path : str
        """

        # get folder path
        file_path = (Path(access_config.phangs_config_dict['muse_data_path']) /
                     access_config.phangs_config_dict['muse_data_ver'] / res / data_prod)
        if (res == 'copt') & (data_prod == 'MUSEDAP'):
            file_path /= ssp_model

        # get file name
        if data_prod == 'MUSEDAP':
            if res == 'copt':
                file_name = '%s-%.2fasec_MAPS.fits' % (self.spec_target_name.upper(),
                                                       obs_info.muse_obs_res_dict[self.spec_target_name]['copt_res'])
            elif res == 'native':
                file_name = '%s_MAPS.fits' % self.spec_target_name.upper()
            elif res == '150pc':
                file_name = '%s-150pc_MAPS.fits' % self.spec_target_name.upper()
            elif res == '15asec':
                file_name = '%s-15asec_MAPS.fits' % self.spec_target_name.upper()
            else:
                raise KeyError(res, ' must be copt, native, 150pc or 15asec')
        elif data_prod == 'datacubes':
            if res == 'copt':
                file_name = '%s-%.2fasec.fits' % (self.spec_target_name.upper(),
                                                       obs_info.muse_obs_res_dict[self.spec_target_name]['copt_res'])
            elif res == 'native':
                file_name = '%s.fits' % self.spec_target_name.upper()
            elif res == '150pc':
                file_name = '%s-150pc.fits' % self.spec_target_name.upper()
            elif res == '15asec':
                file_name = '%s-15asec.fits' % self.spec_target_name.upper()
            else:
                raise KeyError(res, ' must be copt, native, 150pc or 15asec')
        else:
            raise KeyError(data_prod, ' must be either DAPMAP or datacubes')

        return file_path / file_name

    def get_kcwi_data_file_name(self):
        """

        Parameters
        ----------

        Return
        ---------
        fiel_path : str
        """

        # get folder path
        file_path = (Path(access_config.phangs_config_dict['kcwi_data_path']) / 'Datacubes')

        file_name = helper_func.FileTools.target_names_no_zeros(target=self.spec_target_name).upper() + '.fits'

        return file_path / file_name

    def get_nirspec_data_file_name(self, region_name, grating):
        """
        Function to get the file name and path of the JWST nirspec observations
        Parameters
        ----------

        Return
        ---------
        file_path : str
        """

        file_path = (Path(access_config.phangs_config_dict['nirspec_data_path']) / self.nirspec_data_ver)

        file_name = ('%s_nirspec_%s_s3d.fits' %
                     (obs_info.nirspec_obs_target_dict[self.spec_target_name][region_name][self.nirspec_data_ver]
                      ['file_id'],
                      helper_func.FileTools.get_nirspec_grating_file_name_comp(grating=grating)))

        return file_path / file_name

    def get_miri_mrs_data_file_name(self, region_name, channel):
        """
        Function to get the file name and path of the JWST miri_mrs observations
        Parameters
        ----------

        Return
        ---------
        file_path : str
        """

        file_path = (Path(access_config.phangs_config_dict['miri_mrs_data_path']) / self.miri_mrs_data_ver)
        file_name = ('%s_Level3_%s-shortmediumlong_s3d.fits' %
                     (obs_info.miri_mrs_obs_target_dict[self.spec_target_name][region_name][self.miri_mrs_data_ver]
                      ['file_id'], channel.lower()))
        return file_path / file_name

    def load_muse_dap_map(self, res='copt', ssp_model='fiducial', map_type='HA6562_FLUX'):
        """

        Parameters
        ----------
        res : str
        ssp_model : str
            This one is only for the copt DAP results
        map_type : str

        """
        file_path = self.get_muse_data_file_name(res=res, ssp_model=ssp_model)

        # get MUSE data
        muse_hdu = fits.open(file_path)
        muse_map_data = muse_hdu[map_type].data
        muse_map_wcs = WCS(muse_hdu[map_type].header)
        muse_hdu.close()

        data_identifier = helper_func.FileTools.get_dap_data_identifier(res=res, ssp_model=ssp_model)

        self.muse_dap_map_data.update({
            'dap_map_data_%s_%s' % (data_identifier, map_type): muse_map_data,
            'dap_map_wcs_%s_%s' % (data_identifier, map_type): muse_map_wcs
        })

    def check_dap_map_loaded(self, res='copt', ssp_model='fiducial', map_type='HA6562_FLUX'):
        data_identifier = helper_func.FileTools.get_dap_data_identifier(res=res, ssp_model=ssp_model)
        if not 'dap_map_data_%s_%s' % (data_identifier, map_type) in self.muse_dap_map_data.keys():
            self.load_muse_dap_map(res=res, ssp_model=ssp_model, map_type=map_type)

    def get_muse_dap_map_cutout(self, ra_cutout, dec_cutout, cutout_size, map_type_list=None, res='copt', ssp_model='fiducial'):
        if map_type_list is None:
            map_type_list = ['HA6562_FLUX']
        elif isinstance(map_type_list, str):
            map_type_list = [map_type_list]

        cutout_pos = SkyCoord(ra=ra_cutout, dec=dec_cutout, unit=(u.degree, u.degree), frame='icrs')
        cutout_dict = {'cutout_pos': cutout_pos}
        cutout_dict.update({'cutout_size': cutout_size})
        cutout_dict.update({'map_type_list': map_type_list})

        data_identifier = helper_func.FileTools.get_dap_data_identifier(res=res, ssp_model=ssp_model)

        for map_type in map_type_list:
            # make sure that map typ is loaded
            self.check_dap_map_loaded(res=res, ssp_model=ssp_model, map_type=map_type)
            cutout_dict.update({
                '%s_%s_img_cutout' % (data_identifier, map_type):
                    helper_func.CoordTools.get_img_cutout(
                        img=self.muse_dap_map_data['dap_map_data_%s_%s' % (data_identifier, map_type)],
                        wcs=self.muse_dap_map_data['dap_map_wcs_%s_%s' % (data_identifier, map_type)],
                        coord=cutout_pos, cutout_size=cutout_size)})
        return cutout_dict

    def load_muse_cube(self, res='copt'):
        """
        # loads MUSE data cube into the constructor
        Parameters
        ----------
        res : str

        """
        file_path = self.get_muse_data_file_name(data_prod='datacubes', res=res)
        # get MUSE data
        muse_hdu = fits.open(file_path)
        # get header
        hdr_data = muse_hdu['DATA'].header
        hdr_stat = muse_hdu['STAT'].header
        # get wavelength
        wave_muse = (hdr_data['CRVAL3'] + np.arange(hdr_data['NAXIS3']) * hdr_data['CD3_3']) * u.Unit(hdr_data['CUNIT3'])
        # distinguish between vacuum and air wavelength
        if hdr_data['CTYPE3'][:4] == 'AWAV':
            vacuum = False
        else:
            vacuum = True
        # get data and variance cube
        data_cube_muse = muse_hdu['DATA'].data
        var_cube_muse = muse_hdu['STAT'].data
        # get units of data cube and variables
        if hdr_data['BUNIT'] == '':
            data_cube_unit = (1e-20) * u.erg / u.s / u.cm / u.cm / u.Angstrom
            var_cube_unit = ((1e-20) * u.erg / u.s / u.cm / u.cm / u.Angstrom)**2
        else:
            data_cube_unit = u.Unit(hdr_data['BUNIT'])
            var_cube_unit = u.Unit(hdr_stat['BUNIT'])
        # get WCS
        wcs_3d_muse = WCS(hdr_data)
        wcs_2d_muse = wcs_3d_muse.celestial
        # close header
        muse_hdu.close()
        # add data to the class attributes
        self.muse_datacube_data.update({
            'wave_%s' % res: wave_muse,
            'vacuum_%s' % res: vacuum,
            'data_cube_%s' % res: data_cube_muse,
            'var_cube_%s' % res: var_cube_muse,
            'data_cube_unit_%s' % res: data_cube_unit,
            'var_cube_unit_%s' % res: var_cube_unit,
            'hdr_%s' % res: hdr_data,
            'wcs_3d_%s' % res: wcs_3d_muse,
            'wcs_2d_%s' % res: wcs_2d_muse
        })

    def load_kcwi_cube(self):
        """
        load kcwi cube into

        """
        file_path = self.get_kcwi_data_file_name()
        # get MUSE data
        kcwi_hdu = fits.open(file_path)

        hdr_data = kcwi_hdu[0].header
        # hdr_stat = kcwi_hdu[1].header

        # get wavelength
        wave_kcwi = (hdr_data['CRVAL3'] + np.arange(hdr_data['NAXIS3']) * hdr_data['CDELT3']) * u.Angstrom

        # distinguish between vacuum and air wavelength
        # if hdr_data['CTYPE3'][:4] == 'AWAV':
        #     vacuum = False
        # else:
        #     vacuum = True
        vacuum = False

        # get data and variance cube
        data_cube_kcwi = kcwi_hdu[0].data
        var_cube_kcwi = kcwi_hdu[1].data
        mask_cube_kcwi = kcwi_hdu[2].data

        # get units of data cube and variables
        # if hdr_data['BUNIT'] == '':
        #     data_cube_unit = (1e-20) * u.erg / u.s / u.cm / u.cm / u.Angstrom
        #     var_cube_unit = ((1e-20) * u.erg / u.s / u.cm / u.cm / u.Angstrom)**2
        # else:
        #     data_cube_unit = u.Unit(hdr_data['BUNIT'])
        #     var_cube_unit = u.Unit(hdr_stat['BUNIT'])
        data_cube_unit = (1e-16) * u.erg / u.s / u.cm / u.cm / u.Angstrom
        var_cube_unit = ((1e-16) * u.erg / u.s / u.cm / u.cm / u.Angstrom)**2

        # get WCS
        wcs_3d_kcwi = WCS(hdr_data)
        wcs_2d_kcwi = wcs_3d_kcwi.celestial
        # close header
        kcwi_hdu.close()
        # add data to the class attributes
        self.kcwi_datacube_data.update({
            'wave': wave_kcwi,
            'vacuum': vacuum,
            'data_cube': data_cube_kcwi,
            'var_cube': var_cube_kcwi,
            'mask_cube': mask_cube_kcwi,
            'data_cube_unit': data_cube_unit,
            'var_cube_unit': var_cube_unit,
            'hdr': hdr_data,
            'wcs_3d': wcs_3d_kcwi,
            'wcs_2d': wcs_2d_kcwi
        })
        #
        # exit()
        #
        # # get header
        # hdr = kcwi_hdu[0].header
        # # get WCS
        # wcs_3d_kcwi = WCS(hdr)
        # wcs_2d_kcwi = wcs_3d_kcwi.celestial
        #
        # # get wavelength
        # wave_kcwi = hdr['CRVAL3'] + np.arange(hdr['NAXIS3']) * hdr['CDELT3']
        #
        # # get data, variance and mask cube
        # data_cube_kcwi = kcwi_hdu[0].data
        # var_cube_kcwi = kcwi_hdu[1].data
        # mask_cube_kcwi = kcwi_hdu[2].data
        #
        # kcwi_hdu.close()
        #
        # self.kcwi_datacube_data.update({
        #     'wave': wave_kcwi,
        #     'data_cube': data_cube_kcwi,
        #     'var_cube': var_cube_kcwi,
        #     'mask_cube': mask_cube_kcwi,
        #     'hdr': hdr,
        #     'wcs_3d': wcs_3d_kcwi,
        #     'wcs_2d': wcs_2d_kcwi
        #
        # })

    def load_nirspec_cube(self, region_name, grating):
        """
        loads nirspec data cube into the constructor
        Parameters
        ----------
        region_name : str
        grating : str

        """
        file_path = self.get_nirspec_data_file_name(region_name=region_name, grating=grating)
        # get MUSE data
        nirspec_hdu = fits.open(file_path)
        # get header
        hdr_data = nirspec_hdu['SCI'].header
        hdr_stat = nirspec_hdu['ERR'].header
        # get wavelength
        wave_nirspec = ((hdr_data['CRVAL3'] + np.arange(hdr_data['NAXIS3']) * hdr_data['CDELT3']) *
                        u.Unit(hdr_data['CUNIT3']))
        # get data and uncertainty cube
        data_cube_nirspec = nirspec_hdu['SCI'].data
        err_cube_nirspec = nirspec_hdu['ERR'].data
        # get units of data cube and variables
        data_cube_unit = u.Unit(hdr_data['BUNIT'])
        err_cube_unit = u.Unit(hdr_stat['BUNIT'])
        # get WCS
        wcs_3d_nirspec = WCS(hdr_data)
        wcs_2d_nirspec = wcs_3d_nirspec.celestial
        # close header
        nirspec_hdu.close()
        # add data to the class attributes
        self.nirspec_datacube_data.update({
            region_name: {
                'wave_%s' % grating: wave_nirspec,
                'vacuum_%s' % grating: True,
                'data_cube_%s' % grating: data_cube_nirspec,
                'err_cube_%s' % grating: err_cube_nirspec,
                'data_cube_unit_%s' % grating: data_cube_unit,
                'err_cube_unit_%s' % grating: err_cube_unit,
                'hdr_%s' % grating: hdr_data,
                'wcs_3d_%s' % grating: wcs_3d_nirspec,
                'wcs_2d_%s' % grating: wcs_2d_nirspec
            }
        })

    def load_miri_mrs_cube(self, region_name, channel):
        """
        loads miri_mrs data cube into the constructor
        Parameters
        ----------
        region_name : str
        channel : str

        """
        file_path = self.get_miri_mrs_data_file_name(region_name=region_name, channel=channel)
        # get MUSE data
        miri_mrs_hdu = fits.open(file_path)
        # get header
        hdr_data = miri_mrs_hdu['SCI'].header
        hdr_stat = miri_mrs_hdu['ERR'].header
        # get wavelength
        wave_miri_mrs = ((hdr_data['CRVAL3'] + np.arange(hdr_data['NAXIS3']) * hdr_data['CDELT3']) *
                        u.Unit(hdr_data['CUNIT3']))
        # get data and uncertainty cube
        data_cube_miri_mrs = miri_mrs_hdu['SCI'].data
        err_cube_miri_mrs = miri_mrs_hdu['ERR'].data
        # get units of data cube and variables
        data_cube_unit = u.Unit(hdr_data['BUNIT'])
        err_cube_unit = u.Unit(hdr_stat['BUNIT'])
        # get WCS
        wcs_3d_miri_mrs = WCS(hdr_data)
        wcs_2d_miri_mrs = wcs_3d_miri_mrs.celestial
        # close header
        miri_mrs_hdu.close()
        # add data to the class attributes
        self.miri_mrs_datacube_data.update({
            region_name: {
                'wave_%s' % channel: wave_miri_mrs,
                'vacuum_%s' % channel: True,
                'data_cube_%s' % channel: data_cube_miri_mrs,
                'err_cube_%s' % channel: err_cube_miri_mrs,
                'data_cube_unit_%s' % channel: data_cube_unit,
                'err_cube_unit_%s' % channel: err_cube_unit,
                'hdr_%s' % channel: hdr_data,
                'wcs_3d_%s' % channel: wcs_3d_miri_mrs,
                'wcs_2d_%s' % channel: wcs_2d_miri_mrs
            }
        })

    def get_muse_obs_coverage_hull_dict(self):
        """
        Function to load the coverage dict of MUSE observations

        Returns
        -------
        coverage_mask : dict
        """
        with open(self.path2obs_cover_hull / ('%s_muse_obs_hull_dict.pickle' % self.spec_target_name), 'rb') as file_name:
            return pickle.load(file_name)

    def get_kcwi_obs_coverage_hull_dict(self):
        """
        Function to load the coverage dict of KCWI observations

        Returns
        -------
        coverage_mask : dict
        """
        with (open(self.path2obs_cover_hull / ('%s_kcwi_obs_hull_dict.pickle' % self.spec_target_name), 'rb') as
              file_name):
            return pickle.load(file_name)

    def get_nirspec_obs_coverage_hull_dict(self):
        """
        Function to load the coverage dict of NIRSPEC observations

        Returns
        -------
        coverage_mask : dict
        """
        with open(self.path2obs_cover_hull / ('%s_nirspec_obs_hull_dict_%s.pickle' %
                                              (self.spec_target_name, self.nirspec_data_ver)),
                  'rb') as file_name:
            return pickle.load(file_name)

    def get_miri_mrs_obs_coverage_hull_dict(self):
        """
        Function to load the coverage dict of MIRI MRS observations

        Returns
        -------
        coverage_mask : dict
        """
        with open(self.path2obs_cover_hull / ('%s_miri_mrs_obs_hull_dict_%s.pickle' %
                                              (self.spec_target_name, self.miri_mrs_data_ver)),
                  'rb') as file_name:
            return pickle.load(file_name)

    def check_coords_covered_by_muse(self, ra, dec, res='copt', max_dist_dist2hull_arcsec=0.5):
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

        band_hull_dict = self.get_muse_obs_coverage_hull_dict()[res]
        coverage_mask = np.zeros(len(ra), dtype=bool)
        hull_data_ra = np.array([])
        hull_data_dec = np.array([])

        for hull_idx in band_hull_dict.keys():
            ra_hull = band_hull_dict[hull_idx]['ra']
            dec_hull = band_hull_dict[hull_idx]['dec']
            hull_data_ra = np.concatenate([hull_data_ra, ra_hull])
            hull_data_dec = np.concatenate([hull_data_dec, dec_hull])
            coverage_mask += helper_func.GeometryTools.check_points_in_polygon(x_point=ra, y_point=dec,
                                                                               x_data_hull=ra_hull,
                                                                               y_data_hull=dec_hull)

        coverage_mask *= helper_func.GeometryTools.flag_close_points2ensemble(
            x_data=ra, y_data=dec, x_data_ensemble=hull_data_ra,y_data_ensemble=hull_data_dec,
            max_dist2ensemble=max_dist_dist2hull_arcsec/3600)

        return coverage_mask

    def check_coords_covered_by_kcwi(self, ra, dec, max_dist_dist2hull_arcsec=0.5):
        """
        Function to check if coordinate points are inside MUSE observation

        Parameters
        ----------
        ra : float or ``np.ndarray``
        dec : float or ``np.ndarray``
        max_dist_dist2hull_arcsec : float

        Returns
        -------
        coverage_dict : ``np.ndarray``
        """

        if isinstance(ra, float):
            ra = [ra]
            dec = [dec]

        band_hull_dict = self.get_kcwi_obs_coverage_hull_dict()

        coverage_mask = np.zeros(len(ra), dtype=bool)
        hull_data_ra = np.array([])
        hull_data_dec = np.array([])

        for hull_idx in band_hull_dict.keys():
            ra_hull = band_hull_dict[hull_idx]['ra']
            dec_hull = band_hull_dict[hull_idx]['dec']
            hull_data_ra = np.concatenate([hull_data_ra, ra_hull])
            hull_data_dec = np.concatenate([hull_data_dec, dec_hull])
            coverage_mask += helper_func.GeometryTools.check_points_in_polygon(x_point=ra, y_point=dec,
                                                                               x_data_hull=ra_hull,
                                                                               y_data_hull=dec_hull)

        coverage_mask *= helper_func.GeometryTools.flag_close_points2ensemble(
            x_data=ra, y_data=dec, x_data_ensemble=hull_data_ra,y_data_ensemble=hull_data_dec,
            max_dist2ensemble=max_dist_dist2hull_arcsec/3600)

        return coverage_mask

    def check_coords_covered_by_nirspec(self, ra, dec, region_name, grating, max_dist_dist2hull_arcsec=0.5):
        """
        Function to check if coordinate points are inside NIRSpec observation

        Parameters
        ----------
        ra : float or ``np.ndarray``
        dec : float or ``np.ndarray``
        region_name: str
        grating: str
        max_dist_dist2hull_arcsec : float

        Returns
        -------
        coverage_dict : ``np.ndarray``
        """

        if isinstance(ra, float):
            ra = [ra]
            dec = [dec]

        band_hull_dict = self.get_nirspec_obs_coverage_hull_dict()

        coverage_mask = np.zeros(len(ra), dtype=bool)
        hull_data_ra = np.array([])
        hull_data_dec = np.array([])

        for hull_idx in band_hull_dict.keys():
            ra_hull = band_hull_dict[region_name][grating][hull_idx]['ra']
            dec_hull = band_hull_dict[region_name][grating][hull_idx]['dec']
            hull_data_ra = np.concatenate([hull_data_ra, ra_hull])
            hull_data_dec = np.concatenate([hull_data_dec, dec_hull])
            coverage_mask += helper_func.GeometryTools.check_points_in_polygon(x_point=ra, y_point=dec,
                                                                               x_data_hull=ra_hull,
                                                                               y_data_hull=dec_hull)
        coverage_mask *= helper_func.GeometryTools.flag_close_points2ensemble(
            x_data=ra, y_data=dec, x_data_ensemble=hull_data_ra, y_data_ensemble=hull_data_dec,
            max_dist2ensemble=max_dist_dist2hull_arcsec/3600)

        return coverage_mask

    def check_coords_covered_by_miri_mrs(self, ra, dec, region_name, channel, max_dist_dist2hull_arcsec=0.5):
        """
        Function to check if coordinate points are inside NIRSpec observation

        Parameters
        ----------
        ra : float or ``np.ndarray``
        dec : float or ``np.ndarray``
        region_name: str
        channel: str
        max_dist_dist2hull_arcsec : float

        Returns
        -------
        coverage_dict : ``np.ndarray``
        """

        if isinstance(ra, float):
            ra = [ra]
            dec = [dec]

        band_hull_dict = self.get_miri_mrs_obs_coverage_hull_dict()

        coverage_mask = np.zeros(len(ra), dtype=bool)
        hull_data_ra = np.array([])
        hull_data_dec = np.array([])

        for hull_idx in band_hull_dict.keys():
            ra_hull = band_hull_dict[region_name][channel][hull_idx]['ra']
            dec_hull = band_hull_dict[region_name][channel][hull_idx]['dec']
            hull_data_ra = np.concatenate([hull_data_ra, ra_hull])
            hull_data_dec = np.concatenate([hull_data_dec, dec_hull])
            coverage_mask += helper_func.GeometryTools.check_points_in_polygon(x_point=ra, y_point=dec,
                                                                               x_data_hull=ra_hull,
                                                                               y_data_hull=dec_hull)
        coverage_mask *= helper_func.GeometryTools.flag_close_points2ensemble(
            x_data=ra, y_data=dec, x_data_ensemble=hull_data_ra, y_data_ensemble=hull_data_dec,
            max_dist2ensemble=max_dist_dist2hull_arcsec/3600)

        return coverage_mask

    def extract_muse_spec_circ_app(self, ra, dec, rad_arcsec, wave_range=None, res='copt'):
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

        # make sure muse cube is loaded
        if 'data_cube_%s' % res not in self.muse_datacube_data.keys():
            self.load_muse_cube(res=res)

        # get select spectra from coordinates
        obj_coords_world = SkyCoord(ra=ra * u.deg, dec=dec * u.deg)
        obj_coords_muse_pix = self.muse_datacube_data['wcs_2d_%s' % res].world_to_pixel(obj_coords_world)
        selection_radius_pix = helper_func.CoordTools.transform_world2pix_scale(
            length_in_arcsec=rad_arcsec, wcs=self.muse_datacube_data['wcs_2d_%s' % res])
        # select spaxels which are inside the selected radius
        x_lin_muse = np.linspace(1, self.muse_datacube_data['data_cube_%s' % res].shape[2],
                                 self.muse_datacube_data['data_cube_%s' % res].shape[2])
        y_lin_muse = np.linspace(1, self.muse_datacube_data['data_cube_%s' % res].shape[1],
                                 self.muse_datacube_data['data_cube_%s' % res].shape[1])
        x_data_muse, y_data_muse = np.meshgrid(x_lin_muse, y_lin_muse)
        mask_spectrum = (np.sqrt((x_data_muse - obj_coords_muse_pix[0]) ** 2 +
                                 (y_data_muse - obj_coords_muse_pix[1]) ** 2) < selection_radius_pix)

        # extract fluxes
        native_spec_flx = (np.sum(self.muse_datacube_data['data_cube_%s' % res][:, mask_spectrum], axis=1) *
                           self.muse_datacube_data['data_cube_unit_%s' % res])

        native_spec_flx_err = np.sqrt(np.sum(self.muse_datacube_data['var_cube_%s' % res][:, mask_spectrum], axis=1) *
                                      self.muse_datacube_data['var_cube_unit_%s' % res])

        # get wavelengths
        native_wave = self.muse_datacube_data['wave_%s' % res]
        # getting line spread function
        lsf_fwhm = spec_tools.SpecHelper.get_muse_lsf_fwhm(wave = native_wave)

        if sum(np.invert(np.isnan(native_spec_flx))) == 0:
            return None

        # get wavelength range
        if wave_range is None:
            wave_range = [np.nanmin(self.muse_datacube_data['wave_%s' % res][np.invert(np.isnan(native_spec_flx))]),
                          np.nanmax(self.muse_datacube_data['wave_%s' % res][np.invert(np.isnan(native_spec_flx))])]
        else:
            wave_range = wave_range

        # make sure wave length range is applied to all arrays
        mask_wave_range = (native_wave >= wave_range[0]) & (native_wave <= wave_range[1])
        native_spec_flx = native_spec_flx[mask_wave_range]
        native_spec_flx_err = native_spec_flx_err[mask_wave_range]
        # hot fix for nan errors
        mask_problematic_err = np.isnan(native_spec_flx) + np.isinf(native_spec_flx) + np.isnan(native_spec_flx_err) + np.isinf(native_spec_flx_err)
        native_spec_flx_err[mask_problematic_err] = native_spec_flx[mask_problematic_err]
        native_wave = native_wave[mask_wave_range]
        lsf_fwhm = lsf_fwhm[mask_wave_range]
        native_good_pixel_mask = (np.invert(np.isnan(native_spec_flx) + np.isinf(native_spec_flx) +
                                            np.isnan(native_spec_flx_err) + np.isinf(native_spec_flx_err)))

        spec_dict = {
            # general description of spectrum
            'rad_arcsec': rad_arcsec,
            'lsf_fwhm': lsf_fwhm,
            # native spectrum
            'native_wave_range': wave_range,
            'nativ_wave_vaccum': self.muse_datacube_data['vacuum_%s' % res],
            'native_spec_flx': native_spec_flx,
            'native_spec_flx_err': native_spec_flx_err,
            'native_wave': native_wave,
            'native_good_pixel_mask': native_good_pixel_mask,
        }


        # get redshift if possible
        redshift = spec_tools.SpecHelper.get_target_ned_redshift(target=self.spec_target_name)
        sys_vel = spec_tools.SpecHelper.get_target_sys_vel(target=self.spec_target_name, redshift=redshift)
        spec_dict.update({
            'redshift': redshift,
            'sys_vel': sys_vel,
        })

        return spec_dict

    def extract_kcwi_spec_circ_app(self, ra, dec, rad_arcsec, wave_range=None):
        """

        Parameters
        ----------
        ra, dec : float
        rad_arcsec : float
            radius in arcseconds
        wave_range : tuple or None

        Return
        ---------
        spec_dict : dict
            Dictionary with all spectra related values
            fluxes are in units of ``10$^{-16}$ erg cm$^{-2}$ s$^{-1}$ ${\rm \AA^{-1}}$``
        """

        # make sure kcwi cube is loaded
        if 'data_cube' not in self.kcwi_datacube_data.keys():
            self.load_kcwi_cube()

        # get select spectra from coordinates
        obj_coords_world = SkyCoord(ra=ra * u.deg, dec=dec * u.deg)
        obj_coords_kcwi_pix = self.kcwi_datacube_data['wcs_2d'].world_to_pixel(obj_coords_world)
        selection_radius_pix = helper_func.CoordTools.transform_world2pix_scale(
            length_in_arcsec=rad_arcsec, wcs=self.kcwi_datacube_data['wcs_2d'], dim=0)

        # select spaxels which are inside the selected radius
        x_lin_kcwi = np.linspace(1, self.kcwi_datacube_data['data_cube'].shape[2],
                                 self.kcwi_datacube_data['data_cube'].shape[2])
        y_lin_kcwi = np.linspace(1, self.kcwi_datacube_data['data_cube'].shape[1],
                                 self.kcwi_datacube_data['data_cube'].shape[1])
        x_data_kcwi, y_data_kcwi = np.meshgrid(x_lin_kcwi, y_lin_kcwi)
        mask_spectrum = (np.sqrt((x_data_kcwi - obj_coords_kcwi_pix[0]) ** 2 +
                                 (y_data_kcwi - obj_coords_kcwi_pix[1]) ** 2) < selection_radius_pix)

        # extract fluxes
        native_spec_flx = (np.sum(self.kcwi_datacube_data['data_cube'][:, mask_spectrum], axis=1) *
                           self.kcwi_datacube_data['data_cube_unit'])

        native_spec_flx_err = np.sqrt(np.sum(self.kcwi_datacube_data['var_cube'][:, mask_spectrum], axis=1) *
                                      self.kcwi_datacube_data['var_cube_unit'])

        # get wavelengths
        native_wave = self.kcwi_datacube_data['wave']
        # getting line spread function
        lsf_fwhm = spec_tools.SpecHelper.get_kcwi_lsf_fwhm(wave=native_wave)

        if sum(np.invert(np.isnan(native_spec_flx))) == 0:
            return None
        # get wavelength range
        if wave_range is None:
            wave_range = [np.nanmin(self.kcwi_datacube_data['wave'][np.invert(np.isnan(native_spec_flx))]),
                          np.nanmax(self.kcwi_datacube_data['wave'][np.invert(np.isnan(native_spec_flx))])]
        else:
            wave_range = wave_range

        # make sure wave length range is applied to all arrays
        mask_wave_range = (native_wave >= wave_range[0]) & (native_wave <= wave_range[1])
        native_spec_flx = native_spec_flx[mask_wave_range]
        native_spec_flx_err = native_spec_flx_err[mask_wave_range]
        # hot fix for nan errors
        mask_problematic_err = np.isnan(native_spec_flx) + np.isinf(native_spec_flx) + np.isnan(native_spec_flx_err) + np.isinf(native_spec_flx_err)
        native_spec_flx_err[mask_problematic_err] = native_spec_flx[mask_problematic_err]
        native_wave = native_wave[mask_wave_range]
        lsf_fwhm = lsf_fwhm[mask_wave_range]
        native_good_pixel_mask = (np.invert(np.isnan(native_spec_flx) + np.isinf(native_spec_flx) +
                                            np.isnan(native_spec_flx_err) + np.isinf(native_spec_flx_err)))

        spec_dict = {
            # general description of spectrum
            'rad_arcsec': rad_arcsec,
            'lsf_fwhm': lsf_fwhm,
            # native spectrum
            'native_wave_range': wave_range,
            'nativ_wave_vaccum': self.kcwi_datacube_data['vacuum'],
            'native_spec_flx': native_spec_flx,
            'native_spec_flx_err': native_spec_flx_err,
            'native_wave': native_wave,
            'native_good_pixel_mask': native_good_pixel_mask,
        }


        # get redshift if possible
        redshift = spec_tools.SpecHelper.get_target_ned_redshift(target=self.spec_target_name)
        sys_vel = spec_tools.SpecHelper.get_target_sys_vel(target=self.spec_target_name, redshift=redshift)
        spec_dict.update({
            'redshift': redshift,
            'sys_vel': sys_vel,
        })

        return spec_dict
        #
        # # # extract fluxes
        # # spec_flux = np.sum(self.kcwi_datacube_data['data_cube'][:, mask_spectrum], axis=1)
        # # spec_flux_err = np.sqrt(np.sum(self.kcwi_datacube_data['var_cube'][:, mask_spectrum], axis=1))
        # # # rescale to the unit of 10$^{-16}$ erg cm$^{-2}$ s$^{-1}$ ${\rm \AA^{-1}}$
        # # spec_flux *= 1e-16
        # # spec_flux_err *= 1e-16
        # # # get wavelengths
        # # lam = self.kcwi_datacube_data['wave']
        # # # getting line spread function
        # # lsf = get_MUSE_polyFWHM(self.kcwi_datacube_data['wave'], version="udf10")
        #
        # # get wavelength range
        # if wave_range is None:
        #     lam_range = [np.nanmin(self.kcwi_datacube_data['wave'][np.invert(np.isnan(spec_flux))]),
        #                  np.nanmax(self.kcwi_datacube_data['wave'][np.invert(np.isnan(spec_flux))])]
        # else:
        #     lam_range = wave_range
        # # make sure wave length range is applied to all arrays
        # mask_wave_range = (lam > lam_range[0]) & (lam < lam_range[1])
        # spec_flux = spec_flux[mask_wave_range]
        # spec_flux_err = spec_flux_err[mask_wave_range]
        # lam = lam[mask_wave_range]
        # lsf = lsf[mask_wave_range]
        # good_pixel_mask = np.invert(np.isnan(spec_flux) + np.isinf(spec_flux))
        #
        # return {'lam_range': lam_range, 'spec_flux': spec_flux, 'spec_flux_err': spec_flux_err, 'lam': lam,
        #         'lsf': lsf, 'good_pixel_mask': good_pixel_mask, 'rad_arcsec': rad_arcsec}

    def extract_kcwi_sub_cube(self, ra, dec, cutout_size):
        """

        Parameters
        ----------
        ra, dec : float
        cutout_size : tuple

        Return
        ---------
        spec_dict : dict
            Dictionary with all spectra related values
            fluxes are in units of ``10$^{-16}$ erg cm$^{-2}$ s$^{-1}$ ${\rm \AA^{-1}}$``
        """

        # make sure kcwi cube is loaded
        if 'data_cube' not in self.kcwi_datacube_data.keys():
            self.load_kcwi_cube()

        # get select spectra from coordinates
        obj_coords_world = SkyCoord(ra=ra * u.deg, dec=dec * u.deg)
        obj_coords_kcwi_pix = self.kcwi_datacube_data['wcs_2d'].world_to_pixel(obj_coords_world)
        width_pix = helper_func.CoordTools.transform_world2pix_scale(
            length_in_arcsec=cutout_size[0], wcs=self.kcwi_datacube_data['wcs_2d'], dim=0)
        length_pix = helper_func.CoordTools.transform_world2pix_scale(
            length_in_arcsec=cutout_size[1], wcs=self.kcwi_datacube_data['wcs_2d'], dim=0)

        print(width_pix, length_pix)

        #get a sub cutout


        exit()


        selection_radius_pix = helper_func.CoordTools.transform_world2pix_scale(
            length_in_arcsec=rad_arcsec, wcs=self.kcwi_datacube_data['wcs_2d'], dim=0)

        # select spaxels which are inside the selected radius
        x_lin_kcwi = np.linspace(1, self.kcwi_datacube_data['data_cube'].shape[2],
                                 self.kcwi_datacube_data['data_cube'].shape[2])
        y_lin_kcwi = np.linspace(1, self.kcwi_datacube_data['data_cube'].shape[1],
                                 self.kcwi_datacube_data['data_cube'].shape[1])
        x_data_kcwi, y_data_kcwi = np.meshgrid(x_lin_kcwi, y_lin_kcwi)
        mask_spectrum = (np.sqrt((x_data_kcwi - obj_coords_kcwi_pix[0]) ** 2 +
                                 (y_data_kcwi - obj_coords_kcwi_pix[1]) ** 2) < selection_radius_pix)

        # extract fluxes
        spec_flux = np.sum(self.kcwi_datacube_data['data_cube'][:, mask_spectrum], axis=1)
        spec_flux_err = np.sqrt(np.sum(self.kcwi_datacube_data['var_cube'][:, mask_spectrum], axis=1))
        # rescale to the unit of 10$^{-16}$ erg cm$^{-2}$ s$^{-1}$ ${\rm \AA^{-1}}$
        spec_flux *= 1e-16
        spec_flux_err *= 1e-16
        # get wavelengths
        lam = self.kcwi_datacube_data['wave']
        # getting line spread function
        lsf = get_MUSE_polyFWHM(self.kcwi_datacube_data['wave'], version="udf10")

        # get wavelength range
        if wave_range is None:
            lam_range = [np.nanmin(self.kcwi_datacube_data['wave'][np.invert(np.isnan(spec_flux))]),
                         np.nanmax(self.kcwi_datacube_data['wave'][np.invert(np.isnan(spec_flux))])]
        else:
            lam_range = wave_range
        # make sure wave length range is applied to all arrays
        mask_wave_range = (lam > lam_range[0]) & (lam < lam_range[1])
        spec_flux = spec_flux[mask_wave_range]
        spec_flux_err = spec_flux_err[mask_wave_range]
        lam = lam[mask_wave_range]
        lsf = lsf[mask_wave_range]
        good_pixel_mask = np.invert(np.isnan(spec_flux) + np.isinf(spec_flux))

        return {'lam_range': lam_range, 'spec_flux': spec_flux, 'spec_flux_err': spec_flux_err, 'lam': lam,
                'lsf': lsf, 'good_pixel_mask': good_pixel_mask, 'rad_arcsec': rad_arcsec}

    def extract_muse_spec_gauss_app(self, ra, dec, fwhm_arcsec, wave_range=None, res='copt'):
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
            Dictionary with all spectra related values
            fluxes are in units of ``10$^{-16}$ erg cm$^{-2}$ s$^{-1}$ ${\rm \AA^{-1}}$``
        """

        # make sure muse cube is loaded
        if 'data_cube_%s' % res not in self.muse_datacube_data.keys():
            self.load_muse_cube(res=res)

        # get select spectra from coordinates
        obj_coords_world = SkyCoord(ra=ra * u.deg, dec=dec * u.deg)
        obj_coords_muse_pix = self.muse_datacube_data['wcs_2d_%s' % res].world_to_pixel(obj_coords_world)
        selection_fwhm_pix = helper_func.CoordTools.transform_world2pix_scale(
            length_in_arcsec=fwhm_arcsec, wcs=self.muse_datacube_data['wcs_2d_%s' % res])
        selection_sig_pix = helper_func.TransTools.gauss_fwhm2sig(fwhm=selection_fwhm_pix)

        # select spaxels which are inside the selected radius
        x_lin_muse = np.linspace(1, self.muse_datacube_data['data_cube_%s' % res].shape[2],
                                 self.muse_datacube_data['data_cube_%s' % res].shape[2])
        y_lin_muse = np.linspace(1, self.muse_datacube_data['data_cube_%s' % res].shape[1],
                                 self.muse_datacube_data['data_cube_%s' % res].shape[1])
        x_data_muse, y_data_muse = np.meshgrid(x_lin_muse, y_lin_muse)
        # create 2d Gauss
        gauss = helper_func.FuncAndModels.gauss2d_rot(x=x_data_muse, y=y_data_muse, amp=1, x0=obj_coords_muse_pix[0],
                                                      y0=obj_coords_muse_pix[1], sig_x=selection_sig_pix, sig_y=selection_sig_pix, theta=0)

        scaled_spec = self.muse_datacube_data['data_cube_%s' % res] * gauss
        scaled_var = self.muse_datacube_data['var_cube_%s' % res] * gauss

        # import matplotlib.pyplot as plt
        #
        # plt.imshow(np.nansum(scaled_spec, axis=0))
        # plt.show()
        # exit()

        # now to accelerate this we select only those spaxels which are somwhere around this
        mask_spectrum = (np.sqrt((x_data_muse - obj_coords_muse_pix[0]) ** 2 +
                                 (y_data_muse - obj_coords_muse_pix[1]) ** 2) < 10 * selection_fwhm_pix)
        # extract fluxes
        spec_flux = np.nansum(scaled_spec[:, mask_spectrum], axis=1)
        spec_flux_err = np.sqrt(np.nansum(scaled_var[:, mask_spectrum], axis=1))





        # spec_flux = np.nansum(self.muse_datacube_data['data_cube_%s' % res] * gauss, axis=(1, 2))
        # spec_flux_err = np.sqrt(np.nansum(self.muse_datacube_data['var_cube_%s' % res] * gauss, axis=(1, 2)))

        # rescale to the unit of 10$^{-16}$ erg cm$^{-2}$ s$^{-1}$ ${\rm \AA^{-1}}$
        spec_flux *= 1e-20
        spec_flux_err *= 1e-20
        # get wavelengths
        lam = self.muse_datacube_data['wave_%s' % res]
        # getting line spread function
        lsf = get_MUSE_polyFWHM(self.muse_datacube_data['wave_%s' % res], version="udf10")

        # get wavelength range
        if wave_range is None:
            lam_range = [np.nanmin(self.muse_datacube_data['wave_%s' % res][np.invert(np.isnan(spec_flux))]),
                         np.nanmax(self.muse_datacube_data['wave_%s' % res][np.invert(np.isnan(spec_flux))])]
        else:
            lam_range = wave_range
        # make sure wave length range is applied to all arrays
        mask_wave_range = (lam > lam_range[0]) & (lam < lam_range[1])
        spec_flux = spec_flux[mask_wave_range]
        spec_flux_err = spec_flux_err[mask_wave_range]
        lam = lam[mask_wave_range]
        lsf = lsf[mask_wave_range]
        good_pixel_mask = np.invert(np.isnan(spec_flux) + np.isinf(spec_flux))

        return {'lam_range': lam_range, 'spec_flux': np.array(spec_flux), 'spec_flux_err': np.array(spec_flux_err), 'lam': np.array(lam),
                'lsf': np.array(lsf), 'good_pixel_mask': np.array(good_pixel_mask), 'fwhm_arcsec': fwhm_arcsec}

