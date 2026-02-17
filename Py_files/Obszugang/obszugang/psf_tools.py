"""
Tools to help organize observations
"""
import numpy as np
from pathlib import Path
import pickle
from scipy.interpolate import interp1d
from obszugang import obs_info
from werkzeugkiste import phys_params
from werkzeugkiste import helper_func
from obszugang import obs_tools


class PSFTools:
    """
    class to access PSF information of different telescopes
    Note that this is work in progress and should so far not be used as a reference.
    Especially for JWST and HST the output PSF will be either a mean PSF or a psf at a specific detector position.
    """

    @staticmethod
    def get_closest_available_hst_psf_filter(band, instrument):
        if instrument == 'acs':
            available_band_list = obs_info.acs_wfc_psf_band_list
            wave_list = []
            for available_band in available_band_list:
                wave_list.append(obs_tools.ObsTools.get_hst_band_wave(available_band, instrument='acs',
                                                                        wave_estimator='mean_wave', unit='mu'))
            wave = obs_tools.ObsTools.get_hst_band_wave(band, instrument='acs', wave_estimator='mean_wave', unit='mu')
        elif instrument == 'uvis':
            available_band_list = obs_info.wfc3_uv_psf_band_list
            wave_list = []
            for available_band in available_band_list:
                wave_list.append(obs_tools.ObsTools.get_hst_band_wave(available_band, instrument='uvis',
                                                                        wave_estimator='mean_wave', unit='mu'))
            wave = obs_tools.ObsTools.get_hst_band_wave(band, instrument='uvis', wave_estimator='mean_wave',
                                                          unit='mu')
        elif instrument == 'ir':
            available_band_list = obs_info.wfc3_ir_psf_band_list
            wave_list = []
            for available_band in available_band_list:
                wave_list.append(obs_tools.ObsTools.get_hst_band_wave(available_band, instrument='ir',
                                                                        wave_estimator='mean_wave', unit='mu'))
            wave = obs_tools.ObsTools.get_hst_band_wave(band, instrument='ir', wave_estimator='mean_wave', unit='mu')
        else:
            raise KeyError('instrument not understood')
        if band in available_band_list:
            return band
        else:
            # get the closest band:
            min_diff = np.min(np.abs(np.array(wave_list) - wave))
            idx_closest_wave = np.where((np.abs(np.array(wave_list) - wave)) == min_diff)[0][0]
            return available_band_list[idx_closest_wave]

    @staticmethod
    def load_hst_psf_dict(band, instrument):
        """
        Parameters
        ----------
        band : str
            HST band
        instrument : str
            must be acs, uvis or ir

        Returns
        -------
        psf_dict : dict

        """
        # get psf dict path
        path2psf_dict = Path(__file__).parent.parent.resolve() / 'meta_data' / 'psf_data' / 'data_output'
        if instrument == 'acs':
            instrument_str = 'acs_wfc'
        elif instrument == 'uvis':
            instrument_str = 'wfc3_uv'
        elif instrument == 'ir':
            instrument_str = 'wfc3_ir'
        else:
            raise KeyError('instrument not understood')

        psf_dict_filename = 'hst_%s_psf_dict.pickle' % instrument_str
        # now there is not an estimated PSF for every filter at the moment hence we use the closest filter available
        with open(path2psf_dict / psf_dict_filename, 'rb') as file_name:
            psf_dict = pickle.load(file_name)
        # print(psf_dict)
        return psf_dict[PSFTools.get_closest_available_hst_psf_filter(band=band, instrument=instrument)]

    @staticmethod
    def get_hst_psf_rad_profile(band, instrument):
        """
        Parameters
        ----------
        band : str
            HST band
        instrument : str
            must be acs, uvis or ir

        Returns
        -------
        rad, profile : array-like

        """
        psf_dict = PSFTools.load_hst_psf_dict(band=band, instrument=instrument)
        return psf_dict['radius_arcsec'], psf_dict['psf_profile']

    @staticmethod
    def get_hst_psf_fwhm(band, instrument):
        """
        Parameters
        ----------
        band : str
            HST band
        instrument : str
            must be acs, uvis or ir

        Returns
        -------
        rad, profile : array-like

        """
        psf_dict = PSFTools.load_hst_psf_dict(band=band, instrument=instrument)
        return psf_dict['gaussian_fwhm']

    @staticmethod
    def get_hst_psf_gauss_approx(band, instrument, rad_arcsec, amp=1):
        """
        Parameters
        ----------
        band : str
            HST band
        instrument : str
            must be acs, uvis or ir

        Returns
        -------
        rad, profile : array-like

        """
        psf_dict = PSFTools.load_hst_psf_dict(band=band, instrument=instrument)
        mu = psf_dict['gaussian_mean']
        sig = psf_dict['gaussian_std']
        return amp * np.exp(-(rad_arcsec - mu) ** 2 / (2 * sig ** 2))

    @staticmethod
    # get psf dict path
    def load_jwst_psf_dict(band, instrument):
        """
        Parameters
        ----------
        band : str
            HST band
        instrument : str
            must be nircam or miri

        Returns
        -------
        psf_dict : dict

        """
        # get psf dict path
        assert(instrument in ['nircam', 'miri'])

        path2psf_dict = Path(__file__).parent.parent.resolve() / 'meta_data' / 'psf_data' / 'data_output'

        psf_dict_filename = '%s_psf_dict.pickle' % instrument
        # now there is not an estimated PSF for every filter at the moment hence we use the closest filter available
        with open(path2psf_dict / psf_dict_filename, 'rb') as file_name:
            psf_dict = pickle.load(file_name)
        return psf_dict[band]

    @staticmethod
    def get_jwst_psf_rad_profile(band, instrument):
        """
        Parameters
        ----------
        band : str
            NIRCAM or MIRI band
        instrument : str
            must nircam or miri

        Returns
        -------
        rad, profile : array-like

        """
        psf_dict = PSFTools.load_jwst_psf_dict(band=band, instrument=instrument)
        return psf_dict['radius_arcsec'], psf_dict['psf_profile']

    @staticmethod
    def get_jwst_psf_fwhm(band, instrument):
        """
        Parameters
        ----------
        band : str
            NIRCAM or MIRI band
        instrument : str
            must nircam or miri

        Returns
        -------
        rad, profile : array-like

        """
        psf_dict = PSFTools.load_jwst_psf_dict(band=band, instrument=instrument)
        return psf_dict['gaussian_fwhm']

    @staticmethod
    def get_jwst_psf_gauss_approx(band, instrument, rad_arcsec, amp=1):
        """
        Parameters
        ----------
        band : str
            NIRCAM or MIRI band
        instrument : str
            must nircam or miri

        Returns
        -------
        rad, profile : array-like

        """
        psf_dict = PSFTools.load_jwst_psf_dict(band=band, instrument=instrument)
        mu = psf_dict['gaussian_mean']
        sig = psf_dict['gaussian_std']
        return amp * np.exp(-(rad_arcsec - mu) ** 2 / (2 * sig ** 2))

    @staticmethod
    def get_astrosat_psf_fwhm(band):
        """
        Parameters
        ----------
        band : str
            ASTROSAT band


        Returns
        -------
        rad, profile : array-like

        """
        if band in ['F148W', 'F148W_old', 'F148Wa', 'F154W', 'F154W_old', 'F169M', 'F169M_old', 'F172M', 'F172M_old']:
            return phys_params.astrosat_psf_fwhm_arcsec['fuv']
        elif band in ['N219M', 'N219M_old', 'N242W', 'N242W_old', 'N245M', 'N245M_old', 'N263M', 'N263M_old', 'N279N',
                      'N279N_old']:
            return phys_params.astrosat_psf_fwhm_arcsec['nuv']
        else:
            raise KeyError('this band is not covered for PSF estimation or not an ASTROSAT band')

    @staticmethod
    def load_custom_gauss_psf_dict(band):
        """
        Parameters
        ----------
        band : str
            the custom gaussian band which is in this case a pseudo band


        Returns
        -------
        psf_dict : dict

        """
        # get psf dict path
        path2psf_dict = Path(__file__).parent.parent.resolve() / 'meta_data' / 'psf_data' / 'data_output'

        psf_dict_filename = 'custom_gauss_psf_dict.pickle'
        # now there is not an estimated PSF for every filter at the moment hence we use the closest filter available
        with open(path2psf_dict / psf_dict_filename, 'rb') as file_name:
            psf_dict = pickle.load(file_name)

        return psf_dict[band]

    @staticmethod
    def get_custom_gauss_psf_fwhm(band):
        """
        Parameters
        ----------
        band : str

        Returns
        -------

        """
        psf_dict = PSFTools.load_custom_gauss_psf_dict(band=band)
        return psf_dict['gaussian_fwhm']


    @staticmethod
    def load_obs_psf_dict(band, instrument):
        """
        Parameters
        ----------
        band : str
        instrument : str
            must be acs, uvis, ir, nircam or miri

        Returns
        -------
        rad, profile : array-like

        """
        assert(instrument in ['acs', 'uvis', 'ir', 'nircam', 'miri', 'custom_gauss'])
        if instrument in ['acs', 'uvis', 'ir']:
            return PSFTools.load_hst_psf_dict(band=band, instrument=instrument)
        elif instrument in ['nircam', 'miri']:
            return PSFTools.load_jwst_psf_dict(band=band, instrument=instrument)
        elif instrument == 'custom_gauss':
            return PSFTools.load_custom_gauss_psf_dict(band=band)

    @staticmethod
    def get_obs_psf_rad_profile(band, instrument):
        psf_dict = PSFTools.load_obs_psf_dict(band=band, instrument=instrument)
        return psf_dict['radius_arcsec'], psf_dict['psf_profile']

    @staticmethod
    def get_obs_psf_fwhm(band, instrument):
        assert(instrument in ['acs', 'uvis', 'ir', 'nircam', 'miri', 'astrosat', 'custom_gauss'])
        if instrument in ['acs', 'uvis', 'ir']:
            return PSFTools.get_hst_psf_fwhm(band=band, instrument=instrument)
        elif instrument in ['nircam', 'miri']:
            return PSFTools.get_jwst_psf_fwhm(band=band, instrument=instrument)
        elif instrument == 'custom_gauss':
            return PSFTools.get_custom_gauss_psf_fwhm(band=band)
        elif instrument == 'astrosat':
            return PSFTools.get_astrosat_psf_fwhm(band=band)

        # psf_dict = PSFTools.load_obs_psf_dict(band=band, instrument=instrument)
        # return psf_dict['gaussian_fwhm']

    @staticmethod
    def get_obs_psf_std(band, instrument):
        psf_dict = PSFTools.load_obs_psf_dict(band=band, instrument=instrument)
        return psf_dict['gaussian_std']

    @staticmethod
    def get_obs_psf_gauss_approx(band, instrument, rad_arcsec, amp=1):
        psf_dict = PSFTools.load_obs_psf_dict(band=band, instrument=instrument)
        mu = psf_dict['gaussian_mean']
        sig = psf_dict['gaussian_std']
        return amp * np.exp(-(rad_arcsec - mu) ** 2 / (2 * sig ** 2))

    @staticmethod
    def get_psf_gauss_corr_fact(band, instrument, std):
        # load the EE dict
        path2psf_data = Path(__file__).parent.parent.resolve() / 'meta_data' / 'psf_data' / 'data_output'
        if instrument == 'acs':
            ee_corr_file_name = 'acs_wfc_psf_correction_dict.pickle'
            used_band = PSFTools.get_closest_available_hst_psf_filter(band=band, instrument=instrument)
        elif instrument == 'uvis':
            ee_corr_file_name = 'wfc3_uv_psf_correction_dict.pickle'
            used_band = PSFTools.get_closest_available_hst_psf_filter(band=band, instrument=instrument)
        elif instrument == 'ir':
            ee_corr_file_name = 'wfc3_ir_psf_correction_dict.pickle'
            used_band = PSFTools.get_closest_available_hst_psf_filter(band=band, instrument=instrument)
        elif instrument == 'nircam':
            ee_corr_file_name = 'nircam_psf_correction_dict.pickle'
            used_band = band
        elif instrument == 'miri':
            ee_corr_file_name = 'miri_psf_correction_dict.pickle'
            used_band = band
        else:
            raise KeyError(instrument, ' not understood as instrument keyword')

        # psf_correction_dict = pickle.load(ee_corr_file_name)
        with open(path2psf_data / ee_corr_file_name, 'rb') as file_name:
            psf_correction_dict = pickle.load(file_name)

        if std < psf_correction_dict[used_band]['measured_std_values'][0]:
            return psf_correction_dict[used_band]['measured_ee_value'][0]
        elif std > psf_correction_dict[used_band]['measured_std_values'][-1]:
            return psf_correction_dict[used_band]['measured_ee_value'][-1]

        # interpolate the
        interp_func = interp1d(psf_correction_dict[used_band]['measured_std_values'], psf_correction_dict[used_band]['measured_ee_value'])
        return interp_func(std)

    @staticmethod
    def get_apert_gauss_corr_fact(band, instrument, apert_rad, std):
        # load the EE dict
        path2psf_data = (Path(__file__).parent.parent.resolve() / 'meta_data' / 'psf_data' / 'data_output' /
                         'apert_gauss_corr')

        if instrument == 'acs':
            apert_gauss_corr_file_name = 'acs_wfc_apert_gauss_corr_dict.pickle'
            used_band = PSFTools.get_closest_available_hst_psf_filter(band=band, instrument=instrument)
        elif instrument == 'uvis':
            apert_gauss_corr_file_name = 'wfc3_uv_apert_gauss_corr_dict.pickle'
            used_band = PSFTools.get_closest_available_hst_psf_filter(band=band, instrument=instrument)
        elif instrument == 'ir':
            apert_gauss_corr_file_name = 'wfc3_ir_apert_gauss_corr_dict.pickle'
            used_band = PSFTools.get_closest_available_hst_psf_filter(band=band, instrument=instrument)
        elif instrument == 'nircam':
            apert_gauss_corr_file_name = 'nircam_apert_gauss_corr_dict.pickle'
            used_band = band
        elif instrument == 'miri':
            apert_gauss_corr_file_name = 'miri_apert_gauss_corr_dict.pickle'
            used_band = band
        elif instrument == 'custom_gauss':
            apert_gauss_corr_file_name = 'custom_gauss_apert_gauss_corr_dict.pickle'
            used_band = band
        else:
            raise KeyError(instrument, ' not understood as instrument keyword')

        # apert_gauss_corr_dict = pickle.load(apert_gauss_corr_file_name)
        with open(path2psf_data / apert_gauss_corr_file_name, 'rb') as file_name:
            apert_gauss_corr_dict = pickle.load(file_name)

        # get the interpolation_function
        interp_func = helper_func.InterpTools.interp2dgrid(
            x_bins=apert_gauss_corr_dict[used_band]['apert_radii_arcsec_list'],
            y_bins=apert_gauss_corr_dict[used_band]['measured_std_values'],
            func_values=apert_gauss_corr_dict[used_band]['measured_flux_frac_in_apert'],
            method='cubic')

        print(apert_gauss_corr_dict[used_band]['apert_radii_arcsec_list'], apert_gauss_corr_dict[used_band]['measured_std_values'])


        if std < np.min(apert_gauss_corr_dict[used_band]['measured_std_values']):
            std = np.min(apert_gauss_corr_dict[used_band]['measured_std_values'])
        if std > np.max(apert_gauss_corr_dict[used_band]['measured_std_values']):
            std = np.max(apert_gauss_corr_dict[used_band]['measured_std_values']) * 0.999
            too_extended_flag = True
        else:
            too_extended_flag = False

        corr_fact = 1 / helper_func.InterpTools.get2dinterp_value(interp_func=interp_func, x_val=apert_rad, y_val=std)

        return corr_fact, too_extended_flag


