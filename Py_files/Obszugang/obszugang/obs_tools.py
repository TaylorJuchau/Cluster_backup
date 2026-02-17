"""
Tools to help organize observations
"""
import numpy as np
from obszugang import obs_info
from werkzeugkiste import phys_params
from werkzeugkiste.helper_func import UnitTools


class ObsTools:
    """
    Class to sort band names and identify instruments and telescopes
    """

    @staticmethod
    def get_hst_instrument(target, band):
        """
        get the corresponding instrument for hst observations

        Parameters
        ----------
        target :  str
        band :  str

        Returns
        -------
        target_name : str
        """
        if band in obs_info.hst_obs_band_dict[target]['acs']:
            return 'acs'
        elif band in obs_info.hst_obs_band_dict[target]['uvis']:
            return 'uvis'
        elif band in obs_info.hst_obs_band_dict[target]['ir']:
            return 'ir'
        else:
            print(target, ' has no HST observation for the Band ', band)
            return None

    @staticmethod
    def get_hst_ha_instrument(target):
        """
        get the corresponding instrument for hst H-alpha observations

        Parameters
        ----------
        target :  str

        Returns
        -------
        target_name : str
        """
        if (('F657N' in obs_info.hst_obs_band_dict[target]['uvis']) |
                ('F658N' in obs_info.hst_obs_band_dict[target]['uvis'])):
            return 'uvis'
        elif (('F657N' in obs_info.hst_obs_band_dict[target]['acs']) |
              ('F658N' in obs_info.hst_obs_band_dict[target]['acs'])):
            return 'acs'
        else:
            raise KeyError(target, ' has no H-alpha observation ')

    @staticmethod
    def get_hst_band_wave(band, instrument='acs', wave_estimator='mean_wave', unit='mu'):
        """
        Returns mean wavelength of an HST specific band
        Parameters
        ----------
        band : str
        instrument : str
        wave_estimator: str
            can be mean_wave, min_wave or max_wave
        unit : str

        Returns
        -------
        wavelength : float
        """
        if instrument == 'acs':
            return UnitTools.angstrom2unit(wave=phys_params.hst_acs_wfc1_bands_wave[band][wave_estimator], unit=unit)
        elif (instrument == 'uvis') | (instrument == 'uvis1'):
            return UnitTools.angstrom2unit(wave=phys_params.hst_wfc3_uvis1_bands_wave[band][wave_estimator], unit=unit)
        elif instrument == 'uvis2':
            return UnitTools.angstrom2unit(wave=phys_params.hst_wfc3_uvis2_bands_wave[band][wave_estimator], unit=unit)
        elif instrument == 'ir':
            return UnitTools.angstrom2unit(wave=phys_params.hst_wfc3_ir_bands_wave[band][wave_estimator], unit=unit)
        else:
            raise KeyError(instrument, ' is not a HST instrument')

    @staticmethod
    def get_roman_band_wave(band, wave_estimator='mean_wave', unit='mu'):
        """
        Returns mean wavelength of an HST specific band
        Parameters
        ----------
        band : str
        instrument : str
        wave_estimator: str
            can be mean_wave, min_wave or max_wave
        unit : str

        Returns
        -------
        wavelength : float
        """
        return UnitTools.angstrom2unit(wave=phys_params.roman_bands_wave[band][wave_estimator], unit=unit)

    @staticmethod
    def get_jwst_band_wave(band, instrument='nircam', wave_estimator='mean_wave', unit='mu'):
        """
        Returns mean wavelength of an JWST specific band
        Parameters
        ----------
        band : str
        instrument : str
        wave_estimator: str
            can be mean_wave, min_wave or max_wave
        unit : str

        Returns
        -------
        wavelength : float
        """
        if instrument == 'nircam':
            return UnitTools.angstrom2unit(wave=phys_params.nircam_bands_wave[band][wave_estimator], unit=unit)
        elif instrument == 'miri':
            return UnitTools.angstrom2unit(wave=phys_params.miri_bands_wave[band][wave_estimator], unit=unit)
        else:
            raise KeyError(instrument, ' is not a JWST instrument')

    @staticmethod
    def get_astrosat_band_wave(band, wave_estimator='mean_wave', unit='mu'):
        """
        Returns mean wavelength of an JWST specific band
        Parameters
        ----------
        band : str
        wave_estimator: str
            can be mean_wave, min_wave or max_wave
        unit : str

        Returns
        -------
        wavelength : float
        """
        return UnitTools.angstrom2unit(wave=phys_params.astrosat_bands_wave[band][wave_estimator], unit=unit)

    @staticmethod
    def get_obs_wave(band, obs, target=None, instrument=None, wave_estimator='mean_wave', unit='mu'):
        if obs == 'hst':
            if (target is None) & (instrument is None):
                raise KeyError('Either target or instrument needs to be given!')
            return ObsTools.get_hst_band_wave(band=band, instrument=instrument, wave_estimator=wave_estimator,
                                              unit=unit)
        elif obs == 'jwst':
            if instrument is None:
                raise KeyError('instrument is needed to distinguish between nircam and miri')
            return ObsTools.get_jwst_band_wave(band=band, instrument=instrument, wave_estimator=wave_estimator,
                                               unit=unit)
        elif obs == 'astrosat':
            return ObsTools.get_astrosat_band_wave(band=band, wave_estimator=wave_estimator, unit=unit)
        elif obs == 'roman':
            return ObsTools.get_roman_band_wave(band=band, wave_estimator=wave_estimator, unit=unit)
        else:
            raise KeyError(obs, ' is not understand')

    @staticmethod
    def get_hst_obs_band_list(target):
        """
        gets list of bands of HST
        Parameters
        ----------
        target : str

        Returns
        -------
        band_list : list
        """
        acs_band_list = obs_info.hst_obs_band_dict[target]['acs']
        uvis_band_list = obs_info.hst_obs_band_dict[target]['uvis']
        acs_uvis_band_list = obs_info.hst_obs_band_dict[target]['acs_uvis']
        ir_band_list = obs_info.hst_obs_band_dict[target]['ir']
        band_list = acs_band_list + uvis_band_list + acs_uvis_band_list + ir_band_list
        wave_list = []
        for band in acs_band_list:
            wave_list.append(ObsTools.get_hst_band_wave(band=band))
        for band in uvis_band_list:
            wave_list.append(ObsTools.get_hst_band_wave(band=band, instrument='uvis'))
        for band in acs_uvis_band_list:
            wave_list.append(ObsTools.get_hst_band_wave(band=band, instrument='uvis'))
        for band in ir_band_list:
            wave_list.append(ObsTools.get_hst_band_wave(band=band, instrument='ir'))


        return ObsTools.sort_band_list(band_list=band_list, wave_list=wave_list)

    @staticmethod
    def get_hst_obs_broad_band_list(target, allow_medium_bands=True):
        """
        gets list of bands of HST
        Parameters
        ----------
        target : str

        Returns
        -------
        band_list : list
        """
        acs_band_list = obs_info.hst_obs_band_dict[target]['acs']
        uvis_band_list = obs_info.hst_obs_band_dict[target]['uvis']
        acs_uvis_band_list = obs_info.hst_obs_band_dict[target]['acs_uvis']
        # ir_band_list = obs_info.hst_obs_band_dict[target]['ir']
        band_list = acs_band_list + uvis_band_list + acs_uvis_band_list #+ ir_band_list
        wave_list = []
        for band in acs_band_list:
            wave_list.append(ObsTools.get_hst_band_wave(band=band))
        for band in uvis_band_list:
            wave_list.append(ObsTools.get_hst_band_wave(band=band, instrument='uvis'))
        for band in acs_uvis_band_list:
            wave_list.append(ObsTools.get_hst_band_wave(band=band, instrument='uvis'))
        # for band in ir_band_list:
        #     wave_list.append(ObsTools.get_hst_band_wave(band=band, instrument='ir'))

        # kick out bands which are not broad bands
        for band, wave in zip(band_list, wave_list):
            if allow_medium_bands:
                if (band[-1] != 'W') & (band[-1] != 'M'):
                    band_list.remove(band)
                    wave_list.remove(wave)
            else:
                if band[-1] != 'W':
                    band_list.remove(band)
                    wave_list.remove(wave)

        return ObsTools.sort_band_list(band_list=band_list, wave_list=wave_list)

    @staticmethod
    def get_hst_ha_band(target):
        """
        get the corresponding H-alpha band for a target

        Parameters
        ----------
        target :  str

        Returns
        -------
        target_name : str
        """
        if (('F657N' in obs_info.hst_obs_band_dict[target]['uvis']) |
                ('F657N' in obs_info.hst_obs_band_dict[target]['acs'])):
            return 'F657N'
        elif (('F658N' in obs_info.hst_obs_band_dict[target]['uvis']) |
              ('F658N' in obs_info.hst_obs_band_dict[target]['acs'])):
            return 'F658N'
        else:
            raise KeyError(target, ' has no H-alpha observation ')

    @staticmethod
    def get_nircam_obs_band_list(target, version=None):
        """
        gets list of bands of NIRCAM bands
        Parameters
        ----------
        target : str
        version : str
        Returns
        -------
        band_list : list
        """

        if version is None: version = obs_info.nircam_available_data_versions[-1]
        nircam_band_list = getattr(obs_info, 'jwst_obs_band_dict_%s' % version)[target]['nircam_observed_bands']
        wave_list = []
        for band in nircam_band_list:
            wave_list.append(ObsTools.get_jwst_band_wave(band=band))
        return ObsTools.sort_band_list(band_list=nircam_band_list, wave_list=wave_list)

    @staticmethod
    def get_miri_obs_band_list(target, version=None):
        """
        gets list of bands of MIRI bands
        Parameters
        ----------
        target : str
        Returns
        -------
        band_list : list
        """
        if version is None: version = obs_info.miri_available_data_versions[-1]
        miri_band_list = getattr(obs_info, 'jwst_obs_band_dict_%s' % version)[target]['miri_observed_bands']
        wave_list = []
        for band in miri_band_list:
            wave_list.append(ObsTools.get_jwst_band_wave(band=band, instrument='miri'))
        return ObsTools.sort_band_list(band_list=miri_band_list, wave_list=wave_list)

    @staticmethod
    def get_astrosat_obs_band_list(target, version=None):
        """
        gets list of bands of HST
        Parameters
        ----------
        target : str
        Returns
        -------
        band_list : list
        """
        if version is None: version = obs_info.astrosat_available_data_versions[-1]
        astrosat_band_list = getattr(obs_info, 'astrosat_obs_band_dict_%s' % version)[target]['observed_bands']
        wave_list = []
        for band in astrosat_band_list:
            wave_list.append(ObsTools.get_astrosat_band_wave(band=band))
        return ObsTools.sort_band_list(band_list=astrosat_band_list, wave_list=wave_list)

    @staticmethod
    def get_band_list2obs_list(band_list, target, hst_target=None, nircam_version=None, miri_version=None):

        if hst_target is None:
            hst_target = target

        hst_band_list = ObsTools.get_hst_obs_band_list(target=hst_target)
        nircam_band_list = ObsTools.get_nircam_obs_band_list(target=target, version=nircam_version)
        miri_band_list = ObsTools.get_miri_obs_band_list(target=target, version=miri_version)

        obs_list = []
        for band in band_list:
            # get observations
            if band in hst_band_list:
                obs_list.append('hst')
            elif band in nircam_band_list:
                obs_list.append('jwst')
            elif band in miri_band_list:
                obs_list.append('jwst')
            else:
                raise RuntimeError(band, ' is in no observation present.')
        return obs_list

    @staticmethod
    def get_band_list2instrument_list(band_list, target, hst_target=None, nircam_version=None, miri_version=None):

        if hst_target is None:
            hst_target = target

        hst_band_list = ObsTools.get_hst_obs_band_list(target=hst_target)
        nircam_band_list = ObsTools.get_nircam_obs_band_list(target=target, version=nircam_version)
        miri_band_list = ObsTools.get_miri_obs_band_list(target=target, version=miri_version)

        instrument_list = []
        for band in band_list:
            # get observations
            if band in hst_band_list:
                instrument_list.append(
                    ObsTools.get_hst_instrument(target=hst_target, band=band))
            elif band in nircam_band_list:
                instrument_list.append('nircam')
            elif band in miri_band_list:
                instrument_list.append('miri')
            else:
                raise RuntimeError(band, ' is in no observation present.')
        return instrument_list

    @staticmethod
    def sort_band_list(band_list, wave_list):
        """
        sorts a band list with increasing wavelength
        Parameters
        ----------
        band_list : list
        wave_list : list
        Returns
        -------
        sorted_band_list : list
        """
        # sort wavelength bands
        sort = np.argsort(wave_list)
        return list(np.array(band_list)[sort])

    @staticmethod
    def split_obs_band_list(band_list=None,
                            target_name=None, hst_target_name=None, nircam_target_name=None, miri_target_name=None,
                            astrosat_target_name=None,
                            nircam_data_ver=None, miri_data_ver=None, astrosat_data_ver=None,
                            include_hst=True, include_nircam=True, include_miri=True, include_astrosat=True):
        """
        Function to get band lists that are split up by observations

        Parameters
        ----------
        band_list :  None, list or array-like
        target_name : None
        hst_target_name :  str
        nircam_target_name :  str
        miri_target_name :  str
        astrosat_target_name :  str
        nircam_data_ver :  str
        miri_data_ver :  str
        astrosat_data_ver :  str
        include_hst : bool
        include_nircam : bool
        include_miri : bool
        include_astrosat :  bool

        Returns
        -------
        band_list, hst_band_list, nircam_band_list, miri_band_list, astrosat_band_list : list
        """

        # organize the target names
        if (hst_target_name is None) & include_hst:
            if target_name is None:
                raise KeyError('If you want to include hst you either need to provide target_name or hst_target_name')
            hst_target_name = target_name
        if (nircam_target_name is None) & include_nircam:
            if target_name is None:
                raise KeyError('If you want to include nircam you either need to provide target_name or '
                               'nircam_target_name')
            nircam_target_name = target_name
        if (miri_target_name is None) & include_miri:
            if target_name is None:
                raise KeyError('If you want to include miri you either need to provide target_name or miri_target_name')
            miri_target_name = target_name
        if (astrosat_target_name is None) & include_astrosat:
            if target_name is None:
                raise KeyError('If you want to include astrosat you either need to provide target_name or '
                               'astrosat_target_name')
            astrosat_target_name = target_name

        # get entire band list to fit
        if include_hst:
            complete_hst_band_list = ObsTools.get_hst_obs_band_list(target=hst_target_name)
        else:
            complete_hst_band_list = []
        if include_nircam:
            complete_nircam_band_list = ObsTools.get_nircam_obs_band_list(
                target=nircam_target_name, version=nircam_data_ver)
        else:
            complete_nircam_band_list = []

        if include_miri:
            complete_miri_band_list = ObsTools.get_miri_obs_band_list(
                target=miri_target_name, version=miri_data_ver)
        else:
            complete_miri_band_list = []

        if include_astrosat:
            complete_astrosat_band_list = ObsTools.get_astrosat_obs_band_list(
                target=astrosat_target_name, version=astrosat_data_ver)
        else:
            complete_astrosat_band_list = []

        # put together band list for the fitting
        if band_list is None:
            hst_band_list = complete_hst_band_list
            nircam_band_list = complete_nircam_band_list
            miri_band_list = complete_miri_band_list
            astrosat_band_list = complete_astrosat_band_list
            band_list = []
            if include_hst:
                band_list += hst_band_list
            if include_nircam:
                band_list += nircam_band_list
            if include_miri:
                band_list += miri_band_list
            if include_astrosat:
                band_list += astrosat_band_list
        else:
            hst_band_list = []
            nircam_band_list = []
            miri_band_list = []
            astrosat_band_list = []
            for band in band_list:
                if band in complete_hst_band_list:
                    hst_band_list.append(band)
                elif band in complete_nircam_band_list:
                    nircam_band_list.append(band)
                elif band in complete_miri_band_list:
                    miri_band_list.append(band)
                elif band in complete_astrosat_band_list:
                    astrosat_band_list.append(band)

        return band_list, hst_band_list, nircam_band_list, miri_band_list, astrosat_band_list

    @staticmethod
    def get_obs_list(band_list, hst_band_list, nircam_band_list, miri_band_list, astrosat_band_list,
                     instrument_list=None,
                     target_name=None, hst_target_name=None, nircam_target_name=None, miri_target_name=None,
                     astrosat_target_name=None,):
        """
        Function to get list of observations instruments and target names based on band lists

        Parameters
        ----------
        band_list :  None, list or array-like
        hst_band_list : list or array-like
        nircam_band_list : list or array-like
        miri_band_list : list or array-like
        astrosat_band_list : list or array-like
        instrument_list :  None, list or array-like
        target_name : None
        hst_target_name : str
        nircam_target_name : str
        miri_target_name : str
        astrosat_target_name : str

        Returns
        -------
        obs_list, instrument_list, target_name_list : list
        """
        if hst_target_name is None:
            if target_name is None:
                raise KeyError('If you want to include hst you either need to provide target_name or hst_target_name')
            hst_target_name = target_name
        if nircam_target_name is None:
            if target_name is None:
                raise KeyError('If you want to include nircam you either need to provide target_name or '
                               'nircam_target_name')
            nircam_target_name = target_name
        if miri_target_name is None:
            if target_name is None:
                raise KeyError('If you want to include miri you either need to provide target_name or miri_target_name')
            miri_target_name = target_name
        if astrosat_target_name is None:
            if target_name is None:
                raise KeyError('If you want to include astrosat you either need to provide target_name or '
                               'astrosat_target_name')
            astrosat_target_name = target_name

        obs_list = []
        if instrument_list is None:
            get_instrument_flag = True
            instrument_list = []
        else:
            get_instrument_flag = False

        target_name_list = []

        for band in band_list:
            # get observations
            if band in hst_band_list:
                obs_list.append('hst')
                if get_instrument_flag:
                    instrument_list.append(ObsTools.get_hst_instrument(target=hst_target_name, band=band))
                target_name_list.append(hst_target_name)
            elif band in nircam_band_list:
                obs_list.append('jwst')
                if get_instrument_flag:
                    instrument_list.append('nircam')
                target_name_list.append(nircam_target_name)
            elif band in miri_band_list:
                obs_list.append('jwst')
                if get_instrument_flag:
                    instrument_list.append('miri')
                target_name_list.append(miri_target_name)
            elif band in astrosat_band_list:
                obs_list.append('astrosat')
                if get_instrument_flag:
                    instrument_list.append('astrosat')
                target_name_list.append(astrosat_target_name)
            else:
                raise RuntimeError(band, ' is in no observation present.')

        return obs_list, instrument_list, target_name_list

    @staticmethod
    def get_obs_band_dict(band_list=None, target_name=None, hst_target_name=None, nircam_target_name=None,
                          miri_target_name=None, astrosat_target_name=None, nircam_data_ver=None, miri_data_ver=None,
                          astrosat_data_ver=None, include_hst=True, include_nircam=True, include_miri=True,
                          include_astrosat=True):
        """
        Function to create a dict for the

        Parameters
        ----------
        band_list :  None, list or array-like
        target_name : None
        hst_target_name :  str
        nircam_target_name :  str
        miri_target_name :  str
        astrosat_target_name :  str
        nircam_data_ver :  str
        miri_data_ver :  str
        astrosat_data_ver :  str
        include_hst : bool
        include_nircam : bool
        include_miri : bool
        include_astrosat :  bool

        Returns
        -------
        obs_dict : dict
        """

        band_list, hst_band_list, nircam_band_list, miri_band_list, astrosat_band_list = (
            ObsTools.split_obs_band_list(band_list=band_list, target_name=target_name, hst_target_name=hst_target_name,
                                         nircam_target_name=nircam_target_name, miri_target_name=miri_target_name,
                                         astrosat_target_name=astrosat_target_name, nircam_data_ver=nircam_data_ver,
                                         miri_data_ver=miri_data_ver, astrosat_data_ver=astrosat_data_ver,
                                         include_hst=include_hst, include_nircam=include_nircam,
                                         include_miri=include_miri, include_astrosat=include_astrosat))

        obs_list, instrument_list, target_name_list = (
            ObsTools.get_obs_list(band_list=band_list, hst_band_list=hst_band_list, nircam_band_list=nircam_band_list,
                                  miri_band_list=miri_band_list, astrosat_band_list=astrosat_band_list,
                                  target_name=target_name, hst_target_name=hst_target_name,
                                  nircam_target_name=nircam_target_name, miri_target_name=miri_target_name,
                                  astrosat_target_name=astrosat_target_name,))

        obs_dict = {}

        # create dictionary
        if include_hst & (len(hst_band_list) > 0):

            mask_acs = np.array(instrument_list) == 'acs'
            mask_uvis = np.array(instrument_list) == 'uvis'
            mask_ir = np.array(instrument_list) == 'ir'

            if sum(mask_acs) > 0:
                obs_dict.update({'hst': {'acs': list(np.array(band_list)[mask_acs])}})
            if sum(mask_uvis) > 0:
                obs_dict.update({'hst': {'uvis': list(np.array(band_list)[mask_uvis])}})
            if sum(mask_ir) > 0:
                obs_dict.update({'hst': {'ir': list(np.array(band_list)[mask_ir])}})

        if include_nircam & (len(nircam_band_list) > 0):
            obs_dict.update({'jwst': {'nircam': nircam_band_list}})

        if include_miri & (len(miri_band_list) > 0):
            if 'jwst' in obs_dict.keys():
                obs_dict['jwst'].update({'miri': miri_band_list})
            else:
                obs_dict.update({'jwst': {'miri': miri_band_list}})

        if include_astrosat & (len(astrosat_band_list) > 0):
            obs_dict.update({'astrosat': astrosat_band_list})

        return obs_dict

    @staticmethod
    def filter_name2hst_band(target, filter_name):
        """
        Method to get from band-pass filter names to the HST filter names used for this observation.
        """
        if filter_name == 'NUV': return 'F275W'
        elif filter_name == 'U': return 'F336W'
        elif filter_name == 'B':
            if 'F438W' in obs_info.hst_obs_band_dict[target]['uvis']: return 'F438W'
            else: return 'F435W'
        elif filter_name == 'V': return 'F555W'
        elif filter_name == 'I': return 'F814W'
        elif filter_name == 'Ha': return ObsTools.get_hst_ha_band(target=target)
        else:
            raise KeyError(filter_name, ' is not available ')

    @staticmethod
    def hst_band2filter_name(band, latex=True):
        """
        Method to get from HST filter names used for this observation to the band pass name.
        """
        if band == 'F275W': return 'NUV'
        elif band == 'F336W': return 'U'
        elif (band == 'F438W') | (band == 'F435W'): return 'B'
        elif band == 'F555W': return 'V'
        elif band == 'F814W': return 'I'
        elif (band == 'F657N') | (band == 'F658N'):
            if latex: return r'H$\alpha$'
            else: return 'Ha'
        else: raise KeyError(band, ' is not available ')

    @staticmethod
    def check_hst_obs(target):
        """
        check if HST has any observation for target
        Parameters
        ----------
        target :  str
        Returns
        -------
        observation flag : bool
        """

        if target in obs_info.hst_obs_band_dict.keys():
            return True
        else:
            return False

    @staticmethod
    def check_hst_broad_band_obs(target):
        """
        check if HST has broad band observation available for target
        Parameters
        ----------
        target :  str
        Returns
        -------
        observation flag : bool
        """

        if ObsTools.check_hst_obs(target=target):
            band_list = ObsTools.get_hst_obs_broad_band_list(target=target)
            for band in band_list:
                if band[-1] == 'W':
                    return True
            return False
        else:
            return False

    @staticmethod
    def check_hst_ha_obs(target):
        """
        check if HST has broad band observation available for target
        Parameters
        ----------
        target :  str
        Returns
        -------
        observation flag : bool
        """

        if (('F657N' in ObsTools.get_hst_obs_band_list(target=target)) |
                ('F658N' in ObsTools.get_hst_obs_band_list(target=target))):
            return True
        else: return False

    @staticmethod
    def check_hst_ha_cont_sub_obs(target):
        """
        boolean function checking if H-alpha band for a target exists
        """
        if target in list(obs_info.hst_ha_cont_sub_dict.keys()):
            return True
        else: return False

    @staticmethod
    def check_nircam_obs(target, version=None):
        """
        check if NIRCAM observed
        """
        if version is None: version = obs_info.nircam_available_data_versions[-1]
        jwst_obs_band_dict = getattr(obs_info, 'jwst_obs_band_dict_%s' % version)

        if target in jwst_obs_band_dict.keys():
            if jwst_obs_band_dict[target]['nircam_observed_bands']: return True
            else: return False
        else: return False

    @staticmethod
    def check_miri_obs(target, version=None):
        """
        check if MIRI observed
        """
        if version is None: version = obs_info.miri_available_data_versions[-1]
        jwst_obs_band_dict = getattr(obs_info, 'jwst_obs_band_dict_%s' % version)

        if target in jwst_obs_band_dict.keys():
            if jwst_obs_band_dict[target]['miri_observed_bands']: return True
            else: return False
        else: return False

    @staticmethod
    def check_astrosat_obs(target, version=None):
        """
        Check for astrosat obs
        """
        if version is None: version = obs_info.astrosat_available_data_versions[-1]
        astrosat_obs_band_dict = getattr(obs_info, 'astrosat_obs_band_dict_%s' % version)

        if target in astrosat_obs_band_dict.keys(): return True
        else: return False

    @staticmethod
    def check_muse_obs(target):
        """
        check if MUSE observation is available for target
        Parameters
        ----------
        target :  str
        Returns
        -------
        observation flag : bool
        """
        if target in obs_info.phangs_muse_galaxy_list: return True
        else: return False

    @staticmethod
    def check_alma_obs(target):
        """
        check if ALMA observation is available for target
        Parameters
        ----------
        target :  str
        Returns
        -------
        observation flag : bool
        """
        if target in obs_info.phangs_alma_galaxy_list: return True
        else: return False


