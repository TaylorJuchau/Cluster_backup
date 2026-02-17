"""
Combination of all available observations
"""

from obszugang import phot_access, gas_access, spec_access, x_access, radio_access


class ObsAccess(phot_access.PhotAccess, gas_access.GasAccess, spec_access.SpecAccess, x_access.XAccess, radio_access.RadioAccess):
    """
    Class to plot cutouts in multiple bands
    """

    def __init__(self, target_name=None, phot_hst_target_name=None, phot_hst_ha_cont_sub_target_name=None,
                 phot_nircam_target_name=None, phot_miri_target_name=None, phot_astrosat_target_name=None,
                 x_target_name=None, radio_target_name=None,
                 nircam_data_ver='v1p1p1', miri_data_ver='v1p1p1', astrosat_data_ver='v1p0',
                 nirspec_data_ver=None, miri_mrs_data_ver=None):
        phot_access.PhotAccess.__init__(self,
                                        phot_target_name=target_name,
                                        phot_hst_target_name=phot_hst_target_name,
                                        phot_hst_ha_cont_sub_target_name=phot_hst_ha_cont_sub_target_name,
                                        phot_nircam_target_name=phot_nircam_target_name,
                                        phot_miri_target_name=phot_miri_target_name,
                                        phot_astrosat_target_name=phot_astrosat_target_name,
                                        nircam_data_ver=nircam_data_ver,
                                        miri_data_ver=miri_data_ver,
                                        astrosat_data_ver=astrosat_data_ver)
        gas_access.GasAccess.__init__(self, gas_target_name=target_name)
        spec_access.SpecAccess.__init__(self, spec_target_name=target_name,
                                        nirspec_data_ver=nirspec_data_ver,
                                        miri_mrs_data_ver=miri_mrs_data_ver)
        x_access.XAccess.__init__(self, x_target_name=x_target_name)
        radio_access.RadioAccess.__init__(self, radio_target_name=radio_target_name)
