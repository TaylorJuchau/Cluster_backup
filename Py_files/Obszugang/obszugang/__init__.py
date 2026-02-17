# noinspection PyPep8Naming

__version__ = 0.1
__author__ = 'Daniel Maschmann'

__all__ = [
    'ObsAccess',
    'ObsTools',
    'ClusterCatAccess',
    'PhangsSampleAccess',
    'obs_info',
    'phot_access',
    'gas_access',
    'spec_access',
    'x_access',
    'radio_access',

]

from obszugang.obs_tools import ObsTools
from obszugang.obs_access import ObsAccess
from obszugang.cluster_cat_access import ClusterCatAccess
from obszugang.gal_access import PhangsSampleAccess

import obszugang.obs_info
import obszugang.phot_access
import obszugang.gas_access
import obszugang.spec_access
import obszugang.x_access
import obszugang.radio_access
