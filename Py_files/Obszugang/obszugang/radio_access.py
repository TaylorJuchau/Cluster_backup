"""
Construct a data access structure for all kind of X-ray observation
"""
from pathlib import Path
from werkzeugkiste import helper_func, phys_params
from obszugang import access_config


class RadioAccess:
    """
    Access class to organize data structure of chnadra observations
    """

    def __init__(self, radio_target_name=None):
        """

        Parameters
        ----------
        radio_target_name : str
            Default None. Target name
        radio_target_name : str
            Default None. Target name used for Hs observation
        """

        # get target specifications
        # check if the target names are compatible

        if radio_target_name is not None:
            radio_target_name = helper_func.FileTools.target_name_no_directions(target=radio_target_name)

        self.radio_target_name = radio_target_name

        # loaded data dictionaries
        self.radio_data = {}

        # get path to observation coverage hulls
        self.path2obs_cover_hull = (Path(__file__).parent.parent.absolute() / 'meta_data' / 'obs_coverage' /
                                    'data_output')

        super().__init__()

    def get_vla_file_name(self, band='L', map='tt0'):
        """

        Parameters
        ----------
        band : str
        map : str
        Returns
        -------
        data_file_path : ``Path``
        """

        file_path = Path(access_config.phangs_config_dict['radio_data_path'])
        if self.radio_target_name == 'ngc5194':
            file_name = 'm51_A_%s_briggs_0p5.pbcor.image.%s.fits' % (band, map)
        else:
            raise KeyError('Only M51 target has been added.')
        return file_path / file_name

    def load_radio_data(self, band='L', map='tt0', reload=False):
        """

        Parameters
        ----------
        band : str
        map : str
        reload : bool
        """
        # check if data is already loaded
        if (('%s_%s_data_img' % (band, map)) in self.radio_data.keys()) & (not reload):
            return None

        file_name = self.get_vla_file_name(band=band, map=map)

        img_data, img_header, img_wcs = helper_func.FileTools.load_img(file_name=file_name)

        self.radio_data.update({'%s_%s_data_img' % (band, map): img_data,
                                '%s_%s_header_img' % (band, map): img_header,
                                '%s_%s_wcs_img' % (band, map): img_wcs,
                                '%s_%s_unit_img' % (band, map): 'counts',
                                '%s_%s_pixel_area_size_sr_img' % (band, map):
                                    img_wcs.proj_plane_pixel_area().value * phys_params.sr_per_square_deg})
