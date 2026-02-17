"""
script to gather all local access configuration in a dictionary
"""

phangs_config_dict = {

    # sample table
    'phangs_sample_table_path': '/home/benutzer/data/PHANGS_products/sample_table',
    'phangs_sample_table_ver': 'v1p6',
    # gas access
    'alma_data_path': '/media/benutzer/derka_derka/data/alma',
    'alma_data_ver': 'v4p0',
    'alma_conv_map_data_path': '/media/benutzer/derka_derka/data/phangs_data_products/conversion_factor_maps',
    'alma_conv_map_data_ver': 'v2p0',
    'alma_cloud_cat_data_path': '/media/benutzer/derka_derka/data/phangs_data_products/cloud_catalogs/',
    'alma_cloud_cat_data_ver': 'v4p0',
    'alma_cloud_cat_data_release_ver': 'ST1p6',
    # photometric access
    'hst_data_path': '/media/benutzer/derka_derka/data/hst',
    'hst_ha_cont_sub_data_path': '/media/benutzer/derka_derka/data/hst/HST_Halpha_contsub_images',
    'hst_obs_hdr_file_path': '/home/benutzer/data/PHANGS_products/tables',
    'nircam_data_path': '/media/benutzer/derka_derka/data/jwst/',
    'miri_data_path': '/media/benutzer/derka_derka/data/jwst',
    'astrosat_data_path': '/media/benutzer/derka_derka/data/astrosat',
    'hst_ha_cont_sub_ver': 'v1p1',
    # X-ray access
    'chandra_data_path': '/media/benutzer/derka_derka/data/chandra/MAIN_IMAGE_PRODUCTS',
    # X-ray access
    'radio_data_path': '/media/benutzer/derka_derka/data/radio/VLA/',
    # JWST spectra
    'nirspec_data_path': '/media/benutzer/derka_derka/data/jwst/',
    'miri_mrs_data_path': '/media/benutzer/derka_derka/data/jwst/',
    # HST cluster catalog
    'phangs_hst_cluster_cat_data_path': '/home/benutzer/data/PHANGS_products/HST_catalogs',
    'phangs_hst_cluster_cat_release': 'phangs_hst_cc_dr4_cr3_public',
    # 'phangs_hst_cluster_cat_ver': 'v1',
    'phangs_hst_cluster_cat_quick_access_path':
        '/media/benutzer/derka_derka/data/phangs_data_products/cluster_cat_quick_access',
    # internal releases
    'phangs_hst_cluster_cat_extend_photo_path':
        '/media/benutzer/derka_derka/data/phangs_data_products/jimane_photometry',
    'phangs_hst_cluster_cat_extend_photo_ver': 'v5p4',
    'phangs_hst_cluster_cat_extend_sed_fit_path':
        '/media/benutzer/derka_derka/data/phangs_data_products/kiana_sed_fit/NUV_optical_NIR_final_fits/',
    # ALL SOURCE CATALOGS
    # 'phangs_hst_all_src_cat': '/home/benutzer/data/PHANGS_products/HST_catalogs/',
    'hst_all_src_cat_drive_path': 'scratch/HUBBLE_AND_CLUSTER_TECHNICAL_WORK/HST_catalogs/Cluster_pipeline_intermediate_stage_tables/',
    'hst_all_src_cat_local_path': '/media/benutzer/derka_derka/data/phangs_data_products/hst_all_src_cat/',

    # scale_decompositions
    'scale_decomposition_path': '/media/benutzer/derka_derka/data/phangs_data_products/scale_decompositions/',

    # MUSE data access
    'muse_data_path': '/media/benutzer/derka_derka/data/muse',
    'muse_data_ver': 'DR2.2',
    # KCWI data access
    'kcwi_data_path': '/media/benutzer/derka_derka/data/kcwi',

    # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    # !!!!! build data access !!!!!
    # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    # rclone name
    'rclone_name': 'drive',
    # data destinations for local storage
    'hst_obs_data_local_path': '/media/benutzer/derka_derka/data/hst/HST_reduced_images',
    'jwst_obs_data_local_path_v2p0': '/media/benutzer/derka_derka/data/jwst/v2p0',
    # path where to find HST data on Drive
    'hst_obs_data_drive_path': 'scratch/HUBBLE_AND_CLUSTER_TECHNICAL_WORK/HST_image_products/HST_reduced_images/',
    # path where to find JWST data on Drive
    'jwst_obs_data_drive_path_v2p0': 'Archive/JWST/v2p0/',
    # target list or use `all` to get all observed targets
    'hst_obs_target_list': 'all',


    # same for ha continuum subtracted images
    'hst_ha_cont_sub_data_local_path': '/media/benutzer/derka_derka/data/hst',
    'hst_ha_cont_sub_target_list': 'all',

}
