'''
Created on Dec 13, 2023

@author: Karl
'''

import os
from main import vis_ncdf
from functions import list
# 
out_ncdf = r"C:\GIS\ArcPro_Projects\Visibility_Testing\Output\corsica_albedo_5km.nc"
pr_gdb = r'C:\GIS\ArcPro_Projects\Visibility_Testing\Visibility_Testing.gdb'
pt_mask_fc = r'C:\GIS\ArcPro_Projects\Visibility_Testing\Features\case_study_corsica\pt_mask_corsica_01.shp'
in_dem = r'C:\GIS\ArcPro_Projects\Visibility_Testing\Rasters\DEM\gmted_mea150_corsica.tif'
land_fc = r'C:\GIS\ArcPro_Projects\Visibility_Testing\Features\case_study_corsica\obstructions_corsica_01.shp'

xy_spacing = 5000
dist_range = [0,120000,1000]
out_dist_list = list.wmo_dist_list_2(120000)
z_range = [0,0,1]
#sample_raster = r'C:\GIS\ArcPro_Projects\PM_Med_Survey\Rasters\gmted_mea150_clip_02.tif'
sample_raster = r'C:\GIS\ArcPro_Projects\Visibility_Testing\Rasters\sample\corsica_albedo_test.tif'
#sample_raster = ''
lmark_fc = r'C:\GIS\ArcPro_Projects\Visibility_Testing\Features\case_study_corsica\coastal_prominences_corsica_01.shp'
lmark_fc = ''

densify_dist = 500
obs_z_offset = 1.5
write_output = True

pt_mask_json = r'C:\GIS\ArcPro_Projects\Visibility_Testing\json\{}.json'.format(os.path.basename(out_ncdf).split('.')[0])
pt_mask_json = ''

vis_ncdf(out_ncdf, pr_gdb, in_dem, pt_mask_fc, land_fc, xy_spacing, dist_range, densify_dist, z_range, obs_z_offset, sample_raster, lmark_fc, out_dist_list, pt_mask_json, write_output)