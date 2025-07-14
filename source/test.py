'''
Created on Dec 13, 2023

@author: Karl
'''

import os
from main import vis_ncdf
from functions import list

# Set parameters.
out_ncdf = r"C:\GIS\ArcPro_Projects\Visibility_Testing\Output\med_only_5000m_01.nc"
pr_gdb = r'C:\GIS\ArcPro_Projects\Visibility_Testing\Visibility_Testing.gdb'
pt_mask_fc = r'C:\GIS\ArcPro_Projects\Visibility_Testing\Features\pt_mask_medonly.shp'
in_dem = r'C:\GIS\ArcPro_Projects\Visibility_Testing\Rasters\DEM\gmted_mea150_medonly.tif'
land_fc = r'C:\GIS\ArcPro_Projects\Visibility_Testing\Features\obstruction_medonly.shp'

xy_spacing = 5000
dist_range = [0,150000,1000]
z_range = [-2,2,0.5]

sample_raster = ""
lmark_fc = ""
out_dist_list = list.wmo_dist_list_2(150000)
densify_dist = 500
obs_z_offset = 1.5
write_output = True
#pt_mask_json = r'C:\GIS\ArcPro_Projects\Visibility_Testing\json\{}.json'.format(os.path.basename(out_ncdf).split('.')[0])
pt_mask_json = ''


vis_ncdf(out_ncdf, pr_gdb, in_dem, pt_mask_fc, land_fc, xy_spacing, dist_range, densify_dist, z_range, obs_z_offset, sample_raster, lmark_fc, out_dist_list, pt_mask_json, write_output)

