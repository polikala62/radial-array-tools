'''
Created on Dec 13, 2023

@author: Karl
'''

import arcpy
from main import vis_ncdf
from functions import list

out_ncdf = r"C:\GIS\ArcPro_Projects\Visibility_Testing\Output\med_test_01.nc"

pr_gdb = r'C:\GIS\ArcPro_Projects\Visibility_Testing\Visibility_Testing.gdb'
pt_mask_fc = r'C:\GIS\ArcPro_Projects\PM_Med_Survey\Features\med_vis\inputs\vis_buffer_04.shp'
in_dem = r'C:\GIS\DEM\GMTED2010_15_arcsec\gmted_mea150_rc_nulls.tif'
land_fc = r'C:\GIS\ArcPro_Projects\Visibility_Testing\Features\obstruction_med.shp'

xy_spacing = 10000
dist_range = [0,200000,1000]
z_range = [1,1,1]
#sample_raster = r'C:\GIS\ArcPro_Projects\PM_Med_Survey\Rasters\gmted_mea150_clip_02.tif'
sample_raster = ''
lmark_fc = ""
out_dist_list = list.wmo_dist_list_2(200000)
densify_dist = 500
obs_z_offset = 1.5
        
vis_ncdf(out_ncdf, pr_gdb, in_dem, pt_mask_fc, land_fc, xy_spacing, dist_range, densify_dist, z_range, obs_z_offset, sample_raster, lmark_fc, out_dist_list)

print("Script finished!")