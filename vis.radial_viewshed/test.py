'''
Created on Dec 13, 2023

@author: Karl
'''

import arcpy
from main import vis_ncdf

out_ncdf = r"C:\GIS\ArcPro_Projects\Visibility_Testing\Output\test_01.nc"

pr_gdb = r'C:\GIS\ArcPro_Projects\Visibility_Testing\Visibility_Testing.gdb'
pt_mask_fc = r'C:\GIS\ArcPro_Projects\Visibility_Testing\Features\pt_mask_01.shp'
in_dem = r'C:\GIS\ArcPro_Projects\PM_Med_Survey\Rasters\gmted_mea150_clip_02.tif'
land_fc = r'C:\GIS\ArcPro_Projects\PM_Med_Survey\Features\grid_test\obstruction_features_etrs.shp'

xy_spacing = 10000
dist_range = [0,100000,1000]
z_range = [0,1,1]
sample_raster = r'C:\GIS\ArcPro_Projects\PM_Med_Survey\Rasters\gmted_mea150_clip_02.tif'
lmark_fc = ""
        
vis_ncdf(out_ncdf, pr_gdb, in_dem, pt_mask_fc, land_fc, xy_spacing, dist_range, z_range, sample_raster, lmark_fc)

print("Script finished!")