'''
Created on Dec 13, 2023

@author: Karl
'''

import os
from main import vis_ncdf
from functions import list

# Set parameters.
out_ncdf = r"C:\GIS\ArcPro_Projects\PM_Med_Visibility\Output\cyclades_500m_01.nc"
pr_gdb = r'C:\GIS\ArcPro_Projects\Visibility_Testing\Visibility_Testing.gdb'
pt_mask_fc = r'C:\GIS\ArcPro_Projects\PM_Med_Visibility\Feature\cyclades_vis_inputs\cyclades_pt_mask.shp'
in_dem = r'C:\GIS\ArcPro_Projects\PM_Med_Visibility\Raster\cyclades_vis_dem\cyclades_vis_dem_100km.tif'
land_fc = r'C:\GIS\ArcPro_Projects\PM_Med_Visibility\Feature\cyclades_vis_inputs\cyclades_obstructions_detailed_100km.shp'

xy_spacing = 500
dist_range = [0,70000,500]
z_range = [-2,2,0.5]

sample_raster = ""
lmark_fc = ""
out_dist_list = list.wmo_dist_list_2(70000)
densify_dist = 500
obs_z_offset = 1.5
write_output = True
#pt_mask_json = r'C:\GIS\ArcPro_Projects\PM_Med_Visibility\Json\{}.json'.format(os.path.basename(out_ncdf).split('.')[0])
pt_mask_json = ''


vis_ncdf(out_ncdf, pr_gdb, in_dem, pt_mask_fc, land_fc, xy_spacing, dist_range, densify_dist, z_range, obs_z_offset, sample_raster, lmark_fc, out_dist_list, pt_mask_json, write_output)

