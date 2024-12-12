'''
Created on Dec 13, 2023

@author: Karl
'''

import arcpy
#import funcs_visibility
'''
pr_gdb = r'C:\GIS\vis_radial_viewshed\Radial_Viewshed\Radial_Viewshed.gdb'
obs_pt_fc = r'C:\GIS\vis_radial_viewshed\Radial_Viewshed\Radial_Viewshed.gdb\tmp_pts'
#in_dem = r'C:\GIS\vis_radial_viewshed\Radial_Viewshed\Radial_Viewshed.gdb\SRTM_clip'
in_dem = r'C:\GIS\vis_radial_viewshed\Radial_Viewshed\Datasets\SRTM_GEBCO_m.tif'
coast_fc = r'C:\GIS\vis_radial_viewshed\Radial_Viewshed\Radial_Viewshed.gdb\coast_fc'


out_pts = r'C:\GIS\vis_radial_viewshed\Radial_Viewshed\Radial_Viewshed.gdb\vis_pts'
out_skyline = r'C:\GIS\vis_radial_viewshed\Radial_Viewshed\Radial_Viewshed.gdb\skyline'

obs_pt_list = []
with arcpy.da.SearchCursor(obs_pt_fc, ["SHAPE@X", "SHAPE@Y", "SHAPE@Z"]) as cursor: #@UndefinedVariableFromImport
    for row in cursor:
        obs_pt_list.append(row)
        
funcs_visibility.radial_viewshed(obs_pt_list[0], in_dem, pr_gdb, obs_z_offset=0, max_dist=50000, lmark_geom_list="", coast_fc=coast_fc, out_pts=out_pts, out_skyline=out_skyline)

print("Script finished!")
'''

sr = arcpy.SpatialReference(3035)
print(sr.name)
print(sr.linearUnitName)