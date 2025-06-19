'''
Created on 4 Feb 2025

@author: Karl Smith
'''
import arcpy
from functions import clip, list, viewshed

obs_z_list = [1]
dist_list = list.wmo_dist_list_2(200000)
in_dem_ras = r'C:\GIS\ArcPro_Projects\Visibility_Testing\DEM\malta_test_dem_2.tif'
in_dem_res = 482.3567932855411
pr_gdb = r'C:\GIS\ArcPro_Projects\Visibility_Testing\Visibility_Testing.gdb'
pt_crs = arcpy.Describe(r'C:\GIS\ArcPro_Projects\Visibility_Testing\Features\pt_mask_malta_03.shp').spatialReference
densify_dist = 500

land_fc = r'C:\GIS\ArcPro_Projects\Visibility_Testing\Features\obstruction_med.shp'

pt_list = [[4687310,1390380], [4696360, 1396740]]
pt_list = [[4710700, 1424880]]

#------------------------------------------------------------------------------ 

# Get extent and resolution from DEM.
dem_extent = [float(arcpy.management.GetRasterProperties(in_dem_ras, "LEFT").getOutput(0)),
              float(arcpy.management.GetRasterProperties(in_dem_ras, "BOTTOM").getOutput(0)),
              float(arcpy.management.GetRasterProperties(in_dem_ras, "RIGHT").getOutput(0)),
              float(arcpy.management.GetRasterProperties(in_dem_ras, "TOP").getOutput(0))]

dem_resolution = max([float(arcpy.management.GetRasterProperties(in_dem_ras, "CELLSIZEX").getOutput(0)),
                      float(arcpy.management.GetRasterProperties(in_dem_ras, "CELLSIZEY").getOutput(0))])

# Read coast polylines to list of geometry objects.
if land_fc == "":
    land_poly = None
    
else:
    
    # Clip input land mask features to fix raster.
    clip_land_mask = clip.clip_polygon(dem_extent, land_fc, pt_crs, clip_poly_path=r"memory\clip_land_mask")
    
    land_poly = []
    
    with arcpy.da.SearchCursor(clip_land_mask, ["SHAPE@"]) as cursor: #@UndefinedVariableFromImport
        for row in cursor:
            for part in row:
                land_poly.append(part)
    
    # Delete search cursor.
    del cursor
    
    # Delete clipped features.
    arcpy.Delete_management(clip_land_mask)

#------------------------------------------------------------------------------ 

for pt in pt_list:
    
    obs_x, obs_y = pt
    
    obs_dict, benchmark_dict = viewshed.radial_viewshed(obs_x, obs_y, obs_z_list, dist_list, in_dem_ras, in_dem_res, pr_gdb, pt_crs, densify_dist,
                                                        obs_z_offset=0, lmark_geom_list="", land_fc=land_poly, override_min_dist="")
    
    print(obs_dict)
    
    # Obs_dict is dictionary where key is observer height and value is distance dictionary.
    for obs_z in obs_dict.keys():
        
        # Find index for z value.
        z_idx = list.list_index_from_val(obs_z_list, obs_z)
        
        # Dist_dict is dictionary where key is visible distance and value is tuple.
        vdist_dict = obs_dict[obs_z]
        
        # Tuple is format [v_angle_sum, v_angle_range, h_angle_sum, *sample_max, *sample_min, *sample_avg, *landmark_sum]
        for dist_val in vdist_dict.keys():
            
            # Find index for distances.
            d_idx = list.list_index_from_val(dist_list, dist_val)
            
            v_angle_sum, v_angle_max, h_angle_sum = vdist_dict[dist_val][0:3]
            print(obs_z, d_idx, v_angle_sum, v_angle_max, h_angle_sum)
    
print('Script finished!')