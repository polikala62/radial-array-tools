'''
Created on Mar 12, 2024

@author: Karl
'''

from main import vis_ncdf

# (out_ncdf, pr_gdb, in_dem, lmark_fc, land_fc, max_dist, x_range, y_range, z_range, coordinate_system)
vis_ncdf(r"C:\GIS\CAA_2024\Test_Data\test_main_01.nc", 
         r"C:\GIS\CAA_2024\script_testing.gdb", 
         r"C:\GIS\CAA_2023\Data\Script_Datasets\DEM\SRTM_GEBCO_m.tif", 
         r"C:\GIS\CAA_2024\Test_Data\test_landmarks.shp", 
         r"C:\GIS\CAA_2024\Test_Data\test_coast_poly.shp", 
         10000,
         [3245000, 3260000, 2500], [3105000, 3125000, 2500], [0, 1, 1]
         )

#[3258000, 3261000, 1000], [3109000, 3111000, 1000], [0, 3, 1]