'''
Created on Mar 12, 2024

@author: Karl
'''

from main import vis_ncdf

vis_ncdf(r"C:\GIS\Shaw\Waldron_NetCDF\netCDF\TEST_03.nc", 
         r"C:\GIS\Shaw\Waldron_NetCDF\Waldron_NetCDF.gdb", 
         r"C:\GIS\Shaw\Waldron_NetCDF\Waldron_NetCDF.gdb\waldron_dem", 
         r"C:\GIS\Shaw\Waldron_NetCDF\Waldron_NetCDF.gdb\waldron_buffer_5kft_erase",
         r"C:\GIS\Shaw\Waldron_NetCDF\Waldron_NetCDF.gdb\waldron_poly",
         500,
         200,
         [1, 1, 1]
         )


'''

vis_ncdf(r"C:\GIS\Shaw\Waldron_NetCDF\netCDF\TEST_01.nc", 
         r"C:\GIS\Shaw\Waldron_NetCDF\Waldron_NetCDF.gdb", 
         r"C:\GIS\Shaw\Waldron_NetCDF\Waldron_NetCDF.gdb\waldron_dem", 
         r"C:\GIS\Shaw\Waldron_NetCDF\Waldron_NetCDF.gdb\waldron_buffer_5kft_erase",
         r"C:\GIS\Shaw\Waldron_NetCDF\Waldron_NetCDF.gdb\waldron_poly",
         18000,
         200,
         [1, 1, 1]
         )


vis_ncdf(r"C:\GIS\Maritime_Encounters\OrmeSim\Output\TEST_04.nc", 
         r"C:\GIS\CAA_2024\script_testing.gdb", 
         r"C:\GIS\Maritime_Encounters\OrmeSim\GMTED20_dem_rc.tif", 
         r"C:\GIS\Maritime_Encounters\OrmeSim\GMTED20_sea_mask_02.shp", 
         r"", 
         r"C:\GIS\Maritime_Encounters\OrmeSim\GMTED20_land_mask_02.shp", 
         50000,
         50000,
         [2924500, 4300000, 20000], [2800000, 3799900, 20000], [-8, 8, 1]
         )




vis_ncdf(r"C:\GIS\CAA_2024\Test_Data\test_main_01.nc", 
         r"C:\GIS\CAA_2024\script_testing.gdb", 
         r"C:\GIS\CAA_2023\Data\Script_Datasets\DEM\SRTM_GEBCO_m.tif",
         r"C:\GIS\CAA_2024\Test_Data\test_coast_poly.shp", 
         50000,
         [3245000, 3260000, 10000], [3105000, 3125000, 10000], [0, 1, 1]
         )


vis_ncdf(r"C:\GIS\CAA_2024\Test_Data\test_main_01.nc", 
         r"C:\GIS\CAA_2024\script_testing.gdb", 
         r"C:\GIS\CAA_2023\Data\Script_Datasets\DEM\SRTM_GEBCO_m.tif",
         r"C:\GIS\CAA_2024\Test_Data\test_coast_poly.shp", 
         10000,
         [3245000, 3260000, 2500], [3105000, 3125000, 2500], [0, 1, 1],
         lmark_fc=r"C:\GIS\CAA_2024\Test_Data\test_landmarks.shp"
         )
'''