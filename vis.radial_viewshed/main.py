'''
Created on Mar 7, 2024

@author: Karl
'''


import netCDF4 as nc
import numpy as np
import arcpy, os, random, datetime

import console
import funcs_intersect as fi
import funcs_netcdf as fn
import funcs_visibility as fv

def vis_ncdf(out_ncdf, pr_gdb, in_dem, lmark_fc, land_fc, max_dist, x_range, y_range, z_range):
    
    start_time = datetime.datetime.now()
    
    #===========================================================================
    # GET CRS.
    #===========================================================================
    
    dem_crs = arcpy.Describe(in_dem).spatialReference
    lmark_crs = arcpy.Describe(lmark_fc).spatialReference
    coast_crs = arcpy.Describe(land_fc).spatialReference
    
    '''
    print(dem_crs.name)
    print(lmark_crs.name)
    print(coast_crs.name)
    '''
    
    if dem_crs.name == lmark_crs.name and dem_crs.name == coast_crs.name and lmark_crs.name == coast_crs.name: # Yes, it's clumsy.
        
        out_crs = dem_crs
        out_crs_name = out_crs.name
        out_crs_units = out_crs.linearUnitName
    
    else:
        
        raise Exception("Input DEM and feature classes have different coordinate systems.")
    
    #===========================================================================
    # GENERATE ARRAYS.
    #===========================================================================
    console.console("Creating arrays for analysis...")
    
    # Get range values from input (they should all be lists in format: [min, max, increment]).
    if len(x_range) == 3 and len(y_range) == 3 and len(z_range) == 3: #@TODO: better validation.
        
        # Get range values as variables.
        x_min, x_max, x_step = x_range
        y_min, y_max, y_step = y_range
        z_min, z_max, z_step = z_range
        
    else:
        
        raise Exception("Could not read input ranges. Ranges for x, y, and z should be lists with format [minimum_value, maximum_value, increment].")
    
    # Create ranges for x, y, and z.
    x_range = np.arange(x_min, (x_max + x_step), x_step)
    y_range = np.arange(y_min, (y_max + y_step), y_step)
    z_range = np.arange(z_min, (z_max + z_step), z_step)
    
    # Create 2D arrays for x, and y.
    x_array = np.array([x_range for i in y_range])
    y_array = np.rot90(np.array([y_range for i in x_range]))
    
    # Create 3D array for z.
    z_slice = np.empty(x_array.shape)
    z_array = [z_slice.fill(i) for i in z_range]
    
    # Use ranges to create a shape for data outputs.
    out_array_shape = (len(z_range), len(y_range), len(x_range))
    
    # Create arrays for data outputs.
    ts_vis_array = np.zeros(out_array_shape, dtype=float, order='C')
    ms_vis_array = np.zeros(out_array_shape, dtype=float, order='C')
    lk_vis_array = np.zeros(out_array_shape, dtype=float, order='C')
    
    #===========================================================================
    # CALCULATE VISIBILITY VALUES.
    #===========================================================================
    console.console("Calculating visibility values...")
    
    # Prepare common datasets for use.
    
    # Get variables from DEM: extent, CRS.
    dem_crs = arcpy.Describe(in_dem).spatialReference
    # Get raster resolution (use the minimum).
    raster_resx = arcpy.management.GetRasterProperties(in_dem, "CELLSIZEX").getOutput(0)
    raster_resy = arcpy.management.GetRasterProperties(in_dem, "CELLSIZEY").getOutput(0)
    raster_res = int(min([raster_resx, raster_resy])) * 5
    
    in_dem_ras = arcpy.Raster(in_dem)
    
    # Read coast polylines to list of geometry objects.
    if land_fc == "":
        land_poly = None
        
    else:
        land_poly = []
        with arcpy.da.SearchCursor(land_fc, ["SHAPE@"]) as cursor: #@UndefinedVariableFromImport
            for row in cursor:
                for part in row:
                    land_poly.append(part)
    
    # Read landmark polygons to list of geometry objects.
    #@TODO: Scan for duplicate IDs and print warning?
    if lmark_fc == "":
        lmark_geom_list = None
        
    else:
        lmark_geom_list = []
        with arcpy.da.SearchCursor(lmark_fc, ["SHAPE@"]) as cursor: #@UndefinedVariableFromImport
            for row in cursor:
                lmark_geom_list.append(row[0])
    
    #------------------------------------------------------------------------------ 
    
    # Create counter for points processed.
    pr_count = 0
    pr_total = int(len(x_range) * len(y_range) * len(z_range))
    
    console.console("{} points found in array.".format(str(pr_total)),2)
    
    # Iterate through indices.
    for index, x in np.ndenumerate(ts_vis_array): #@UnusedVariable
        
        # Get indeces from iterator.
        z_idx, y_idx, x_idx = index
        
        # Get x, y, z values for iterated data point.
        z_val = float(z_range[z_idx])
        x_val = float(x_range[x_idx])
        y_val = float(y_range[y_idx])
        
        #print("Processing point at {} / {} / {}...".format(x_val, y_val, z_val))
        
        # Check if point is within mask.
        pr_pt = arcpy.PointGeometry(arcpy.Point(x_val, y_val))
        
        # Check that the point does not intersect the land mask.
        if fi.check_disjoint(pr_pt, land_poly):
            
            # Calculate the minimum distance from the observer point to land mask.
            shore_dist = int(min([pr_pt.distanceTo(i) for i in land_poly]))
            
            # Check that the distance to shore is less than the maximum allowed distance.
            if shore_dist <= max_dist:
                
                # Get visibility values.
                ts_vis_val, ms_vis_val, lk_vis_val = fv.radial_viewshed([x_val, y_val, z_val], in_dem_ras, pr_gdb, shore_dist, max_dist, out_crs, lmark_geom_list, land_poly, out_pts="", out_skyline="")
                
                # TEST: ASSIGN RANDOM INTEGERS SO THAT I CAN FINISH THE NETCDF EXPORT PART.
                #ts_vis_val, ms_vis_val, lk_vis_val = [random.randint(0,101) for i in range(0,3,1)]
                
                # Update arrays.
                ts_vis_array[z_idx, y_idx, x_idx] = ts_vis_val
                ms_vis_array[z_idx, y_idx, x_idx] = ms_vis_val
                lk_vis_array[z_idx, y_idx, x_idx] = lk_vis_val
        
        
        pr_count += 1
        
        console.prcnt_complete(pr_count, pr_total, 5, start_time, leading_spaces=2, leading_text="")
    #===========================================================================
    # GENERATE NETCDF FILE.
    #===========================================================================
    console.console("Writing NetCDF to file...")
    
    '''
    print(ts_vis_array)
    print(ms_vis_array)
    print(lk_vis_array)
    '''
    
    # Define dictionaries.
    file_dict = {'title':'maritime visibility parameters',
                 'source':os.path.basename(out_ncdf),
                 'coordinate_system':out_crs_name,
                 'date_created':datetime.datetime.now().strftime('%Y-%m-%d')}
    
    dim_dict = {"x":len(x_range),
                "y":len(y_range),
                "z":len(z_range)
                }
    
    var_dict = {"x":{'datatype':np.float64,
                     'dimensions':('x')},
                "y":{'datatype':np.float64,
                     'dimensions':('y')},
                "z":{'datatype':np.float64,
                     'dimensions':('z')},
                "ts":{'datatype':np.float64,
                      'dimensions':('z','x','y')},
                "ms":{'datatype':np.float64,
                      'dimensions':('z','x','y')},
                "lk":{'datatype':np.float64,
                      'dimensions':('z','x','y')}
             }
    
    var_atts_dict = {"x":{'long_name':'easting',
                          'units':out_crs_units},
                     "y":{'long_name':'northing',
                          'units':out_crs_units},
                     "z":{'long_name':'elevation',
                          'units':'{} above DEM zero point'.format(out_crs_units)},
                     "ts":{'long_name':'total_coastal_subtended_area',
                           'units':'square_degrees_of_arc'},
                     "ms":{'long_name':'maximum_coastal_subtended_area',
                           'units':'square_degrees_of_arc'},
                     "lk":{'long_name':'landmark_count',
                           'units':'count_of_landmarks_sighted'}
                     }
    
    var_arrays_dict = {"x":x_array,
                     "y":y_array,
                     "z":z_array,
                     "ts":ts_vis_array,
                     "ms":ms_vis_array,
                     "lk":lk_vis_array
                     }
    # Create NetCDF.
    #fn.create_netcdf(out_ncdf, file_dict, dim_dict, var_dict, var_atts_dict, var_arrays_dict)
    
    with nc.Dataset(out_ncdf, 'w', format="NETCDF4") as output:
        
        # Set global attributes for file.
        output.setncatts(file_dict)
        
        # Create dimensions.
        output.createDimension("x", len(x_range))
        output.createDimension("y", len(y_range))
        output.createDimension("z", len(z_range))
        
        # Create the variable.
        output.createVariable('easting', np.float64, ('y','x'))
        output.createVariable('northing', np.float64, ('y','x'))
        
        output.createVariable('ts', np.float64, ('z','y','x'))
        output.createVariable('ms', np.float64, ('z','y','x'))
        output.createVariable('lk', np.float64, ('z','y','x'))
        
        output['easting'][:] = x_array
        output['northing'][:] = y_array
        #output['z'][:] = z_array
        
        output['ts'][:] = ts_vis_array
        output['ms'][:] = ms_vis_array
        output['lk'][:] = lk_vis_array
    
    # Print message to console.
    console.console("Script finished.")
    
# (out_ncdf, pr_gdb, in_dem, lmark_fc, land_fc, max_dist, x_range, y_range, z_range, coordinate_system)
vis_ncdf(r"C:\GIS\CAA_2024\Test_Data\test_main_01.nc", 
         r"C:\GIS\CAA_2024\script_testing.gdb", 
         r"C:\GIS\CAA_2023\Data\Script_Datasets\DEM\SRTM_GEBCO_m.tif", 
         r"C:\GIS\CAA_2024\Test_Data\test_landmarks.shp", 
         r"C:\GIS\CAA_2024\Test_Data\test_coast_poly.shp", 
         10000,
         [3245000, 3261000, 1000], [3109000, 3126000, 1000], [0, 2, 1]
         )

#[3258000, 3261000, 1000], [3109000, 3111000, 1000], [0, 3, 1]