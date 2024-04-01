'''
Created on Mar 7, 2024

@author: Karl
'''


import numpy as np
import arcpy, os, random, datetime

# Import local functions.
from functions import clip, console, intersect, list, netcdf, viewshed
from functions.benchmark import benchmark
from functions.benchmark import print_benchmark_to_console

def vis_ncdf(out_ncdf, pr_gdb, in_dem, pt_mask_fc, land_fc, max_dist, xy_spacing, z_range, lmark_fc=""):
    
    # Print starting message to console.
    console.console("Starting script...")
    
    # Get start time.
    start_time = datetime.datetime.now()
    
    benchmark_dict = {}
    
    #===========================================================================
    # FIND COORDINATE SYSTEM.
    #===========================================================================
    
    if lmark_fc != "":
    
        dem_crs = arcpy.Describe(in_dem).spatialReference
        lmark_crs = arcpy.Describe(lmark_fc).spatialReference
        coast_crs = arcpy.Describe(land_fc).spatialReference
        
        if dem_crs.name == lmark_crs.name and dem_crs.name == coast_crs.name and lmark_crs.name == coast_crs.name: # Yes, it's clumsy.
            
            out_crs = dem_crs
            out_crs_name = out_crs.name
            out_crs_units = out_crs.linearUnitName
        
        else:
            
            raise Exception("Input DEM and feature classes have different coordinate systems.")
        
    else:
        
        dem_crs = arcpy.Describe(in_dem).spatialReference
        coast_crs = arcpy.Describe(land_fc).spatialReference
        
        if dem_crs.name == coast_crs.name:
            
            out_crs = dem_crs
            out_crs_name = out_crs.name
            out_crs_units = out_crs.linearUnitName
        
        else:
            
            raise Exception("Input DEM and feature classes have different coordinate systems.")
    
    #===========================================================================
    # GET X, Y RANGES FROM INPUT MASK SHAPE.
    #===========================================================================
    
    land_fc_extent = arcpy.Describe(pt_mask_fc).extent
    # Round west and south up to the nearest spacing value.
    land_fc_west = land_fc_extent.XMin + (xy_spacing - (land_fc_extent.XMin % xy_spacing))
    land_fc_south = land_fc_extent.YMin + (xy_spacing - (land_fc_extent.YMin % xy_spacing))
    
    # Round east and north down to the nearest spacing value.
    land_fc_east = land_fc_extent.XMax - (land_fc_extent.XMax % xy_spacing)
    land_fc_north = land_fc_extent.YMax - (land_fc_extent.YMax % xy_spacing)
    
    # Calculate ranges.
    x_range = [land_fc_west, land_fc_east, xy_spacing]
    y_range = [land_fc_south, land_fc_north, xy_spacing]
    
    del land_fc_extent
    
    #===========================================================================
    # HANDLE SOURCE DATA
    #===========================================================================
    
    # Get extent and resolution from DEM.
    dem_extent = [float(arcpy.management.GetRasterProperties(in_dem, "LEFT").getOutput(0)),
                  float(arcpy.management.GetRasterProperties(in_dem, "BOTTOM").getOutput(0)),
                  float(arcpy.management.GetRasterProperties(in_dem, "RIGHT").getOutput(0)),
                  float(arcpy.management.GetRasterProperties(in_dem, "TOP").getOutput(0))]
    
    dem_resolution = max([float(arcpy.management.GetRasterProperties(in_dem, "CELLSIZEX").getOutput(0)),
                          float(arcpy.management.GetRasterProperties(in_dem, "CELLSIZEY").getOutput(0))])
    
    # Read coast polylines to list of geometry objects.
    if land_fc == "":
        land_poly = None
        
    else:
        
        # Clip input land mask features to fix raster.
        clip_land_mask = clip.clip_polygon(dem_extent, land_fc, out_crs, clip_poly_path=r"memory\clip_land_mask")
        
        land_poly = []
        
        with arcpy.da.SearchCursor(clip_land_mask, ["SHAPE@"]) as cursor: #@UndefinedVariableFromImport
            for row in cursor:
                for part in row:
                    land_poly.append(part)
        
        # Delete search cursor.
        del cursor
        
        # Delete clipped features.
        arcpy.Delete_management(clip_land_mask)
    
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
    x_array = np.array([x_range for i in y_range]) #@UnusedVariable
    y_array = np.rot90(np.array([y_range for i in x_range])) #@UnusedVariable
    
    # Use ranges to create a shape for data outputs.
    out_array_shape = (len(z_range), len(y_range), len(x_range))
    
    # Create arrays for data outputs.
    sub_area_array = np.zeros(out_array_shape, dtype=float, order='C')
    sub_sum_x_array = np.zeros(out_array_shape, dtype=float, order='C')
    sub_max_y_array = np.zeros(out_array_shape, dtype=float, order='C')
    landmark_count_array = np.zeros(out_array_shape, dtype=float, order='C')
    
    # Create array for mask intersection.
    mask_int_array_shape = (len(y_range), len(x_range))
    mask_int_array = np.zeros(mask_int_array_shape, dtype=float, order='C')
    
    # Read landmark polygons to list of geometry objects.
    #@TODO: Scan for duplicate IDs and print warning?
    if lmark_fc == "":
        
        # Return empty list.
        lmark_geom_list = []
        
        # Print message to console.
        console.console("WARNING: No landmark dataset specified, all landmark values will be zero.", 2)
        
    else:
        lmark_geom_list = []
        with arcpy.da.SearchCursor(lmark_fc, ["SHAPE@"]) as cursor: #@UndefinedVariableFromImport
            for row in cursor:
                lmark_geom_list.append(row[0])
        
        del cursor
        
    #------------------------------------------------------------------------------ 
    
    # Create counter for points processed.
    pr_count = 0
    pr_total = int(len(x_range) * len(y_range) * len(z_range))
    
    # Print dimension and point count to console.
    console.console("X dimension has length: {}.".format(str(len(x_range))),2)
    console.console("Y dimension has length: {}.".format(str(len(y_range))),2)
    console.console("Z dimension has length: {}.".format(str(len(z_range))),2)
    console.console("3D array has {} data points.".format(str(pr_total)),2)
    
    #===========================================================================
    # CHECK POINTS AGAINST POINT MASK.
    #===========================================================================
    console.console("Checking points against point mask...")
    
    # Set time for benchmarking.
    point_mask_start_time = datetime.datetime.now() #@TODO: RE-MERGE THESE TWO, GET RID OF TIME.
    
    checkpoints_pr_count = 0
    checkpoints_pr_total = int(len(x_range) * len(y_range))
    
    mask_poly = []
    
    with arcpy.da.SearchCursor(pt_mask_fc, ["SHAPE@"]) as cursor: #@UndefinedVariableFromImport
            for row in cursor:
                for part in row:
                    mask_poly.append(part)
    
    del cursor
    
    # Iterate through indices.
    for index, x in np.ndenumerate(mask_int_array): #@UnusedVariable
        
        # Get indices from iterator.
        y_idx, x_idx = index
        
        # Get x, y values for iterated data point.
        x_val = float(x_range[x_idx])
        y_val = float(y_range[y_idx])
        
        pr_pt = arcpy.PointGeometry(arcpy.Point(x_val, y_val))
        
        if intersect.check_disjoint(pr_pt, mask_poly) == False:
            
            mask_int_array[y_idx, x_idx] = 1
            
        checkpoints_pr_count += 1
        
        console.prcnt_complete(checkpoints_pr_count, checkpoints_pr_total, 5, point_mask_start_time, leading_spaces=2, leading_text="")
    
    del mask_poly
    
    # Add to benchmark dictionary.
    benchmark_dict = benchmark(point_mask_start_time, benchmark_dict, "chek_pts_against_pt_mask")
    
    pr_total = int(np.sum(mask_int_array))
    
    console.console("2D array has {} data points that intersect land mask.".format(str(int(pr_total))),2)
    
    #===========================================================================
    # CALCULATE VISIBILITY VALUES.
    #===========================================================================
    console.console("Calculating visibility values...")
    
    iter_vis_start_time = datetime.datetime.now()
    
    # Iterate through indices in mask array.
    for index, x in np.ndenumerate(mask_int_array): #@UnusedVariable
        
        # Get indices from iterator.
        y_idx, x_idx = index
        
        # Get x, y values for iterated data point.
        x_val = float(x_range[x_idx])
        y_val = float(y_range[y_idx])
        
        #print("Processing point at {} / {} / {}...".format(x_val, y_val, z_val))
        if mask_int_array[y_idx, x_idx] == 1:
        
            # Set time for benchmarking.
            function_start_time = datetime.datetime.now()
            
            # Check if point is within mask.
            pr_pt = arcpy.PointGeometry(arcpy.Point(x_val, y_val))
            
            # Check the distance from the point to polygons in the input mask (returns empty list if point intersects mask).
            pr_pt_dist_list = intersect.check_distance_to(pr_pt, land_poly)
            
            # Add to benchmark dictionary.
            benchmark_dict = benchmark(function_start_time, benchmark_dict, "point_land_mask_distanceTo")
            
            # Only collect values if point does not intersect the mask.
            if len(pr_pt_dist_list) > 0:
                
                # Calculate the minimum distance from the observer point to land mask.
                shore_dist = int(min(pr_pt_dist_list))
                
                # Check that the distance to shore is less than the maximum allowed distance.
                if shore_dist <= max_dist:
                    
                    # Get visibility values.
                    vis_dict, benchmark_dict = viewshed.radial_viewshed(x_val, y_val, z_range, in_dem, dem_resolution, pr_gdb, shore_dist, max_dist, out_crs, lmark_geom_list, 
                                                                          land_poly, out_pts="", out_skyline="", benchmark_dict=benchmark_dict)
                    
                    # TEST: ASSIGN RANDOM INTEGERS SO THAT I CAN FINISH THE NETCDF EXPORT PART.
                    #ts_vis_val, ms_vis_val, lk_vis_val = [random.randint(0,101) for i in range(0,3,1)]
                    
                    for z_val in vis_dict.keys():
                        
                        sub_area, sub_max_y, sub_sum_x, lmark_int = vis_dict[z_val]
                        
                        z_idx = list.list_index_from_val(z_range, z_val)
                        
                        # Update arrays.
                        sub_area_array[z_idx, y_idx, x_idx] = sub_area
                        sub_sum_x_array[z_idx, y_idx, x_idx] = sub_sum_x
                        sub_max_y_array[z_idx, y_idx, x_idx] = sub_max_y
                        landmark_count_array[z_idx, y_idx, x_idx] = lmark_int
            
            
            pr_count += 1
            
            console.prcnt_complete(pr_count, pr_total, 1, iter_vis_start_time, leading_spaces=2, leading_text="")
    #===========================================================================
    # GENERATE NETCDF FILE.
    #===========================================================================
    console.console("Writing NetCDF to file...")
    
    # Define dictionaries.
    file_dict = {'title':'maritime_visibility_parameters',
                 'output_filename':os.path.basename(out_ncdf),
                 'coordinate_system':out_crs_name,
                 'date_created':datetime.datetime.now().strftime('%Y-%m-%d')}
    
    dim_dict = {"x":len(x_range),
                "y":len(y_range),
                "z":None} # Make z an unlimited dimension (so more data can be added later).
    
    var_dict = {"easting":{'datatype':np.float64,
                           'dimensions':('y','x')},
                "northing":{'datatype':np.float64,
                            'dimensions':('y','x')},
                "sub_area":{'datatype':np.float64,
                            'dimensions':('z','y','x')},
                "sub_sum_x":{'datatype':np.float64,
                             'dimensions':('z','y','x')},
                "sub_max_y":{'datatype':np.float64,
                             'dimensions':('z','y','x')},
                "lmrk_count":{'datatype':np.float64,
                              'dimensions':('z','y','x')}}
    
    var_atts_dict = {"easting":{'long_name':'easting',
                                'units':out_crs_units},
                     "northing":{'long_name':'northing',
                                 'units':out_crs_units},
                     "sub_area":{'long_name':'total_subtended_area',
                                 'units':'square_degrees_of_arc'},
                     "sub_sum_x":{'long_name':'sum_of_horizontal_subtended_area',
                                  'units':'degrees_of_arc'},
                     "sub_max_y":{'long_name':'maximum_vertical_subtended_area',
                                  'units':'degrees_of_arc'},
                     "lmrk_count":{'long_name':'landmark_count',
                                   'units':'number_of_landmarks'}}
    
    var_arrays_dict = {"easting":x_array,
                       "northing":y_array,
                       "sub_area":sub_area_array,
                       "sub_sum_x":sub_sum_x_array,
                       "sub_max_y":sub_max_y_array,
                       "lmrk_count":landmark_count_array}
    
    # Create NetCDF.
    netcdf.create_netcdf(out_ncdf, file_dict, dim_dict, var_dict, var_atts_dict, var_arrays_dict)
    '''
    #===========================================================================
    # GENERATE RASTER.
    #===========================================================================
    console.console("Writing Raster(s) to file...")
    
    # Set environment crs.
    arcpy.env.outputCoordinateSystem = out_crs
    arcpy.env.overwriteOutput = True
    
    raster_path = r"{}\subarea_{}".format(pr_gdb, os.path.basename(out_ncdf).split(".")[0])
    mod_sub_array = np.flip(sub_area_array[0],0)
    
    out_ras = arcpy.NumPyArrayToRaster(mod_sub_array,
                                 lower_left_corner=arcpy.Point((land_fc_west-(xy_spacing*0.5)), (land_fc_south-(xy_spacing*0.5))), 
                                 x_cell_size=xy_spacing, 
                                 y_cell_size=xy_spacing,
                                 value_to_nodata=None)
    
    out_ras.save(raster_path)
    '''
    # Print message to console.
    console.console("Script finished.")
    
    # Print benchmarking messages.
    print()
    print_benchmark_to_console(start_time, benchmark_dict)
    