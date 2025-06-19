'''
Created on Mar 7, 2024

@author: Karl
'''

# Import libraries.
import numpy as np
import arcpy, json, os, datetime, tqdm

# Disable log history.
if arcpy.GetLogHistory():
    arcpy.SetLogHistory(False)

# Import local functions.
from functions import clip, console, intersect, list, netcdf, validate, viewshed
from functions.benchmark import benchmark
from functions.benchmark import print_benchmark_to_console
from functions import json as pr_json

# Define function for main loop.
def vis_ncdf(out_ncdf, pr_gdb, in_dem, pt_mask_fc, land_fc, xy_spacing, dist_range, densify_dist, z_range, obs_z_offset=1.5, 
             sample_raster="", lmark_fc="", out_dist_list="", pt_mask_json='', write_log=False):
    
    # Print starting message to console.
    console.console("Starting script...")
    
    #===========================================================================
    # CONFIRM OVERWRITE
    #===========================================================================
    
    # Never again.
    if os.path.exists(out_ncdf):
        
        console.console("Output '{}' already exists. Enter 'y' to confirm overwrite. Press 'enter' to exit.".format(os.path.basename(out_ncdf)))
        overwrite_response = input()
        
        if overwrite_response.lower() != 'y':
            
            raise SystemExit
    
    #===========================================================================
    # SET INITIALISATION VARIABLES
    #===========================================================================
    
    # Get start time.
    start_time = datetime.datetime.now()
    
    # Create dictionary to hold benchmarking values.
    benchmark_dict = {}
    
    #===========================================================================
    # FIND COORDINATE SYSTEM.
    #===========================================================================
    
    # Get coordinate system by comparing input features (raises exception if CRSsses are different).
    out_crs = validate.validate_crs([in_dem, pt_mask_fc, land_fc, lmark_fc])
    
    # Get coordinate system name as variable.
    out_crs_name = out_crs.name
    
    # Get coordinate system linear unit as variable.
    out_crs_units = out_crs.linearUnitName
    
    #===========================================================================
    # GET X, Y RANGES FROM INPUT MASK SHAPE.
    #===========================================================================
    
    # Create an extent object for the land mask.
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
    
    # Delete describe object.
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
        
        raise Exception("Could not read input ranges. Ranges for x, y, z should be lists with format [minimum_value, maximum_value, increment].")
    
    # Create ranges for x, y, z.
    x_range = np.arange(x_min, (x_max + x_step), x_step)
    y_range = np.arange(y_min, (y_max + y_step), y_step)
    z_range = np.arange(z_min, (z_max + z_step), z_step)
    
    #------------------------------------------------------------------------------ 
    
    if out_dist_list == "":
        
        d_min, d_max, d_step = dist_range
        d_range = np.arange(d_min, (d_max + d_step), d_step)
        
    else:
        
        d_min = min(out_dist_list)
        d_max = max(out_dist_list)
        d_range = np.array(out_dist_list)
    
    # Create 2D arrays for x, y, z, and distance.
    x_array = np.array([x_range for i in y_range]) #@UnusedVariable
    y_array = np.rot90(np.array([y_range for i in x_range])) #@UnusedVariable
    
    # Use ranges to create a shape for data outputs.
    out_array_shape = (len(d_range), len(z_range), len(y_range), len(x_range))
    
    # Create default arrays for data outputs.
    sub_area_array = np.zeros(out_array_shape, dtype=float, order='C')
    sub_sum_x_array = np.zeros(out_array_shape, dtype=float, order='C')
    sub_max_y_array = np.zeros(out_array_shape, dtype=float, order='C')
    
    # Create arrays for sample values, if enabled.
    if sample_raster != "":
        sample_min_array = np.zeros(out_array_shape, dtype=float, order='C')
        sample_max_array = np.zeros(out_array_shape, dtype=float, order='C')
        sample_avg_array = np.zeros(out_array_shape, dtype=float, order='C')
    
    # Create arrays for landmark values, if enabled.   
    if lmark_fc != "":
        landmark_count_array = np.zeros(out_array_shape, dtype=float, order='C')
    
    # Create array for mask intersection.
    mask_int_array_shape = (len(y_range), len(x_range))
    mask_int_array = np.zeros(mask_int_array_shape, dtype=float, order='C')
    
    
    # Read landmark polygons to list of geometry objects.
    #@TODO: Scan for duplicate IDs and print warning?
    if lmark_fc == "":
        
        # Return empty string.
        lmark_geom_list = ""
        
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
    console.console("D dimension has length: {}.".format(str(len(d_range))),2)
    console.console("4D array has {} data points.".format(str(pr_total)),2)
    
    #===========================================================================
    # CHECK POINTS AGAINST POINT MASK.
    #===========================================================================
    
    if os.path.exists(pt_mask_json):
        
        console.console("Importing checked points dictionary from '{}'...".format(os.path.basename(pt_mask_json)))
        
        mask_int_array = pr_json.load_json(pt_mask_json)
        
    else:
    
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
        
        # Generate list of mask indices.
        mask_int_indices = [i for i, j in np.ndenumerate(mask_int_array)] #@UnusedVariable
        
        # Iterate through indices.
        for index in tqdm.tqdm(mask_int_indices, disable=True):
            
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
        
        console.console("2D array has {} data points that intersect point mask.".format(str(int(pr_total))),2)
        
        # Write points to json, if enabled.
        if pt_mask_json != '':
            
            console.console("Writing checked points dictionary to '{}'...".format(os.path.basename(pt_mask_json)))
            
            pr_json.write_json(mask_int_array, pt_mask_json)
    
    #===========================================================================
    # CALCULATE VISIBILITY VALUES.
    #===========================================================================
    
    console.console("Calculating visibility values...")
    
    iter_vis_start_time = datetime.datetime.now()
    
    # (Re)generate list of mask indices.
    mask_int_indices = [i for i, j in np.ndenumerate(mask_int_array)] #@UnusedVariable
    
    # Iterate through indices in mask array.
    for index in tqdm.tqdm(mask_int_indices):
        
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
            
            
            if len(pr_pt_dist_list) > 0:
                pt_pt_min_land_dist = min(pr_pt_dist_list)
            else:
                pt_pt_min_land_dist = 0
            
            # Add to benchmark dictionary.
            benchmark_dict = benchmark(function_start_time, benchmark_dict, "point_land_mask_distanceTo")
            
            # Only collect values if point does not intersect the mask.
            if len(pr_pt_dist_list) > 0:
                    
                # Get visibility values.
                obs_dict, benchmark_dict = viewshed.radial_viewshed(x_val, y_val, z_range, d_range, in_dem, dem_resolution, pr_gdb, out_crs, densify_dist, obs_z_offset, lmark_geom_list, 
                                                                      land_poly, override_min_dist=pt_pt_min_land_dist, sample_ras=sample_raster, benchmark_dict=benchmark_dict)
                
                # Obs_dict is dictionary where key is observer height and value is distance dictionary.
                for obs_z in obs_dict.keys():
                    
                    # Find index for z value.
                    z_idx = list.list_index_from_val(z_range, obs_z)
                    
                    # Dist_dict is dictionary where key is visible distance and value is tuple.
                    vdist_dict = obs_dict[obs_z]
                    
                    # Tuple is format [v_angle_sum, v_angle_range, h_angle_sum, *sample_max, *sample_min, *sample_avg, *landmark_sum]
                    for dist_val in vdist_dict.keys():
                        
                        # Find index for distances.
                        d_idx = list.list_index_from_val(d_range, dist_val)
                        
                        v_angle_sum, v_angle_max, h_angle_sum = vdist_dict[dist_val][0:3]
                        #print(obs_z, d_idx, v_angle_sum, v_angle_max, h_angle_sum)
                        # Update default arrays.
                        sub_area_array[d_idx, z_idx, y_idx, x_idx] = v_angle_sum
                        sub_max_y_array[d_idx, z_idx, y_idx, x_idx] = v_angle_max
                        sub_sum_x_array[d_idx, z_idx, y_idx, x_idx] = h_angle_sum
                        
                        if sample_raster != "":
                            
                            # Retrieve values from viewshed output.
                            sample_max, sample_min, sample_avg = vdist_dict[dist_val][3:6]
                            
                            # Update arrays.
                            sample_min_array[d_idx, z_idx, y_idx, x_idx] = sample_min
                            sample_max_array[d_idx, z_idx, y_idx, x_idx] = sample_max
                            sample_avg_array[d_idx, z_idx, y_idx, x_idx] = sample_avg
                            
                        if lmark_fc != "":
                            
                            # Retrieve values from viewshed output.
                            landmark_sum = vdist_dict[dist_val][-1]
                            
                            # Update arrays.
                            landmark_count_array[d_idx, z_idx, y_idx, x_idx] = landmark_sum            
            
            pr_count += 1
            
            #console.prcnt_complete(pr_count, pr_total, 1, iter_vis_start_time, leading_spaces=2, leading_text="")
    
    #===========================================================================
    # GENERATE NETCDF FILE.
    #===========================================================================
    
    console.console('')
    console.console("Writing NetCDF to file...")
    
    # Define dictionaries.
    file_dict = {'title':'maritime_visibility_parameters',
                 'output_filename':os.path.basename(out_ncdf),
                 'coordinate_system':out_crs_name,
                 'date_created':datetime.datetime.now().strftime('%Y-%m-%d')}
    
    dim_dict = {"x":len(x_range),
                "y":len(y_range),
                "z":len(z_range),
                "d":len(d_range)}
    
    #------------------------------------------------------------------------------ 
    
    # Create variable dictionary for default variables.
    var_dict = {"x":{'datatype':np.float64, 'dimensions':('x')},
                "y":{'datatype':np.float64, 'dimensions':('y')},
                "easting":{'datatype':np.float64, 'dimensions':('y','x')},
                "northing":{'datatype':np.float64, 'dimensions':('y','x')},
                "z":{'datatype':np.float64, 'dimensions':('z')},
                "d":{'datatype':np.float64, 'dimensions':('d')},
                "sub_area":{'datatype':np.float64, 'dimensions':('d','z','y','x')},
                "sub_sum_x":{'datatype':np.float64, 'dimensions':('d','z','y','x')},
                "sub_max_y":{'datatype':np.float64, 'dimensions':('d','z','y','x')}}
    
    # Create variable dictionary for sample data, if enabled.
    if sample_raster != "":
        var_dict['sample_max'] = {'datatype':np.float64, 'dimensions':('d','z','y','x')}
        var_dict['sample_min'] = {'datatype':np.float64, 'dimensions':('d','z','y','x')}
        var_dict['sample_avg'] = {'datatype':np.float64, 'dimensions':('d','z','y','x')}
    
    # Create variable dictionary for landmark data, if enabled.  
    if lmark_fc != "":
        var_dict["lmrk_count"] = {'datatype':np.float64, 'dimensions':('d','z','y','x')}
    
    #------------------------------------------------------------------------------ 
    
    # Create attributes dictionary for default variables.
    var_atts_dict = {"x":{'long_name':'x', 'units':out_crs_units, 'coordinates':'easting northing', 'grid_mapping':'spatial_ref'},
                     "y":{'long_name':'y', 'units':out_crs_units, 'coordinates':'easting northing', 'grid_mapping':'spatial_ref'},
                     "easting":{'long_name':'easting', 'units':out_crs_units, 'coordinates':'easting northing', 'grid_mapping':'spatial_ref'},
                     "northing":{'long_name':'northing', 'units':out_crs_units, 'coordinates':'easting northing', 'grid_mapping':'spatial_ref'},
                     "z":{'long_name':'observer_height', 'units':'{}_relative_to_vertical_datum'.format(out_crs_units), 'coordinates':'easting northing', 'grid_mapping':'spatial_ref'},
                     "d":{'long_name':'visible_distance_limit', 'units':out_crs_units, 'coordinates':'easting northing', 'grid_mapping':'spatial_ref'},
                     "sub_area":{'long_name':'total_subtended_area', 'units':'square_degrees_of_arc', 'coordinates':'easting northing', 'grid_mapping':'spatial_ref'},
                     "sub_sum_x":{'long_name':'sum_of_horizontal_subtended_area', 'units':'degrees_of_arc', 'coordinates':'easting northing', 'grid_mapping':'spatial_ref'},
                     "sub_max_y":{'long_name':'maximum_vertical_subtended_area', 'units':'degrees_of_arc', 'coordinates':'easting northing', 'grid_mapping':'spatial_ref'}}
    
    # Create variable dictionary for sample data, if enabled.
    if sample_raster != "":
        var_atts_dict['sample_max'] = {'long_name':'sample_raster_max_value', 'units':'cell_value'}
        var_atts_dict['sample_min'] = {'long_name':'sample_raster_min_value', 'units':'cell_value'}
        var_atts_dict['sample_avg'] = {'long_name':'sample_raster_avg_value', 'units':'cell_value'}
    
    # Create variable dictionary for landmark data, if enabled.
    if lmark_fc != "":
        var_atts_dict["lmrk_count"] = {'long_name':'landmark_count', 'units':'number_of_landmarks'}
    
    #------------------------------------------------------------------------------ 
    
    # Create array dictionary for default variables.
    var_arrays_dict = {"x":x_range,
                       "y":y_range,
                       "easting":x_array,
                       "northing":y_array,
                       "z":z_range,
                       "d":d_range,
                       "sub_area":sub_area_array,
                       "sub_sum_x":sub_sum_x_array,
                       "sub_max_y":sub_max_y_array}

    # Create array dictionary for sample data, if enabled.
    if sample_raster != "":
        var_arrays_dict['sample_max'] = sample_max_array
        var_arrays_dict['sample_min'] = sample_min_array
        var_arrays_dict['sample_avg'] = sample_avg_array
    
    # Create array dictionary for landmark data, if enabled.
    if lmark_fc != "":
        var_arrays_dict["lmrk_count"] = landmark_count_array
        
    #------------------------------------------------------------------------------ 
    
    # Create NetCDF.
    netcdf.create_netcdf(out_ncdf, out_crs, file_dict, dim_dict, var_dict, var_atts_dict, var_arrays_dict)
    
    console.console("NetCDF written to '{}'.".format(out_ncdf))
    
    #===========================================================================
    # WRITE BENCHMARKING & PROCESSING LOGS
    #===========================================================================
    
    param_log_dict = {'out_ncdf':str(out_ncdf),'pr_gdb':str(pr_gdb),'in_dem':str(in_dem),'pt_mask_fc':str(pt_mask_fc),'land_fc':str(land_fc),
                      'xy_spacing':str(xy_spacing),'dist_range':str(dist_range),'densify_dist':str(densify_dist),'z_range':str(z_range),'obs_z_offset':str(obs_z_offset),
                      'sample_raster':str(sample_raster),'lmark_fc':str(lmark_fc),'out_dist_list':str(out_dist_list)}
    
    # Print benchmarking messages.
    console.console('')
    benchmark_log_dict = print_benchmark_to_console(start_time, benchmark_dict)
    console.console('')
    
    # Write log, if enabled.
    if write_log:
        
        log_fp = os.path.join(os.path.dirname(out_ncdf), "{}_log.json".format(os.path.splitext(os.path.basename(out_ncdf))[0]))
        with open(log_fp, 'w') as fp:
            json.dump({'INPUT_PARAMETERS':param_log_dict, 'SCRIPT_TIME': str(datetime.datetime.now()-start_time),'BENCHMARK_TABLE':benchmark_log_dict}, fp, indent=4)
            
        console.console("Log file written to '{}'.".format(log_fp))
        console.console('')
        
    # Print message to console.
    console.console("Script finished!")