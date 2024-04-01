'''
Created on Dec 13, 2023

@author: Karl
'''

import arcpy, datetime
from functions import array, clip, heading
import output
from functions.benchmark import benchmark

def radial_viewshed(obs_x, obs_y, obs_z_list, in_dem_ras, in_dem_res, pr_gdb, min_dist, max_dist, pt_crs, lmark_geom_list="", land_fc="", land_poly_list="", out_pts="", out_skyline="", benchmark_dict=""):
    
    #mod_obs_pt = [obs_x, obs_y, obs_z]
    
    #mod_obs_pt_xy = arcpy.PointGeometry(arcpy.Point(obs_x, obs_y))
    
    # Set start time.
    start_time = datetime.datetime.now()
    
    # Read check_features to polygons if they aren't already.
    
    # Check whether features are already polygon geometries.
    '''
    # Create variables for output.
    array_sub_area = 0
    array_max_v_angle = 0
    array_sum_h_angle = 0
    array_lmark_int = 0
    '''
    # Create radial array.
    
    # Generate 2d radial array.
    function_start_time = datetime.datetime.now()
    # Generate 2D radial array.
    check_coordinates_list, heading_polylines_dict = array.generate_2d_radial_array(obs_x, obs_y, min_dist, max_dist, pr_gdb, degree_interval=1)
    benchmark_dict = benchmark(function_start_time, benchmark_dict, "generate_2d_radial_array")
    
    # Clip raster to array extent.
    function_start_time = datetime.datetime.now()
    # Calculate clipping extent.
    clip_extent_list = clip.clip_extent(check_coordinates_list, (in_dem_res*2))
    # Clip the raster.
    clip_ras = clip.clip_raster(clip_extent_list, in_dem_ras)
    benchmark_dict = benchmark(function_start_time, benchmark_dict, "clip_raster")
    
    # Clip land polygons to array extent.
    function_start_time = datetime.datetime.now()
    reclip_land_mask = clip.clip_polygon(clip_extent_list, land_fc, pt_crs, clip_poly_path=r"memory\reclip_land_mask")
    
    reclip_land_poly_list = []
    
    with arcpy.da.SearchCursor(reclip_land_mask, ["SHAPE@"]) as cursor: #@UndefinedVariableFromImport
        for row in cursor:
            for part in row:
                reclip_land_poly_list.append(part)
    
    # Delete search cursor.
    del cursor
    benchmark_dict = benchmark(function_start_time, benchmark_dict, "reclip_land_mask")
    
    # Delete clipped features.
    arcpy.Delete_management(reclip_land_mask)
    
    # Filter rays in radial array.
    function_start_time = datetime.datetime.now()
    insert_cursor_list = array.filter_2d_radial_array(heading_polylines_dict, reclip_land_poly_list)
    benchmark_dict = benchmark(function_start_time, benchmark_dict, "filter_2d_radial_array")
    
    # Interpolate 2d array to 3d array, and return list of points.
    function_start_time = datetime.datetime.now()
    radial_array_dict = array.interpolate_2d_radial_array(obs_x, obs_y, obs_z_list, insert_cursor_list, clip_ras, pr_gdb, pt_crs)
    benchmark_dict = benchmark(function_start_time, benchmark_dict, "interpolate_2d_radial_array")
    
    # Delete clip raster.
    arcpy.Delete_management(clip_ras)
    
    # Create list of skyline points.
    skyline_dict = {}
    
    # Create list of sub angle points.
    out_pt_dict = {}
    
    pr_count = 0
    
    # Create dictionary for output.
    out_dict = {}
    
    # Loop through radial arrays in radial array dictionary.
    for obs_z in radial_array_dict.keys():
        
        radial_array = radial_array_dict[obs_z]
        
        obs_pt = [float(i) for i in [obs_x, obs_y, obs_z]]
        
        # Create variables for output.
        array_sub_area = 0
        array_max_v_angle = 0
        array_sum_h_angle = 0
        array_lmark_int = 0
    
        # Loop through array.
        for heading_dict in radial_array:
            
            if len(heading_dict.keys()) > 0:
                
                # If heading intersects the land, increment sum of horizontal angles.
                array_sum_h_angle += 1
                
                function_start_time = datetime.datetime.now()
                
                # For heading, create subtended angle dict.
                heading_sub_dict = heading.heading_sub_angle(obs_pt, heading_dict)
                
                benchmark_dict = benchmark(function_start_time, benchmark_dict, "heading_sub_angle")
                
                # Calculate maximum subtended angle.
                heading_max_sub_angle = max([heading_sub_dict[key][-1] for key in heading_sub_dict.keys()])
                
                # Update max sub angle.
                if heading_max_sub_angle > array_max_v_angle:
                    array_max_v_angle = heading_max_sub_angle
                
                # Update sub angle value by adding maxiumum subtended angle for points in heading.
                array_sub_area += heading_max_sub_angle
                
                function_start_time = datetime.datetime.now()
                
                # For heading, create vis dict.
                heading_vis_dict = heading.heading_vis_pts(heading_sub_dict, obs_pt[2])
                
                benchmark_dict = benchmark(function_start_time, benchmark_dict, "heading_vis_pts")
                
                # Populate heading_vis_dict.
                for key in heading_dict.keys():
                    
                    out_pt_dict[key] = heading_vis_dict[key]
                
                # Populate vertex list if output skyline is being generated.
                if out_skyline != "":
                    
                    # Create dict for export.
                    slct_list = []
                    
                    # Loop through heading_sub_angle_list to find skyline pt.
                    for key in heading_vis_dict.keys():
                        
                        if heading_vis_dict[key][-1] == 1:
                            
                            slct_list.append(key)
                            
                    if len(slct_list) > 0:
                        
                        skyline_dict[slct_list[-1]] = heading_vis_dict[slct_list[-1]]
                        
                    else:
                        
                        skyline_dict[list(heading_vis_dict.keys())[-1]] = heading_vis_dict[list(heading_vis_dict.keys())[-1]]
                
                pr_count += 1
                
                #console.prcnt_complete(pr_count, len(radial_array), 5, start_time, leading_spaces=0, leading_text="Headings ")
              
        # If there are input landmarks features, loop through them and check for intersecting visible points.
        if len(lmark_geom_list) > 0:
            
            # Select from geometry list where geometry is within distance of observer point.
            lmark_mod_geom_list = [i for i in lmark_geom_list if i.distanceTo(obs_pt) <= max_dist]
            
            # If there are any landmark geometries left, loop through them.
            if len(lmark_mod_geom_list) > 0:
                
                function_start_time = datetime.datetime.now()
                
                # Create multipoint geometry consisting of visible points within the radial array.
                vis_pts_list = [out_pt_dict[i] for i in out_pt_dict.keys()]
                
                # Check that visible points exist.
                if len(vis_pts_list) > 0:
                    
                    # Create multipoint from visible points.
                    vis_pts_multipoint = arcpy.Multipoint(arcpy.Array([arcpy.Point(*coords) for coords in vis_pts_list]))
                    
                    # Iterate through landmarks features.
                    for lmark_geom in lmark_mod_geom_list:
                        
                        # Increment visible landmarks count if landmark and visible point geometries are not disjoint (i.e. intersect).
                        if lmark_geom.disjoint(vis_pts_multipoint) == False:
                            
                            # Add to landmark count.
                            array_lmark_int += 1
                
                benchmark_dict = benchmark(function_start_time, benchmark_dict, "check_landmark_disjoint")
                
        # Export features (if enabled).
        if out_pts != "":
            output.output_pts(out_pt_dict, out_pts, pt_crs)
        
        if out_skyline != "":
            
            # Add first skyline entry to the end of the dictionary to complete the feature.
            skyline_dict["{}-1".format(list(skyline_dict.keys())[0])] = skyline_dict[list(skyline_dict.keys())[0]]
            
            output.output_skyline(skyline_dict, out_skyline, pt_crs)
        
        # Return: max subtended angle, total subtended angle, visible landmarks.
        iter_row = [round(i, 4) for i in [array_sub_area, array_max_v_angle, array_sum_h_angle, array_lmark_int]]
        
        out_dict[obs_z] = iter_row
        
    return out_dict, benchmark_dict