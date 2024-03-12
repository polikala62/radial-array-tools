'''
Created on Dec 13, 2023

@author: Karl
'''

import arcpy, datetime, time
from functions import array, clip, heading
import output
from functions.benchmark import benchmark

def radial_viewshed(obs_pt, in_dem_ras, in_dem_res, pr_gdb, min_dist, max_dist, pt_crs, lmark_geom_list="", land_poly="", out_pts="", out_skyline="", benchmark_dict=""):
    
    # Set start time.
    start_time = datetime.datetime.now()
    
    # Get coordinates from input, check that they are within DEM extent.
    obs_x = float(obs_pt[0])
    obs_y = float(obs_pt[1])
    obs_z = float(obs_pt[2])
    
    mod_obs_pt = [obs_x, obs_y, obs_z]
    
    mod_obs_pt_xy = arcpy.PointGeometry(arcpy.Point(obs_x, obs_y))
    
    # Read check_features to polygons if they aren't already.
    
    # Check whether features are already polygon geometries.
    
    # Create variables for output.
    array_sub_area = 0
    array_max_v_angle = 0
    array_sum_h_angle = 0
    array_lmark_int = 0
    
    # Create radial array.
    
    # Generate 2d radial array.
    function_start_time = time.time()
    check_coordinates_list, insert_cursor_list = array.generate_2d_radial_array(mod_obs_pt, min_dist, max_dist, pr_gdb, coast_polylines=land_poly, degree_interval=1)
    benchmark_dict = benchmark(function_start_time, benchmark_dict, "generate_2d_radial_array")
    
    # Clip raster to array extent.
    function_start_time = time.time()
    # Calculate clipping extent.
    clip_extent_list = clip.clip_extent(check_coordinates_list, (in_dem_res*2))
    clip_ras = clip.clip_raster(clip_extent_list, in_dem_ras)
    benchmark_dict = benchmark(function_start_time, benchmark_dict, "clip_datasets")
                
    # Interpolate 2d array to 3d array, and return list of points.
    function_start_time = time.time()
    radial_array = array.interpolate_2d_radial_array(obs_pt, insert_cursor_list, clip_ras, pr_gdb, pt_crs)
    benchmark_dict = benchmark(function_start_time, benchmark_dict, "interpolate_2d_radial_array")
    
    # Delete clip raster.
    arcpy.Delete_management(clip_ras)
    
    #radial_array = funcs_array.generate_radial_array(obs_pt, in_dem_ras, max_dist, raster_res=raster_res, degree_interval=2)
    #radial_array = array.interpolate3d_radial_array(mod_obs_pt, in_dem_ras, min_dist, max_dist, pr_gdb, pt_crs, coast_polylines=land_poly, raster_res=5, degree_interval=1)
    
    # Create list of skyline points.
    skyline_dict = {}
    
    # Create list of sub angle points.
    out_pt_dict = {}
    
    pr_count = 0
    
    # Loop through array.
    for heading_dict in radial_array:
        
        if len(heading_dict.keys()) > 0:
            
            # If heading intersects the land, increment sum of horizontal angles.
            array_sum_h_angle += 1
            
            function_start_time = time.time()
            
            # For heading, create subtended angle dict.
            heading_sub_dict = heading.heading_sub_angle(mod_obs_pt, heading_dict)
            
            benchmark_dict = benchmark(function_start_time, benchmark_dict, "heading_sub_angle")
            
            # Calculate maximum subtended angle.
            heading_max_sub_angle = max([heading_sub_dict[key][-1] for key in heading_sub_dict.keys()])
            
            # Update max sub angle.
            if heading_max_sub_angle > array_max_v_angle:
                array_max_v_angle = heading_max_sub_angle
            
            # Update sub angle value by adding maxiumum subtended angle for points in heading.
            array_sub_area += heading_max_sub_angle
            
            function_start_time = time.time()
            
            # For heading, create vis dict.
            heading_vis_dict = heading.heading_vis_pts(heading_sub_dict, mod_obs_pt[2])
            
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
    if lmark_geom_list != "":
        
        # Select from geometry list where geometry is within distance of observer point.
        lmark_mod_geom_list = [i for i in lmark_geom_list if i.distanceTo(mod_obs_pt_xy) <= max_dist]
        
        # If there are any landmark geometries left, loop through them.
        if len(lmark_mod_geom_list) > 0:
            
            function_start_time = time.time()
            
            # Create multipoint geometry consisting of visible points within the radial array.
            vis_pts_list = [out_pt_dict[i] for i in out_pt_dict.keys()]
            vis_pts_multipoint = arcpy.Multipoint(arcpy.Array([arcpy.Point(*coords) for coords in vis_pts_list]))
            
            # Iterate through landmarks features.
            for lmark_geom in lmark_mod_geom_list:
                
                #@TODO: CODE SPATIAL INTERESECT (IN SEPARATE FUNCTION?)
                if lmark_geom.disjoint(vis_pts_multipoint) == False:
                    
                    # Increment visible landmarks count if landmark and visible point geometries intersect.
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
    return [array_sub_area, array_max_v_angle, array_sum_h_angle, array_lmark_int, benchmark_dict]