'''
Created on Dec 13, 2023

@author: Karl
'''

import arcpy, datetime
from functions import array, clip, heading, ray, validate
import output
from functions.benchmark import benchmark

def radial_viewshed(obs_x, obs_y, obs_z_list, d_range, in_dem_ras, in_dem_res, pr_gdb, pt_crs, 
                    lmark_geom_list="", land_fc="", override_min_dist="", sample_ras="", benchmark_dict=""):
    
    # Get distance min, max, and increment from variable.
    min_dist, max_dist, dist_inc = d_range
    
    # Set minimum distance to land distance, if override value is supplied.
    if override_min_dist != "" and validate.floatable(override_min_dist):
        min_dist = override_min_dist
        
        #TODO: ADD BETTER VALIDATION.
    
    #------------------------------------------------------------------------------ 
    
    function_start_time = datetime.datetime.now()
    
    # Generate 2D radial array.
    check_coordinates_list, heading_polylines_dict = array.generate_2d_radial_array(obs_x, obs_y, min_dist, max_dist, pr_gdb, degree_interval=1)
    
    benchmark_dict = benchmark(function_start_time, benchmark_dict, "array.generate_2d_radial_array")
    
    #------------------------------------------------------------------------------ 
    
    function_start_time = datetime.datetime.now()
    
    # Calculate clipping extent.
    clip_extent_list = clip.clip_extent(check_coordinates_list, (in_dem_res*2))
    
    # Clip raster to array extent.
    clip_ras = clip.clip_raster(clip_extent_list, in_dem_ras)
    
    benchmark_dict = benchmark(function_start_time, benchmark_dict, "clip.clip_raster")
    
    #------------------------------------------------------------------------------ 
    
    function_start_time = datetime.datetime.now()
    
    # Clip land polygons to array extent.
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
    
    #------------------------------------------------------------------------------ 
    
    function_start_time = datetime.datetime.now()
    
    # Filter rays in radial array.
    function_start_time = datetime.datetime.now()
    insert_cursor_list = array.filter_2d_radial_array(heading_polylines_dict, reclip_land_poly_list)
    
    benchmark_dict = benchmark(function_start_time, benchmark_dict, "array.filter_2d_radial_arra")
    
    #------------------------------------------------------------------------------ 
    
    function_start_time = datetime.datetime.now()
    
    # Interpolate 2d array to 3d array, and return nested list where: [[vertices in line1], [vertices in line2]...].
    function_start_time = datetime.datetime.now()
    radial_array_list = array.interpolate_2d_radial_array_2(obs_x, obs_y, obs_z_list, insert_cursor_list, clip_ras, pr_gdb, pt_crs)
    
    benchmark_dict = benchmark(function_start_time, benchmark_dict, "array.interpolate_2d_radial_array_2")
    
    # Delete clip raster.
    arcpy.Delete_management(clip_ras)
    
    #------------------------------------------------------------------------------ 
    
    # Create list to hold updated rays.
    array_pt_list = []
    array_sample_list = []
    array_landmark_list = []
    
    # Loop through rays.
    for ray_vertex_list in radial_array_list:
        
        # Check that input list has points - ignore rays that have no geometry.
        if len(ray_vertex_list) > 0:
            
            function_start_time = datetime.datetime.now()
            
            # Get 2D distance list for vertices in ray.
            sorted_ray_pt_list, ray_dist_list = ray.sort_vertices_2d(ray_vertex_list)
            
            # Get z values for vertices as a list.
            sorted_ray_z_list = [i[2] for i in sorted_ray_pt_list]
            
            # Densify ray, returning new point list.
            densified_ray_pt_list, densified_ray_dist_list = ray.densify_3d_ray(sorted_ray_pt_list[0], sorted_ray_pt_list[-1], ray_dist_list, sorted_ray_z_list, dist_inc)
            
            benchmark_dict = benchmark(function_start_time, benchmark_dict, "ray.densify_3d_ray")
            
            #------------------------------------------------------------------------------ 
            
            function_start_time = datetime.datetime.now()
            
            # Adjust z values in densified point list to account for curvature.
            pr_ray_pt_list = ray.adjust_curvature(obs_x, obs_y, densified_ray_pt_list)
            
            # Add modified array to list.
            array_pt_list.append(pr_ray_pt_list)
            
            benchmark_dict = benchmark(function_start_time, benchmark_dict, "ray.adjust_curvature")
            
            #------------------------------------------------------------------------------ 
            
            # Create list of sample values, if sampling is enabled, and add to list.
            if sample_ras != "":
                
                function_start_time = datetime.datetime.now()
                
                array_sample_list.append(ray.sample_raster(sample_ras, pr_ray_pt_list, pt_crs))
                
                benchmark_dict = benchmark(function_start_time, benchmark_dict, "ray.sample_raster")
            
            #------------------------------------------------------------------------------ 
                
            # Create list of landmark values, if enabled, and add to list.
            if lmark_geom_list != "":
                
                function_start_time = datetime.datetime.now()
                
                array_landmark_list.append(ray.count_landmarks(obs_x, obs_y, max_dist, lmark_geom_list, pr_ray_pt_list, pt_crs))
                
                benchmark_dict = benchmark(function_start_time, benchmark_dict, "ray.count_landmarks")
            
    #------------------------------------------------------------------------------ 
    
    # Create dictionary for output.
    obs_dict = {}
    
    # Loop through observer points, calculate visibility for all points.
    for obs_z in obs_z_list:
        
        obs_v_angle_list = []
        obs_vis_list = []
        
        # Loop through rays in radial array list.
        for iter_ray in array_pt_list:
            
            # Calculate angle for vertices in list, add to observer point list.
            obs_v_angle_list.append(ray.angle_list(obs_x, obs_y, obs_z, iter_ray))
            
            # Calculate visibility for vertices in list, add to observer point list.
            obs_vis_list.append(ray.visibility_list(iter_ray))
        
        dist_dict = {}
           
        # Loop through indices in radial array list.
        for dist_idx, dist_val in enumerate(densified_ray_dist_list):
            
            function_start_time = datetime.datetime.now()
            
            # Get stats for vertical angle.
            v_angle_sum = array.array_stats(obs_v_angle_list, obs_vis_list, dist_idx, "MAX", "SUM")
            v_angle_max = array.array_stats(obs_v_angle_list, obs_vis_list, dist_idx, "MAX", "MAX")
            v_angle_min = array.array_stats(obs_v_angle_list, obs_vis_list, dist_idx, "MIN", "MIN")
            v_angle_range = abs(v_angle_max - v_angle_min)
            
            # Get stats for horizontal angle.
            h_angle_sum = len(array_pt_list)
            
            # Create output list.
            out_row = [v_angle_sum, v_angle_range, h_angle_sum]
            
            # Get stats for sample raster, if enabled.
            if sample_ras != "":
                sample_max = array.array_stats(array_sample_list, obs_vis_list, dist_idx, "MAX", "MAX")
                sample_min = array.array_stats(array_sample_list, obs_vis_list, dist_idx, "MIN", "MIN")
                sample_avg = array.array_stats(array_sample_list, obs_vis_list, dist_idx, "AVG", "AVG")
            
                for i in [sample_max, sample_min, sample_avg]:
                    out_row.append(i)
            
            # Get stats for landmarks, if enabled.
            if lmark_geom_list != "":
                landmark_sum = array.array_stats(array_landmark_list, obs_vis_list, dist_idx, "SUM", "SUM")
            
                out_row.append(landmark_sum)
            
            benchmark_dict = benchmark(function_start_time, benchmark_dict, "array.array_stats")
            
            #------------------------------------------------------------------------------ 
            
            # Update distance dictionary with rounded values.
            dist_dict[dist_val] = [round(i, 4) for i in [out_row]]
        
        # Add distance dictionary to observer dictionary as subdictionary.
        obs_dict[obs_z] = dist_dict
        
    # Return observer dictionary.
    return obs_dict, benchmark_dict
