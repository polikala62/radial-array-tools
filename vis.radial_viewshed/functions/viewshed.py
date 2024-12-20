'''
Created on Dec 13, 2023

@author: Karl
'''

import arcpy, datetime
from functions import array, clip, heading, ray, validate, list
import output
from functions.benchmark import benchmark

def radial_viewshed(obs_x, obs_y, obs_z_list, dist_list, in_dem_ras, in_dem_res, pr_gdb, pt_crs, densify_dist, obs_z_offset=0,
                    lmark_geom_list="", land_fc="", override_min_dist="", sample_ras="", benchmark_dict=""):
    
    # Get distance min, max, and increment from variable.
    min_dist = min(dist_list)
    max_dist = max(dist_list)
    
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
    array_2d_dist_list = []
    array_null_list = []
    array_sample_list = []
    array_landmark_list = []
    
    # Create counter to track processed rays.
    pr_rays = 0
    
    # Loop through rays.
    for ray_vertex_list in radial_array_list:
        
        # Check that input list has points - ignore rays that have no geometry.
        if len(ray_vertex_list) > 0:
            
            # Increment processed rays.
            pr_rays += 1
            
            function_start_time = datetime.datetime.now()
            
            # Get 2D distance list for vertices in ray.
            sorted_ray_pt_list, ray_dist_list = ray.sort_vertices_2d(obs_x, obs_y, ray_vertex_list)
            
            # Get z values for vertices as a list.
            sorted_ray_z_list = [i[2] for i in sorted_ray_pt_list]
            
            # Densify ray, returning new point list.
            densified_ray_pt_list, densified_ray_dist_list, densified_ray_null_list = ray.densify_3d_ray([obs_x, obs_y], ray_vertex_list[-1], ray_dist_list, sorted_ray_z_list, densify_dist, max_dist)
            
            array_2d_dist_list.append(densified_ray_dist_list)
            array_null_list.append(densified_ray_null_list)
            
            benchmark_dict = benchmark(function_start_time, benchmark_dict, "ray.densify_3d_ray")
            
            #------------------------------------------------------------------------------ 
            
            function_start_time = datetime.datetime.now()
            
            # Adjust z values in densified point list to account for curvature.
            pr_ray_pt_list = ray.adjust_curvature_2(obs_x, obs_y, densified_ray_pt_list)
            
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
    
    # Check that there are rays in the radial array.
    if pr_rays > 0:
    
        # Create dictionary for output.
        obs_dict = {}
         
        # Loop through observer points, calculate visibility for all points.
        for obs_z in obs_z_list:
            
            obs_v_angle_list = []
            obs_vis_list = []
            
            # Loop through rays in radial array list.
            for iter_ray_idx, iter_ray in enumerate(array_pt_list):
                
                # Get null list for ray.
                iter_ray_nulls = array_null_list[iter_ray_idx]
                
                # Calculate angle for vertices in list, add to observer point list.
                ray_angle_list = ray.angle_list(obs_x, obs_y, obs_z, obs_z_offset, iter_ray, iter_ray_nulls)
                obs_v_angle_list.append(ray_angle_list)
                
                function_start_time = datetime.datetime.now()
                
                # Get z values from list.
                ray_z_vals = [i[2] for i in iter_ray]
                
                # Zip distances and z vals to get list for visibility.
                vis_list = [[array_2d_dist_list[iter_ray_idx][i], ray_z_vals[i]] for i in range(0, len(ray_z_vals))]
                
                #simp_vis_list = list.simplify_vis_list(vis_list)
                
                # Calculate visibility for vertices in list, add to observer point list.
                #obs_vis_list.append(ray.visibility_list(ray_z_vals, obs_z))
                obs_vis_list.append(ray.visibility_list_2(obs_z, obs_z_offset, vis_list, iter_ray_nulls))
                
                benchmark_dict = benchmark(function_start_time, benchmark_dict, "ray.visibility_list")
                
            dist_dict = {}
            
            # Loop through indices in radial array list.
            for dist_idx, dist_val in enumerate(densified_ray_dist_list):
                
                # Only add to output if it's in the input list.
                if dist_val in dist_list:
                    
                    function_start_time = datetime.datetime.now()
                    #print('v_angle stats', dist_val)
                    # Get stats for vertical angle.
                    v_angle_sum = array.array_stats(obs_v_angle_list, obs_vis_list, dist_idx, "MAX", "SUM")
                    v_angle_max = array.array_stats(obs_v_angle_list, obs_vis_list, dist_idx, "MAX", "MAX")
                    #print('v angle stats end')
                    # Get stats for horizontal angle.
                    h_angle_sum = array.array_stats(obs_vis_list, obs_vis_list, dist_idx, "MAX", "SUM")
                    
                    # Create output list.
                    out_row = [v_angle_sum, v_angle_max, h_angle_sum]
                    
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
                    dist_dict[dist_val] = [round(float(i), 4) for i in out_row]
            
            # Add distance dictionary to observer dictionary as subdictionary.
            obs_dict[obs_z] = dist_dict
            
        # Return observer dictionary.
        return obs_dict, benchmark_dict
    
    else:
        
        out_fieldcount = 4
        if sample_ras != "":
            out_fieldcount += 3
        if lmark_geom_list != "":
            out_fieldcount += 1
        
        # Create dictionary for output.
        obs_dict = {}
         
        # Loop through observer points, calculate visibility for all points.
        for obs_z in obs_z_list:
            
            dist_dict = {}
            
            # Loop through rays in radial array list.
            for iter_ray in array_pt_list:
                
                # Update distance dictionary with zero values.
                dist_dict[dist_val] = [0 for i in range(0, out_fieldcount)]
                
            # Add distance dictionary to observer dictionary as subdictionary.
            obs_dict[obs_z] = dist_dict
            
        # Return observer dictionary.
        return obs_dict, benchmark_dict
