'''
Created on Dec 13, 2023

@author: Karl
'''

import arcpy, math

# Disable log history.
if arcpy.GetLogHistory():
    arcpy.SetLogHistory(False)

from functions.intersect import check_disjoint

def sort_pts_by_dist_2d(origin_pt, pt_list, multiplier=10000):
    
    out_list = []
    
    pr_dict = {}
    
    # Get first 2 coordinates from input.
    mod_origin_pt = origin_pt[:2]
    
    for idx, pt_coords in enumerate(pt_list):
        
        mod_pt_coords = pt_coords[:2]
        
        dist_val = int(math.dist(mod_origin_pt, mod_pt_coords) * multiplier)
        
        pr_dict[dist_val] = idx
        
    for idx in sorted(pr_dict.keys()):
        
        out_list.append(pt_list[pr_dict[idx]])
    
    return out_list


#------------------------------------------------------------------------------ 
def generate_2d_radial_array(obs_x, obs_y, min_dist, max_dist, pr_gdb, degree_interval=1):
    
    # Set the geoprocessing workspace
    arcpy.env.workspace = pr_gdb
    
    # Create list for raster clipping.
    raster_clip_list = []
    
    out_dict = {}
    
    # Loop through degree intervals.
    for deg_heading in range(0,360,degree_interval):
        
        # Convert heading to radians.
        heading = math.radians(deg_heading)
        
        # Calculate start point for polyline.
        start_x = (math.sin(heading)*min_dist) + obs_x
        start_y = (math.cos(heading)*min_dist) + obs_y
        
        # Calculate endpoint for polyline.
        pt_x = (math.sin(heading)*max_dist) + obs_x
        pt_y = (math.cos(heading)*max_dist) + obs_y
        
        # Add start and end coordinates to raster clipping list.
        raster_clip_list.append([start_x, start_y])
        raster_clip_list.append([pt_x, pt_y])
        
        # Create a polyline geometry
        heading_array = arcpy.Array([arcpy.Point(start_x, start_y),arcpy.Point(pt_x, pt_y)])
        heading_polyline = arcpy.Polyline(heading_array)
        
        out_dict[deg_heading] = heading_polyline
        
    # Return list for export cursor.
    return raster_clip_list, out_dict

def filter_2d_radial_array(heading_polylines_dict, coast_poly_list):
    
    insert_cursor_list = []
    
    for deg_heading in heading_polylines_dict.keys():
    
        heading_polyline = heading_polylines_dict[deg_heading]
        
        # If there are no land mask features, add polyline to list.
        if coast_poly_list == None:
            
            insert_cursor_list.append([deg_heading, heading_polyline])
        
        # If there are land mask features, only add polyline to list if it is not disjoint (i.e. intersects).
        else:
            
            if check_disjoint(heading_polyline, coast_poly_list) == False:
                
                insert_cursor_list.append([deg_heading, heading_polyline])
                
    return insert_cursor_list

#------------------------------------------------------------------------------ 
def interpolate_2d_radial_array_2(obs_x, obs_y, obs_z_list, insert_cursor_list, in_ras, pr_gdb, pr_crs):
    
    # Set the geoprocessing workspace
    arcpy.env.workspace = pr_gdb
    
    # Create dict for output.
    out_list = []
    
    if len(insert_cursor_list) > 0:
    
        # Save polyline features to memory.
        arcpy.management.CreateFeatureclass(r"memory", "array_2d", geometry_type="POLYLINE", has_m="DISABLED", has_z="DISABLED", spatial_reference=pr_crs)
        
        with arcpy.da.InsertCursor(r"memory\array_2d", ["OID@", "SHAPE@"]) as cursor: #@UndefinedVariableFromImport
            
            for row in insert_cursor_list:
                
                # Insert row.
                cursor.insertRow(row)
        
        del cursor
        
        # Interpolate polyline features.
        arcpy.ddd.InterpolateShape(in_ras, r"memory\array_2d", r"memory\array_3d")
        
                
        # Read feature vertices to dictionary, with values [X, Y, Z, HEADING].
        with arcpy.da.SearchCursor(r"memory\array_3d", ["SHAPE@"]) as cursor: #@UndefinedVariableFromImport
        
            # Loop through rows in cursor.
            for row in cursor:
                
                # Create dictionary to hold vertices.
                vertex_list = []
                
                # Check that row has a shape.
                if row[0] != None:
                    
                    for part in row[0]:
                        
                        for pnt in part:
                            
                            # Add point to list.
                            vertex_list.append([pnt.X, pnt.Y, pnt.Z])
                
                # Add list to dictionary.
                out_list.append(vertex_list)
                
        # Delete cursor.
        del cursor
        
        # Delete processing datasets.
        arcpy.Delete_management(r"memory\array_2d")
        arcpy.Delete_management(r"memory\array_3d")
        
    # Return output dictionary.
    return out_list

#------------------------------------------------------------------------------ 

def array_stats(in_array, vis_list, stop_idx, ray_method, array_method, verbose=False):    
    
    pr_list = []
    
    for ray_idx, ray_list in enumerate(in_array):
        
        clip_list = ray_list[0:stop_idx]
        
        if verbose:
            print('stop_idx:', stop_idx)
            print('ray_idx:', ray_idx)
            print(len(clip_list), 'in c', [round(i, 4) for i in clip_list])
            print(len(vis_list), 'in v:', vis_list)
            print(ray_idx, len(in_array))
            print('v_rounded: ', [round(i, 4) for i in vis_list[ray_idx]])
            
        iter_list = [clip_list[i] for i in range(0, len(clip_list)) if vis_list[ray_idx][i] == 1]
        
        if verbose:
            print('i', [round(i, 4) for i in iter_list])
            print()
        #iter_list = filter_vis(ray_list[0:ray_idx+1], vis_list, ray_idx)
        
        mod_iter_list = [i for i in iter_list if i is not None]
        
        if len(mod_iter_list) > 0:
            
            # Summarise ray according to method.
            if ray_method == "MIN":
                pr_list.append(min(mod_iter_list))
            elif ray_method == "MAX":
                pr_list.append(max(mod_iter_list))
            elif ray_method == "SUM":
                pr_list.append(sum(mod_iter_list))
            elif ray_method == "AVG":
                pr_list.append(sum(mod_iter_list)/len(iter_list))
                
        del iter_list
        
    if len(pr_list) > 0:
        
        #print(stop_idx, pr_list)
        
        # Summarise array according to method.
        if array_method == "MIN":
            return min(pr_list)
        elif array_method == "MAX":
            return max(pr_list)
        elif array_method == "SUM":
            return sum(pr_list)
        elif array_method == "AVG":
            return sum(pr_list)/len(pr_list)
    
    #@TODO: WORK OUT WHAT TO DO WITH THESE CASES. ZERO PROBABLY ISN'T THE RIGHT ANSWER.
    else:
        return 0

#------------------------------------------------------------------------------ 
