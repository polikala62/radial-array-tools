'''
Created on Dec 13, 2023

@author: Karl
'''

import arcpy, math

from functions.clip import clip_raster
from functions.intersect import check_disjoint
'''
def curvature_drop(z_surface, dist): # Assumes a spherical earth.
    
    # Adapted from a confusing formula at: https://pro.arcgis.com/en/pro-app/latest/tool-reference/3d-analyst/using-viewshed-and-observer-points-for-visibility.htm
    
    diam_earth = 12740000
    
    r_refr = 0.13
    
    z_actual = z_surface - (dist^2)/diam_earth + r_refr * (dist^2)/diam_earth
    
    return z_actual
'''
def curvature_drop(z_surface, dist):
    
    rad_earth = 6370000
    
    drop = math.sqrt((rad_earth**2)+(dist**2)) - rad_earth
    
    return z_surface - drop

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

# Creates a radial array, and selects points from array based on the raster GetCellValue tool.
def getcellvalue_radial_array(obs_pt, in_ras, max_dist, raster_res=5, degree_interval=1): # Returns dictionary objects.
    
    # Get observation point coordinates.
    obs_x, obs_y, obs_z = [float(i) for i in obs_pt]
    
    # Create list for output.
    out_list = []
    
    # Loop through degree intervals.
    for deg_heading in range(0,360,degree_interval):
        
        heading = math.radians(deg_heading)
        
        # Create heading list.
        heading_dict = {}
        
        # Create counter for points along heading.
        heading_pt_count = 0
        
        # Loop through distance intervals.
        for vector_distance in range(0,max_dist,raster_res):
            
            # Increment point counter.
            heading_pt_count += 1
            
            # Create point key.
            pt_key = "{}-{}".format(deg_heading, heading_pt_count)
            
            # Calculate endpoint for vector x value.
            pt_x = (math.sin(heading)*vector_distance) + obs_x
            
            # Calculate endpoint for vector y value.
            pt_y = (math.cos(heading)*vector_distance) + obs_y
            
            # Sample raster at point.
            pt_z_result = arcpy.management.GetCellValue(in_ras, "{} {}".format(pt_x, pt_y))
            # Create list of sample values (1 per band).
            pt_z_list = pt_z_result.getOutput(0).split("\n")
            # Get z point from input band (1 is default).
            pt_z = float(pt_z_list[0])
            
            # Discard z point if it is less than obsever z.
            if pt_z < obs_z:
                
                heading_dict[pt_key] = [pt_x, pt_y, obs_z]
                
            # Check that point has z value.
            elif pt_z != "NoData":
                
                # Add point to heading list.
                heading_dict[pt_key] = [pt_x, pt_y, pt_z]
        
        # Add heading list to output list.
        out_list.append(heading_dict)
    
    # Return list of all coordinates for all headings.
    return out_list

#------------------------------------------------------------------------------ 

# Creates a radial array, and selects points from the array based on the Extract Multivalue tool.
def extractmultivalue_radial_array(obs_pt, in_ras, max_dist, pr_gdb, coast_polygon=None, raster_res=5, degree_interval=1):
    pass
    # Loop through degrees, create point dictionary with values [X, Y, HEADING, DIST_TO_OBS].
    
    # Save point features to memory.
    
    # Extract raster values to point features.
    
    # Read points to dictionary, with values [X, Y, Z, HEADING, DIST_TO_OBS].
    
    # Sort feature vertices by HEADING, DIST_TO_OBS.

#------------------------------------------------------------------------------ 

def interpolate_3d_radial_array(obs_pt, in_ras, min_dist, max_dist, pr_gdb, pr_crs, coast_polylines=None, raster_res=5, degree_interval=1):
    
    # Set the geoprocessing workspace
    arcpy.env.workspace = pr_gdb
    
    # Loop through degrees, create polyline objects if polylines intersect coast.
    
    # Get observation point coordinates.
    obs_x, obs_y, obs_z = [float(i) for i in obs_pt] #@UnusedVariable
    
    # Create list for output.
    out_list = []
    
    # Create list for raster clipping.
    raster_clip_list = []
    
    # Create list for insert cursor.
    insert_cursor_list = []
    
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
        
        if coast_polylines == None:
            
            insert_cursor_list.append([deg_heading, heading_polyline])
            
        else:
            
            if check_disjoint(heading_polyline, coast_polylines) == False:
            #if True in [heading_polyline.crosses(i) for i in coast_polylines]:
                
                insert_cursor_list.append([deg_heading, heading_polyline])
    
    # Clip raster to fit polyline features.
    clip_ras = clip_raster(raster_clip_list, (raster_res*2), in_ras)
    
    if len(insert_cursor_list) > 0:
    
        # Save polyline features to memory.
        arcpy.management.CreateFeatureclass(r"memory", "array_2d", geometry_type="POLYLINE", has_m="DISABLED", has_z="ENABLED", spatial_reference=pr_crs)
        
        with arcpy.da.InsertCursor(r"memory\array_2d", ["OID@", "SHAPE@"]) as cursor: #@UndefinedVariableFromImport
            
            for row in insert_cursor_list:
                
                # Insert row.
                cursor.insertRow(row)
        
        del cursor
        
        # Interpolate polyline features.
        arcpy.ddd.InterpolateShape(clip_ras, r"memory\array_2d", r"memory\array_3d")
        
        # Delete raster.
        arcpy.Delete_management(clip_ras)
        
        # Read feature vertices to dictionary, with values [X, Y, Z, HEADING].
        with arcpy.da.SearchCursor(r"memory\array_3d", ["OID@", "SHAPE@"]) as cursor: #@UndefinedVariableFromImport
            
            for row in cursor:
                
                heading_dict = {}
                
                pr_list = []
                
                for part in row[1]:
                    
                    for pnt in part:
                        
                        # Find the drop in distance due to curvature.
                        pt_obs_dist = math.dist([obs_x, obs_y], [pnt.X, pnt.Y])
                        rad_earth = 6370000 # Earth's radius, cribbed from: https://pro.arcgis.com/en/pro-app/latest/tool-reference/3d-analyst/using-viewshed-and-observer-points-for-visibility.htm
                        curv_offset = (math.sqrt((rad_earth**2)+(pt_obs_dist**2)) - rad_earth)
                        
                        # TODO: Try to factor in refraction? But that's also contingent on temperature.
                        
                        # Modify point z to take curvature into account.
                        mod_pnt_z = pnt.Z - curv_offset
                        
                        # Discard z point if point + offset is less than observer z.
                        if mod_pnt_z >= obs_z:
                            
                            # Add point to output.
                            pr_list.append([pnt.X, pnt.Y, mod_pnt_z])
                
                # Sort list by distance.
                for idx, pt in enumerate(sort_pts_by_dist_2d([obs_x, obs_y], pr_list)):
                    
                    # Create point key and add to dictionary.
                    pt_key = "{}-{}".format(row[0], idx)
                    heading_dict[pt_key] = pt
                
                # Add dictionary to output.
                out_list.append(heading_dict)
        
        # Delete cursor.
        del cursor
        
        # Delete processing datasets.
        arcpy.Delete_management(r"memory\array_2d")
        arcpy.Delete_management(r"memory\array_3d")
        
    # Return list of all coordinates for all headings.
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

def interpolate_2d_radial_array(obs_x, obs_y, obs_z_list, insert_cursor_list, in_ras, pr_gdb, pr_crs):
    
    # Set the geoprocessing workspace
    arcpy.env.workspace = pr_gdb
    
    # Create dictionary for output.
    out_dict = {}
    
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
        
        # Loop through z values.
        for obs_z in obs_z_list:
            
            # Create list to hold heading dictionaries.
            iter_heading_dict_list = []
                
            # Read feature vertices to dictionary, with values [X, Y, Z, HEADING].
            with arcpy.da.SearchCursor(r"memory\array_3d", ["OID@", "SHAPE@"]) as cursor: #@UndefinedVariableFromImport
            
                # Loop through rows in cursor.
                for row in cursor:
                    
                    # Create dictionary to hold visibility values per heading.
                    heading_dict = {}
                    
                    # Create processing list to hold filtered values.
                    pr_list = []
                    
                    # Check that row has a shape.
                    if row[1] != None:
                        
                        for part in row[1]:
                            
                            for pnt in part:
                                
                                # Find the drop in distance due to curvature.
                                pt_obs_dist = math.dist([obs_x, obs_y], [pnt.X, pnt.Y])
                                rad_earth = 6370000 # Earth's radius, cribbed from: https://pro.arcgis.com/en/pro-app/latest/tool-reference/3d-analyst/using-viewshed-and-observer-points-for-visibility.htm
                                curv_offset = (math.sqrt((rad_earth**2)+(pt_obs_dist**2)) - rad_earth)
                                
                                # TODO: Try to factor in refraction? But that's also contingent on temperature.
                                
                                # Modify point z to take curvature into account.
                                mod_pnt_z = pnt.Z - curv_offset
                                
                                # Discard z point if point + offset is less than observer z.
                                if mod_pnt_z >= obs_z:
                                    
                                    # Add point to output.
                                    pr_list.append([pnt.X, pnt.Y, mod_pnt_z])
                        
                        # Sort list by distance.
                        for idx, pt in enumerate(sort_pts_by_dist_2d([obs_x, obs_y], pr_list)):
                            
                            # Create point key and add to dictionary.
                            pt_key = "{}-{}".format(row[0], idx)
                            heading_dict[pt_key] = pt
                    
                    # Add dictionary to output.
                    if len(heading_dict.keys()) > 0:
                        iter_heading_dict_list.append(heading_dict)
                
                # Add list of heading dicts to output dictionary, with obs_z as key.
                out_dict[obs_z] = iter_heading_dict_list
                
            # Delete cursor.
            del cursor
        
        # Delete processing datasets.
        arcpy.Delete_management(r"memory\array_2d")
        arcpy.Delete_management(r"memory\array_3d")
        
    # Return output dictionary.
    return out_dict

#------------------------------------------------------------------------------ 
