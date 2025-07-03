'''
Created on Dec 9, 2024

@author: Karl
'''
import arcpy, math
from arcpy.sa import ExtractValuesToPoints

#===============================================================================
# CLOCKWISE_ANGLE
# Calculates the clockwise angle between two vectors.
#===============================================================================

def segment_angle(in_segment): # Segment must be list with projected XY coordinate-pair (i.e. [x,y]). Values are in degrees from 0=N.
    
    # Translate segment to origin.
    segment_pt_1, segment_pt_2 = in_segment
    mod_vector = [segment_pt_2[0] - segment_pt_1[0], segment_pt_2[1] - segment_pt_1[1]]
    
    # Define vector where y=1.
    #ref_vector = [0, 1]
    
    # Calculate the dot product of input vectors.
    #dot = sum([i*j for (i, j) in zip(ref_vector, mod_vector)]) # Calculate without numpy.
    dot = mod_vector[1] # Calculate without numpy.
    
    # Calculate the magnitudes of input vector.
    vector2_mag = math.sqrt((mod_vector[0]**2) + (mod_vector[1]**2))
    
    # Calculate the determinant.
    #determinant = ref_vector[0]*mod_vector[1]-ref_vector[1]*mod_vector[0]
    determinant = 0-mod_vector[0]
    
    # Attempt to caclulate the angle between the two vectors.
    try:
        
        # Calculate the angle between the two vectors, in degrees.
        degree_angle = math.degrees(math.acos(dot/vector2_mag))
    
    except ZeroDivisionError:
        
        raise Exception("Could not calculate angle between vectors with magnitude 1 and {} - magnitude of input angle is 0.".format(vector2_mag))
        
    except:
        
        # If there is another error, file report.
        raise Exception("Could not calculate angle between vectors with magnitude 1 and {} - unknown error occured.".format(vector2_mag))
        
    if determinant < 0:
        
        return degree_angle

    else:
    
        return 360 - degree_angle

#===============================================================================
# DISTANCE_BEARING_TO_VECTOR
#===============================================================================

def distance_bearing_to_vector(bearing, distance): # Assumes ref angle of 0 = NORTH and clockwise orientation.
    
    x_component = math.sin(math.radians(bearing))*distance
    y_component = math.cos(math.radians(bearing))*distance
    
    # Change values based on the quadrant the vector is in.
    if bearing >= 90 and bearing < 180:
        
        y_component *= -1
        
    elif bearing >= 180 and bearing < 270:
        
        x_component *= -1
        y_component *= -1
        
    elif bearing >= 270:
        
        x_component *= -1
        
    return x_component, y_component

#===============================================================================
# 
#===============================================================================

def find_segment_from_list(in_dist, dist_list, z_list):
    
    for idx in range(1, len(dist_list)):
        # If distance is less than first distance value in list or greater than last value, return None.
        if in_dist < dist_list[0] or in_dist > dist_list[-1]:
            return None
        elif in_dist >= dist_list[idx-1] and in_dist <= dist_list[idx]:
            return [[dist_list[idx-1], z_list[idx-1]], [dist_list[idx], z_list[idx]]]

#===============================================================================
# 
#===============================================================================

def sort_vertices_2d(obs_x, obs_y, ray_pts_list):
    
    ray_dist_list = []
    
    # Loop through vertices in ray.
    for pt in ray_pts_list:
        
        pt_x = pt[0]
        pt_y = pt[1]
        
        # Calculate the 2D distance from the observer point to the target point.
        point_dist_cartesian = math.sqrt((pt_x-obs_x)**2 + (pt_y-obs_y)**2)
        
        # Add values to list.
        ray_dist_list.append(point_dist_cartesian)
        
    # Sort lists by 2D distance.
    zip_list = sorted(zip(ray_dist_list, ray_pts_list))
    sorted_dist_list = [i[0] for i in zip_list]
    sorted_pts_list = [i[1] for i in zip_list]
    
    return [sorted_pts_list, sorted_dist_list]

#===============================================================================
# 
#===============================================================================

def densify_3d_ray(start_pt, end_pt, ray_dist_list, ray_z_list, densify_dist, min_dist, max_dist):
    
    # Create lists to hold output.
    out_pt_list = []
    out_dist_list = []
    out_null_list = []
    
    # Find length of ray.
    #ray_length = math.sqrt((start_pt[1]-end_pt[1])**2 + (start_pt[0]-end_pt[0])**2)
    ray_length = max_dist
    
    # Find angle for ray.
    ray_angle_2d = segment_angle([start_pt, end_pt])
    
    # Loop through distance increments.
    iter_dist = densify_dist
    while iter_dist <= (ray_length):
    
        # Find point x and y values.
        iter_x, iter_y = distance_bearing_to_vector(ray_angle_2d, iter_dist)
        
        if iter_dist > min_dist:
        
            # Use the distance list to find start and endpoints for the line segment that intersects the point.
            z_segment = find_segment_from_list(iter_dist, ray_dist_list, ray_z_list)
            
            # If point is not coincident with line, return elevation zero (i.e. sea level).
            if z_segment == None:
                
                iter_z = 0
                out_null_list.append(True)
        
            else:
                segment_angle_3d = segment_angle(z_segment)
                
                # The z segment's x units are 2d distance from origin, its y units are elevation.
                z_segment_start_x, z_segment_start_y = z_segment[0]
                
                # If z val is the adjacent side of a right triangle in profile view, then 2d distance is opposite side, and segment angle is the angle of the triangle.
                iter_z = ((iter_dist - z_segment_start_x) / math.tan(math.radians(segment_angle_3d))) + z_segment_start_y
                
                out_null_list.append(False)
        
        else:
            
            iter_z = 0
            out_null_list.append(True)
            
        out_pt_list.append([iter_x + start_pt[0], iter_y + start_pt[1], iter_z])
        out_dist_list.append(iter_dist)
        
        # Increment iter_dist.
        iter_dist += densify_dist
        
    # Return list.
    return out_pt_list, out_dist_list, out_null_list

#===============================================================================
# 
#===============================================================================
'''
def adjust_curvature(obs_x, obs_y, ray_pts_list):
    
    out_list = []
    
    # Loop through vertices in ray.
    for iter_pt in ray_pts_list:
        
        pt_x, pt_y, pt_z = iter_pt
        
        # Calculate the 2D distance from the observer point to the target point.
        point_dist_cartesian = math.sqrt((pt_x-obs_x)**2 + (pt_y-obs_y)**2)
        
        # Calculate earth curvature offset.
        rad_earth = 6370000 # Earth's radius, cribbed from: https://pro.arcgis.com/en/pro-app/latest/tool-reference/3d-analyst/using-viewshed-and-observer-points-for-visibility.htm
        curvature_offset = (math.sqrt((rad_earth**2)+(point_dist_cartesian**2)) - rad_earth)
        
        # Subtract curvature offset from z value.
        mod_pt_z = pt_z - curvature_offset
        
        out_list.append([pt_x, pt_y, mod_pt_z])
        
    return out_list

'''
def calc_curvature_offset(pt_dist):
        
    # Calculate the angle (in km) from the centre of the earth to the end of the line.
    a = 0.0089932161 * pt_dist
    
    # Get the radius of the earth. (see https://earthcurvature.com/)
    r = 6371
    
    # Calculate earth curvature offset, convert back to metres.
    curvature_offset = (r * (1 - math.cos(math.radians(a)))) * 1000
    
    # Return curvature drop.
    return float(curvature_offset)

def adjust_curvature_2(obs_x, obs_y, ray_pts_list):
    
    out_list = []
    
    # Loop through vertices in ray.
    for iter_pt in ray_pts_list:
        
        pt_x, pt_y, pt_z = iter_pt
        
        # Calculate the 2D distance from the observer point to the target point, convert to km.
        point_dist_cartesian = math.sqrt((pt_x-obs_x)**2 + (pt_y-obs_y)**2) / 1000
        
        # Calculate curvature drop.
        curvature_offset = calc_curvature_offset(point_dist_cartesian)
        
        # Subtract curvature offset from z value.
        mod_pt_z = float(pt_z - curvature_offset)
        
        out_list.append([pt_x, pt_y, mod_pt_z])
        
    return out_list

#===============================================================================
# ANGLE LIST
# This function iterates the 'angle_above_horizon' function for elements of an input list, eliminating entries that are in a list of null values.
#===============================================================================

def angle_list(obs_x, obs_y, obs_surface_z, obs_z_offset, pt_list, null_list):
    
    # Create list for ouptut.
    out_angle_list = []
    
    # Loop through points, indices in input point list.
    for iter_idx, iter_pt in enumerate(pt_list):
        
        # Get point coordinates.
        pt_x, pt_y, pt_z = iter_pt
        
        # Check if point index is 'True' in null_list, and append None if so.
        if null_list[iter_idx]:
            out_angle_list.append(None)
        
        # Compute angle above horizon and add it to list.
        else:
            pt_angle = angle_above_horizon(obs_x, obs_y, obs_surface_z, obs_z_offset, pt_x, pt_y, pt_z)
            out_angle_list.append(pt_angle)
    
    # Return output list.
    return out_angle_list

#===============================================================================
# ANGLE ABOVE HORIZON
# This function measures the angle defined by an observer point, the lowest visibile point of land, and the highest visible point of land. It accounts
# for observer height by calculating angles below and above the viewing horizon, if applicable.
#===============================================================================

def angle_above_horizon(obs_x, obs_y, obs_surface_z, obs_z_offset, pt_x, pt_y, pt_z):
    
    # Calculate the elevation of the observer.
    obs_z = obs_surface_z + obs_z_offset
    
    # Calculate difference between point elevation and observer elevation. If it's negative, the point is below the horizon.
    obs_pt_z_diff = pt_z - obs_z
    
    # Calculate the 2D and 3D distances from observer to target.
    pt_dist_2d = math.sqrt((pt_x-obs_x)**2 + (pt_y-obs_y)**2)
    pt_dist_3d_above_offset = math.sqrt((pt_x-obs_x)**2 + (pt_y-obs_y)**2 + ((pt_z-(obs_surface_z+obs_z_offset)))**2)
    
    # Assume that the 2D and 3D distances are equivalent to adjacent side / hypotenuse of a right triangle in profile view, and get the result in degrees.
    pt_angle_above_offset = math.degrees(math.acos(pt_dist_2d/pt_dist_3d_above_offset))
    
    # If observer is above sea level and point is above the horizon, calculate offset by finding angle where hypotenuse is defined by obsever point and sea level at target.
    if obs_z_offset > 0 and obs_pt_z_diff > 0:
        
        pt_distance_3d_below_offset = math.sqrt((pt_x-obs_x)**2 + (pt_y-obs_y)**2 + obs_z_offset**2)
        
        # Assume that the 2D and 3D distances are equivalent to adjacent side / hypotenuse of a right triangle in profile view, and get the result in degrees.
        pt_angle_below_offset = math.degrees(math.acos(pt_dist_2d/pt_distance_3d_below_offset))
        
        # Return sum of angles above and below the line defined by obs_z.
        return (pt_angle_above_offset + pt_angle_below_offset)
    
    # If point is below the horizon, return a negative angle to indicate the point does not rise above the horizon.
    elif obs_pt_z_diff <= 0:
        return pt_angle_above_offset * -1
    
    # Otherwise, return the angle above the horizon.
    else:
        return (pt_angle_above_offset)

#===============================================================================
# 
#===============================================================================

# An updated version of the original algorithm.
def visibility_list_3(obs_z, obs_offset, check_pt_list, null_list): # check_pt_list is a tuple in format [distance-from-origin, elevation].
    
    # Generate output list.
    out_list = []
    
    # Set maximum z value to observer height.
    max_val = obs_z
    
    # Loop through points.
    for iter_idx, check_pt in enumerate(check_pt_list):
        '''
        if null_list[iter_idx]:
            check_pt[1] += obs_z
        '''
        vis = 0
        
        # Skip the first point (assume it's in the null list, and therefore invisible).
        if iter_idx > 0:
            
            if null_list[iter_idx] == False:
                
                if check_pt[1] > max_val:
                    vis = 1
                    max_val = check_pt[1]
        
        out_list.append(vis)
        #print(check_pt, vis)
    
    return out_list
    
#===============================================================================
# 
#===============================================================================

def sample_raster(in_ras, pt_list, vis_list, coordinate_system):
    
    # Create list for output.
    out_list = []
    
    # Generate temporary feature class.
    arcpy.CreateFeatureclass_management(r"memory", "sr1", geometry_type="POINT", has_m="DISABLED", has_z="ENABLED", spatial_reference=coordinate_system)
    
    # Use insert cursor to update feature class.
    with arcpy.da.InsertCursor(r"memory\sr1", ["OID@", "SHAPE@X", "SHAPE@Y", "SHAPE@Z"]) as cursor: #@UndefinedVariableFromImport
        
        for oid, row in enumerate(pt_list):
            
            if vis_list[oid] > 0:
            
                # Insert row.
                cursor.insertRow([oid] + row)
    
    del cursor
    
    # Get values from raster.
    ExtractValuesToPoints(r"memory\sr1", in_ras, r"memory\sr2", "NONE", "VALUE_ONLY")
    
    # Read values from points.
    with arcpy.da.SearchCursor(r"memory\sr2", ["RASTERVALU"]) as cursor: #@UndefinedVariableFromImport
        
        for row in cursor:
            
            if row[0] not in ["", None]:
            
                # Add row value to output list.
                out_list.append(row[0])
                
            else:
                
                out_list.append(None)
    
    del cursor
    
    # Delete processing datasets.
    for pr_dataset in [r"memory\sr1", r"memory\sr2"]:
        arcpy.Delete_management(pr_dataset)
    
    return out_list

#===============================================================================
# 
#===============================================================================

def count_landmarks(obs_x, obs_y, max_dist, lmark_geom_list, check_pts_list, vis_list, in_crs):
    
    out_list = []
    
    arcpy.env.outputCoordinateSystem = in_crs
    
    # Select from geometry list where geometry is within distance of observer point.
    lmark_mod_geom_list = [i for i in lmark_geom_list if i.distanceTo(arcpy.Point(obs_x, obs_y)) <= max_dist]
    
    # If there are any landmark geometries left, loop through them.
    if len(lmark_mod_geom_list) > 0:
        
        for check_idx, check_pt in enumerate(check_pts_list):
            
            landmark_int = 0
            
            if vis_list[check_idx] > 0:
                
                # Check point may have three coordinates, use only the first two (x and y).
                mod_check_pt = arcpy.Point(check_pt[0], check_pt[1])
                
                for check_geometry in lmark_mod_geom_list:
                    
                    if check_geometry.disjoint(mod_check_pt) == False:
                        
                        landmark_int += 1
                    
            out_list.append(landmark_int)
            
        return out_list
    
    else:
        
        return [0 for i in check_pts_list]
    