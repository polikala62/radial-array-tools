'''
Created on Dec 13, 2023

@author: Karl
'''

#===============================================================================
# Checks for disjoint between an input point and a list of polygons.
# Function escapes as soon as an intersection is detected.
#===============================================================================

def check_disjoint(in_pt, in_poly_list):
    
    # Loop through polygons in list.
    for poly in in_poly_list:
        
        # If the polygon was not disjoint to the point, return false.
        if in_pt.disjoint(poly) == False:
            
            return False
    
    # If function has not been escaped, then all polygons in list were disjoint, return True.
    return True

#===============================================================================
# Returns list of distances from point to polygons, and returns empty list if 
# point is inside polygons.
#===============================================================================

def check_distance_to(in_pt, in_poly_list):
    
    # Create list for output.
    out_list = []
    
    # Loop through polygons in list.
    for poly in in_poly_list:
        
        # Calculate distance to.
        pt_dist = in_pt.distanceTo(poly)
        
        # If the point is within the polygon, return empty list.
        if pt_dist == 0:
            
            return []
        
        # If not, add point distance to the list.
        else:
            
            out_list.append(pt_dist)
    
    # If function has not been escaped, then all polygons in list were disjoint, return output list.
    return out_list