'''
Created on Dec 13, 2023

@author: Karl
'''

def vis_pts_int(pt_list, check_geom):
    pass

# Check for disjoint slightly more efficiently by escaping the function as soon as an intersection is detected.
def check_disjoint(in_pt, in_poly_list):
    
    # Loop through polygons in list.
    for poly in in_poly_list:
        
        # If the polygon was not disjoint to the point, return false.
        if in_pt.disjoint(poly) == False:
            
            return False
    
    # If function has not been escaped, then all polygons in list were disjoint, return True.
    return True