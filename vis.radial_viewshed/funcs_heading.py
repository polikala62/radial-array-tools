'''
Created on Dec 13, 2023

@author: Karl
'''

import math

#------------------------------------------------------------------------------ 

def heading_sub_angle(obs_pt, heading_dict): # Returns maximum subtended angle.
    
    # Get observer coordinates as separate variables.
    obs_x, obs_y, obs_z = [float(i) for i in obs_pt]
    
    # Create dictionary for output.
    out_dict = {}
    
    # Loop through points in input heading.
    for key in heading_dict.keys():
        
        # Get coordinates as separate variables.
        pt_x, pt_y, pt_z = [float(i) for i in heading_dict[key][0:3]]
        
        # Create entry in dictionary.
        out_dict[key] = heading_dict[key]
        
        # Get distance from observer to target.
        e_dist = math.dist([obs_x, obs_y], [pt_x, pt_y])
        
        # Check if observer elevation and point elevation are the same.
        if pt_z <= obs_z or e_dist == 0:
            
            # Set sub angle to zero.
            out_dict[key].append(0)
            
        else:
            
            # Get elevation difference.
            z_diff = abs(pt_z - obs_z)
            
            # Find subtended angle, where e_dist is adjacent side and z_diff is opposite side.
            sub_angle = math.degrees(math.atan((z_diff/e_dist)))
            
            # Add to output list.
            out_dict[key].append(sub_angle)
    
    # Return output list.
    return out_dict

#------------------------------------------------------------------------------ 

def heading_vis_pts(heading_dict, min_z): # Returns list of visible coordinates.
    
    # Create dict for output.
    out_dict = {}
    
    # Set max_z to input value.
    max_z = min_z
    
    # Loop through points in heading. Points should already be sorted from closest to -> furthest from obs pt.
    for key in heading_dict.keys():
        
        out_dict[key] = heading_dict[key]
        
        # Get z value for point in heading list.
        iter_z = out_dict[key][2]
        
        if iter_z >= max_z:
            # Add coordinates to output.
            out_dict[key].append(1)
            # Update max_z.
            max_z = iter_z
            
        else:
            # Add coordinates to output.
            out_dict[key].append(0)
    
    # Return list of visible points in heading.
    return out_dict

                