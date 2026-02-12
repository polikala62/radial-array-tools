'''
Created on Mar 18, 2024

@author: Karl
'''

#===============================================================================
# Searches an input list for a value and returns the index.
#===============================================================================

def list_index_from_val(in_list, in_val):
    
    # Get list of indices for entries in list that match value.
    check_list = [i[0] for i in enumerate(in_list) if i[1] == in_val]
    
    # Check result of check_list - if there's one value, then return it.
    if len(check_list) == 1:
        
        return check_list[0]
    
    # If there are none, there were no matches.
    elif check_list == 0:
        
        raise Exception("Value '{}' not found in input list '{}'.".format(in_val, ",".join(in_list)))
    
    # If list has more than one entry, there are multiple instances of in_val in it.
    else:
        
        raise Exception("Found multiple instances of value '{}' in list '{}'.".format(in_val, ",".join([str(i) for i in in_list])))

#===============================================================================
# Finds the nearest value in a list to an input value. Returns None if no values
# in list are within threshold distance to input value.
#===============================================================================

def nearest_dist_val(in_val, in_list, threshold):
    
    out_val = None
    out_diff = None
    
    # Loop through values in input.
    for list_val in in_list:
        
        # Find distance between search value and iterated value.
        check_diff = in_val-list_val
        
        # Check if distance is within limits.
        if abs(check_diff) <= threshold:
            
            # Update out_val and out_diff if first iteration.
            if out_val == None:
                
                out_val = list_val
                out_diff = abs(in_val - out_val)
            
            # Update out_val and out_diff if iterated difference is lower (and also less than zero).
            elif abs(in_val-out_val) < out_diff and in_val-out_val <= 0:
                
                out_val = list_val
                out_diff = in_val - out_val
                
    return out_val