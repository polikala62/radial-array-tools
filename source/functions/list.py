'''
Created on Mar 18, 2024

@author: Karl
'''

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
# 
#===============================================================================

def match_val_to_list(in_val, in_list, threshold):
    
    diff_list = [abs(in_val-i) for i in in_list if (in_val-i) <= 0]
    
    print(diff_list)
    if min(diff_list) <= threshold:
        
        return True
    
    else:
        
        return False
    
def nearest_dist_val(in_val, in_list, threshold):
    
    out_val = None
    out_diff = None
    
    for list_val in in_list:
        
        check_diff = in_val-list_val
        
        if abs(check_diff) <= threshold:
            
            if out_val == None:
                
                out_val = list_val
                out_diff = abs(in_val - out_val)
                
            elif abs(in_val-out_val) < out_diff and in_val-out_val <= 0:
                
                out_val = list_val
                out_diff = in_val - out_val
                
    return out_val