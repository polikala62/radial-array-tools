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
        
        raise Exception("Found multiple instances of value '{}' in list '{}'.".format(in_val, ",".join(in_list)))