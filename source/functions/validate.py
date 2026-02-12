'''
Created on Apr 14, 2024

@author: Karl
'''

import arcpy, os
from functions import console

#===============================================================================
# Check if value can be cast to float.
#===============================================================================

def floatable(in_val):
    try:
        float(in_val)
        return True
    except:
        raise Exception("Could not cast value '{}' to float.".format(in_val))

#===============================================================================
# Checks that input datasets have the same coordinate system.
#===============================================================================

def validate_crs(fc_list):
    
    # Get list of spatial reference object names for features (ignore empty strings).
    sr_list = [arcpy.Describe(i).spatialReference.name for i in fc_list if i != ""]
    
    # Check that all spatial reference objects are the same.
    if len(set(sr_list)) == 1:
        
        # Return spatial reference object.
        return arcpy.Describe(fc_list[0]).spatialReference
    
    # If not (or there are no datasets), raise an exception.
    else:
        print(fc_list)
        print(sr_list)
        raise Exception("Input datasets do not have the same CRS.")
    
#===============================================================================
# Checks if strings are valid paths.
#===============================================================================

def check_paths(path_list):
    
    path_msglist = []
    script_abort = False
    
    for path in path_list:
        
        # If path exists, it's valid.
        if os.path.exists(path) == False:
            
            try:
                
                # Create and delete dummy file to check path.
                open(path, 'x')
                os.remove(path)
                
            except:
                
                # If OS couldn't create dummy file, path is not valid.
                path_msglist.append("    {}\n".format(path))
                script_abort = True
    
    if script_abort:
        raise Exception("The following paths are invalid:\n{}".format("".join(path_msglist)))
    
#===============================================================================
# Checks that 'densify distance' is not greater than the 'distance range', and 
# prints a warning if the distance is too large.
#===============================================================================

def check_densify_dist(densify_dist, dist_list):
    
    # Find the minimum distance between elements in the dist_list.
    min_dist_dist = min([abs(dist_list[i]-dist_list[i-1]) for i, j in enumerate(dist_list[1:])]) #@UnusedVariable
    
    if densify_dist >= min_dist_dist:
        
        raise Exception("The parameter 'densify_distance' cannot exceed the minimum distance between bands in the 'dist_range' (or override values).\ndensify_distance={}\nminimum_distance_range={}".format(str(densify_dist), str(min_dist_dist)))
            
    elif densify_dist >= (min_dist_dist * 0.5):
        
        console.console("WARNING: The parameter 'densify_dist' is greater than half the minimum distance between visibility bands.", 2)
        console.console("A 'densify_dist' roughly equivalent to the resolution of the input DEM is recommended.", 11)
        console.console("Visibility indices may not be accurate.", 11)
            
            