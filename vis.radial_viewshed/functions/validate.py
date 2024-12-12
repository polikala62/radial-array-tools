'''
Created on Apr 14, 2024

@author: Karl
'''

import arcpy

def floatable(in_val):
    try:
        float(in_val)
        return True
    except:
        raise Exception("Could not cast value '{}' to float.".format(in_val))

# Create function to check that input datasets have the same coordinate system.
def validate_crs(fc_list):
    
    # Get list of spatial reference object names for features (ignore empty strings).
    sr_list = [arcpy.Describe(i).spatialReference.name for i in fc_list if i != ""]
    
    # Check that all spatial reference objects are the same.
    if len(set(sr_list)) == 1:
        
        # Return spatial reference object.
        return arcpy.Describe(fc_list[0]).spatialReference
    
    # If not (or there are no datasets), raise an exception.
    else:
        
        raise Exception("Input datasets do not have the same CRS.")