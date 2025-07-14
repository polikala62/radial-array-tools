'''
Created on Feb 15, 2024

@author: Karl

These functions retrieve data from, create, and append to NetCDF files. Written for CAA 2024.
'''

import netCDF4 as nc
import numpy as np

def create_netcdf(out_ncdf_path, file_dict, dim_dict, var_dict, var_atts_dict, var_arrays_dict):
    
    with nc.Dataset(out_ncdf_path, 'w', format="NETCDF4") as output:
        
        # Set global attributes for file.
        output.setncatts(file_dict)
        
        # Set dimensions. Dimension dictionary should look like {Y: <class 'netCDF4._netCDF4.Dimension'>: name = 'Y', size = 471}
        for dim_name in dim_dict.keys():
            dimension = dim_dict[dim_name]
            output.createDimension(dim_name, dimension)
            #output.createDimension(dim_name, (len(dimension) if not dimension.isunlimited() else None))
            
        for var_name in var_dict.keys():
            variable = var_dict[var_name]
            variable_attributes = var_atts_dict[var_name]
            variable_array = var_arrays_dict[var_name]
            
            # Create the variable.
            output.createVariable(var_name, variable['datatype'], variable['dimensions'])
            
            # Add attributes to variable.
            output[var_name].setncatts(variable_attributes)
            
            # Add data to variable.
            output[var_name][:] = variable_array
            