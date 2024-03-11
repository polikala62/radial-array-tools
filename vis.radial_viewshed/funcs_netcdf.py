'''
Created on Feb 15, 2024

@author: Karl

These functions retrieve data from, create, and append to NetCDF files. Written for CAA 2024.
'''

import netCDF4 as nc
import numpy as np

def read_nc(in_path):
    
    # Create dataset.
    rootgrp = nc.Dataset(in_path, "r", format="NETCDF4")
    
    # Get rid of masked values (!!!)
    rootgrp.set_auto_mask(False)
    
    print("DATA MODEL =", rootgrp.data_model)
    print("GROUPS =", rootgrp.groups)
    print("DIMENSIONS =", rootgrp.dimensions)
    print("VARIABLES =", rootgrp.variables)
    var_dict = rootgrp.variables
    
    dimvals = []
    
    for var in var_dict.keys():
        
        dimvals.append(rootgrp.variables[var][:])
    
    print(dimvals)
    
    # Convert vals list to numpy array.
    
    rootgrp.close()
    
def netcdf_var_to_array(ncdf_path, var_name, preserve_masked_vals=False):
    
    # Create dataset.
    rootgrp = nc.Dataset(ncdf_path, "r", format="NETCDF4")
    
    # Get rid of masked values (!!!)
    rootgrp.set_auto_mask(preserve_masked_vals)
    
    var_dict = rootgrp.variables
    
    if var_name in var_dict.keys():
        
        out_array = rootgrp.variables[var_name][:]
        
        rootgrp.close()
        
        return out_array
    
    else:
        
        raise ValueError("Variable '{}' not found in dataset '{}'.".format(var_name, ncdf_path))

def array_idx_list(array_list):
    
    check_list = []
    
    for array in array_list:
        
        if array.shape not in check_list:
            
            check_list.append(array.shape)
    
    if len(check_list) == 0:
        
        raise Exception("No arrays in input list.")
    
    elif len(check_list) > 1:
        
        raise Exception("Array shapes are not identical - identified {} unique shapes.".format(len(check_list)))
    # 
    else:
        
        return [i for i, j in np.ndenumerate(array_list[0])] #@UnusedVariable

def create_empty_array(array_list):
    
    check_list = []
    
    for array in array_list:
        
        if array.shape not in check_list:
            
            check_list.append(array.shape)
    
    if len(check_list) == 0:
        
        raise Exception("No arrays in input list.")
    
    elif len(check_list) > 1:
        
        raise Exception("Array shapes are not identical - identified {} unique shapes.".format(len(check_list)))
    # 
    else:
        print(check_list[0])
        return np.empty(check_list[0])
    
def add_empty_variable(nc_path, null_val, sample_array, var_name, var_datatype, var_dimensions, var_att_dict):
    
    with nc.Dataset(nc_path, 'w', format="NETCDF4") as dst:
        
        # Create variable.
        dst.createVariable(var_name, var_datatype, var_dimensions)
        
        # Copy variable attributes all at once via dictionary
        dst[var_name].setncatts(var_att_dict)
        
        # Get shape from sample array.
        array_shape = sample_array.shape
        
        # Create empty numpy array.
        array_empty = np.empty(array_shape)
        
        # Fill array with null value.
        array_fill = array_empty.fill(null_val)
        
        # Write array to netcdf.
        dst[var_name][:] = array_fill

def copy_netcdf(in_ncdf_path, out_ncdf_path, var_list): # Field dict should be in format [field_name]:[upper bound, lower bound]
    
    with nc.Dataset(in_ncdf_path, 'r', format="NETCDF4") as src, nc.Dataset(out_ncdf_path, 'w', format="NETCDF4") as dst:
        # copy global attributes all at once via dictionary
        dst.setncatts(src.__dict__)
        
        # copy dimensions
        for name, dimension in src.dimensions.items():
            #print(name)
            #print(dimension)
            #print()
            dst.createDimension(
                name, (len(dimension) if not dimension.isunlimited() else None))
        # copy all file data except for the excluded
        for name, variable in src.variables.items():
            if name in var_list:
                #print(name)
                #print(variable)
                print(variable.datatype)
                
                x = dst.createVariable(name, variable.datatype, variable.dimensions)
                
                # copy variable attributes all at once via dictionary
                dst[name].setncatts(src[name].__dict__)
                #print(src[name].__dict__)
                dst[name][:] = src[name][:]
                
def netcdf_dicts(in_ncdf_path): # Field dict should be in format [field_name]:[upper bound, lower bound]
    
    with nc.Dataset(in_ncdf_path, 'r', format="NETCDF4") as src:
        # copy global attributes all at once via dictionary
        print("GLOBAL DICT:")
        print(src.__dict__)
        print("")
        print("")
        
        # copy all file data except for the excluded
        for name, variable in src.variables.items():
            print(name)
            print(variable)
            '''
            print("VAR ATTS DICT:")
            print(src[name].__dict__)
            print("")
            '''
            
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
            
    
test_path = r"C:\GIS\CAA_2024\Test_Data\HYDRODYN-SURF_HYCOM3D-SURF_R1000_MANGASC60_20180522.nc"
test_out_path = r"C:\GIS\CAA_2024\Test_Data\test_copy.nc"

#read_nc(test_out_path)

#netcdf_dicts(test_path)

"""
copy_netcdf(test_path, test_out_path, ['lat','lon','salinity'])


#read_nc(test_path)
test_array_1 = netcdf_var_to_array(test_path, "lon")
test_array_2 = netcdf_var_to_array(test_path, "temperature")
#print(test_array_1.shape)
#test_array_3 = create_empty_array([test_array_1,test_array_2])
#print(test_array_3)



idx_list = [i for i, j in np.ndenumerate(test_array_1)]

for i in array_idx_list([test_array_1, test_array_2]):
    print(i)


for idx in idx_list:
    print(test_array_1[idx], test_array_2[idx])
"""
