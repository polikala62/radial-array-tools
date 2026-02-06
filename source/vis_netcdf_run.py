'''
Created on 21 Jan 2026

@author: Karl Smith
'''

# Import libraries.
import os, json
from vis_netcdf_main import vis_ncdf

# Change as necessary.
default_json_path = r'C:\Users\Karl Smith\git\radial-array-tools\demo\vis_netcdf\troubleshooting\vis_netcdf_demo_troubleshooting_parameters_02.json'

# Get filepath from input.
input_path = input("Enter path to json file, then press ENTER.\nPressing ENTER without defining a path will use default path at '{}'.\n".format(default_json_path))

if input_path == '':
    input_path = default_json_path
    
# Check if path exists.
if os.path.exists(input_path):
    
    with open(input_path) as in_json:
        run_dict = json.load(in_json)
    
    # Get list of mandatory arguments.
    arg_list = [run_dict[i] for i in run_dict.keys() if str(i) in ['out_ncdf', 'pr_gdb', 'in_dem', 'pt_mask_fc', 'land_fc', 'xy_spacing', 'dist_range', 'densify_dist', 'obs_z_range']]
    
    # Get dictionary of optional arguments.
    kwarg_dict = {}
    for key in run_dict.keys():
        if key in ['obs_z_offset', 'sample_raster', 'landmark_fc', 'override_dist_list', 'pt_mask_json', 'write_log']:
            kwarg_dict[key] = run_dict[key]
    
    # Run script with mandatory and optional arguments.
    vis_ncdf(*arg_list, **kwarg_dict)
        
else:
    
    print("Path '{}' does not exist. Please check the path and try again.".format(input_path))
    input("Press any key to exit.")