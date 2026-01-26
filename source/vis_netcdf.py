'''
Created on 21 Jan 2026

@author: Karl Smith
'''

import os, json
from main import vis_ncdf

# Change as necessary.
default_json_path = r'C:\GIS\ArcPro_Projects\Visibility_Testing\json\test_cyclades.json'

# Get filepath from input.
input_path = input("Enter path to json file, then press ENTER.\nPressing ENTER without defining a path will use default path at '{}'.\n".format(default_json_path))

if input_path == '':
    input_path = default_json_path
    
# Check if path exists.
if os.path.exists(input_path):
    
    with open(input_path) as in_json:
        run_dict = json.load(in_json)
    
    run_vals = [run_dict[i] for i in run_dict.keys()]
    
    vis_ncdf(*run_vals)
        
else:
    
    print("Path '{}' does not exist. Please try again.".format(input_path))