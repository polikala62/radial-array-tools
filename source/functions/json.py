'''
Created on 24 Mar 2025

@author: Karl Smith
'''
import os, json

def write_json(in_dict, in_path):
    
    if os.path.exists(in_path):
        os.remove(in_path) 
        
    out_file = open(in_path, "w")
    json.dump(in_dict, out_file)
    out_file.close()
    
def load_json(in_path):
    out_file = open(in_path)
    out_dict = json.load(out_file)
    out_file.close()
    return out_dict