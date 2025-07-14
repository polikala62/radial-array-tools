'''
Created on Dec 13, 2023

@author: Karl
'''

import arcpy, os

def output_pts(out_pt_dict, out_fc, out_crs):
    
    # Set processing variables.
    arcpy.env.overwriteOutput = True
    
    #------------------------------------------------------------------------------ 
    
        
    #------------------------------------------------------------------------------ 
    # Create output and add fields.
    
    out_dir, out_filename = os.path.split(out_fc)
    
    try:
        
        # Create output features.
        arcpy.management.CreateFeatureclass(out_dir, out_filename, geometry_type="POINT", has_m="DISABLED", has_z="ENABLED", spatial_reference=out_crs)
        arcpy.management.AddField(out_fc, "PT_ID", "TEXT")
        arcpy.management.AddField(out_fc, "PT_X", "DOUBLE")
        arcpy.management.AddField(out_fc, "PT_Y", "DOUBLE")
        arcpy.management.AddField(out_fc, "PT_Z", "DOUBLE")
        arcpy.management.AddField(out_fc, "SUB_VAL", "DOUBLE")
        arcpy.management.AddField(out_fc, "VIS_VAL", "LONG")
        
    except:
        
        print("Could not create point features!!")
        
    # Use cursor to add values.
    with arcpy.da.InsertCursor(out_fc, ["SHAPE@", "PT_ID", "PT_X", "PT_Y", "PT_Z", "SUB_VAL", "VIS_VAL"]) as cursor: #@UndefinedVariableFromImport
        
        for key in out_pt_dict.keys():
            
            # Get row from dictionary.
            iter_row = [key] + out_pt_dict[key]
            
            # Create single list of points for export.
            out_pt = arcpy.Point(iter_row[1], iter_row[2], iter_row[3])
            out_pt_geometry = arcpy.PointGeometry(out_pt, spatial_reference=out_crs)
            
            # Insert row.
            cursor.insertRow([out_pt_geometry] + iter_row)
            
#------------------------------------------------------------------------------ 

def output_skyline(pt_dict, out_fc, out_crs):
    
    # Set processing variables.
    arcpy.env.overwriteOutput = True
    
    out_dir, out_filename = os.path.split(out_fc)
    
    try:
        
        # Create output features.
        arcpy.management.CreateFeatureclass(out_dir, out_filename, geometry_type="POLYLINE", has_m="DISABLED", has_z="ENABLED", spatial_reference=out_crs)
        
    except:
        
        print("Could not create skyline features!!")
    
    pass
    
    # Use cursor to add values.
    with arcpy.da.InsertCursor(out_fc, ['SHAPE@']) as cursor: #@UndefinedVariableFromImport
        
        # Create list of point objects.
        pt_list = []
        
        for key in pt_dict.keys():
            
            pt_x, pt_y, pt_z = pt_dict[key][0:3]
            
            pt_list.append(arcpy.Point(pt_x, pt_y, pt_z))
    
        # Create a polyline geometry
        array = arcpy.Array(pt_list)
        polyline = arcpy.Polyline(array)
        
        cursor.insertRow([polyline])