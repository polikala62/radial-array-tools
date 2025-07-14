'''
Created on Mar 12, 2024

@author: Karl
'''
import arcpy

land_fc = r"C:\GIS\CAA_2024\Test_Data\test_coast_poly_slct.shp"
land_poly = []
with arcpy.da.SearchCursor(land_fc, ["SHAPE@"]) as cursor: #@UndefinedVariableFromImport
    for row in cursor:
        for part in row:
            land_poly.append(part)
            
check_coord = arcpy.PointGeometry(arcpy.Point(3555000,3286000,1000))

for i in land_poly:
    print(check_coord.distanceTo(i))