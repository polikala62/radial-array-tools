'''
Created on Mar 8, 2024

@author: Karl
'''

import xarray
import netCDF4 as nc

test_path = r"C:\GIS\CAA_2024\Test_Data\HYDRODYN-SURF_HYCOM3D-SURF_R1000_MANGASC60_20180522.nc"
test_out_path = r"C:\GIS\CAA_2024\Test_Data\test_copy.nc"

dataset = xarray.open_dataset(test_out_path)

#data_array = dataset.to_dataarray(dim='salinity', name=None)

print(dataset["salinity"])
#DS.to_dataframe().to_csv(r"C:\GIS\CAA_2024\Test_Data\test_csv.csv")