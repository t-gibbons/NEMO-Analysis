"""
author: Tahya Weiss-Gibbons

Take the average of the runoff from 2001 to 2018, and then scale that average by a set percentage every year afterwards
A rough way to continue to have increasing and varying runoff without continued data output
Scaling amoung comes from trend analysis in Stadnyk 2021
"""
import glob
import xarray as xr

runoff_path = "/project/6007519/weissgib/ANHA4-I/RUNOFF/HydroGFD_temp/"
runoff_files = glob.glob(runoff_path+'*.nc')

df = xr.open_mfdataset(runoff_files)

print(df)

df.close()
