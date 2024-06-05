"""
author: Tahya Weiss-Gibbons

Take the runoff from 2018, and then scale every year by a set amount
A rough way to continue to have increasing and varying runoff without continued data output
Scaling amoung comes from trend analysis in Stadnyk 2021
"""

import xarray as xr
import pandas as pd

hype_path = '/project/6007519/ANHA4-I/RUNOFF/HydroGFD/ANHA4_ReNat_HydroGFD_HBC_runoff_monthly_y2018.nc'

ds = xr.open_mfdataset(hype_path, decode_times=False)
reference_date = '1/1/2018'
runoff = ds['runoff']
ds['time_counter'] = pd.date_range(start=reference_date, periods=ds.sizes['time_counter'], freq='MS')

print(ds)

start_year = 2019
end_year = 2022

scaling_factor = 0.03 #from stadnyk (2021), scaling over the entire area by 3% per year

for y in range(start_year, end_year+1):

   ref_date = '1/1/'+str(y)
   times = pd.date_range(start=ref_date, periods=ds.sizes['time_counter'], freq='MS')
   runoff = runoff*0.03
   print(runoff)
   ns = ds.copy(deep=True)
   ns.assign_coords(time_counter=times)
   ns.time_counter.encoding['units'] = "days since Jan-1-0000 00:00:00"
   ns.time_counter.encoding['calendar'] = 'noleap'
   ns['runoff'] = runoff

   print(ns.time_counter)
   ns.to_netcdf('ANHA4_ReNat_HydroGFD_HBC_runoff_monthly_y'+str(y)+'.nc')
