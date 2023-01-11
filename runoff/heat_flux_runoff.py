import numpy as np
import xarray as xr
import pandas as pd

path = '/project/6007519/ANHA4-I/RUNOFF/HydroGFD_temp/'

#constants
cp = 4.184
rho = 10**6
days_in_month = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]

#runoff files
start_year = 2002
end_year = 2017
data = []

for y in range(start_year, end_year+1):
    files = path+'ANHA4_ReNat_HydroGFD_HBC_runoff_monthly_y'+str(y)+'.nc'
    ds = xr.open_mfdataset(files, decode_times=False)
    reference_date = '1/1/'+str(y)
    ds['time_counter'] = pd.date_range(start=reference_date, periods=ds.sizes['time_counter'], freq='MS')
    data.append(ds)

new_runoff = xr.concat(data, dim='time_counter')

#model grid information for converting units
grid_file = '/project/6007519/weissgib/plotting/ANHA4_mesh_mask.nc'
mesh = xr.open_mfdataset(grid_file)

e1v = mesh['e1v']
e2u = mesh['e2u']

#and the mask files for getting coastal regions
mask_path = '/project/6007519/weissgib/plotting/regions_mask.nc'

mask_data = xr.open_mfdataset(mask_path)

#masks = {'ss_mask': 'Siberian Shelf', 'cs_mask': 'Canadian Shelf'}
masks = {'hb_mask': 'Hudson Bay and Greenland', 'caa_mask': 'Canadian Arctic Archipelago', 'bs_mask': 'Bering Strait', 'ec_mask': 'Eastern Coast', 'nc_mask': 'Norway and British Isles'}

#lets make some time series over these regions
for m in masks:
    md = mask_data[m][0,0]
    masked_new_runoff = new_runoff['runoff'].where(md ==2)

    e1v_mask = e1v.where(md==2)[0]
    e2u_mask = e2u.where(md==2)[0]

    new_convert = masked_new_runoff*e1v_mask*e2u_mask*0.001
    new_timeseries = new_convert.sum(('x', 'y'))

    new_temp = new_runoff['rotemper'].where(md==2)
    new_ts_temp = new_temp.sum(('x', 'y'))

    #now want the average for each month
    runoff_mean = new_timeseries.groupby('time_counter.month').mean()
    temp_mean = new_ts_temp.groupby('time_counter.month').mean()

    #and now calculate the heat flux for each month
    heat_flux = []
    for m in range(12):
        n = days_in_month[m]
        q = runoff_mean[m]
        wt = temp_mean[m]

        hf = 86400*cp*rho*q*wt*(n/(10**12))
        heat_flux.append(hf)
    print(heat_flux)

new_runoff.close()
mask_data.close()
mesh.close()
