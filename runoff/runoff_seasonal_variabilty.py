import numpy as np
import xarray as xr
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.path as mpath
import cartopy.crs as ccrs
import cartopy.feature as feature

#new runoff
path = '/project/6007519/ANHA4-I/RUNOFF/HydroGFD/'

start_year = 2002
end_year = 2019
data = []

for y in range(start_year, end_year+1):
    files = path+'ANHA4_ReNat_HydroGFD_HBC_runoff_monthly_y'+str(y)+'.nc'
    ds = xr.open_mfdataset(files, decode_times=False)
    reference_date = '1/1/'+str(y)
    v = ds['runoff'].values
    ds['time_counter'] = pd.date_range(start=reference_date, periods=ds.sizes['time_counter'], freq='MS')
    data.append(ds)

new_runoff = xr.concat(data, dim='time_counter')

#old runoff
old_path = '/project/6007519/pmyers/ANHA4-I/RUNOFF/Bamber2012/'
files = []

for y in range(start_year, end_year+1):
    files.append(old_path+'ANHA4_runoff_monthly_combined_Dai_Trenberth_Bamber_y'+str(y)+'.nc')

old_runoff = xr.open_mfdataset(files)

#and the mask files for getting coastal regions
mask_path = '/project/6007519/weissgib/plotting/data_files/anha4_files/runoff_temp_regions_mask.nc'
#mask_path = '/project/6007519/weissgib/plotting/regions_mask.nc'

mask_data = xr.open_mfdataset(mask_path)

masks = {'bs_mask': 'Mackenzie Region', 'kara_mask': 'Kara Sea', 'laptev_mask': 'Laptev Sea', 'bs_east_mask': 'Eastern Bering Strait'}
#masks = {'hb_mask': 'Hudson Bay and Greenland', 'caa_mask': 'Canadian Arctic Archipelago', 'bs_mask': 'Bering Strait', 'ec_mask': 'Eastern Coast', 'nc_mask': 'Norway and British Isles'}

#lets make time series of the average monthly runoff
for m in masks:
    md = mask_data[m][0,0]
    masked_new_runoff = new_runoff['runoff'].where(md ==2)
    masked_old_runoff = old_runoff['runoff'].where(md ==2)
    new_timeseries = masked_new_runoff.sum(('x', 'y'))
    old_timeseries = masked_old_runoff.sum(('x', 'y'))

    dai_mean = old_timeseries.groupby('time_counter.month').mean()
    hype_mean = new_timeseries.groupby('time_counter.month').mean()

    print(dai_mean)
    print(hype_mean)

    nt = hype_mean.values
    ot = dai_mean.values

    dn = hype_mean['month'].values
    
    plt.plot(dn, nt, label='HYPE')
    plt.plot(dn, ot, label='Dai and Trenberth')
    plt.legend()
    plt.title(masks[m])
    plt.ylabel('runoff (kg/s/m^2)')
    #plt.show()
    plt.savefig('/project/6007519/weissgib/plotting/figs/runoff/'+m+'_monthly_average_runoff_comp.png')
    plt.clf()

exit()
#and lets make the same plot but for the whole region
full_new_timeseries = new_runoff['runoff'].sum(('x', 'y'))
full_old_timeseries = old_runoff['runoff'].sum(('x', 'y'))

full_new_timeseries = full_new_timeseries.groupby('time_counter.month').mean()
full_old_timeseries = full_old_timeseries.groupby('time_counter.month').mean()

nt = full_new_timeseries.values
ot = full_old_timeseries.values

dn = full_new_timeseries['month'].values

plt.plot(dn, nt, label='HYPE')
plt.plot(dn, ot, label='Dai and Trenberth')
plt.legend()
plt.ylabel('runoff (kg/s/m^2)')

plt.savefig('/project/6007519/weissgib/plotting/figs/runoff/total_runoff_seasonal_comp.png')
plt.clf()

new_runoff.close()
old_runoff.close()
mask_data.close()
