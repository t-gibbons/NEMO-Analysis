import numpy as np
import xarray as xr
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.path as mpath
import cartopy.crs as ccrs
import cartopy.feature as feature

#new runoff
path = '/project/6007519/ANHA4-I/RUNOFF/HydroGFD_temp/'
fig_path = '/project/6007519/weissgib/plotting/runoff_figs/'

start_year = 2002
end_year = 2019
data = []
months = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'June', 'July', 'Aug', 'Sept', 'Oct', 'Nov', 'Dec']

for y in range(start_year, end_year):
    #if y == 2011: continue
    files = path+'ANHA4_ReNat_HydroGFD_HBC_runoff_monthly_y'+str(y)+'.nc'
    print(files)
    ds = xr.open_mfdataset(files, decode_times=False)
    reference_date = '1/1/'+str(y)
    v = ds['rotemper'].values
    ds['time_counter'] = pd.date_range(start=reference_date, periods=ds.sizes['time_counter'], freq='MS')
    data.append(ds)

new_runoff = xr.concat(data, dim='time_counter')
new_runoff = new_runoff.where(new_runoff['rotemper'] > -999)
new_runoff = new_runoff.where(new_runoff['runoff'] != 0)

#and the mask files for getting coastal regions
mask_path = '/project/6000276/weissgib/model_files/runoff_temp_regions_mask.nc'

mask_data = xr.open_mfdataset(mask_path)

#masks = {'bs_mask': 'Mackenzie River Region'}
masks = {'hb_mask': 'Hudson Bay', 'caa_mask': 'Canadian Arctic Archipelago', 'bs_mask': 'Mackenzie River Region', 'bs_east_mask': 'Eastern Bering Strait', 'kara_mask': 'Kara Sea', 'laptev_mask': 'Laptev Sea'}

#lets make time series of the average monthly runoff
for m in masks:
    md = mask_data[m][0,0]
    masked_runoff = new_runoff['rotemper'].where(md ==2)
    timeseries = masked_runoff.mean(('x', 'y'))

    hype_mean = timeseries.groupby('time_counter.month').mean()

    print(hype_mean)

    nt = hype_mean.values
    
    plt.plot(months, nt)
    x1,x2,y1,y2 = plt.axis()  
    plt.axis((x1,x2,-0.2,10))

    plt.title(masks[m])
    plt.ylabel('temperature (C)')
    #plt.show()
    plt.savefig(fig_path+m+'_monthly_average_runoff_temp.png')
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
