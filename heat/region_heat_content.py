import numpy as np
import xarray as xr
import matplotlib.pyplot as plt

###this is for running script backend###
import matplotlib
matplotlib.use('Agg')
###----------------------------------###

figs_path = '/project/6007519/weissgib/plotting/heat/'

#old run
path_old = '/project/6007519/weissgib/plotting/heat/EPM151_heat_content.nc'
d_nt = xr.open_mfdataset(path_old)
datetimeindex = d_nt.indexes['time_counter'].to_datetimeindex()
d_nt['time_counter'] = datetimeindex
times_old = datetimeindex.values

#new run 
path_new = '/project/6007519/weissgib/plotting/heat/ETW161_heat_content.nc'
d_t = xr.open_mfdataset(path_new)
datetimeindex = d_t.indexes['time_counter'].to_datetimeindex()
d_t['time_counter'] = datetimeindex
times_new = datetimeindex.values

#and the mask files for getting coastal regions
mask_path = '/project/6007519/weissgib/plotting/data_files/anha4_files/runoff_comp_regions.nc'

mask_data = xr.open_mfdataset(mask_path)

mask_data = mask_data.rename({'x': 'x_grid_T', 'y': 'y_grid_T'})

#masks = {'full_arctic': 'Arctic'}
masks = {'hb_mask': 'Hudson Bay', 'caa_mask': 'Canadian Arctic Archipelago', 'bs_mask': 'Bering Strait', 'bs_east_mask': 'McKenzie River Region', 'laptev_mask': 'Laptev Sea', 'kara_mask': 'Kara Sea', 'nc_mask': 'Northern Coast'}

#lets make some time series over these regions
for m in masks:
    print(m)
    md = mask_data[m][0,0]
    masked_new = d_t['votemper'].where(md ==2)
    masked_old = d_nt['votemper'].where(md ==2)

    new_ts = masked_new.sum(('x_grid_T', 'y_grid_T'))
    old_ts = masked_old.sum(('x_grid_T', 'y_grid_T'))

    l = new_ts.shape[0]
    print(new_ts)
    print(old_ts)

    old_ts = old_ts[:l]
    times_old = times_old[:l]
    print(new_ts.shape)
    print(old_ts.shape)

    #and plot
    diff = (new_ts-old_ts)
    print(diff)
    #plt.plot(times_new, new_ts, label='ETW161')
    #plt.plot(times_old, old_ts, label='EPM151')

    print(diff.shape)
    print(times_old.shape)
    
    diff.plot(x='time_counter') 
    #plt.legend()
    plt.title(masks[m])
    plt.ylabel('Difference in Heat Content')
    plt.tight_layout()
    plt.savefig(figs_path+m+'_heat_content_change_lim2.png')
    plt.clf()

d_nt.close()
d_t.close()
mask_data.close()
    
