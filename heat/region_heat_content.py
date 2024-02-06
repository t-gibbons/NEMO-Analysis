import numpy as np
import xarray as xr
import matplotlib.pyplot as plt

#old run
path_old = '/mnt/storage6/tahya/model_files/EPM161_heat_content.nc'
d_nt = xr.open_mfdataset(path_old)
datetimeindex = d_nt.indexes['time_counter'].to_datetimeindex()
times_old = datetimeindex.values

#new run 
path_new = '/mnt/storage6/tahya/model_files/ETW162_heat_content.nc'
d_t = xr.open_mfdataset(path_new)
datetimeindex = d_t.indexes['time_counter'].to_datetimeindex()
times_new = datetimeindex.values

#and the mask files for getting coastal regions
mask_path = '/mnt/storage4/tahya/runoff/runoff_comp_regions.nc'

mask_data = xr.open_mfdataset(mask_path)

mask_data = mask_data.rename({'x': 'x_grid_T', 'y': 'y_grid_T'})

masks = {'full_arctic': 'Arctic'}
#masks = {'hb_mask': 'Hudson Bay', 'caa_mask': 'Canadian Arctic Archipelago', 'bs_mask': 'Bering Strait', 'bs_east_mask': 'McKenzie River Region', 'laptev_mask': 'Laptev Sea', 'kara_mask': 'Kara Sea', 'nc_mask': 'Northern Coast'}

#lets make some time series over these regions
for m in masks:
    print(m)
    md = mask_data[m][0,0]
    masked_new = d_t['votemper'].where(md ==2)
    masked_old = d_nt['votemper'].where(md ==2)

    new_ts = masked_new.sum(('x_grid_T', 'y_grid_T'))
    old_ts = masked_old.sum(('x_grid_T', 'y_grid_T'))

    print(new_ts.shape)
    print(old_ts.shape)

    old_ts = old_ts[:-1]

    #and plot
    diff = ((new_ts-old_ts)/old_ts)*100
    #plt.plot(times_new, new_ts, label='ETW161')
    #plt.plot(times_old, old_ts, label='EPM151')

    print(diff.shape)
    print(times_old.shape)
    
    plt.plot(times_old, diff)
    #plt.legend()
    plt.title(masks[m])
    plt.savefig(m+'_heat_content_change_lim3.png')
    plt.clf()

d_nt.close()
d_t.close()
mask_data.close()
    
