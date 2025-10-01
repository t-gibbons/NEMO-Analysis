import numpy as np
import xarray as xr
import matplotlib.pyplot as plt

###this is for running script backend###
#import matplotlib
#matplotlib.use('Agg')
###----------------------------------###

figs_path = '/project/6007519/weissgib/plotting/heat/'

#old run
path_old = '/project/6007519/weissgib/plotting/heat/EPM151_heat_content.nc'
d_nt = xr.open_mfdataset(path_old, chunks={'time_counter': 100})
datetimeindex = d_nt.indexes['time_counter'].to_datetimeindex()
d_nt['time_counter'] = datetimeindex
times_old = datetimeindex.values
print(d_nt)

#new run 
path_new = '/project/6007519/weissgib/plotting/heat/ETW161_heat_content.nc'
d_t = xr.open_mfdataset(path_new, chunks={'time_counter': 100})
datetimeindex = d_t.indexes['time_counter'].to_datetimeindex()
d_t['time_counter'] = datetimeindex
times_new = datetimeindex.values
print(d_t)

#and the mask files for getting coastal regions
mask_path = '/project/6007519/weissgib/plotting/data_files/anha4_files/runoff_comp_regions.nc'

mask_data = xr.open_mfdataset(mask_path)

mask_data = mask_data.rename({'x': 'x_grid_T', 'y': 'y_grid_T'})

e1t = mask_data['e1t'].sel(z=0)
e2t = mask_data['e2t'].sel(z=0)
print(e2t)

masks = {'full_arctic': 'Arctic'}
#masks = {'hb_mask': 'Hudson Bay', 'caa_mask': 'Canadian Arctic Archipelago', 'bs_mask': 'Bering Strait', 'bs_east_mask': 'McKenzie River Region', 'laptev_mask': 'Laptev Sea', 'kara_mask': 'Kara Sea', 'nc_mask': 'Northern Coast'}

#lets make some time series over these regions
for m in masks:
    print(m)
    md = mask_data[m][0,0]
    print(md)
    masked_new = d_t['votemper'].where(md ==2)
    masked_old = d_nt['votemper'].where(md ==2)

    #need to take into account the size of each cell for the sum
    e1t_mask = e1t.where(md == 2)[0]
    e2t_mask = e2t.where(md == 2)[0]
    new_convert = masked_new*e1t_mask*e2t_mask
    old_convert = masked_old*e1t_mask*e2t_mask

    new_ts = new_convert.sum(('x_grid_T', 'y_grid_T'))
    old_ts = old_convert.sum(('x_grid_T', 'y_grid_T'))

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

    diff.to_netcdf(figs_path+m+'_diff_heat_lim2.nc')
    exit()
    
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
    
