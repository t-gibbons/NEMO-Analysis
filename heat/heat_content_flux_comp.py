import numpy as np
import xarray as xr
import pandas as pd
import matplotlib.pyplot as plt

###this is for running script backend###
import matplotlib
matplotlib.use('Agg')
###----------------------------------###

#constants
cp = 4.184
rho = 10**6
days_in_month = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
months = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'June', 'July', 'Aug', 'Sept', 'Oct', 'Nov', 'Dec']

figs_path = '/project/6007519/weissgib/plotting/heat/'

#old run
path_old = '/project/6007519/weissgib/plotting/heat/EPM161_heat_content_new.nc'
d_nt = xr.open_mfdataset(path_old, chunks={'time_counter': 100})
datetimeindex = d_nt.indexes['time_counter'].to_datetimeindex()
d_nt['time_counter'] = datetimeindex
times_old = datetimeindex.values
print(d_nt)

#new run 
path_new = '/project/6007519/weissgib/plotting/heat/ETW162_heat_content_new.nc'
d_t = xr.open_mfdataset(path_new, chunks={'time_counter': 100})
datetimeindex = d_t.indexes['time_counter'].to_datetimeindex()
d_t['time_counter'] = datetimeindex
times_new = datetimeindex.values
print(d_t)

#runnoff files
runoff_path = '/project/6007519/ANHA4-I/RUNOFF/HydroGFD_temp/'

start_year = 2002
end_year = 2018
data = []

for y in range(start_year, end_year+1):
    #if y == 2011: continue
    files = runoff_path+'ANHA4_ReNat_HydroGFD_HBC_runoff_monthly_y'+str(y)+'.nc'
    ds = xr.open_mfdataset(files, decode_times=False)
    reference_date = '1/1/'+str(y)
    ds['time_counter'] = pd.date_range(start=reference_date, periods=ds.sizes['time_counter'], freq='MS')
    data.append(ds)

runoff = xr.concat(data, dim='time_counter')
runoff = runoff.rename({'x': 'x_grid_T', 'y': 'y_grid_T'})

#and the mask files for getting coastal regions
mask_path = '/project/6007519/weissgib/plotting/data_files/anha4_files/runoff_comp_regions.nc'

mask_data = xr.open_mfdataset(mask_path)

mask_data = mask_data.rename({'x': 'x_grid_T', 'y': 'y_grid_T'})

e1t = mask_data['e1t'].sel(z=0)
e2t = mask_data['e2t'].sel(z=0)

#masks = {'full_arctic': 'Arctic'}
#masks = {'hb_mask': 'Hudson Bay', 'caa_mask': 'Canadian Arctic Archipelago', 'bs_mask': 'Bering Strait', 'bs_east_mask': 'McKenzie River Region', 'laptev_mask': 'Laptev Sea', 'kara_mask': 'Kara Sea', 'nc_mask': 'Northern Coast'}

masks = {'hb_mask': 'Hudson Bay', 'kara_mask': 'Kara Sea'}

#lets make some time series over these regions
for m in masks:
    print(m)
    md = mask_data[m][0,0]
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

    old_ts = old_ts[:l]
    times_old = times_old[:l]

    #and convert the runoff data
    masked_runoff = runoff['runoff'].where(md == 2)
    runoff_convert = masked_runoff*e1t_mask*e2t_mask*0.001

    temp_mask = runoff['rotemper'].where(runoff['rotemper'] != -999)
    temp = temp_mask.where(md == 2)
    runoff_temp = temp.where(runoff_convert != 0)

    #now we need to calculate the heat flux for each month
    heat_flux = []
    for t in runoff_convert['time_counter']:
        
        month = t.dt.month.values
        n = days_in_month[month-1]

        q = runoff_convert.sel(time_counter=t)
        wt = runoff_temp.sel(time_counter=t)
        
        hf = 86400*cp*rho*q*wt*n

        hf = hf.sum(('x_grid_T', 'y_grid_T')).values
        heat_flux.append(hf)

    diff = (new_ts-old_ts)

    #find the correlation between the two curves
    diff_month = diff.resample(time_counter='M').mean().values
    
    corr = np.corrcoef(diff_month, heat_flux)
    print(corr)
    continue

    fig, ax1 = plt.subplots()

    #and plot the difference heat content
    
    diff.plot(x='time_counter', color='tab:blue')
    ax1.set_ylabel('Difference in Heat Content (J)', color='tab:blue')

    #same plot, other axis heat flux
    ax2 = ax1.twinx()
    times = runoff_convert['time_counter'].values
    ax2.plot(times, heat_flux, label='heat flux', color='tab:green')
    ax2.set_ylabel('Heat Flux from Runoff (J)', color='tab:green')

    fig.tight_layout()
    fig.suptitle(masks[m])
    #plt.show()
    plt.savefig(figs_path+m+'_heat_content_flux_comp_lim3.png')
    plt.clf()

d_nt.close()
d_t.close()
mask_data.close()
    
