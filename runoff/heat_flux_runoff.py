import numpy as np
import xarray as xr
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt

def yang_heat_flux():
    #use the yang definition to calculate the heat flux for different river regions
    #and plot them all on one figure

    path = '/project/6007519/ANHA4-I/RUNOFF/HydroGFD_temp/'
    fig_path = '/project/6007519/weissgib/plotting/runoff_figs/'

    #constants
    cp = 4.184
    rho = 10**6
    days_in_month = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
    months = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'June', 'July', 'Aug', 'Sept', 'Oct', 'Nov', 'Dec']

    #runoff files
    start_year = 2002
    end_year = 2016
    data = []

    for y in range(start_year, end_year+1):
        #if y == 2011: continue
        files = path+'ANHA4_ReNat_HydroGFD_HBC_runoff_monthly_y'+str(y)+'.nc'
        ds = xr.open_mfdataset(files, decode_times=False)
        reference_date = '1/1/'+str(y)
        ds['time_counter'] = pd.date_range(start=reference_date, periods=ds.sizes['time_counter'], freq='MS')
        data.append(ds)

    new_runoff = xr.concat(data, dim='time_counter')

    #model grid information for converting units
    grid_file = '/project/6007519/weissgib/plotting/data_files/anha4_files/ANHA4_mesh_mask.nc'
    mesh = xr.open_mfdataset(grid_file)

    e1v = mesh['e1v']
    e2u = mesh['e2u']

    #and the mask files for getting coastal regions
    mask_path = '/project/6000276/weissgib/model_files/runoff_temp_regions_mask.nc'

    mask_data = xr.open_mfdataset(mask_path)

    #masks = {'hb_mask': 'Hudson Bay', 'bs_mask': 'Mackenzie River Region', 'bs_east_mask': 'Eastern Bering Strait', 'laptev_mask': 'Laptev Sea'}
    #masks = {'caa_mask': 'Canadian Arctic Archipelago','kara_mask': 'Kara Sea'}
    masks = {'hb_mask': 'Hudson Bay', 'bs_mask': 'Mackenzie River Region', 'bs_east_mask': 'Eastern Bering Strait', 'laptev_mask': 'Laptev Sea', 'caa_mask': 'Canadian Arctic Archipelago','kara_mask': 'Kara Sea'}
    #masks = {'all_masks': 'All Mask Area'}

    total_heat_flux = []
    #lets make some time series over these regions
    for m in masks:
        print(m)
        md = mask_data[m][0,0]
        masked_new_runoff = new_runoff['runoff'].where(md ==2)

        e1v_mask = e1v.where(md==2)[0]
        e2u_mask = e2u.where(md==2)[0]

        new_convert = masked_new_runoff*e1v_mask*e2u_mask*0.001
        new_timeseries = new_convert.sum(('x', 'y'))
        print(new_timeseries.groupby('time_counter.month').mean().values)

        temp_mask = new_runoff.where(new_runoff['rotemper'] != -999)
        new_temp = temp_mask['rotemper'].where(md==2)
        new_temp = new_temp.where(new_convert != 0)
        #new_ts_temp = new_temp.mean(('x', 'y'))

        #now want the average for each month
        runoff_mean = new_convert.groupby('time_counter.month').mean()
        temp_mean = new_temp.groupby('time_counter.month').mean()

        #and now calculate the heat flux for each month
        heat_flux = []
        Q = []
        WT = []
        for t in range(12):
            n = days_in_month[t]
            q = runoff_mean[t]
            wt = temp_mean[t]

            hf = 86400*cp*rho*q*wt*(n/(10**12))

            #and take the sum of the heat flux over the region
            hf = hf.sum(('x', 'y')).values
            heat_flux.append(hf)
            Q.append(q)
            WT.append(wt)
        total_heat_flux.append(np.mean(heat_flux))
        plt.plot(months, heat_flux, label=masks[m])

    #print(total_heat_flux)
    #print(np.sum(total_heat_flux))
    plt.title('Average Monthly Heat Flux')
    plt.ylabel('heat flux (10^6 MJ)')
    plt.legend()
    ax = plt.gca()
    ax.yaxis.set_major_formatter(mpl.ticker.StrMethodFormatter('{x:,.0f}'))
    plt.tight_layout()
    plt.savefig(fig_path+'all_region_heat_flux.png')
    #plt.show()

    new_runoff.close()
    mask_data.close()
    mesh.close()

def whitefield_total_heat():
    #use the whitefield heat flux definition to calculate total heat flux from arctic rivers
    #does this include hudson bay rivers??

    path = '/project/6007519/ANHA4-I/RUNOFF/HydroGFD_temp/'
    fig_path = '/project/6007519/weissgib/plotting/runoff_figs/'

    #constants
    cp = 4.184
    rho = 10**6
    months = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'June', 'July', 'Aug', 'Sept', 'Oct', 'Nov', 'Dec']

    #runoff files
    start_year = 2002
    end_year = 2018
    data = []

    for y in range(start_year, end_year+1):
       files = path+'ANHA4_ReNat_HydroGFD_HBC_runoff_monthly_y'+str(y)+'.nc'
       ds = xr.open_mfdataset(files, decode_times=False)
       reference_date = '1/1/'+str(y)
       ds['time_counter'] = pd.date_range(start=reference_date, periods=ds.sizes['time_counter'], freq='MS')
       data.append(ds)

    new_runoff = xr.concat(data, dim='time_counter')

    #model grid information for converting units
    grid_file = '/project/6007519/weissgib/plotting/data_files/anha4_files/ANHA4_mesh_mask.nc'
    mesh = xr.open_mfdataset(grid_file)

    e1v = mesh['e1v']
    e2u = mesh['e2u']

    mesh.close()

    #lets calculate the total heat flux into the model from all rivers
    #including hudson bay complex for now
   
    #mask out the river temp
    temp_mask = new_runoff.where(new_runoff['rotemper'] != -999)

    #need the monthly flow amounts at the same locations
    river_flow = temp_mask['runoff']*e1v*e2u*0.001
    
    water_temp = temp_mask['rotemper'].groupby('time_counter.month').mean()
    river_flow = river_flow.groupby('time_counter.month').mean()

    #and calculate the heat flux
    heat_flux = rho*cp*water_temp*river_flow

    #and then find the total heat flux over all rivers for each month
    total_hf = heat_flux.sum(('x','y')).values

    #convert to terrawatts
    total_hf = total_hf*10**-12
    print(total_hf)

    #and lets plot it
    plt.plot(months,total_hf)
    plt.ylabel('Heat Flux (TW)')
    ax = plt.gca()
    ax.yaxis.set_major_formatter(mpl.ticker.StrMethodFormatter('{x:,.0f}'))
    plt.tight_layout()
    #plt.show()
    plt.savefig(fig_path+'whitefield_heat_flux.png')


if __name__ == "__main__":
    whitefield_total_heat()



