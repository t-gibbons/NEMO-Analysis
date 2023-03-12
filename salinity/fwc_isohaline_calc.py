"""
fwc_isohaline_calc.py
author: Tahya Weiss-Gibbons

Calculates the freshwater content in ANHA4 down to the 34.8 isohaline, relative to 34.8

"""
import glob
import numpy as np
import netCDF4 as nc
import datetime
import xarray as xr
import matplotlib.path as mpath
import matplotlib.pyplot as plt
import pandas as pd
import cartopy.crs as ccrs
import cartopy.feature as feature

###this is for running script backend###
#import matplotlib
#matplotlib.use('Agg')
###----------------------------------###


def fwc_isohaline_calc(runid, endyear, endmonth, endday, startyear=2002, startmonth=1, startday=5):
    grid_file = '/project/6007519/weissgib/plotting/data_files/anha4_files/ANHA4_mesh_mask.nc'
    fig_path = '/project/6000276/weissgib/freshwater_content/'
    
    #lets start by reading in the model files

    #path to the folder with the model runs
    path = "/project/6000276/weissgib/ANHA4/ANHA4-"+runid+"-S/"

    mdl_files = glob.glob(path+'ANHA4-'+runid+'_gridT_*.nc')

    last_file = path+'ANHA4-ETW101_gridT_20140605-20140609.nc' #want to figure this out automatically ideally
    mdl_files.remove(last_file)
    
    d = xr.open_mfdataset(mdl_files, concat_dim='time_counter', data_vars='minimal', coords='minimal', compat='override')
    print(d)

    
    #also want to read in the mesh grid info
    mesh = nc.Dataset(grid_file)

    #lons = np.array(mesh.variables['nav_lon'])
    #lats = np.array(mesh.variables['nav_lat'])
    mask = np.array(mesh.variables['tmask'])

    mesh.close()

    #mask the data

    d.coords['mask'] = (('deptht', 'y_grid_T', 'x_grid_T'), mask[0,:,:,:])
    d = d.where(d.mask == 1)

    #annual averages
    annual_average = d
    #annual_average = d.groupby('time_counter.year').mean('time_counter')
    #annual_average = d['vosaline'].resample(time_counter='M').mean()
    full_depth = list(d['deptht'].values)

    #want to exclude the top 65m to compare with observations
    #annual_average = annual_average.where(annual_average['deptht'] > 65.0, drop=True)

    #and now we want to look at above the 34.8 isohaline
    d = xr.where(annual_average['vosaline'] < 34.8, annual_average['vosaline'], np.nan)

    #calculate the weights of the depth levels
    n = len(d['deptht'])
    zero_start = True
    weight = np.zeros(n)
    dz = np.zeros(n)
    dd = d['deptht'][n-1]
    for i in range(n):
        if zero_start:
            if i == 0:
                weight[i] = d['deptht'][i]/dd
                dz[i] = d['deptht'][i]
            else:
                weight[i] = (d['deptht'][i] - d['deptht'][i-1])/dd
                dz[i] = d['deptht'][i] - d['deptht'][i-1]
        else:
            if i == 0:
                k = full_depth.index(d['deptht'][i])
                dz[i] = d['deptht'][i] - full_depth[k-1]
            else:
                dz[i] = d['deptht'][i] - d['deptht'][i-1]

    """
    weights = xr.DataArray(weight, coords=[d['deptht']], dims=['deptht'])

    #and take the average 
    d_weighted = d.weighted(weights)
    surface_salinity = d_weighted.mean(dim='deptht', skipna=True)
    """

    #lets calculate the fwc
    #using equation fwc = sum[dz*(s-34.8)/34.8]
    #s = d['vosaline']
    tmp = (34.8-d)/34.8

    x = d.sizes['x_grid_T']
    y = d.sizes['y_grid_T']

    #get dz on the grid
    dz_grid = np.tile(dz[:,np.newaxis,np.newaxis], (1,y,x))

    t1 = tmp*dz_grid

    fwc = t1.sum(dim='deptht', skipna=True)
    
    #lets just output the fwc to a netcdf
    #fwc.to_netcdf(fig_path+runid+'_fwc_34.8_isohaline_0m.nc')

if __name__ == "__main__":
    fwc_isohaline_calc(runid='ETW101', endyear=2003, endmonth=4, endday=5)
