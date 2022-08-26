"""
fwc_depth_calc.py
author: Tahya Weiss-Gibbons

Calculate the freshwater content for ANHA4 down to a set depth level, relative to 34.8
Output calculated values as a netCDF file
"""
import numpy as np
import netCDF4 as nc
import datetime
import xarray as xr
import pandas as pd

def fwc_depth_calc(runid, endyear, endmonth, endday, startyear=2002, startmonth=1, startday=5):
    output_path = '/project/6007519/weissgib/plotting/data_files/freshwater_content/'
    #path to the folder with the model runs
    path = "/project/6007519/pmyers/ANHA4/ANHA4-"+runid+"-S/"

    start_time = datetime.date(startyear, startmonth, startday)

    end_time = datetime.date(endyear, endmonth, endday)

    #figure out all the dates we have model files
    delta = end_time - start_time
    times = []

    i = 0
    while i < delta.days+1:
        t = start_time + datetime.timedelta(days=i)
        if t.month == 2 and t.day == 29:
            t = datetime.date(t.year, 3, 1)
            i = i+6
        else:
            i = i+5
        times.append(t)

    #and now make a list of model files to read
    mdl_files = []
    for t in times:
        mdl_files.append(path+"ANHA4-"+runid+"_y"+str(t.year)+"m"+str(t.month).zfill(2)+"d"+str(t.day).zfill(2)+"_gridT.nc")

    d = xr.open_mfdataset(mdl_files, concat_dim='time_counter', data_vars='minimal', coords='minimal', compat='override')

    #also want to read in the mesh grid info
    grid_file = '/project/6007519/weissgib/plotting/ANHA4_mesh_mask.nc'
    mesh = nc.Dataset(grid_file)

    mask = np.array(mesh.variables['tmask'])

    mesh.close()

    #mask the data

    d.coords['mask'] = (('deptht', 'y_grid_T', 'x_grid_T'), mask[0,:,:,:])
    d = d.where(d.mask == 1)

    #want to do annual averages
    #annual_average = d.groupby('time_counter.month').mean('time_counter')
    annual_average = d['vosaline'].resample(time_counter='M').mean()

    #and now we want to average over the top 200m of the water column
    d = annual_average.where(annual_average['deptht'] < 200.0, drop=True)
    print(d)

    #calculate the weights
    n = len(d['deptht'])
    weight = np.zeros(n)
    dz = np.zeros(n)
    dd = d['deptht'][n-1]
    for i in range(n):
        if i == 0:
            weight[i] = d['deptht'][i]/dd
            dz[i] = d['deptht'][i]
        else:
            weight[i] = (d['deptht'][i] - d['deptht'][i-1])/dd
            dz[i] = d['deptht'][i] - d['deptht'][i-1]

    weights = xr.DataArray(weight, coords=[d['deptht']], dims=['deptht'])

    #and take the average
    d_weighted = d.weighted(weights)
    surface_salinity = d_weighted.mean(dim='deptht', skipna=True)

    #lets also try calculating the fwc
    #using equation fwc = sum[dz*(s-34.8)/34.8]
    #s = d['vosaline']
    tmp = (34.8-d)/34.8

    #get dz on the grid
    dz_grid = np.tile(dz[:,np.newaxis,np.newaxis], (1,800,544))

    t1 = tmp*dz_grid

    fwc = t1.sum(dim='deptht', skipna=True)

    #lets just output the fwc to a netcdf
    fwc.to_netcdf(output_path+runid+'_monthly_avg_fwc_200m.nc')

if __name__ == "__main__":
    fwc_depth_calc(runid='EPM101', endyear=2019, endmonth=4, endday=5)
    fwc_depth_calc(runid='EPM102', endyear=2019, endmonth=6, endday=9)
    fwc_depth_calc(runid='EPM014', endyear=2019, endmonth=8, endday=23)
    fwc_depth_calc(runid='EPM015', endyear=2019, endmonth=12, endday=31)
