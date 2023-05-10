"""
February 2022

calculate and output to netcdf files the sea ice volume
"""
import datetime
import numpy as np
import xarray as xr
import netCDF4 as nc

def sea_ice_volume(runid, endyear, endmonth, endday, startyear=2002, startmonth=1, startday=5):
    figs_path = '/project/6007519/weissgib/plotting/sea_ice/'
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
        mdl_files.append(path+"ANHA4-"+runid+"_y"+str(t.year)+"m"+str(t.month).zfill(2)+"d"+str(t.day).zfill(2)+"_icemod.nc")
    
    #also want to read in the mesh grid info
    grid_file = '/project/6007519/weissgib/plotting/data_files/anha4_files/ANHA4_mesh_mask.nc'
    mesh = nc.Dataset(grid_file)

    lons = np.array(mesh.variables['nav_lon'])
    lats = np.array(mesh.variables['nav_lat'])
    e1v = np.array(mesh.variables['e1v'])[0,:,:]
    e2u = np.array(mesh.variables['e2u'])[0,:,:]

    mesh.close()
    
    d = xr.open_mfdataset(mdl_files, concat_dim='time_counter', data_vars='minimal', coords='minimal', compat='override')
    
    x = d.dims['x']
    y = d.dims['y']
    t = d.dims['time_counter']
    
    #calculate the sea ice volume at each grid cell
    volume = np.zeros((t,y,x))
    ice_thic = d['iicethic'].values
    
    for j in range(y):
        for i in range(x):
            v = ice_thic[:,j,i]*e1v[j,i]*e2u[j,i]
            volume[:,j,i] = v

    #output to a nc file
    
    time = d.coords['time_counter'].values
    x = np.arange(x)
    y = np.arange(y)
    
    dv = xr.DataArray(volume, coords=[time,y,x], dims=['time_counter', 'y', 'x'])
    dv.to_netcdf(figs_path+'sea_ice_volume_'+runid+'.nc')
    
    dv.close()
    d.close()

if __name__ == "__main__":
    sea_ice_volume(runid='EPM151', endyear=2019, endmonth=4, endday=5)
    #sea_ice_volume(runid='EPM102', endyear=2019, endmonth=6, endday=9)
    #sea_ice_volume(runid='EPM014', endyear=2019, endmonth=8, endday=23)
    #sea_ice_volume(runid='EPM015', endyear=2019, endmonth=12, endday=31)
