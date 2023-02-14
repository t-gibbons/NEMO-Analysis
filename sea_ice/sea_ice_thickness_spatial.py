"""
February 2022
Updated July 2022

plot the yearly average sea ice thickness from the model files
"""
import glob
import datetime
import numpy as np
import xarray as xr
import netCDF4 as nc
import matplotlib.pyplot as plt
import matplotlib.path as mpath
import cartopy.crs as ccrs
import cartopy.feature as feature

###this is for running script backend###
import matplotlib
matplotlib.use('Agg')
###----------------------------------###

def sea_ice_thickness_spatial(runid, endyear, endmonth, endday, startyear=2002, startmonth=1, startday=5):
    figs_path = '/project/6000276/weissgib/sea_ice/'
    path = "/project/6000276/weissgib/ANHA4/ANHA4-"+runid+"-S/"

    mdl_files = glob.glob(path+'ANHA4-'+runid+'_icemod_*.nc')

    #also want to read in the mesh grid info
    grid_file = '/project/6007519/weissgib/plotting/data_files/anha4_files/ANHA4_mesh_mask.nc'
    mesh = nc.Dataset(grid_file)

    lons = np.array(mesh.variables['nav_lon'])
    lats = np.array(mesh.variables['nav_lat'])

    mesh.close()
    
    d = xr.open_mfdataset(mdl_files, concat_dim='time_counter', data_vars='minimal', coords='minimal', compat='override')
    
    #start with an annual average
    annual_avg = d.groupby('time_counter.year').mean('time_counter')
    yrs = annual_avg['year'].values
    
    for y in range(annual_avg.dims['year']):
        sea_ice = annual_avg['iicethic'].isel(year=y).values
        yr = yrs[y]

        #north pole stero projection
        land_50m = feature.NaturalEarthFeature('physical', 'land', '50m', edgecolor='black', facecolor='gray', linewidth=0.5)
        projection=ccrs.NorthPolarStereo()

        fig = plt.figure(figsize=(10, 9))
        ax = plt.subplot(1, 1, 1, projection=projection)

        ax.set_extent([-280, 80, 80, 35], crs=ccrs.PlateCarree())
        ax.add_feature(land_50m, color=[0.8, 0.8, 0.8])
        ax.coastlines(resolution='50m')

        # Compute a circle in axes coordinates, which we can use as a boundary
        # for the map. We can pan/zoom as much as we like - the boundary will be
        # permanently circular.
        theta = np.linspace(0, 2*np.pi, 100)
        center, radius = [0.5, 0.5], 0.5
        verts = np.vstack([np.sin(theta), np.cos(theta)]).T
        circle = mpath.Path(verts * radius + center)
        ax.set_boundary(circle, transform=ax.transAxes)
        p1 = ax.pcolormesh(lons, lats, sea_ice, transform=ccrs.PlateCarree(), cmap='gist_ncar', vmin=0, vmax=5)
        ax_cb = plt.axes([0.92, 0.25, 0.015, 0.5])
        cb = plt.colorbar(p1,cax=ax_cb, orientation='vertical')
        cb.ax.set_ylabel('Sea Ice Thickness (m)')
        ax.gridlines()
        plt.savefig(figs_path+'sea_ice_thickness_'+runid+'_'+str(yr)+'.png')
        #plt.show()
        plt.clf()

    d.close()

def sea_ice_thickness_diff(endyear, endmonth, endday, startyear=2002, startmonth=1, startday=5):

    figs_path = '/project/6007519/weissgib/plotting/sea_ice/'

    grid_file = '/project/6007519/weissgib/plotting/data_files/anha4_files/ANHA4_mesh_mask.nc'
    mesh = nc.Dataset(grid_file)

    lons = np.array(mesh.variables['nav_lon'])
    lats = np.array(mesh.variables['nav_lat'])

    mesh.close()

    #read temp and no temp run
    path_temp = "/project/6007519/weissgib/ANHA4/ANHA4-ETW151-S/"

    path = "/project/6007519/pmyers/ANHA4/ANHA4-EPM151-S/"
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
    mdl_files_temp = []
    mdl_files_notemp = []
    for t in times:
        mdl_files_temp.append(path_temp+"ANHA4-ETW151_y"+str(t.year)+"m"+str(t.month).zfill(2)+"d"+str(t.day).zfill(2)+"_icemod.nc")
        mdl_files_notemp.append(path+"ANHA4-EPM151_y"+str(t.year)+"m"+str(t.month).zfill(2)+"d"+str(t.day).zfill(2)+"_icemod.nc")

    d_temp = xr.open_mfdataset(mdl_files_temp, concat_dim='time_counter', data_vars='minimal', coords='minimal', compat='override')
    d_no_temp = xr.open_mfdataset(mdl_files_notemp, concat_dim='time_counter', data_vars='minimal', coords='minimal', compat='override')

    no_temp_avg = d_no_temp['iicethic'].resample(time_counter='Q-NOV').mean()
    temp_avg = d_temp['iicethic'].resample(time_counter='Q-NOV').mean()

    print(temp_avg)
    print(no_temp_avg)

    times = temp_avg['time_counter'].values

    for t in times:
        st = t.strftime('%Y-%m-%d')
        print(st)
        if t.month == 2: season = 'DJF'
        if t.month == 5: season = 'MAM'
        if t.month == 8: season = 'JJA'
        if t.month == 11: season = 'SON'

        d1 = temp_avg.sel(time_counter=st).values
        d2 = no_temp_avg.sel(time_counter=st).values
        diff = d1-d2
        diff = diff[0,:,:]

        #north pole stero projection
        land_50m = feature.NaturalEarthFeature('physical', 'land', '50m', edgecolor='black', facecolor='gray', linewidth=0.5)
        projection=ccrs.NorthPolarStereo()

        fig = plt.figure(figsize=(10, 9))
        ax = plt.subplot(1, 1, 1, projection=projection)

        #ax.set_extent([-280, 80, 80, 35], crs=ccrs.PlateCarree())
        ax.set_extent([-180,180,90,60], crs=ccrs.PlateCarree())
        ax.add_feature(land_50m, color=[0.8, 0.8, 0.8])
        ax.coastlines(resolution='50m')

        # Compute a circle in axes coordinates, which we can use as a boundary
        # for the map. We can pan/zoom as much as we like - the boundary will be
        # permanently circular.
        theta = np.linspace(0, 2*np.pi, 100)
        center, radius = [0.5, 0.5], 0.5
        verts = np.vstack([np.sin(theta), np.cos(theta)]).T
        circle = mpath.Path(verts * radius + center)
        ax.set_boundary(circle, transform=ax.transAxes)
        p1 = ax.pcolormesh(lons, lats, diff, transform=ccrs.PlateCarree(), cmap='bwr', vmin=-0.5, vmax=0.5)
        ax_cb = plt.axes([0.92, 0.25, 0.015, 0.5])
        cb = plt.colorbar(p1, cax=ax_cb, orientation='vertical')
        cb.ax.set_ylabel('Sea Ice Thickness (m)')
        ax.gridlines()
        plt.savefig(figs_path+'sea_ice_thickness_diff_'+season+'_'+st+'.png')
        #plt.show()
        plt.clf()

    d_temp.close()
    d_no_temp.close()


if __name__ == "__main__":
    sea_ice_thickness_diff(2012, 12 ,31)
