"""
April 2023

plot the difference in the heat content for each season
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
#import matplotlib
#matplotlib.use('Agg')
###----------------------------------###

def heat_cont_diff(endyear, endmonth, endday, temp_id, notemp_id, startyear=2004, startmonth=1, startday=5):

    figs_path = '/project/6007519/weissgib/plotting/heat/'

    grid_file = '/project/6007519/weissgib/plotting/data_files/anha4_files/ANHA4_mesh_mask.nc'
    mesh = nc.Dataset(grid_file)

    lons = np.array(mesh.variables['nav_lon'])
    lats = np.array(mesh.variables['nav_lat'])

    mesh.close()

    #read temp and no temp run
    path = "/project/6007519/weissgib/plotting/heat/"

    mdl_files_temp = figs_path+temp_id+'_heat_content.nc'
    mdl_files_notemp= figs_path+notemp_id+'_heat_content.nc'

    print('trying to read files')
    d_temp = xr.open_mfdataset(mdl_files_temp)
    d_no_temp = xr.open_mfdataset(mdl_files_notemp)

    #no_temp_avg = d_no_temp['votemper'].resample(time_counter='Q-NOV').mean()
    #temp_avg = d_temp['votemper'].resample(time_counter='Q-NOV').mean()

    no_temp_avg = d_no_temp.groupby('time_counter.year').mean('time_counter')
    temp_avg = d_temp.groupby('time_counter.year').mean('time_counter')

    print(temp_avg)
    print(no_temp_avg)

    times = temp_avg['year'].values

    #for t in times:
    for y in range(temp_avg.dims['year']):

        #st = t.strftime('%Y-%m-%d')
        #print(st)
        t = times[y]
        print(t)
        #if t.month == 2: season = 'DJF'
        #if t.month == 5: season = 'MAM'
        #if t.month == 8: season = 'JJA'
        #if t.month == 11: season = 'SON'

        #d1 = temp_avg.sel(time_counter=st).values
        #d2 = no_temp_avg.sel(time_counter=st).values
        d1 = temp_avg.isel(year=y).values
        d2 = no_temp_avg.isel(year=y).values
        print(d1)
        diff = ((d2-d1)/d1)*100
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
        
        """
        #hudson bay projection
        land_50m = feature.NaturalEarthFeature('physical', 'land', '50m', edgecolor='black', facecolor='gray', linewidth=0.5)
        projection=projection=ccrs.Mercator(central_longitude=-80)

        fig = plt.figure(figsize=(10, 9))
        ax = plt.subplot(1, 1, 1, projection=projection)

        ax.set_extent([-96,-68,50,67], crs=ccrs.PlateCarree())
        ax.add_feature(land_50m, color=[0.8, 0.8, 0.8])
        ax.coastlines(resolution='50m')
        """

        p1 = ax.pcolormesh(lons, lats, diff, transform=ccrs.PlateCarree(), cmap='bwr', vmin=-15, vmax=15)
        ax_cb = plt.axes([0.92, 0.25, 0.015, 0.5])
        cb = plt.colorbar(p1, cax=ax_cb, orientation='vertical')
        cb.ax.set_ylabel('Percentage Difference in Heat Content')
        ax.gridlines()
        #plt.savefig(figs_path+'heat_cont_diff_'+temp_id+'_'+notemp_id+'_'+season+'_'+st+'.png')
        plt.savefig(figs_path+'heat_cont_diff_'+temp_id+'_'+notemp_id+'_'+str(t)+'.png')
        #plt.show()
        plt.clf()

    d_temp.close()
    d_no_temp.close()


if __name__ == "__main__":
    heat_cont_diff(2017, 12 ,31, temp_id='ETW161', notemp_id='EPM151')
