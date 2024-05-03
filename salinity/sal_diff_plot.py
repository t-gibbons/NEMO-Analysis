import numpy as np
import xesmf as xe
import xarray as xr
import matplotlib.path as mpath
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as feature

#compare the surface salinity from glorys to a ANHA4 run
def glorys_comp(runid="EPM151"):

    #path to the folder with the compute surface salinity
    path = "/mnt/storage6/tahya/model_files/"
    
    #lets start by reading in the model files
    model_data = xr.open_mfdataset(path+runid+'_surface_salinity_50m.nc')

    #take the annual average
    annual_average_model = model_data.groupby('time_counter.year').mean('time_counter')

    lons = model_data['nav_lon_grid_T'].values
    lats = model_data['nav_lat_grid_T'].values

    annual_average_model = annual_average_model['vosaline']
    print(annual_average_model)

    #now read in the glorys data
    data_files = []
    year_slice = [2017, 2019]
    for y in range(year_slice[0], year_slice[1]+1):
        data_files.append(path+'glorys12v1_sal_50m_'+str(y)+'.nc')

    glorys_data = xr.open_mfdataset(data_files)
    print(glorys_data)

    #now we need to regrid glorys to the anha4 grid
    annual_average_model = annual_average_model.rename({'nav_lon_grid_T': 'lon', 'nav_lat_grid_T': 'lat'})
    ds_out = model_data.rename({'nav_lon_grid_T': 'lon', 'nav_lat_grid_T': 'lat'}) 
    ds_out = ds_out[['lon','lat']]

    regridder = xe.Regridder(glorys_data, ds_out, 'bilinear', extrap_method='inverse_dist')

    var_out = regridder(glorys_data['__xarray_dataarray_variable__'])

    print(var_out)

    #take the annual average
    annual_average_glorys = var_out.groupby('time.year').mean('time')
    print(annual_average_glorys)

    #now we can take the difference
    diff = annual_average_glorys-annual_average_model
    print(diff)

    #and finally plot
    yr = str(year_slice[0])+'_'+str(year_slice[1])
    sal2d = diff.sel(year=slice(year_slice[0], year_slice[1])).mean('year')
    print(sal2d)
    plot_diff(lons, lats, sal2d, yr)

    glorys_data.close()
    model_data.close()

#compare the surface salinity for a time slice between two anha4 runs
def model_comp(runids=['EPM015', 'EPM151']):

    #lets start by reading in the model files

    #path to the folder with the model runs
    path = "/mnt/storage6/tahya/model_files/"

    d1 = xr.open_mfdataset(path+runids[0]+'_surface_salinity_50m.nc')
    d2 = xr.open_mfdataset(path+runids[1]+'_surface_salinity_50m.nc')

    annual_average1 = d1.groupby('time_counter.year').mean('time_counter')
    annual_average2 = d2.groupby('time_counter.year').mean('time_counter')

    lons = d1['nav_lon_grid_T'].values
    lats = d2['nav_lat_grid_T'].values

    diff = annual_average1-annual_average2
    yrs = diff['year'].values

    #year_slices = [['2005','2007'],['2008','2010'],['2011','2013'],['2014','2016'],['2017','2019']]
    year_slices = [['2017','2019']]

    #for y in range(diff.dims['year']):
    for y in year_slices:
        yr = y[0]+'_'+y[1]
        sal2d = diff['vosaline'].sel(year=slice(y[0], y[1])).mean('year')
        plot_diff(lons, lats, sal2d, yr)

def plot_diff(lons, lats, diff, yr):

    #north pole stero projection
    land_50m = feature.NaturalEarthFeature('physical', 'land', '50m', edgecolor='black', facecolor='gray', linewidth=0.5)
    projection=ccrs.NorthPolarStereo()

    fig = plt.figure(figsize=(10, 9))
    ax = plt.subplot(1, 1, 1, projection=projection)

    ax.set_extent([-180, 180, 90, 60], crs=ccrs.PlateCarree())
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
    p1 = ax.pcolormesh(lons, lats, diff, transform=ccrs.PlateCarree(), cmap='RdBu', vmin=-2, vmax=2)
    ax_cb = plt.axes([0.92, 0.25, 0.015, 0.5])
    cb = plt.colorbar(p1,cax=ax_cb, orientation='vertical')
    cb.ax.set_ylabel('Difference Surface Salinity')
    ax.gridlines()
    #plt.show()
    plt.savefig('surface_sal_diff_cgrf_'+str(yr)+'.png')
    plt.clf()

if __name__ == "__main__":
    #glorys_comp('EPM015')
    model_comp()
