import numpy as np
import xarray as xr
import netCDF4 as nc
import matplotlib.pyplot as plt
import matplotlib.path as mpath
import cartopy.crs as ccrs
import cartopy.feature as feature

path = '/project/6007519/weissgib/plotting/data_files/freshwater_content/'
#runids = {'CGRF': ['EPM015', 'EPM101'], 'ERA': ['EPM014', 'EPM102'] }
runids = {'CGRF': ['EPM151', 'EPM101'], 'ERA': ['EPM152', 'EPM102'] }

#read in the mesh grid file
grid_file = '/project/6007519/weissgib/plotting/data_files/anha4_files/ANHA4_mesh_mask.nc'
mesh = nc.Dataset(grid_file)

lons = np.array(mesh.variables['nav_lon'])
lats = np.array(mesh.variables['nav_lat'])

mesh.close()

for atmo in runids:
    #read in the model files
    print(atmo)

    dai = path+runids[atmo][0]+'_annual_avg_fwc_200m.nc'
    hype = path+runids[atmo][1]+'_annual_avg_fwc_200m.nc'

    dai_data = xr.open_dataset(dai)
    hype_data = xr.open_dataset(hype)

    #take an overall average
    dai_mean = dai_data.sel(year=slice("2017", "2019")).mean('year')
    hype_mean = hype_data.sel(year=slice("2017", "2019")).mean('year')

    dd = dai_mean['vosaline'].values
    hd = hype_mean['vosaline'].values
    print(dd)
    print(hd)

    diff = dd-hd

    #now plot the difference
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
    p1 = ax.pcolormesh(lons, lats, diff, transform=ccrs.PlateCarree(), cmap='bwr', vmin=-5, vmax=5)
    ax_cb = plt.axes([0.92, 0.25, 0.015, 0.5])
    cb = plt.colorbar(p1,cax=ax_cb, orientation='vertical')
    cb.ax.set_ylabel('Freshwater Content (m)')
    ax.gridlines()
    #plt.show()
    plt.savefig(path+atmo+'_2017_2019_diff_fwc.png')
    plt.clf()

    dai_data.close()
    hype_data.close()
