"""
plot the surface water temperature for a region over a time period
"""
import glob
import numpy as np
import xarray as xr
import cartopy.crs as ccrs
import cartopy.feature as feature
import matplotlib.pyplot as plt
import matplotlib.path as mpath

###this is for running script backend###
import matplotlib
matplotlib.use('Agg')
###----------------------------------###

runid = 'ETW162'

path = "/project/6007519/weissgib/ANHA4/ANHA4-"+runid+"-S/"
figs_path = '/project/6007519/weissgib/plotting/heat/figs/surface_temp/'

mdl_files = glob.glob(path+'ANHA4-'+runid+'_y201*_gridT.nc')

df = xr.open_mfdataset(mdl_files)

print(df)

df = df.where(df['deptht'] < 50.0, drop=True)

df = df['votemper']

#calculate the weights of the depth levels
n = len(df['deptht'])
weight = np.zeros(n)
dz = np.zeros(n)
dd = df['deptht'][n-1]
for i in range(n):
    if i == 0:
        weight[i] = df['deptht'][i]/dd
        dz[i] = df['deptht'][i]
    else:
        weight[i] = (df['deptht'][i] - df['deptht'][i-1])/dd
        dz[i] = df['deptht'][i] - df['deptht'][i-1]

x = df.sizes['x_grid_T']
y = df.sizes['y_grid_T']

weights = xr.DataArray(weight, coords=[df['deptht']], dims=['deptht'])

#and take the average 
d_weighted = df.weighted(weights)
surface_temp = d_weighted.mean(dim='deptht', skipna=True)

print(surface_temp)

times = surface_temp['time_counter'].values

lons = df['nav_lon_grid_T'].values
lats = df['nav_lat_grid_T'].values

for t in times:

    st = t.strftime('%Y-%m-%d')
    print(st)
    
    #mackenzie river region
    land_50m = feature.NaturalEarthFeature('physical', 'land', '50m', edgecolor='black', facecolor='gray', linewidth=0.5)
    projection=projection=ccrs.Mercator(central_longitude=-80)

    fig = plt.figure(figsize=(10, 9))
    ax = plt.subplot(1, 1, 1, projection=projection)

    ax.set_extent([-163,-119,67,73], crs=ccrs.PlateCarree())
    ax.add_feature(land_50m, color=[0.8, 0.8, 0.8])
    ax.coastlines(resolution='50m')

    #temp = df.isel(time_counter=0, deptht=0).values

    #p1 = ax.pcolormesh(lons, lats, temp, transform=ccrs.PlateCarree(), cmap='winter')
    #p1 = df.isel(time_counter=0, deptht=0).plot(subplot_kws=dict(projection=projection), transform=ccrs.PlateCarree(), cmap='viridis')#, vmin=-5, vmax=20)
    p1 = surface_temp.sel(time_counter=t).plot(x='nav_lon_grid_T', y='nav_lat_grid_T', subplot_kws=dict(projection=projection), transform=ccrs.PlateCarree(), cmap='winter', vmin=-2, vmax=10)
    #ax_cb = plt.axes([0.92, 0.25, 0.015, 0.5])
    #cb = plt.colorbar(p1, cax=ax_cb, orientation='vertical')
    #cb.ax.set_ylabel('Surface Temperature')
    ax.gridlines()
    plt.savefig(figs_path+'surface_temp_'+runid+'_'+st+'.png')
    #plt.show()
    plt.clf()

