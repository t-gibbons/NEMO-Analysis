import numpy as np
import xarray as xr
import pandas as pd
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import matplotlib.path as mpath
import cartopy.feature as feature

start_year = 2002
end_year = 2019


old_path = '/project/6007519/pmyers/ANHA4-I/RUNOFF/Bamber2012/'
files = []

for y in range(start_year, end_year+1):
    files.append(old_path+'ANHA4_runoff_monthly_combined_Dai_Trenberth_Bamber_y'+str(y)+'.nc')

old_runoff = xr.open_mfdataset(files)

"""
path = '/project/6007519/ANHA4-I/RUNOFF/HydroGFD/'
data = []

for y in range(start_year, end_year+1):
        files = path+'ANHA4_ReNat_HydroGFD_HBC_runoff_monthly_y'+str(y)+'.nc'
        ds = xr.open_mfdataset(files, decode_times=False)
        reference_date = '1/1/'+str(y)
        v = ds['runoff'].values
        ds['time_counter'] = pd.date_range(start=reference_date, periods=ds.sizes['time_counter'], freq='MS')
        data.append(ds)

new_runoff = xr.concat(data, dim='time_counter')

runoff_average = new_runoff.mean('time_counter')
"""

print(old_runoff)

runoff_average = old_runoff.mean('time_counter')

print(runoff_average)

#fill zero values with nan
runoff_average = runoff_average.where(runoff_average['runoff'] > 0)

nav_lon = runoff_average['nav_lon'].values
nav_lat = runoff_average['nav_lat'].values

#and scattter plot for entire domain
land_50m = feature.NaturalEarthFeature('physical', 'land', '50m', edgecolor='black', facecolor='gray', linewidth=0.5)
projection = ccrs.NorthPolarStereo()

fig = plt.figure(figsize=(10,9))
ax = plt.subplot(1, 1, 1, projection=projection)

ax.set_extent([-180,180,90,60], crs=ccrs.PlateCarree())
#ax.add_feature(land_50m, color=[0.8, 0.8, 0.8])
ax.coastlines(resolution='50m')

theta = np.linspace(0, 2*np.pi, 100)
center, radius = [0.5,0.5], 0.5
verts = np.vstack([np.sin(theta), np.cos(theta)]).T
circle = mpath.Path(verts*radius+center)
ax.set_boundary(circle, transform=ax.transAxes)

p1 = ax.scatter(nav_lon, nav_lat, c=runoff_average['runoff'], transform=ccrs.PlateCarree(), cmap='gist_rainbow')
ax_cb = plt.axes([0.92,0.25,0.015,0.5])
cb = plt.colorbar(p1,cax=ax_cb,orientation='vertical')

ax.gridlines()
plt.savefig('dai_all_runoff_locations.png')
#plt.show()
