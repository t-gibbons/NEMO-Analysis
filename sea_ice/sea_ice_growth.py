"""
sea_ice-growth.py - author: weissgib@ualberta.ca

Compare the melt and growth parameters for sea ice between model runs
We want to compare between a LIM2 vs LIM3 run
"""

import xarray as xr
import cartopy.crs as ccrs
import matplotlib.pyplot as plt

#lim2 run
lim2_path = '/project/6007519/weissgib/ANHA4/ANHA4-ETW161-S/ANHA4-ETW161_y2014m09d27_icemod.nc'
lim2 = xr.open_mfdataset(lim2_path)

#lim3 run
lim3_path = '/project/6007519/weissgib/ANHA4/ANHA4-ETW162-S/ANHA4-ETW162_y2014m09d27_icemod.nc'
lim3 = xr.open_mfdataset(lim3_path)

#convert the units of lim2
lim2_ice = (lim2['iiceprod'].values)*86400
lim3_ice = lim3['itoticegrowmelt'].values
diff = lim2_ice-lim3_ice

lons = lim2['nav_lon'].values
lats = lim2['nav_lat'].values

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
p1 = ax.pcolormesh(lons, lats, diff, transform=ccrs.PlateCarree(), cmap='bwr')
ax_cb = plt.axes([0.92, 0.25, 0.015, 0.5])
cb = plt.colorbar(p1,cax=ax_cb, orientation='vertical')
cb.ax.set_ylabel('Difference in Sea Ice Melt/Growth')
ax.gridlines()
#plt.savefig(figs_path+'sea_ice_thickness_'+runid+'_'+str(yr)+'.png')
plt.show()

lim2.close()
lim3.close()
