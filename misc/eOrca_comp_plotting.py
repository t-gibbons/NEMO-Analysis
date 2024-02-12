import iris
import iris.plot as iplt
import iris.quickplot as qplt
import iris.analysis.cartography

import numpy as np
import cartopy.crs as ccrs
import cartopy.feature as feature
import matplotlib.pyplot as plt
import matplotlib.path as mpath

path = '/mnt/storage6/tahya/model_files/eOrca_diff_fwc_1994.nc'
cube = iris.load_cube(path)

"""
path_36 = '/mnt/storage6/tahya/model_files/eOrca_36_temp_1994.nc'
cube_36 = iris.load_cube(path_36)

print(cube_42)
print(cube_36)

cube = cube_36-cube_42
"""
print(cube) 

projection = ccrs.NorthPolarStereo()

pcarree = ccrs.PlateCarree()

#transform cube to target projection
new_cube, extent = iris.analysis.cartography.project(cube, pcarree, nx=1207, ny=1442)

new_cube = new_cube.extract(iris.Constraint(year=1994))

print(new_cube)

land_50m = feature.NaturalEarthFeature('physical', 'land', '50m', edgecolor='black', facecolor='gray', linewidth=0.5)

fig = plt.figure(figsize=(10,9))
fig.suptitle("Difference in Freshwater Content")
ax = plt.subplot(projection=projection)
ax.set_extent([-180, 180, 90, 60], crs=ccrs.PlateCarree())
ax.add_feature(land_50m, color=[0.8,0.8,0.8])

# Compute a circle in axes coordinates, which we can use as a boundary
# for the map. We can pan/zoom as much as we like - the boundary will be
# permanently circular.
theta = np.linspace(0, 2*np.pi, 100)
center, radius = [0.5, 0.5], 0.5
verts = np.vstack([np.sin(theta), np.cos(theta)]).T
circle = mpath.Path(verts * radius + center)
ax.set_boundary(circle, transform=ax.transAxes)

p1 = qplt.pcolormesh(new_cube, cmap='bwr', vmin=-2, vmax=2)

ax.coastlines(resolution='50m')

ax_cb = plt.axes([0.92,0.25,0.015,0.5])
cb = plt.colorbar(p1, cax = ax_cb, orientation='vertical')

plt.savefig('fwc_diff_1994.png')
