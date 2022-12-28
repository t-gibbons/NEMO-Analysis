"""
Map a transport section on map of the domain
Need the i,j start and end point
Uses same line calculation as transport script
"""
import numpy as np
import netCDF4 as nc
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import matplotlib.path as mpath
import cartopy.feature as feature

#first calculate the section data
x1 = 256
x2 = 282
y1 = 560
y2 = 546

ii = []
jj = []

dx = x2-x1
dy = y2-y1
yi = 1

if dy < 0:
    yi = -1
    dy = -dy

D = (2*dy) - dx
y = y1
slope = []

for x in range(x1,x2):
    ii.append(x)
    jj.append(y)
    if D > 0:
        y = y +yi
        D = D + (2*(dy-dx))
        ii.append(x)
        jj.append(y)
    else:
        D = D + 2*dy

print(ii)
print(jj)

#next read in the mesh grid file 
mask_file = '/project/6007519/weissgib/plotting/data_files/anha4_files/ANHA4_mesh_mask.nc'
mf = nc.Dataset(mask_file)

nav_lon = np.array(mf.variables['nav_lon'])
nav_lat = np.array(mf.variables['nav_lat'])

mf.close()

#get the grid coordinates on map coordinates
ln = []
lt = []
for k in range(len(ii)):
    ln.append(nav_lon[jj[k], ii[k]])
    lt.append(nav_lat[jj[k], ii[k]])

print(ln)
print(lt)

#now we can plot the line on a map 
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

ax.plot(ln, lt, linewidth=3.0, transform=ccrs.PlateCarree())
#plt.show()
plt.savefig('greenland_north_section.png')

