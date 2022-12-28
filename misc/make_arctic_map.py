import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import matplotlib.path as mpath
import cartopy.feature as feature


def get_section_line(x1,x2,y1,y2,nav_lon,nav_lat):
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

    #get the grid coordinates on map coordinates
    ln = []
    lt = []
    for k in range(len(ii)):
        ln.append(nav_lon[jj[k], ii[k]])
        lt.append(nav_lat[jj[k], ii[k]])

    return(ln,lt)

#first read in the bathymetry data

bath_file = '/project/6007519/weissgib/plotting/data_files/anha4_files/ANHA4_bathymetry.nc'

bf = nc.Dataset(bath_file)

bathy = np.array(bf.variables['bathyt'])

bf.close()

#and get the mesh grid information

grid_file = '/project/6007519/weissgib/plotting/data_files/anha4_files/ANHA4_mesh_mask.nc'

mesh = nc.Dataset(grid_file)
    
lons = np.array(mesh.variables['nav_lon'])
lats = np.array(mesh.variables['nav_lat'])
               
mesh.close()

#get the sections for the major gateways

davis_ln,davis_lt = get_section_line(175,214,443,443,lons,lats)

fram_ln,fram_lt = get_section_line(304,360,503,526,lons,lats)

#and plot

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
p1 = ax.pcolormesh(lons, lats, bathy, transform=ccrs.PlateCarree(), cmap='YlGnBu')
ax_cb = plt.axes([0.92, 0.25, 0.015, 0.5])
cb = plt.colorbar(p1,cax=ax_cb, orientation='vertical')
cb.ax.set_ylabel('Depth (m)')
ax.gridlines()
ax.plot(davis_ln, davis_lt, linewidth=4.0, color='k', transform=ccrs.PlateCarree())
ax.plot(fram_ln, fram_lt, linewidth=4.0, color='k', transform=ccrs.PlateCarree())
#plt.show()
plt.savefig('arctic_map.svg', format='svg', dpi=1200)
plt.clf()
