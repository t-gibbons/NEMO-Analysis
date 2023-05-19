import datetime
import numpy as np
import numpy.ma as ma
import netCDF4 as nc
import matplotlib.pyplot as plt
import matplotlib.path as mpath
import cartopy.crs as ccrs
import cartopy.feature as feature

###this is for running script backend###
import matplotlib
matplotlib.use('Agg')
###----------------------------------###

def river_tracers(runid, endyear, endmonth, endday, startyear=2002, startmonth=1, startday=5):
    figs_path = '/project/6007519/weissgib/plotting/figs/tracers/'
    path = "/project/6007519/weissgib/ANHA4/ANHA4-"+runid+"-S/"
    home_path = '/project/6007519/weissgib/plotting/'

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

    #just need the T grid variables
    mdl_files = []
        
    for t in times:
        mdl_files.append(path+"ANHA4-"+runid+"_y"+str(t.year)+"m"+str(t.month).zfill(2)+"d"+str(t.day).zfill(2)+"_gridT.nc")

    #read in the mesh grid info
    mesh_file = home_path+'ANHA4_mesh_mask.nc'
    mf = nc.Dataset(mesh_file)

    nav_lon = np.array(mf.variables['nav_lon'])
    nav_lat = np.array(mf.variables['nav_lat'])

    mf.close()

    #tracers = ['TRC13', 'TRC14', 'TRC15', 'TRC16']
    tracers = ['TRC15']

    #north pole stero projection
    land_50m = feature.NaturalEarthFeature('physical', 'land', '50m', edgecolor='black', facecolor='gray', linewidth=0.5)
    projection = ccrs.NorthPolarStereo()
    fig = plt.figure(figsize=(10,9))

    for i in range(len(mdl_files)):

        mdl_file = nc.Dataset(mdl_files[i])
        tim = times[i]

        #make plots of the tracers
        for t in tracers:

            trac = np.array(mdl_file[t])[0,:,:,:]
            miss_val = mdl_file[t].missing_value

            trac[trac==miss_val] = np.nan
            
            #sum the tracer over the water column
            trac = np.sum(trac, axis=0)

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
            p1 = ax.pcolormesh(nav_lon, nav_lat, trac, transform=ccrs.PlateCarree(), cmap='jet', vmin=0, vmax=0.5)
            ax_cb = plt.axes([0.92, 0.25, 0.015, 0.5])
            cb = plt.colorbar(p1,cax=ax_cb, orientation='vertical')
            ax.gridlines()
            plt.title(str(tim))
            plt.savefig(figs_path+runid+'_'+t+'_'+str(tim)+'.png')
            #plt.show()
            plt.clf()

        mdl_file.close()


if __name__ == "__main__":
        river_tracers(runid='ETW001', startyear=2002, endyear=2019, endmonth=6, endday=14)

