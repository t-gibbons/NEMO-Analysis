import datetime
import numpy as np
import xarray as xr
import netCDF4 as nc
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as feature

###this is for running script backend###
import matplotlib
matplotlib.use('Agg')
###----------------------------------###

def z_masked_overlap(axe, X, Y, Z, source_projection=None):
    """
    for data in projection axe.projection
    find and mask the overlaps (more 1/2 the axe.projection range)

    X, Y either the coordinates in axe.projection or longitudes latitudes
    Z the data
    operation one of 'pcorlor', 'pcolormesh', 'countour', 'countourf'

    if source_projection is a geodetic CRS data is in geodetic coordinates
    and should first be projected in axe.projection

    X, Y are 2D same dimension as Z for contour and contourf
    same dimension as Z or with an extra row and column for pcolor
    and pcolormesh

    return ptx, pty, Z
    """
    if not hasattr(axe, 'projection'):
        return X, Y, Z
    if not isinstance(axe.projection, ccrs.Projection):
        return X, Y, Z

    if len(X.shape) != 2 or len(Y.shape) != 2:
        return X, Y, Z

    if (source_projection is not None and
            isinstance(source_projection, ccrs.Geodetic)):
        transformed_pts = axe.projection.transform_points(
            source_projection, X, Y)
        ptx, pty = transformed_pts[..., 0], transformed_pts[..., 1]
    else:
        ptx, pty = X, Y


    with np.errstate(invalid='ignore'):
        # diagonals have one less row and one less columns
        diagonal0_lengths = np.hypot(
            ptx[1:, 1:] - ptx[:-1, :-1],
            pty[1:, 1:] - pty[:-1, :-1]
        )
        diagonal1_lengths = np.hypot(
            ptx[1:, :-1] - ptx[:-1, 1:],
            pty[1:, :-1] - pty[:-1, 1:]
        )
        to_mask = (
            (diagonal0_lengths > (
                abs(axe.projection.x_limits[1]
                    - axe.projection.x_limits[0])) / 2) |
            np.isnan(diagonal0_lengths) |
            (diagonal1_lengths > (
                abs(axe.projection.x_limits[1]
                    - axe.projection.x_limits[0])) / 2) |
            np.isnan(diagonal1_lengths)
        )

        # TODO check if we need to do something about surrounding vertices

        # add one extra colum and row for contour and contourf
        if (to_mask.shape[0] == Z.shape[0] - 1 and
                to_mask.shape[1] == Z.shape[1] - 1):
            to_mask_extended = np.zeros(Z.shape, dtype=bool)
            to_mask_extended[:-1, :-1] = to_mask
            to_mask_extended[-1, :] = to_mask_extended[-2, :]
            to_mask_extended[:, -1] = to_mask_extended[:, -2]
            to_mask = to_mask_extended
        if np.any(to_mask):

            Z_mask = getattr(Z, 'mask', None)
            to_mask = to_mask if Z_mask is None else to_mask | Z_mask

            Z = ma.masked_where(to_mask, Z)

        return ptx, pty, Z


def salinity_anamoly_plot(runid, endyear, endmonth, endday, startyear=2002, startmonth=1, startday=5):
    output_path = '/project/6007519/weissgib/plotting/figs/salinity_plots/salinity_anamoly/'
    path = "/project/6007519/pmyers/ANHA4/ANHA4-"+runid+"-S/"
    grid_file = '/project/6007519/weissgib/plotting/data_files/anha4_files/ANHA4_mesh_mask.nc'
    data_output = '/project/6007519/weissgib/plotting/data_files/salinity/'

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

    #and now make a list of model files to read
    mdl_files = []
    for t in times:
        mdl_files.append(path+"ANHA4-"+runid+"_y"+str(t.year)+"m"+str(t.month).zfill(2)+"d"+str(t.day).zfill(2)+"_gridT.nc")

    d = xr.open_mfdataset(mdl_files, concat_dim='time_counter', data_vars='minimal', coords='minimal', compat='override')

    #also want to read in the mesh grid info
    mesh = nc.Dataset(grid_file)

    mask = np.array(mesh.variables['tmask'])
    lons = np.array(mesh.variables['nav_lon'])
    lats = np.array(mesh.variables['nav_lat'])

    mesh.close()
    
    #mask the data

    d.coords['mask'] = (('deptht', 'y_grid_T', 'x_grid_T'), mask[0,:,:,:])
    d = d.where(d.mask == 1)

    full_depth = list(d['deptht'].values)

    #want to do annual averages and drop the spinup period (2002-2005)
    annual_average = d.groupby('time_counter.year').mean('time_counter')
    annual_average = annual_average.isel(year=annual_average.year.isin([2006,2007,2008,2009,2010,2011,2012,2013,2014,2015,2016,2017,2018,2019]))
    
    print(annual_average)

    """
    #looking at the 200-1000m level
    d = annual_average.where(annual_average['deptht'] > 200.0, drop=True)
    d = d.where(d['deptht'] < 1000.0, drop=True)
    """

    #0-220m level
    d = annual_average.where(annual_average['deptht'] < 200.0, drop=True)

    print(d['deptht'])

    #calculate the weights
    zero_start = True
    n = len(d['deptht'])
    weight = np.zeros(n)
    dd = d['deptht'][n-1] #SHOULD THIS BE THE FROM FULL DEPTH??
    for i in range(n):
        if zero_start:
            if i == 0:
                weight[i] = d['deptht'][i]/dd
            else:
                weight[i] = (d['deptht'][i] - d['deptht'][i-1])/dd
        else:
            if i == 0:
                k = full_depth.index(d['deptht'][i])
                weight[i] = (d['deptht'][i] - full_depth[k-1])/dd
            else:
                weight[i] = (d['deptht'][i] - d['deptht'][i-1])/dd


    weights = xr.DataArray(weight, coords=[d['deptht']], dims=['deptht'])

    #and take the average
    d_weighted = d.weighted(weights)
    surface_salinity = d_weighted.mean(dim='deptht', skipna=True)
    
    """
    #first want to get the mean salinity to calculate anamolies from
    
    mean_salinity = surface_salinity.mean('year')
    print(mean_salinity)

    #output the mean so can get a combined run mean
    mean_salinity['vosaline'].to_netcdf(data_output+runid+'_mean_salinity_0_200m.nc')

    
    print(mean_salinity)
    mean = mean_salinity['vosaline'].values
    print(mean.shape)
    """

    #read in the mean for all model runs
    mean = np.load(data_output+'total_mean_0-200.npy')

    #now calculate the salinity anamoly for each year from the mean
    years = surface_salinity['year'].values

    for y in range(surface_salinity.dims['year']):
        sal = surface_salinity['vosaline'].isel(year=y).values
        anamoly = sal-mean

        #plot just the Labrador Sea region
        ax = plt.axes(projection=ccrs.Mercator(central_longitude=-45))
        ax.coastlines(resolution='50m')
        land_50m = feature.NaturalEarthFeature('physical', 'land', '50m', edgecolor='black', facecolor='gray', linewidth=0.5)
        ax.add_feature(land_50m,color=[0.8, 0.8, 0.8])

        X,Y,masked_anamoly = z_masked_overlap(ax,lons,lats,anamoly,source_projection=ccrs.Mercator(central_longitude=-45))

        #p1 = ax.pcolormesh(lons,lats,anamoly,transform=ccrs.PlateCarree(), cmap='bwr', vmin=-0.2, vmax=0.2)
        p1 = ax.contourf(X,Y,masked_anamoly,transform=ccrs.PlateCarree(), cmap='bwr', levels=[-0.2,-0.15,-0.1,-0.05,0,0.05,0.1,0.15,0.2])
        ax.set_extent([-100,30,20,80], crs=ccrs.PlateCarree())
        ax_cb = plt.axes([0.92, 0.25, 0.015, 0.5])
        cb = plt.colorbar(p1,cax=ax_cb, orientation='vertical')

        plt.savefig(output_path+'salinity_anamoly_total_0_200_'+runid+'_'+str(years[y])+'.png')
        #plt.show()
        plt.clf()
    
if __name__ == "__main__":
    salinity_anamoly_plot(runid='EPM151', endyear=2019, endmonth=12, endday=31)
    salinity_anamoly_plot(runid='EPM152', endyear=2019, endmonth=8, endday=23)
    salinity_anamoly_plot(runid='EPM014', endyear=2019, endmonth=8, endday=23)
    salinity_anamoly_plot(runid='EPM015', endyear=2019, endmonth=12, endday=31)
    
    
