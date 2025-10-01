"""
Create the ANHA4x bathy and coordinate file

Based off old ANHA4 bathy file, with slice taken from eOrca for extended boundary

Extended to make other supporting files for new ANHA4x configuration
"""
import re
import glob
import pathlib
import numpy as np
import xarray as xr
import netCDF4 as nc
#import cartopy.crs as ccrs
#import matplotlib.pyplot as plt
#import cartopy.feature as feature
#import matplotlib.colors as colors


def make_bathy():
    anha4_file = '/project/6007519/ANHA4-I/ANHA4_bathy_etopo1_gebco1_smoothed_coast_corrected_mar10_Tide.nc'
    eorca_file = '/project/6007519/pennelly/eORCA025-I/eORCA025_Bathymetry.nc'

    af = xr.open_dataset(anha4_file)
    ef = xr.open_dataset(eorca_file)

    #get the boundary from old anha4

    anha_lon = af['nav_lon'].values
    anha_lat = af['nav_lat'].values
    anha_bathy = af['Bathymetry'].values
    print(anha_bathy.shape)

    eorca_bathy = ef['Bathymetry'].values
    eorca_lon = ef['nav_lon'].values
    eorca_lat = ef['nav_lat'].values

    """
    #this is how figured out where to cut eorca to match anha
    lon_bd = anha_lon[799, :]
    lat_bd = anha_lat[799,:]

    #find where that lines up in the eOrca grid

    lon_start = np.where(eorca_lon == lon_bd[0])
    lat_start = np.where(eorca_lat == lat_bd[0])

    lon_end = np.where(eorca_lon == lon_bd[-1])
    lat_end = np.where(eorca_lat == lat_bd[-1])
    """

    #cheating and hardcoding where I think the matching line is in eOrca
    #should be j=1009, i=143:686, found using output above
    #want to take a chunk which extends down to j=898

    eorca_bathy = ef['Bathymetry'].values
    eorca_lon = ef['nav_lon'].values
    eorca_lat = ef['nav_lat'].values

    eorca_chunk = eorca_bathy[898:1009, 143:687]
    lon_chunk = eorca_lon[898:1009, 143:687]
    lat_chunk = eorca_lat[898:1009, 143:687]
    print(eorca_chunk.shape)

    #now need to rotate this chunk so we can properly stitch them together
    #use rot90 twice to rotate 180 degrees

    r1 = np.rot90(eorca_chunk)
    rot_chunk = np.rot90(r1)
    print(rot_chunk.shape)

    l1 = np.rot90(lon_chunk)
    lon_rot = np.rot90(l1)

    l2 = np.rot90(lat_chunk)
    lat_rot = np.rot90(l2)

    anhax_bathy = np.concatenate((anha_bathy, rot_chunk), axis=0)
    print(anhax_bathy.shape)

    anhax_lon = np.concatenate((anha_lon, lon_rot), axis=0)
    anhax_lat = np.concatenate((anha_lat, lat_rot), axis=0)

    #plt.pcolormesh(anhax_bathy, vmin=0, vmax=100)
    #plt.show()


    """
    #lets try plotting the new anhax on a map
    land_50m = feature.NaturalEarthFeature('physical', 'land', '50m', edgecolor='black', facecolor='gray', linewidth=0.5)

    fig = plt.figure(figsize=(10,9))
    ax = plt.subplot(1,1,1, projection=ccrs.NorthPolarStereo())

    #ax.add_feature(land_50m, color=[0.8,0.8,0.8])
    ax.coastlines(resolution='50m')

    p1 = ax.pcolormesh(anhax_lon, anhax_lat, anhax_bathy, transform=ccrs.PlateCarree(), norm=colors.LogNorm(), cmap='ocean')
    ax_cb = plt.axes([0.92, 0.25, 0.015, 0.5])
    cb = plt.colorbar(p1, cax=ax_cb, orientation='vertical')
    #plt.show()
    plt.savefig('anhax_bathy.png')

    #now lets make a new bathy file

    anhax_data = xr.Dataset(data_vars=dict(nav_lon=(["y","x"], anhax_lon), nav_lat=(["y","x"], anhax_lat), Bathymetry=(["y","x"], anhax_bathy)))
    print(anhax_data)

    anhax_data.to_netcdf('ANHA4x_bathy.nc')
    """

    af.close()
    ef.close()

def make_coords():
    anha4_file = '/project/6007519/ANHA4-I/ANHA4_coordinates.nc'
    eorca_file = '/project/6007519/pennelly/eORCA025-I/eORCA025_coordinates.nc'

    af = xr.open_dataset(anha4_file)
    ef = xr.open_dataset(eorca_file)

    #we know the dims will be y=911 and x=544 for the new anhax
    anhax = xr.Dataset(None, coords=dict(y=range(911), x=range(544)))

    for v in list(af.keys()):
        print(v)
        anha_var = af[v].values
        ev = ef[v].values

        #get the variable chunk to add from eorca
        eorca_chunk = ev[898:1009, 143:687]

        #now rotate to match the anha grid

        r1 = np.rot90(eorca_chunk)
        rot_chunk = np.rot90(r1)

        #and concatenate with anha
        anhax_var = np.concatenate((anha_var, rot_chunk), axis=0)
        print(anhax_var.shape)

        anhax[v] = (('y','x'), anhax_var)


    """
    #hardcoding the slice from eorca that matches the anha4 grid we want to take
    eorca_chunk = ef.isel(y=slice(898,1009), x=slice(143,687))
    print(eorca_chunk)

    #need to rotate in order to stitch with anha
    eorca_rot = eorca_chunk.isel(x=slice(None, None, -1))
    print(eorca_rot)

    #and concatenate them together
    print(af)
    print(eorca_rot)

    anhax = xr.concat([af, eorca_rot], dim="y")
    print(anhax)

    anhax.to_netcdf('ANHAx_coordinates.nc')
    """

    print(anhax)
    
    anhax.to_netcdf('ANHAx_coordinates.nc')

    af.close()
    ef.close()


def make_runoff():

    anha_loc = '/project/6007519/ANHA4-I/RUNOFF/HydroGFD/'
    eorca_loc = '/project/6007519/eORCA025-I/Runoff/'

    start_year = 2002
    end_year = 2018

    for y in range(start_year, end_year+1):

        eorca_file = eorca_loc+'eORCA025_runoffICBmelt_monthly_combined_Dai_Trenberth_Bamber2016_y'+str(y)+'.nc'
        anha_file = eorca_loc+'ANHA4_ReNat_HydroGFD_HBC_runoff_monthly_y'+str(y)+'.nc'

        af = xr.open_dataset(anha_file, decode_times=False)
        ef = xr.open_dataset(eorca_file)

        #we know the dims will be y=911 and x=544 for the new anhax
        anhax = xr.Dataset(None, coords=dict(y=range(911), x=range(544)))

        for v in list(af.keys()):
            print(v)
            anha_var = af[v].values
            ev = ef[v].values

            if v == 'runoff':

                eorca_chunk = ev[:, 898:1009, 143:687]
                print(eorca_chunk.shape)
                print(anha_var.shape)

                #rotate
                r1 = np.rot90(eorca_chunk, axes=(1,2))
                rot_chunk = np.rot90(r1, axes=(1,2))

                print(rot_chunk.shape)

                anhax_var = np.concatenate((anha_var, rot_chunk), axis=1)
                print(anhax_var.shape)

                anhax[v] = (('time_counter','y','x'), anhax_var)

                continue

            #get the variable chunk to add from eorca
            eorca_chunk = ev[898:1009, 143:687]

            #now rotate to match the anha grid

            r1 = np.rot90(eorca_chunk)
            rot_chunk = np.rot90(r1)

            #and concatenate with anha
            anhax_var = np.concatenate((anha_var, rot_chunk), axis=0)
            print(anhax_var.shape)

            anhax[v] = (('y','x'), anhax_var)


        anhax = anhax.drop_vars("x")
        anhax = anhax.drop_vars("y")
        print(anhax)

        anhax.to_netcdf('ANHAx_HydroGFD_Dai_runoff_monthly_y'+str(y)+'.nc', unlimited_dims={'time_counter': True})

        af.close()
        ef.close()


def make_tracer_mask():

    #a simple hack to duplicate the tracer mask for anhax using the existing eOrca mask

    anha4_file = '/project/6007519/ANHA4-I/TRC/tracer_mask_ANHA4_Oct2016.nc'
    eorca_file = '/project/6007519/eORCA025-I/TRC/eORCA025_tracer_mask_Dec2023.nc'

    af = xr.open_dataset(anha4_file)
    ef = xr.open_dataset(eorca_file)

    #we know the dims will be y=911 and x=544 for the new anhax
    anhax = xr.Dataset(None, coords=dict(y=range(911), x=range(544)))

    for v in list(af.keys()):
        print(v)
        anha_var = af[v].values
        ev = ef[v].values

        #get the variable chunk to add from eorca
        eorca_chunk = ev[898:1009, 143:687]

        #now rotate to match the anha grid

        r1 = np.rot90(eorca_chunk)
        rot_chunk = np.rot90(r1)

        #and concatenate with anha
        anhax_var = np.concatenate((anha_var, rot_chunk), axis=0)
        print(anhax_var.shape)

        anhax[v] = (('y','x'), anhax_var)

    print(anhax)

    anhax.to_netcdf('ANHA4x_trc_mask_Oct2024.nc')

    af.close()
    ef.close()

def make_runoff_mask():

    #lets just make a runoff mask from the runoff files

    anha4_file = '/project/6007519/weissgib/ANHA4x-I/RUNOFF/ANHAx_HydroGFD_Dai_runoff_monthly_y2015.nc'

    af = xr.open_dataset(anha4_file)

    af = af.drop_vars(['runoff'])

    print(af)

    af.to_netcdf('ANHA4x_runoff_mask_long_simulation.nc')

    af.close()


def make_calving_files():

    #get the navlon/lat info from the mesh file
    ax_mesh = '/project/6007519/weissgib/ANHA4x-I/ANHA4x_mesh_mask.nc'

    ax = xr.open_mfdataset(ax_mesh)

    nav_lon = ax['nav_lon']
    nav_lat = ax['nav_lat']

    cal_path = '/project/6007519/ANHA4-I/CALVING/Bamber2018/'
    save_path = '/project/6007519/weissgib/ANHA4x-I/CALVING/Bamber2018/'

    #make the calving.nc file
    anahx = xr.Dataset(None, coords=dict(y=range(911), x=range(544)))
    anahx['nav_lon'] = nav_lon
    anahx['nav_lat'] = nav_lat
    anahx['maxclass'] = (('y', 'x'), np.full((911,544), 10))

    anahx.to_netcdf(save_path+'calving.nc')
    exit()

    calving_files = glob.glob(cal_path+'ANHA4*.nc')

    for cf in calving_files:

        print(cf)

        old_filename = pathlib.Path(cf).name
        new_filename = re.sub(r'(ANHA4)', r'\1x', old_filename)

        af = xr.open_mfdataset(cf)

        #create the blank template for new file
        anhax = xr.Dataset(None, coords=dict(y=range(911), x=range(544)))
        anhax['nav_lon'] = nav_lon
        anhax['nav_lat'] = nav_lat
        anhax['time_counter'] = af['time_counter']

        #get the data from the anha4 mask and expand for anha4x
        calving_mask = af['calvingmask'].values
        socoefr = af['socoefr'].values

        calving_mask = np.resize(calving_mask, (12,911,544))
        socoefr = np.resize(socoefr, (911,544))

        anhax['calvingmask'] = (('time_counter', 'y', 'x'), calving_mask)
        anhax['socoefr'] = (('y', 'x'), socoefr)

        print(anhax)

        anhax.to_netcdf(save_path+new_filename)

        af.close()


if __name__ == "__main__":
   make_calving_files()

