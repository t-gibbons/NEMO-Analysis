import glob
import numpy as np
import xarray as xr

def mean_kinetic_energy(runid):
    rho = 0.91 #density of sea ice, using average estimate of 0.91 Mg m^-3
    save_path = '/project/6007519/weissgib/plotting/sea_ice/'

    #first get the sea ice volume
    vol_path = '/project/6007519/weissgib/plotting/sea_ice/sea_ice_volume_'+runid+'.nc'

    ice_vol = xr.open_mfdataset(vol_path)
    print(ice_vol)

    #and now the velocity data
    model_path = "/project/6007519/pmyers/ANHA4/ANHA4-"+runid+"-S/"

    mdl_files = glob.glob(model_path+'ANHA4-'+runid+'*_icemod.nc')

    ice_data = xr.open_mfdataset(mdl_files)

    print(ice_data)

    ice_mag = np.sqrt((ice_data['iicevelv']**2) + (ice_data['iicevelu']**2))
    print(ice_mag)

    #now calculate the mean kinetic energy
    mke_ice = 0.5*ice_vol*rho*ice_mag
    print(mke_ice)

    mke_ice.to_netcdf(save_path+'mke_ice_'+runid+'.nc')

    ice_vol.close()
    ice_data.close()

def region_mke(runids = []):

    path = '/project/6007519/weissgib/plotting/sea_ice/'
    mask_file = '/project/6007519/weissgib/plotting/data_files/anha4_files/runoff_comp_regions.nc'

    mask_data = xr.open_mfdataset(mask_file)
    mask_data = mask_data.rename({'x': 'x_grid_T', 'y': 'y_grid_T'})

    regions = {'full_arctic': 'Arctic'}

    experiment = []
    mke = []
    date = []
    region = []
    for r in runids:

        try:
            mke_data = xr.open_mfdataset(path+'mke_ice_'+r+'.nc')
        except:
            print("Opps, need to calculate MKE for "+r+" first...")
            print("Excluding from plot and continuing")
            continue

        #mask the data for each region
        for reg in regions:
            md = mask_data[reg][0,0]
            masked_mke = mke_data.where(md ==2)

            print(masked_mke)

            ts = masked_mke.sum(('x_grid_T', 'y_grid_T'))
            print(ts)

            ts['__xarray_dataarray_variable__'].plot(x='time_counter', label=r)

        plt.title('Arctic')
        plt.legend()
        plt.tight_layout()
        plt.show()

        plt.clf()



if __name__ =="__main__":

    region_mke(runids=['EPM151', 'ETW161'])
