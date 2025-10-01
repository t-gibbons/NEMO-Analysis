"""
Create the ANHA4x runoff files and masks

Based off old ANHA4 runoff file, with slice taken from eOrca for extended boundary
"""
import numpy as np
import xarray as xr
import netCDF4 as nc

def make_runoff():

    anha_loc = '/project/6007519/ANHA4-I/RUNOFF/HydroGFD/'
    eorca_loc = '/project/6007519/eORCA025-I/Runoff/'

    start_year = 2000
    end_year = 2001

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


def make_runoff_mask():

    anha4_file = '/project/6007519/ANHA4-I/RUNOFF/HydroGFD/ANHA4_KitikmeotSea_mask.nc'
    eorca_file = '/project/6007519/weissgib/eORCA025-I/Runoff/eORCA025_KitikmeotSea_mask.nc'

    af = xr.open_dataset(anha4_file)
    ef = xr.open_dataset(eorca_file)

    #we know the dims will be y=911 and x=544 for the new anhax
    anhax = xr.Dataset(None, coords=dict(y=range(911), x=range(544)))

    for v in list(af.keys()):
        print(v)

        anha_var = af[v].values

        if v == 'lon':

            anhax['nav_lon'] = anha_var
            continue

        if v == 'lat':
        
            anhax['nav_lat'] = anha_var
            continue

        ev = ef[v].values

        #get the variable chunk to add from eorca
        eorca_chunk = ev[898:1009, 143:687]

        #now rotate to match the anha grid
        r1 = np.rot90(eorca_chunk)
        rot_chunk = np.rot90(r1)

        print(anha_var.shape)
        print(rot_chunk.shape)

        #and concatenate with anha
        anhax_var = np.concatenate((anha_var, rot_chunk), axis=0)
        print(anhax_var.shape)

        anhax[v] = (('y','x'), anhax_var)

    print(anhax)

    anhax.to_netcdf('ANHA4x_KitikmeotSea_mask.nc')

    af.close()
    ef.close()

if __name__ == "__main__":
   make_runoff()

