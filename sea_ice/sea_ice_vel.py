"""
July 2024
Tahya Weiss-Gibbons (weissgib@ualberta.ca)

calculate the magnitude of the sea ice velocity at each point on the model grid
"""
import glob
import os.path
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt

def sea_ice_vel(runid):
    figs_path = '/project/6007519/weissgib/sea_ice/'
    path = "/project/6007519/weissgib/ANHA4/ANHA4-"+runid+"-S/"

    mdl_files = glob.glob(path+'ANHA4-'+runid+'*_icemod.nc')

    ice_data = xr.open_mfdataset(mdl_files)

    print(ice_data)

    ice_data['magnitude'] = np.sqrt((ice_data['iicevelv']**2) + (ice_data['iicevelu']**2))

    print(ice_data['magnitude'])

    ice_data['magnitude'][0].plot()
    plt.show()

    ice_data.close()

def vel_diff(temp_run, notemp_run):

    figs_path = '/project/6007519/weissgib/plotting/sea_ice/'
    temp_path = "/project/6007519/weissgib/ANHA4/ANHA4-"+temp_run+"-S/"
    bad_path = "/project/6007519/weissgib/plotting/data_files/bad_files/"

    if notemp_run == "EPM151":
        notemp_path = "/project/6007519/pmyers/ANHA4/ANHA4-"+notemp_run+"-S/"
    else:
        notemp_path = "/project/6007519/weissgib/ANHA4/ANHA4-"+notemp_run+"-S/"

    temp_files = glob.glob(temp_path+'ANHA4-'+temp_run+'*_icemod.nc')
    notemp_files = glob.glob(notemp_path+'ANHA4-'+notemp_run+'*_icemod.nc')

    #dont try to open the corrupted ice files
    if os.path.isfile(bad_path+temp_run+"_bad_ice_files.txt"):
        with open (bad_path+temp_run+"_bad_ice_files.txt", "r") as file:
            bad_files = eval(file.readline())

        temp_files = [x for x in temp_files if x not in bad_files]

    if os.path.isfile(bad_path+notemp_run+"_bad_ice_files.txt"):
        with open (bad_path+notemp_run+"_bad_ice_files.txt", "r") as file:
            bad_files = eval(file.readline())

        notemp_files = [x for x in notemp_files if x not in bad_files]

    temp_ice = xr.open_mfdataset(temp_files)
    notemp_ice = xr.open_mfdataset(notemp_files)

    temp_ice['magnitude'] = np.sqrt((temp_ice['iicevelv']**2) + (temp_ice['iicevelu']**2))
    notemp_ice['magnitude'] = np.sqrt((notemp_ice['iicevelv']**2) + (notemp_ice['iicevelu']**2))

    diff_mag = temp_ice['magnitude'] - notemp_ice['magnitude']

    print(diff_mag)

    diff_mag.to_netcdf(figs_path+temp_run+"-"+notemp_run+"_sea_ice_mag.nc")

    temp_ice.close()
    notemp_ice.close()

if __name__ == "__main__":
    vel_diff('ETW161', 'EPM151')

