"""
Make a spatial correlation between the change in the heat content and the MLD
Do this at each point to get a spatial map
"""
import xarray as xr

def heat_mld():

    #first read in the heat content files and take the difference
    heat_dir = '/mnt/storage6/tahya/model_files/'

    noheat_file = heat_dir+'EPM161_heat_content_new.nc'
    heat_file = heat_dir+'ETW162_heat_content_new.nc'

    hf = xr.open_mfdataset(heat_file)
    nf = xr.open_mfdataset(noheat_file)

    print(hf)

    diff_heat = hf-nf

    print(diff_heat)

    #now do the same for the mld
    #heat run
    heat_path = '/mnt/storage6/tahya/model_files/ANHA4-ETW162/'
    heat_files = glob.glob(heat_path+'*gridT.nc')
    heat_mld = xr.open_mfdataset(heat_files)
    datetimeindex = heat_mld.indexes['time_counter'].to_datetimeindex()
    heat_mld['time_counter'] = datetimeindex
    print(heat)

    #noheat run
    noheat_path = '/mnt/storage4/tahya/model_files/ANHA4-EPM161/'
    noheat_files = glob.glob(noheat_path+'ANHA4-EPM161*gridT.nc')
    noheat_mld = xr.open_mfdataset(noheat_files)
    datetimeindex = noheat_mld.indexes['time_counter'].to_datetimeindex()
    noheat_mld['time_counter'] = datetimeindex
    print(noheat_mld)

    diff_mld = heat_mld['somxlts'] - noheat_mld['somxlts']
    print(diff_mld)

    #now lets take the spatial correlation between the two diffs
    #take it for the whole time series to start

    correl = xr.corr(diff_heat.votemper, diff_mld.somxlts, dim="time_counter")
    print(correl)

    hf.close()
    nf.close()

    heat_mld.close()
    noheat_mld.close()

if __name__ == "__main__":
    heat_mld()
