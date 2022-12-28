import datetime
import numpy as np
import xarray as xr
import netCDF4 as nc
import matplotlib.pyplot as plt

def gridT_avg_calc(runid, endyear, endmonth, endday, startyear=2002, startmonth=1, startday=5):
    output_path = '//project/6007519/weissgib/plotting/data_files/surf_bio/'
    path = "/project/6007519/pmyers/ANHA4/ANHA4-"+runid+"-S/"
    grid_file = '/project/6007519/weissgib/plotting/data_files/anha4_files/ANHA4_mesh_mask.nc'
    mask_file = '/project/6007519/weissgib/plotting/data_files/anha4_files/regions_mask.nc'

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

    mesh.close()

    d.coords['mask'] = (('deptht', 'y_grid_T', 'x_grid_T'), mask[0,:,:,:])
    d = d.where(d.mask == 1)

    #just want the top 100m
    d = d.where(d['deptht'] < 100.0, drop=True)

    #calculate the weights
    n = len(d['deptht'])
    weight = np.zeros(n)
    dd = d['deptht'][n-1]
    for i in range(n):
        if i == 0:
            weight[i] = d['deptht'][i]/dd
        else:
            weight[i] = (d['deptht'][i] - d['deptht'][i-1])/dd
    
    weights = xr.DataArray(weight, coords=[d['deptht']], dims=['deptht'])
    d_weighted = d.weighted(weights)
    surface_variables = d_weighted.mean(dim='deptht', skipna=True)
    
    vooxy = surface_variables['vooxy']
    vodic = surface_variables['vodic']
    voalk = surface_variables['voalk']

    mf = nc.Dataset(mask_file)
    
    regions = {'caa_mask': 'Canadian Arctic Archipelago'}
    #regions = {'caa_mask': 'Canadian Arctic Archipelago', 'ca_mask': 'Central Arctic', 'cs_mask': 'Canadian Shelf', 'cb_mask': 'Canadian Basin', 'eb_mask': 'Eurasian Basin', 'ss_mask': 'Siberian Shelf', 'ds_mask': 'Davis Strait', 'hb_mask': 'Hudson Bay', 'bs_mask': 'Bering Strait', 'ls_mask': 'Labrador Sea', 'ns_mask': 'Nares Strait', 'fs_mask': 'Fram Strait', 'lc_mask': 'Labrador Current'}

    for r in regions.keys():
        rmask = mf[r][0,0]
        masked_vooxy = vooxy.where(rmask == 2)
        vooxy_ts = masked_vooxy.mean(('x_grid_T', 'y_grid_T'))
        vooxy_ts.to_netcdf(output_path+runid+'_vooxy_'+r+'.nc')
        
        masked_vodic = vodic.where(rmask == 2)
        vodic_ts = masked_vodic.mean(('x_grid_T', 'y_grid_T'))
        vodic_ts.to_netcdf(output_path+runid+'_vodic_'+r+'.nc')
        
        masked_voalk = voalk.where(rmask == 2)
        voalk_ts = masked_voalk.mean(('x_grid_T', 'y_grid_T'))
        voalk_ts.to_netcdf(output_path+runid+'_voalk_'+r+'.nc')
        
    mf.close()
    d.close()
    
def gridB_avg_calc(runid, endyear, endmonth, endday, startyear=2002, startmonth=1, startday=5):

    output_path = '/project/6007519/weissgib/plotting/data_files/surf_bio/'
    path = "/project/6007519/pmyers/ANHA4/ANHA4-"+runid+"-S/"
    mask_file = '/project/6007519/weissgib/plotting/data_files/anha4_files/regions_mask.nc'

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
        mdl_files.append(path+"ANHA4-"+runid+"_y"+str(t.year)+"m"+str(t.month).zfill(2)+"d"+str(t.day).zfill(2)+"_gridB.nc")

    d = xr.open_mfdataset(mdl_files, concat_dim='time_counter', data_vars='minimal', coords='minimal', compat='override')

    pco2 = d['pco2_surf']

    mf = nc.Dataset(mask_file)

    regions = {'caa_mask': 'Canadian Arctic Archipelago'}
    #regions = {'caa_mask': 'Canadian Arctic Archipelago', 'ca_mask': 'Central Arctic', 'cs_mask': 'Canadian Shelf', 'cb_mask': 'Canadian Basin', 'eb_mask': 'Eurasian Basin', 'ss_mask': 'Siberian Shelf', 'ds_mask': 'Davis Strait', 'hb_mask': 'Hudson Bay', 'bs_mask': 'Bering Strait', 'ls_mask': 'Labrador Sea', 'ns_mask': 'Nares Strait', 'fs_mask': 'Fram Strait', 'lc_mask': 'Labrador Current'}

    for r in regions.keys():
        rmask = mf[r][0,0]
        masked_pco2 = pco2.where(rmask == 2)
        pco2_ts = masked_pco2.mean(('x', 'y'))
        pco2_ts.to_netcdf(output_path+runid+'_pco2_'+r+'.nc')

    d.close()

if __name__ == "__main__":
    gridB_avg_calc(runid='EPM101', endyear=2019, endmonth=4, endday=5)
    gridB_avg_calc(runid='EPM102', endyear=2019, endmonth=6, endday=9)
    gridB_avg_calc(runid='EPM014', endyear=2019, endmonth=8, endday=23)
    gridB_avg_calc(runid='EPM015', endyear=2019, endmonth=12, endday=31)

