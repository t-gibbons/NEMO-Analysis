"""
all_transport_calcuations.py
author: Tahya Weiss-Gibbons

Calculate the volume, freshwater, heat and salt transport across a section for ANHA4
Need to provide the start and end model grid points, and a straight line in the model grid will be calculated from there
Outputs netCDF files with the output

"""
import math
import datetime
import numpy as np
import xarray as xr
import netCDF4 as nc

#gets the list of grid points for a straight section between two points
#uses the bresenham line algorithm
def section_calculation(x1, x2, y1, y2):

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

    return ii, jj

def transport_calculations(runid, endyear, endmonth, endday, startyear=2004, startmonth=1, startday=5):
    figs_path = '/project/6007519/weissgib/plotting/figs/transports/'
    path = "/project/6007519/pmyers/ANHA4/ANHA4-"+runid+"-S/"
    other_path = '/project/6007519/weissgib/plotting/data_files/anha4_files/'

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
	
    #need both the u and v components of velocity
    mdl_files_v = []
    mdl_files_u = []
    mdl_files_t = []
    for t in times:
        mdl_files_v.append(path+"ANHA4-"+runid+"_y"+str(t.year)+"m"+str(t.month).zfill(2)+"d"+str(t.day).zfill(2)+"_gridV.nc")
        mdl_files_u.append(path+"ANHA4-"+runid+"_y"+str(t.year)+"m"+str(t.month).zfill(2)+"d"+str(t.day).zfill(2)+"_gridU.nc")
        mdl_files_t.append(path+"ANHA4-"+runid+"_y"+str(t.year)+"m"+str(t.month).zfill(2)+"d"+str(t.day).zfill(2)+"_gridT.nc")

    dv = xr.open_mfdataset(mdl_files_v, concat_dim='time_counter', data_vars='minimal', coords='minimal', compat='override')
    du = xr.open_mfdataset(mdl_files_u, concat_dim='time_counter', data_vars='minimal', coords='minimal', compat='override')
    dt = xr.open_mfdataset(mdl_files_t, concat_dim='time_counter', data_vars='minimal', coords='minimal', compat='override')

    #read in the mask file
    mask_file = other_path+'ANHA4_mesh_mask.nc'
    mf = nc.Dataset(mask_file)
	
    nav_lon = np.array(mf.variables['nav_lon'])
    nav_lat = np.array(mf.variables['nav_lat'])
    umask = np.array(mf.variables['umask'])
    vmask = np.array(mf.variables['vmask'])
    tmask = np.array(mf.variables['tmask'])
    e1v = np.array(mf.variables['e1v'])[0,:,:]
    e2u = np.array(mf.variables['e2u'])[0,:,:]
    
    mf.close()

    #and finally the calcualted cell thickness
    cell_thickness_u = np.load(other_path+'computed_u_thickness.npy')
    cell_thickness_v = np.load(other_path+'computed_v_thickness.npy')

    #mask the data
    dv.coords['vmask'] = (('depthv', 'y', 'x'), vmask[0,:,:,:])
    du.coords['umask'] = (('depthu', 'y', 'x'), umask[0,:,:,:])
    dt.coords['tmask'] = (('deptht', 'y_grid_T', 'x_grid_T'), tmask[0,:,:,:])

    dv = dv.where(dv.vmask == 1)
    du = du.where(du.umask == 1)
    dt = dt.where(dt.tmask == 1)

    #need the list of i/j points for the section
    #can then use this to calculate the area of the section
    #and to see which velocity components are needed at each grid space
	
    #straight line
    #section = 'fram_strait'
    #ii, jj = section_calculation(304, 360, 503, 526)

    #section = 'davis_strait'
    #ii, jj = section_calculation(175,214,443,443)

    #section = 'bering_strait'
    #ii,jj = section_calculation(222, 237, 783, 791)

    #section = 'nares_strait'
    #ii, jj = section_calculation(197,214,537,522)

    #section = 'labrador_current_2'
    #ii, jj = section_calculation(174, 199, 330, 308)

    #section = 'siberian_shelf'
    #ii, jj = section_calculation(335,391,687,671)

    #section = 'barrow_strait'
    #ii, jj = section_calculation(156,164,550,550)

    section = 'fram_south'
    ii, jj = section_calculation(327,336,508,510)

    t = du.dims['time_counter']
    total_volume = []
    total_freshwater = []
    total_heat = []
    total_salt = []

    fw_trans = []

    lon = []
    lat = []
    comp = []
    d = []

    for n in range(0,len(ii)-1):
        i1 = int(ii[n])
        i2 = int(ii[n+1])
        j1 = int(jj[n])
        j2 = int(jj[n+1])
        lon.append(nav_lon[j1, i1])
        lat.append(nav_lat[j1,i1])

        #are we going south or west?
        negative=False
        if nav_lat[j1,i1] > nav_lat[j2,i2]: negative = True #south
        if nav_lon[j1,i1] > nav_lon[j2,i2]: negative = True #west

        #so we should only either have i change or j change
        if j1 == j2:

            comp.append('v')
            
            vt = dv['vomecrty']
            v = vt.isel(y=j1, x=i1)
            v = v.rename('vel')
            v = v.rename({'depthv': 'depth'})

            cell_thickness = cell_thickness_u[:,j1,i1]
            cell_width = np.ones(v.shape)*e1v[j1,i1]

            #salinity and temperature are on the T grid, need them on v
            sal1 = dt['vosaline'].isel(y_grid_T=j1, x_grid_T=i1)
            sal2 = dt['vosaline'].isel(y_grid_T=j1, x_grid_T=i1+1)
            sal = (sal1+sal2)/2

            temp1 = dt['votemper'].isel(y_grid_T=j1, x_grid_T=i1)
            temp2 = dt['votemper'].isel(y_grid_T=j1, x_grid_T=i1+1)
            temp = (temp1+temp2)/2

        elif i1 == i1:

            #change in x
            comp.append('u')

            vt = du['vozocrtx']
            v = vt.isel(y=j1, x=i1)
            v = v.rename('vel')
            v = v.rename({'depthu': 'depth'})
            
            cell_thickness = cell_thickness_v[:,j1,i1]
            cell_width = np.ones(v.shape)*e2u[j1,i1]

            #salinity and temperature are on the T grid, need them on u
            sal1 = dt['vosaline'].isel(y_grid_T=j1, x_grid_T=i1)
            sal2 = dt['vosaline'].isel(y_grid_T=j1+1, x_grid_T=i1)
            sal = (sal1+sal2)/2

            temp1 = dt['votemper'].isel(y_grid_T=j1, x_grid_T=i1)
            temp2 = dt['votemper'].isel(y_grid_T=j1+1, x_grid_T=i1)
            temp = (temp1+temp2)/2

        sal = sal.rename({'deptht': 'depth'})
        temp = temp.rename({'deptht': 'depth'})

        #using a reference salinity of 34.8
        fwc = (34.8-sal)/34.8

        if negative:
            vol = -v*cell_thickness*cell_width
            fw = -v*fwc*cell_thickness*cell_width
            heat = -v*temp*cell_thickness*cell_width
            salt = -v*sal*cell_thickness*cell_width
        else:
            vol = v*cell_thickness*cell_width
            fw = v*fwc*cell_thickness*cell_width
            heat = v*temp*cell_thickness*cell_width
            salt = v*sal*cell_thickness*cell_width

        vl = vol.sum(dim='depth')
        f = fw.sum(dim='depth', skipna=True)
        h = heat.sum(dim='depth')
        s = salt.sum(dim='depth')

        fw_trans.append(vol)
        total_volume.append(vl)
        total_freshwater.append(f)
        total_heat.append(h)
        total_salt.append(s)

    #now sum the data at each location    
    volume_transport = total_volume[0]
    freshwater_transport = total_freshwater[0]
    heat_transport = total_heat[0]
    salt_transport = total_salt[0]

    for t in range(1, len(total_volume)):
        volume_transport = volume_transport+total_volume[t]
        freshwater_transport = freshwater_transport+total_freshwater[t]
        heat_transport = heat_transport+total_heat[t]
        salt_transport = salt_transport+total_salt[t]

    #convert to sverdrups
    volume_transport = volume_transport*0.000001
    freshwater_transport = freshwater_transport*0.000001
    heat_transport = heat_transport*0.000001
    salt_transport = salt_transport*0.000001

    print(volume_transport.mean(dim='time_counter').values)

    #and output calcualted data to a netcdf file
    volume_transport.to_netcdf(figs_path+section+'_volume_transport_'+runid+'.nc')
    freshwater_transport.to_netcdf(figs_path+section+'_freshwater_transport_'+runid+'.nc')
    heat_transport.to_netcdf(figs_path+section+'_heat_transport_'+runid+'.nc')
    salt_transport.to_netcdf(figs_path+section+'_salt_transport_'+runid+'.nc')
    
    dv.close()
    du.close()

if __name__ == "__main__":
    transport_calculations(runid='EPM151', endyear=2018, endmonth=12, endday=31)
    transport_calculations(runid='EPM152', endyear=2018, endmonth=12, endday=31)
    transport_calculations(runid='EPM014', endyear=2019, endmonth=8, endday=23)
    transport_calculations(runid='EPM015', endyear=2019, endmonth=12, endday=31)
