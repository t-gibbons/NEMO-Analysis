"""
all_transport_calcuations.py
author: Tahya Weiss-Gibbons

Calculate the volume, freshwater, heat and salt transport across a section for ANHA4
Need to provide the start and end model grid points, and a straight line in the model grid will be calculated from there
Outputs netCDF files with the output

NOTE: This is a slightly updated version, which includes some changes to variable names and time decoding to be consistent with new NEMO5 output

"""
import math
import datetime
import numpy as np
import xarray as xr
import netCDF4 as nc

#dont try to open the corrupted files
def remove_bad_files(bad_path, runid, file_type, mdl_files):
    with open (bad_path+runid+"_bad_"+file_type+"_files.txt", "r") as file:
        bad_files = eval(file.readline())

    return [x for x in mdl_files if x not in bad_files]

#gets the list of grid points for a straight section between two points
#uses the bresenham line algorithm

def plotLineLow(x0, x1, y0, y1):

    print("low line!")

    ii = []
    jj = []

    dx = x1 -x0
    dy = y1 - y0
    yi = 1

    if dy < 0:
        yi = -1
        dy = -dy

    D = (2*dy) - dx
    y = y0

    for x in range(x0, x1):
        ii.append(x)
        jj.append(y)
        if D > 0:
            y = y+yi
            D = D+(2*(dy-dx))
            ii.append(x)
            jj.append(y)
        else:
            D = D+2*dy

    return ii, jj


def plotLineHigh(x0, x1, y0, y1):

    print("high line!")

    ii = []
    jj = []

    dx = x1 - x0
    dy = y1 - y0
    xi = 1

    if dx < 0:
        xi = -1
        dx = -dx
    D = (2*dx) - dy
    x = x0

    for y in range(y0, y1):
        ii.append(x)
        jj.append(y)
        if D > 0:
            x = x+xi
            D = D + (2*(dx-dy))
            ii.append(x)
            jj.append(y)
        else:
            D = D + 2*dx

    return ii, jj


def section_calculation(x0, x1, y0, y1):

    if abs(y1-y0) < abs(x1-x0):
        if x0 > x1:
            ii, jj = plotLineLow(x1, x0, y1, y0)
        else:
            ii, jj = plotLineLow(x0, x1, y0, y1)

    else:
            if y0 > y1:
                ii, jj = plotLineHigh(x1, x0, y1, y0)
            else:
                ii, jj = plotLineHigh(x0, x1, y0, y1)

    return ii, jj


def transport_calculations(runid, endyear, endmonth, endday, startyear=2004, startmonth=1, startday=5):
    figs_path = '/project/6007519/weissgib/plotting/transports/'
    #path = "/project/6007519/weissgib/ANHA4x/ANHA4x-"+runid+"-S/"
    path = "/project/6007519/weissgib/eORCA025/eORCA025-"+runid+"-S/"
    other_path = '/project/6007519/weissgib/plotting/data_files/eorca_files/'
    bad_path = '/project/6007519/weissgib/plotting/data_files/bad_files/'

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
        mdl_files_v.append(path+"eORCA025-"+runid+"_y"+str(t.year)+"m"+str(t.month).zfill(2)+"d"+str(t.day).zfill(2)+"_gridV.nc")
        mdl_files_u.append(path+"eORCA025-"+runid+"_y"+str(t.year)+"m"+str(t.month).zfill(2)+"d"+str(t.day).zfill(2)+"_gridU.nc")
        mdl_files_t.append(path+"eORCA025-"+runid+"_y"+str(t.year)+"m"+str(t.month).zfill(2)+"d"+str(t.day).zfill(2)+"_gridT.nc")

    #mdl_files_v = remove_bad_files(bad_path, runid, 'v', mdl_files_v)
    #mdl_files_u = remove_bad_files(bad_path, runid, 'u', mdl_files_u)
    #mdl_files_t = remove_bad_files(bad_path, runid, 't', mdl_files_t)

    dv = xr.open_mfdataset(mdl_files_v, concat_dim='time_counter', data_vars='minimal', coords='minimal', compat='override', decode_times=False)
    du = xr.open_mfdataset(mdl_files_u, concat_dim='time_counter', data_vars='minimal', coords='minimal', compat='override', decode_times=False)
    dt = xr.open_mfdataset(mdl_files_t, concat_dim='time_counter', data_vars='minimal', coords='minimal', compat='override', decode_times=False)

    time_convert = nc.num2date(dv['time_counter'].values, dv['time_counter'].units, dv['time_counter'].calendar)
    dv['time_counter'] = time_convert

    time_convert = nc.num2date(du['time_counter'].values, du['time_counter'].units, du['time_counter'].calendar)
    du['time_counter'] = time_convert

    time_convert = nc.num2date(dt['time_counter'].values, dt['time_counter'].units, dt['time_counter'].calendar)
    dt['time_counter'] = time_convert

    #take the monthly average to deal with missing files
    dv = dv.resample(time_counter='M').mean()
    du = du.resample(time_counter='M').mean()
    dt = dt.resample(time_counter='M').mean()

    print(dv)
    print(du)
    print(dt)

    #read in the mask file
    #mask_file = other_path+'ANHA4_mesh_mask.nc'
    mask_file = other_path+'ANHA4x_mesh_mask.nc'
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
    cell_thickness_u = np.load(other_path+'computed_u_thickness_anha4x.npy')
    cell_thickness_v = np.load(other_path+'computed_v_thickness_anha4x.npy')

    #mask the data
    dv.coords['vmask'] = (('depthv', 'y_grid_V', 'x_grid_V'), vmask[0,:,:,:])
    du.coords['umask'] = (('depthu', 'y_grid_U', 'x_grid_U'), umask[0,:,:,:])
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

    #section = 'davis_strait_south'
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

    #section = 'fram_south'
    #ii, jj = section_calculation(327,336,508,510)
    
    #section = 'barent_sea'
    #ii, jj = section_calculation(371, 413, 512, 480)

    #section = 'mack_coast'
    #ii, jj = section_calculation(140, 154, 716, 709)

    #section = 'sib_coast_2'
    #ii, jj = section_calculation(395, 380, 706, 686)

    section = 'boundary_test'
    ii, jj = section_calculation(33, 429, 900, 900)

    t = du.dims['time_counter']
    total_volume = []
    total_freshwater = []
    total_heat = []
    total_salt = []

    #if you want to look at the section, use these too
    sec_volume = []
    sec_freshwater = []
    sec_heat = []
    sec_salt = []

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

        #if negative == False: continue #just want the southward component for davis strait

        #so we should only either have i change or j change
        if j1 == j2:

            comp.append('v')
            
            vt = dv['vomecrty']
            v = vt.isel(y_grid_V=j1, x_grid_V=i1)
            v = v.rename('vel')
            v = v.rename({'depthv': 'depth', 'nav_lon_grid_V': 'nav_lon', 'nav_lat_grid_V': 'nav_lat', 'vmask': 'mask'})

            cell_thickness = cell_thickness_u[:,j1,i1]
            cell_width = np.ones(v.shape)*e1v[j1,i1]

            #salinity and temperature are on the T grid, need them on v
            sal1 = dt['vosaline'].isel(y_grid_T=j1, x_grid_T=i1)
            sal2 = dt['vosaline'].isel(y_grid_T=j1, x_grid_T=i1+1)
            sal = (sal1+sal2)/2

            temp1 = dt['votemper'].isel(y_grid_T=j1, x_grid_T=i1)
            temp2 = dt['votemper'].isel(y_grid_T=j1, x_grid_T=i1+1)
            temp = (temp1+temp2)/2

        elif i1 == i2:

            #change in x
            comp.append('u')

            vt = du['vozocrtx']
            v = vt.isel(y_grid_U=j1, x_grid_U=i1)
            v = v.rename('vel')
            v = v.rename({'depthu': 'depth', 'nav_lon_grid_U': 'nav_lon', 'nav_lat_grid_U': 'nav_lat', 'umask': 'mask'})
            
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

        #NOTE if you want to look at the transport section, don't sum over depth
        sec_volume.append(vol)
        sec_freshwater.append(fw)
        sec_heat.append(heat)
        sec_salt.append(salt)
        
        vl = vol.sum(dim='depth')
        f = fw.sum(dim='depth', skipna=True)
        h = heat.sum(dim='depth')
        s = salt.sum(dim='depth')

        fw_trans.append(vol)
        total_volume.append(vl)
        total_freshwater.append(f)
        total_heat.append(h)
        total_salt.append(s)


    #also for transport section, don't sum over the section
    #we'll just combine on a new dimension
    volume_section = xr.concat(sec_volume, dim='sec_loc')
    #freshwater_section = xr.concat(sec_freshwater, dim='sec_loc')
    #heat_section = xr.concat(sec_heat, dim='sec_loc')
    #salt_section = xr.concat(sec_salt, dim='sec_loc')

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

    #and output calcualted data to a netcdf file
    volume_transport.to_netcdf(figs_path+section+'_volume_transport_anha4x_'+runid+'.nc')
    #freshwater_transport.to_netcdf(figs_path+section+'_freshwater_transport_anha4x_'+runid+'.nc')
    #heat_transport.to_netcdf(figs_path+section+'_heat_transport_anha4x_'+runid+'.nc')
    #salt_transport.to_netcdf(figs_path+section+'_salt_transport_anha4x_'+runid+'.nc')

    #if you made it, also save the section data
    volume_section.to_netcdf(figs_path+section+'_volume_section_anha4x_'+runid+'.nc')
    #freshwater_section.to_netcdf(figs_path+section+'_freshwater_section_anha4x_'+runid+'.nc')
    #heat_section.to_netcdf(figs_path+section+'_heat_section_anha4x_'+runid+'.nc')
    #salt_section.to_netcdf(figs_path+section+'_salt_section_anha4x_'+runid+'.nc')
    
    dv.close()
    du.close()

if __name__ == "__main__":
    #transport_calculations(runid='ETW502', endyear=2003, endmonth=12, endday=31, startyear=2002)
    transport_calculations(runid='ETW501', endyear=2003, endmonth=12, endday=31, startyear=2002)
