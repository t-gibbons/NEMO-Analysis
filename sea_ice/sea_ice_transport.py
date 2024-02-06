"""
October 2022
Author: Tahya Weiss-Gibbons (weissgib@ualberta.ca)

Calculates the sea ice transport across a section
Need to know the begining and end grid points of your line, as well as time period and model run to calculate for
"""
import math
import glob
import datetime
import numpy as np
import xarray as xr
import netCDF4 as nc
import matplotlib.pyplot as plt
import matplotlib.path as mpath
import cartopy.crs as ccrs
import cartopy.feature as feature

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

def transport_calculations(runid, endyear, endmonth, endday, lim3=False, startyear=2002, startmonth=1, startday=5):
    figs_path = '/project/6007519/weissgib/plotting/transports/'
    
    other_path = '/project/6007519/weissgib/plotting/data_files/anha4_files/'

    start_time = datetime.date(startyear, startmonth, startday)

    end_time = datetime.date(endyear, endmonth, endday)

    #figure out all the dates we have model files
    delta = end_time - start_time
    times = []
    """
    i = 0
    while i < delta.days+1:
        t = start_time + datetime.timedelta(days=i)
        if t.month == 2 and t.day == 29:
            t = datetime.date(t.year, 3, 1)
            i = i+6
        else:
            i = i+5
        times.append(t)
	
    mdl_files_ice = []
    for t in times:
        mdl_files_ice.append(path+"ANHA4-"+runid+"_y"+str(t.year)+"m"+str(t.month).zfill(2)+"d"+str(t.day).zfill(2)+"_icemod.nc")
    """
    path = "/project/6007519/weissgib/ANHA4/ANHA4-"+runid+"-S/"

    mdl_files_ice = glob.glob(path+'ANHA4-'+runid+'*_icemod.nc')
    di = xr.open_mfdataset(mdl_files_ice, concat_dim='time_counter', data_vars='minimal', coords='minimal', compat='override')

    #read in the mask file
    mask_file = other_path+'ANHA4_mesh_mask.nc'
    mf = nc.Dataset(mask_file)
	
    nav_lon = np.array(mf.variables['nav_lon'])
    nav_lat = np.array(mf.variables['nav_lat'])
    e1v = np.array(mf.variables['e1v'])[0,:,:]
    e2u = np.array(mf.variables['e2u'])[0,:,:]
    
    mf.close()

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

    section= 'greenland_north'
    ii, jj = section_calculation(256,282,560,546)

    t = di.dims['time_counter']
    total_ice = []

    fw_trans = []

    lon = []
    lat = []

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
            
            vt = di['iicevelv']
            v = vt.isel(y=j1, x=i1)
            v = v.rename('vel')
            
            cell_width = np.ones(v.shape)*e1v[j1,i1]

        elif i1 == i1:

            vt = di['iicevelu']
            v = vt.isel(y=j1, x=i1)
            v = v.rename('vel')
            
            cell_width = np.ones(v.shape)*e2u[j1,i1]

        #if lim3, get the ice thickness by finding the thickness in each category
        if lim3:
            t_ice = di['iiceethick_cat']*di['ileadfra_cat']

            #and sum over the categories
            ice_thick = np.sum(t_ice, axis=1)
            ice_thick = ice_thick.isel(y=j1, x=i1)

        else:

            ice_thick = di['iicethic'].isel(y=j1, x=i1)

        ice_frac = di['ileadfra'].isel(y=j1, x=i1)

        if negative:
            vol = -v*ice_thick*cell_width*ice_frac
        else:
            vol = v*ice_thick*cell_width*ice_frac

        total_ice.append(vol)

    #now sum the data at each location
    print(total_ice)
    ice_transport = total_ice[0]

    for t in range(1, len(total_ice)):
        ice_transport = ice_transport+total_ice[t]

    #and output calcualted data to a netcdf file
    ice_transport.to_netcdf(figs_path+section+'_seaice_transport_'+runid+'.nc')
    
    di.close()
			
if __name__ == "__main__":
    transport_calculations(runid='EPM161', endyear=2018, endmonth=12, endday=31, lim3=True)
