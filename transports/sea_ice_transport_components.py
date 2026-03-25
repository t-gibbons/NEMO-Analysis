"""
October 2025
Author: Tahya Weiss-Gibbons (weissgib@ualberta.ca)

Calculates the speed and volume components of the sea ice transport across a set line
"""
import math
import glob
import datetime
import numpy as np
import xarray as xr
import netCDF4 as nc
import matplotlib.pyplot as plt
#import matplotlib.path as mpath
#import cartopy.crs as ccrs
#import cartopy.feature as feature

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


def transport_calculations(runid, endyear, endmonth, endday, startyear=2002, startmonth=1, startday=5, lim3=False):
    figs_path = '/project/6007519/weissgib/plotting/sea_ice/'
    if runid == 'EPM151':
        path = "/project/6007519/pmyers/ANHA4/ANHA4-"+runid+"-S/"
    else:
        path = "/project/6007519/weissgib/ANHA4/ANHA4-"+runid+"-S/"
    other_path = '/project/6007519/weissgib/plotting/data_files/anha4_files/'
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
	
    mdl_files_ice = []
    for t in times:
        mdl_files_ice.append(path+"ANHA4-"+runid+"_y"+str(t.year)+"m"+str(t.month).zfill(2)+"d"+str(t.day).zfill(2)+"_icemod.nc")

     #dont try to open the corrupted ice files
    if runid in ('EPM161', 'ETW162'):
        with open (bad_path+runid+"_bad_ice_files.txt", "r") as file:
            bad_files = eval(file.readline())

        mdl_files_ice = [x for x in mdl_files_ice if x not in bad_files]

    di = xr.open_mfdataset(mdl_files_ice, concat_dim='time_counter', data_vars='minimal', coords='minimal', compat='override')

    #read in the mask file
    mask_file = other_path+'ANHA4_mesh_mask.nc'
    mf = nc.Dataset(mask_file)
	
    nav_lon = np.array(mf.variables['nav_lon'])
    nav_lat = np.array(mf.variables['nav_lat'])
    e1t = np.array(mf.variables['e1t'])[0,:,:]
    e2t = np.array(mf.variables['e2t'])[0,:,:]
    
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

    section = 'nares_strait'
    ii, jj = section_calculation(197,214,537,522)

    #section = 'labrador_current_2'
    #ii, jj = section_calculation(174, 199, 330, 308)

    #section = 'siberian_shelf'
    #ii, jj = section_calculation(335,391,687,671)

    #section= 'greenland_north'
    #ii, jj = section_calculation(256,282,560,546)

    t = di.dims['time_counter']
    total_ice = []
    speed_ice = []
    vol_ice = []

    fw_trans = []

    lon = []
    lat = []

    #also need the pre-calculated averages of speed and volume across the section
    speed_file = xr.open_mfdataset(figs_path+section+'_seaice_transport_speed_avg_'+runid+'.nc')
    vol_file = xr.open_mfdataset(figs_path+section+'_seaice_transport_vol_avg_'+runid+'.nc')

    speed_avg = speed_file['vel']
    if lim3:
        vol_avg = vol_file['__xarray_dataarray_variable__']
    else:
        vol_avg = vol_file['iicethic']
    print(speed_avg)
    print(vol_avg)

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
            
            cell_width = np.ones(v.shape)*e1t[j1,i1]

        elif i1 == i1:

            vt = di['iicevelu']
            v = vt.isel(y=j1, x=i1)
            v = v.rename('vel')
            
            cell_width = np.ones(v.shape)*e2t[j1,i1]

        #if lim3, get the ice thickness by finding the thickness in each category
        if lim3:
            t_ice = di['iiceethick_cat']*di['ileadfra_cat']

            #and sum over the categories
            ice_thick = np.sum(t_ice, axis=1)
            ice_thick = ice_thick.isel(y=j1, x=i1)

        else:

            ice_thick = di['iicethic'].isel(y=j1, x=i1)


        if negative:
            total = -v*ice_thick*cell_width
            speed = -v*vol_avg*cell_width
            vol = -ice_thick*speed_avg*cell_width
        else:
            total = v*ice_thick*cell_width
            speed = v*vol_avg*cell_width
            vol = ice_thick*speed_avg*cell_width

        total_ice.append(total)
        speed_ice.append(speed)
        vol_ice.append(vol)

    #now sum the data at each location    
    ice_transport = total_ice[0]
    ice_speed = speed_ice[0]
    ice_vol = vol_ice[0]

    for t in range(1, len(total_ice)):
        ice_transport = ice_transport+total_ice[t]
        ice_speed = ice_speed+speed_ice[t]
        ice_vol = ice_vol+vol_ice[t]

    #and output calcualted data to a netcdf file
    #ice_transport.to_netcdf(figs_path+section+'_seaice_transport_'+runid+'.nc')
    ice_speed.to_netcdf(figs_path+section+'_seaice_transport_speed_comp_'+runid+'.nc')
    ice_vol.to_netcdf(figs_path+section+'_seaice_transport_vol_comp_'+runid+'.nc')
    
    di.close()
			
if __name__ == "__main__":
    transport_calculations(runid='ETW161', endyear=2018, endmonth=12, endday=31, lim3 = False)
    transport_calculations(runid='ETW162', endyear=2018, endmonth=12, endday=31, lim3 = True)
    transport_calculations(runid='EPM151', endyear=2018, endmonth=12, endday=31, lim3 = False)
    transport_calculations(runid='EPM161', endyear=2018, endmonth=12, endday=31, lim3 = True)
