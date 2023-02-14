"""
Calcualtes the montgomery stream function
Does this for indiviual time steps, and then outputs the results for each timestep as npy files

**There is definitely a better way to do this, but works for now**
"""
import os.path
import datetime
import numpy as np
import xarray as xr
import netCDF4 as nc
from scipy.integrate import quad

T_0 = 0 #reference temperature in celcius
S_0 = 34.8 #reference salinity in psu
rho_0 = 998.23 #denisty of freshwater
P_a = 101325 #atmospheric pressure
g = 9.8 #acceleration due to gravity
z_0 = 200 #depth you want to calcualte the stream function over, in m


def rho(T,S,P):
    r = 998.23*((0.15)*(T-T_0)+(7.6*10**(-4))*(S-S_0)+(4.1*10**(-10))*(P-P_a))
    return(r)

def r_b(T,S,z):
        P=(z*(1.025*10**3)*9.8)+P_a
        rho_2 = rho(T,S,P)
        rho_1 = rho(T_0,S_0,P)
        return(rho_2-rho_1+rho_0)

def mont_stream_fn(T,S,eta,dz):
    sum = 0
    n = dz.size
    for i in range(n):
        sum = sum+(r_b(T[i],S[i],dz[i])*dz[i])
    M_b = (P_a+(g*eta+rho_0)+(g*sum)+(r_b(T[n-1], S[n-1],z_0)*g*z_0))/rho_0
    return(M_b)

def calc_stream(runid, endyear, endmonth, endday, startyear=2002, startmonth=1, startday=5):

    path = "/project/6007519/pmyers/ANHA4/ANHA4-"+runid+"-S/"
    ssh_path = '/project/6007519/weissgib/plotting/figs/ssh_plots/'
    other_path = '/project/6007519/weissgib/plotting/data_files/stream_function/'

    #read in the model files
    #will need the ssh anomaly and gridT files

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

    mdl_files_ssh = ssh_path+'ssh_anamoly_'+runid+'.nc'
    ds = nc.Dataset(mdl_files_ssh)    
    eta = np.array(ds.variables['sossheig'])
    ds.close()
    
    for k in range(len(times)):
        
        t = times[k]
        #check if already made this file
        check_file = other_path+'stream_fn_'+runid+'_'+str(t.year)+'-'+str(t.month)+'-'+str(t.day)+'.npy'
        if os.path.isfile(check_file):
            print(check_file)
            continue

        mdl_file_t = path+"ANHA4-"+runid+"_y"+str(t.year)+"m"+str(t.month).zfill(2)+"d"+str(t.day).zfill(2)+"_gridT.nc"

        dt = nc.Dataset(mdl_file_t)
        
        sal = np.array(dt.variables['vosaline'])
        temp = np.array(dt.variables['votemper'])
        depth = np.array(dt.variables['deptht'])

        x = len(dt.dimensions['x_grid_T'])
        y = len(dt.dimensions['y_grid_T'])
        
        dt.close()
    
        #figure out the closest model level to desired z
        idx = (np.abs(depth - z_0)).argmin()+1
        depth = depth[:idx]
    
        #calculate dz
        n = len(depth)
        dz = np.zeros(n)
        for i in range(n):
            if i == 0:
                dz[i] = depth[i]
            else:
                dz[i] = depth[i] - depth[i-1]
 
        stream_fn = np.zeros((y,x))

        for j in range(y):
            for i in range(x):

                s = sal[0,:idx,j,i]
                h = temp[0,:idx,j,i]

                e = eta[k,j,i]

                stream_fn[j,i] = mont_stream_fn(h,s,e,dz)
        
        np.save(other_path+'stream_fn_'+runid+'_'+str(t.year)+'-'+str(t.month)+'-'+str(t.day)+'.npy', stream_fn)
    


if __name__ == "__main__":
    calc_stream(runid='EPM101', endyear=2019, endmonth=4, endday=5)
    calc_stream(runid='EPM102', endyear=2019, endmonth=6, endday=9)
    calc_stream(runid='EPM014', endyear=2019, endmonth=8, endday=23)
    calc_stream(runid='EPM015', endyear=2019, endmonth=12, endday=31)

