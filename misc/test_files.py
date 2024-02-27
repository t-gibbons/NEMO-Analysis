import glob
import datetime
import xarray as xr

runid = 'EPM161'
path = '/project/6007519/weissgib/ANHA4/ANHA4-'+runid+'-S/'

#mdl_files = glob.glob(path+'ANHA4-'+runid+'*.nc')

startyear = 2002
startmonth = 1
startday = 5

endyear = 2018
endmonth = 12
endday = 31

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
    mdl_files.append(path+"ANHA4-"+runid+"_y"+str(t.year)+"m"+str(t.month).zfill(2)+"d"+str(t.day).zfill(2)+"_icemod.nc")

bad_files = 0

for m in mdl_files:
    try:
        d = xr.open_mfdataset(m)
        d.close()
    except:
        bad_files = bad_files+1
        print(m)

if bad_files == 0:
    print('Successfully opened all files!')
else:
    print("Couldn't open "+str(bad_files)+" files")
