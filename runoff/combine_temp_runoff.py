"""
combine the runoff and temperature files in the right format

"""
import numpy as np
import netCDF4 as nc

years = [str(y) for y in range(1981,2019)]

for start_year in years:
    runoff_file = '/project/6007519/ANHA4-I/RUNOFF/HydroGFD/ANHA4_ReNat_HydroGFD_HBC_runoff_monthly_y'+start_year+'.nc'
    temp_file = '/project/6007519/weissgib/plotting/data_files/temp_input/ANHA4_ReNat_HydroGFD_HBC_temp_monthly_y'+start_year+'.nc'
    new_file = 'data_files/temp_input/ANHA4_ReNat_HydroGFD_HBC_runoff_monthly_y'+start_year+'.nc'

    #first copy everything from the runoff file to the new file
    with nc.Dataset(runoff_file, mode='r') as src, nc.Dataset(temp_file, mode='r') as nrc, nc.Dataset(new_file, mode='w') as dst:
        # copy global attributes all at once via dictionary
        dst.setncatts(src.__dict__)
        # copy dimensions
        for name, dimension in src.dimensions.items():
            dst.createDimension(
                name, (len(dimension) if not dimension.isunlimited() else None))
        # copy all file data 
        for name, variable in src.variables.items():
            x = dst.createVariable(name, variable.datatype, variable.dimensions)
            dst[name][:] = src[name][:]
            # copy variable attributes all at once via dictionary
            dst[name].setncatts(src[name].__dict__)

        #now need to add just the temperature information
        temp = np.array(nrc['runoff'][:])
        tm = nrc['runoff']
        temp[temp == 9999] = -999

        t = dst.createVariable("rotemper",tm.datatype,tm.dimensions)
        dst["rotemper"][:] = temp
