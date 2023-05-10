import netCDF4 as nc
import numpy as np
import numpy.ma as ma

#read in the model grid info
mask_file = '/project/6007519/weissgib/plotting/data_files/anha4_files/ANHA4_mesh_mask.nc'
mf = nc.Dataset(mask_file)

x = len(mf.dimensions['x'])
y = len(mf.dimensions['y'])
z = len(mf.dimensions['z'])

e3u = np.array(mf.variables['e3u'])[0,:,:,:]
e3v = np.array(mf.variables['e3v'])[0,:,:,:]

mf.close()

#and we need the bathymetry info
bathy_file = '/project/6007519/weissgib/plotting/data_files/anha4_files/ANHA4_bathymetry.nc'
bf = nc.Dataset(bathy_file)

bathyu = np.array(bf.variables['bathyu'])
bathyv = np.array(bf.variables['bathyv'])

mf.close()

cu_cell_thickness = np.zeros((z,y,x))
cv_cell_thickness = np.zeros((z,y,x))

#now compute the cell thickness
for j in range(y):
    for i in range(x):
        depthu = bathyu[j,i]
        depthv = bathyv[j,i]

        pu_thickness = e3u[:,j,i]
        pv_thickness = e3v[:,j,i]

        tu = 0
        tv = 0

        u_comp = False
        v_comp = False

        for k in range(z):
            tu = tu+pu_thickness[k]
            tv = tv+pv_thickness[k]

            if tu > depthu:
                if u_comp:
                    cu_cell_thickness[k,j,i] = 0
                else:
                    cu_cell_thickness[k,j,i] = tu-depthu
                    u_comp = True
            elif tu < depthu:
                cu_cell_thickness[k,j,i] = e3u[k,j,i]

            if tv > depthv:
                if v_comp:
                    cv_cell_thickness[k,j,i] = 0
                else:
                    cv_cell_thickness[k,j,i] = tv-depthv
                    v_comp = True
            elif tv < depthv:
                cv_cell_thickness[k,j,i] = e3v[k,j,i]

#save calculated cell thickness
np.save('/project/6007519/weissgib/plotting/computed_u_thickness.npy', cu_cell_thickness)
np.save('/project/6007519/weissgib/plotting/computed_v_thickness.npy', cv_cell_thickness)

