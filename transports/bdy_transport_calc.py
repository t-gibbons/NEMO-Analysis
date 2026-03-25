import numpy as np
import xarray as xr
import matplotlib.pyplot as plt

#define your file paths here
bdy_path = '/project/6007519/ANHA4-I/BDY/eORCA025-ECP013/'
bath_file = '/project/6007519/weissgib/ANHA4-I/ANHA4_bathy_etopo1_gebco1_smoothed_coast_corrected_mar10.nc'
grid_file = '/project/6007519/weissgib/ANHA4-I/ANHA4_mesh_mask.nc'
thickness_file = '/project/6007519/weissgib/plotting/data_files/anha4_files/computed_u_thickness.npy'
fig_dir = '/project/6007519/weissgib/plotting/transports/figs/'

#we want to calculate the volume transport in the regional boundary files
def calc_bdy_trans(year, direction):

    if direction == 'north':
        i = 799
        j = 541
    elif direction == 'south':
        i = 0
        j = 542
    else:
        print('Not a supported boundary')
        exit()

    bdy_data = bdy_path+'ANHA4_bdy_'+direction+'_eORCA025-ECP013_y'+year+'.nc'

    bf = xr.open_mfdataset(bdy_data)

    print(bf)

    #will need the grid file

    gf = xr.open_mfdataset(grid_file)

    e2v = gf['e2v'][0, i, :].values

    #get in the right shape
    e2v = np.tile(e2v.reshape(1, 1, 1, 544), (12, 50, 1, 1))

    e2v = e2v[:, :, :, 2:j+2]

    print(e2v.shape)

    gf.close()

    bt = xr.open_mfdataset(bath_file)

    bath = bt['Bathymetry'][i, :].values
    bath = np.where(bath != 0, 1, 0)
    bath = np.tile(bath.reshape(1, 1, 1, 544), (12, 50, 1, 1))

    bath = bath[:, :, :, 2:j+2]

    bt.close()

    print(bath.shape)

    #and final file is the cell thickness

    thickness = np.load(thickness_file)

    thickness = thickness[:,i,2:j+2]

    thickness = np.tile(thickness.reshape(1,50,1,j), (12, 1, 1, 1))

    print(thickness.shape)

    #and now we can calculate the transport

    print(bf['vomecrty'])

    masked_vel = bf['vomecrty']*bath

    vol_trans = masked_vel*thickness*e2v

    print(vol_trans)

    #average the transport over depth and across the line
    trans = vol_trans.sum(dim='z')
    trans = trans.sum(dim='x')

    #convert to sverdrups
    trans = trans*0.000001

    return(trans)

def plot_bdy_transport():

    #get the bdy transport for both north and south
    trans_north = calc_bdy_trans('2002', 'north')
    trans_south = calc_bdy_trans('2003', 'south')

    months = ['Jan', 'Feb', 'March', 'April', 'May', 'June', 'July', 'Aug', 'Sept', 'Oct', 'Nov', 'Dec']

    avg_trans_north = trans_north.mean('time_counter')
    avg_trans_south = trans_south.mean('time_counter')
    print("Average North BDY Transport: "+str(avg_trans_north.values[0]))
    print("Average South BDY Transport: "+str(avg_trans_south.values[0]))

    plt.plot(months, trans_north[:,0], color='green', label='North Transport')
    plt.plot(months, trans_south[:,0], color='blue', label='South Transport')

    plt.axhline(y=avg_trans_north, color='green', linestyle='dashed')
    plt.axhline(y=avg_trans_south, color='blue', linestyle='dashed')

    plt.legend()
    plt.ylabel('transport (Sv)')
    plt.title('Boundary Transport ANHA4 eORCA025-ECP013 2002')

    plt.savefig(fig_dir+'bdy_trans_eORCA025-ECP013_anha4.png')
    #plt.show()

    trans_north.close()
    trans_south.close()

plot_bdy_transport()
