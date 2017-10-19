import matplotlib.pylab as plt
import numpy as np
import os.path
from netCDF4 import Dataset
import glob

# path to ascii data
datdir = '/Users/gabin/Documents/FMI/ABOA_data/ORIG_10m/'
# path to netcdf data
outdir1 = '/Users/gabin/Documents/FMI/ABOA_data/G_NCDF_10m/'
outdir2 = '/Users/gabin/Documents/FMI/ABOA_data/B_NCDF_10m/'

# ascii file names
# One day at a time! 
fnames = np.sort(np.array(glob.glob(datdir+'10m-101226*.dat')))

print("number of data files:" + str(len(fnames)))
print(fnames)
hgt = 10

good_file = 0
bad_file = 0
# Loop over ascii files:
for fname in fnames:
    #
    # Read data
    [time,x,y,z,T] = np.transpose(np.genfromtxt(fname,usecols=(0,3,5,7,9)))
    #
    # Split data into 10 min sequencies:
    u10min = np.unique(np.floor(time[:]/1.0e+3-np.floor(time[:]/1.0e+6)*1.0e+3))
    for i in u10min:
        #
        # Find all data points belonging to each 10 min sequence:
        msk = np.floor(time[:]/1.0e+3-np.floor(time[:]/1.0e+6)*1.0e+3)==i
        #
        # Check that every second has 20 data points.
        sec = np.unique(time[msk]) 
        for ti in sec:
            msk2 = time[msk] == ti
            if len(msk2[msk2]) == 20 and len(sec) == 600:
                Datgood = True
            else: 
                Datgood = False 
                break
        
        # check if the number of data points is correct (10 min of 20 Hz data):
        
        if Datgood:
            print('good data:' + str(i) )
            
            good_file = good_file + 1
            t10 = time[msk]
            x10 = x[msk]
            y10 = y[msk]
            z10 = z[msk]
            T10 = T[msk]
            jno = np.arange(len(t10),dtype=int)+1 # order number of the sample
            #
            tlabel = '20'+np.floor(t10[0]/1.0e+6).astype(int).astype(str)
            if i<10:
                tlabel = tlabel+'00'+str(int(i))+'0'
            elif i<100:
                tlabel = tlabel+'0'+str(int(i))+'0'
            else:
                tlabel = tlabel+str(int(i))+'0'
            #
            #
            # create/open nc-file
            if os.path.isfile(outdir1+'AWS5_sonic_20hz_'+tlabel+'.nc'):
                nc = Dataset(outdir1+'AWS5_sonic_20hz_'+tlabel+'.nc','r+',format='NETCDF4')
            else:
                nc = Dataset(outdir1+'AWS5_sonic_20hz_'+tlabel+'.nc','w',format='NETCDF4')
            # create/open group & variables
            gkeys = nc.groups.keys()
            if (str(hgt)+'m') in gkeys:
                gg = nc.groups[str(hgt)+'m']
                st = gg.variables['time']
                sjno = gg.variables['jno']
                sX = gg.variables['X']
                sY = gg.variables['Y']
                sZ = gg.variables['Z']
                sT = gg.variables['T']
                sS = gg.variables['status']
            else:
                gg = nc.createGroup(str(hgt)+'m')
                gg.createDimension('dim',12000)
                st = gg.createVariable('time', 'i8', ('dim'))
                sjno = gg.createVariable('jno', 'i8', ('dim'))
                sX = gg.createVariable('X', 'f4', ('dim'))
                sY = gg.createVariable('Y', 'f4', ('dim'))
                sZ = gg.createVariable('Z', 'f4', ('dim'))
                sT = gg.createVariable('T', 'f4', ('dim'))
                sS = gg.createVariable('status', 'i4', ('dim'))
            #
            # save data into nc-file:
            st[:] =  t10.copy()
            sjno[:] =  jno.copy()
            sX[:] =  x10.copy()
            sY[:] =  y10.copy()
            sZ[:] =  z10.copy()
            sT[:] =  T10.copy()
            sS[:] =  np.zeros(12000)
            nc.close()
        else: 
            print('bad data:' + str(i))
            
            bad_file = bad_file + 1
            
            t10 = time[msk]
            x10 = x[msk]
            y10 = y[msk]
            z10 = z[msk]
            T10 = T[msk]
            jno = np.arange(len(t10),dtype=int)+1 # order number of the sample
            #
            tlabel = '20'+np.floor(t10[0]/1.0e+6).astype(int).astype(str)
            if i<10:
                tlabel = tlabel+'00'+str(int(i))+'0'
            elif i<100:
                tlabel = tlabel+'0'+str(int(i))+'0'
            else:
                tlabel = tlabel+str(int(i))+'0'
            #
            #
            # create/open nc-file
            if os.path.isfile(outdir2+'AWS5_sonic_20hz_'+ tlabel + '.nc'):
                nc = Dataset(outdir2+'AWS5_sonic_20hz_'+tlabel+'.nc','r+',format='NETCDF4')
            else:
                nc = Dataset(outdir2+'AWS5_sonic_20hz_'+tlabel+'.nc','w',format='NETCDF4')
            # create/open group & variables
            gkeys = nc.groups.keys()
            if (str(hgt)+'m_BAD') in gkeys:
                gg = nc.groups[str(hgt)+'m_BAD']
                st = gg.variables['time']
                sjno = gg.variables['jno']
                sX = gg.variables['X']
                sY = gg.variables['Y']
                sZ = gg.variables['Z']
                sT = gg.variables['T']
                sS = gg.variables['status']
            else:
                gg = nc.createGroup(str(hgt)+'m_BAD')
                gg.createDimension('dim',len(t10))
                st = gg.createVariable('time', 'i8', ('dim'))
                sjno = gg.createVariable('jno', 'i8', ('dim'))
                sX = gg.createVariable('X', 'f4', ('dim'))
                sY = gg.createVariable('Y', 'f4', ('dim'))
                sZ = gg.createVariable('Z', 'f4', ('dim'))
                sT = gg.createVariable('T', 'f4', ('dim'))
                sS = gg.createVariable('status', 'i4', ('dim'))
            #
            # save data into nc-file:
            st[:] =  t10.copy()
            sjno[:] = jno.copy()
            sX[:] =  x10.copy()
            sY[:] =  y10.copy()
            sZ[:] =  z10.copy()
            sT[:] =  T10.copy()
            sS[:] =  np.zeros(len(t10))
            nc.close()
            
            
print( str(np.floor(good_file/(bad_file + good_file)*100)) +'% Success'  )  
            
            
            
            