##---------------------------------------------------------------------------
# R1_conv_ncdf.py 
# G.Urbancic 
#
# Takes raw .dat files: 
#    (1) removes NaN values and replace by interpolation
#    (2) converts flags into decimal format 
#    (3) saves data into NetCDF 
##---------------------------------------------------------------------------
 
import numpy as np
import os.path
from netCDF4 import Dataset
import glob

# path to ascii raw data 
datdir = '/Users/gabin/Documents/ABOA/DATA/CSAT_10m/'

# ascii file names
fnames = np.sort(np.array(glob.glob(datdir+'10m*.dat')))

print("number of data files:" + str(len(fnames)))

for fname in fnames:

    Ddate = fname.replace(datdir,"")
    Ddate = Ddate.replace("_0000.dat","")
    
    print("*************************************")
    print(Ddate)
    print("*************************************")
     
    # Read data
    [time,x,y,z,T] = np.transpose(np.genfromtxt(fname,usecols=(0,3,5,7,9)))
    flags = np.genfromtxt(fname, dtype=str,usecols=(10)) 
    
    #--------------------------------------------------------------------------
    # Reformatting the flags 

    j = 0
    for i in flags :
        flags[j] = int(i.replace("F=",""),16)
        j = j+1
    Flags = flags.astype(int)

    #--------------------------------------------------------------------------
    # Remove NaN and replace by interpolation, count % NaN/min. 
    
    u1min = np.unique(np.floor(time[:]/1.0e+2))         # YYMMHHmm
    nan_count = [[],[],[],[],[]]                        # count (time,x,y,z,T) 
    
    for i in u1min:
        # isolate one minute intervals
        msk1 = np.floor(time/100) == i
        nan_count[0] = nan_count[0]+[i]
        
        x1 = x[msk1]
        y1 = y[msk1]
        z1 = z[msk1]
        T1 = T[msk1]
        flag1 = Flags[msk1] 
        
        j = 0
        for D in [x1,y1,z1,T1]:
        
            # correct for NaN Values 
            msk_nans = np.isnan(D)
            if len(msk_nans[msk_nans])>0:
                xdata = np.arange(len(D))
                D[msk_nans] = np.interp(xdata[msk_nans], xdata[~msk_nans],D[~msk_nans])
                flag1[msk_nans] = [-1]*int(len(msk_nans[msk_nans]))
                
            # compute percentage NaN
            perc_nan = (len(msk_nans[msk_nans])/len(D)) * 100
            nan_count[j+1] = nan_count[j+1] + [perc_nan]
            j = j+1 
            
        # update the full data set for corrected one minute interpolations   
        x[msk1] = x1
        y[msk1] = y1 
        z[msk1] = z1
        T[msk1] = T1 
        Flags[msk1] = flag1
    
    maxx = max(nan_count[1])
    maxy = max(nan_count[2])
    maxz = max(nan_count[3])
    maxT = max(nan_count[4])

    print("WORST MINUTE STATS:")
    print("max % NaN: " + str(max([maxx,maxy,maxz,maxT] )))

    #--------------------------------------------------------------------------
    # create/open nc-file
        
    if os.path.isfile(datdir +'24h_'+ Ddate + '.nc'):
        nc = Dataset(datdir +'24h_'+ Ddate + '.nc','r+',format='NETCDF4')
    else:
        nc = Dataset(datdir +'24h_'+ Ddate + '.nc','w',format='NETCDF4')
            
    # create/open group & variables
    gkeys = nc.groups.keys()
    if 'RAW' in gkeys:
        gg = nc.groups['RAW']
        st = gg.variables['time']
        sX = gg.variables['X']
        sY = gg.variables['Y']
        sZ = gg.variables['Z']
        sT = gg.variables['T']
        sF = gg.variables['F']

    else:
        gg = nc.createGroup('RAW')
        gg.createDimension('dim',len(time))
        st = gg.createVariable('time', 'i8', 'dim')
        sX = gg.createVariable('X', 'f4', 'dim')
        sY = gg.createVariable('Y', 'f4', 'dim')
        sZ = gg.createVariable('Z', 'f4', 'dim')
        sT = gg.createVariable('T', 'f4', 'dim')
        sF = gg.createVariable('F', 'i8', 'dim')
        
    if 'NaN_min' in gkeys: 
        gnan = nc.groups['NaN_min']
        t1 = gnan.variables['t']
        Xn = gnan.variables['X']
        Yn = gnan.variables['Y']
        Zn = gnan.variables['Z']
        Tn = gnan.variables['T'] 
        
    else:
        gnan = nc.createGroup('NaN_min')
        gnan.createDimension('dim2',len(nan_count[0]))
        t1 = gnan.createVariable('t', 'i8', 'dim2')
        Xn = gnan.createVariable('X', 'f4', 'dim2')
        Yn = gnan.createVariable('Y', 'f4', 'dim2')
        Zn = gnan.createVariable('Z', 'f4', 'dim2')
        Tn = gnan.createVariable('T', 'f4', 'dim2')
            
    # save data into nc-file:
    st[:] =  time
    sX[:] =  x
    sY[:] =  y
    sZ[:] =  z
    sT[:] =  T
    sF[:] =  Flags
    t1[:] = nan_count[0]
    Xn[:] = nan_count[1]
    Yn[:] = nan_count[2]
    Zn[:] = nan_count[3]
    Tn[:] = nan_count[4]
    
    nc.close()
    
    


