import matplotlib.pylab as plt
import numpy as np
import os.path
from netCDF4 import Dataset
import glob

# Read the raw data to determine valid 10-minute intervals where values are realistic.
#change
# path to ascii data
datdir = '/Users/gabin/Documents/FMI/ABOA_data/ORIG_10m/'

# ascii file names
fnames = np.sort(np.array(glob.glob(datdir+'10m-*.dat')))

print("number of data files:" + str(len(fnames)))
print("                    ")

size = 20

xcrit = 30
ycrit = 30
zcrit = 30
Tcrit = 10

# Loop over ascii files:
for fname in fnames:
    print(fname)
    samp = [[],[],[]]
    #
    # Read data
    [time,x,y,z,T] = np.transpose(np.genfromtxt(fname,usecols=(0,3,5,7,9)))
    day = np.floor(time[10]/1.0e+6)
    #
    # Split data into 1h sequencies:x
    u1hr = np.unique(np.floor(time[:]/1.0e+4)*1.0e+4)           # YYMMHH0000
    u10min = np.unique(np.floor(time[:]/1.0e+3)*1.0e+3)         # YYMMHHm000
    
    
    ## HERE COMPUTE MEANS for x,y,z,T 
    j = 0
    mn=[[],[],[],[]]
    for D in [x,y,z,T]:
        msk_nans = np.isnan(D)
        if len(msk_nans[msk_nans])>0:
            xdata = np.arange(len(D))
            D[msk_nans] = np.interp(xdata[msk_nans], xdata[~msk_nans],D[~msk_nans])
            
        dat = np.array([ D[x:x+size] for x in range( len(D) - size + 1 ) ])    
        tim = np.array([ time[i] for i in range( len(D) - size + 1 ) ])
        #calculate statistical parameters
        mn[j] = np.mean(dat, axis = -1)  
        j = j + 1
    
    for i in u1hr:
        msk1 = np.floor(u10min/1.0e+4)*1.0e+4==i
        u10min_h = u10min[msk1]
        
        good = 0
        bad = 0
        
        for m in u10min_h :
            
            
            msk2 = np.floor(tim[:]/1.0e+3) == m/1.0e+3

            xm = mn[0][msk2]
            ym = mn[1][msk2]
            zm = mn[2][msk2]
            Tm = mn[3][msk2]
    
            for x in xm:
                if x > xcrit:
                    datgood = False 
                    break
                else: 
                    datgood = True
            for y in ym:
                if y > ycrit or datgood == False:
                    datgood = False
                    break
                else: 
                    datgood = True 
            for z in zm:
                if z > zcrit or datgood == False:
                    datgood = False
                    break
                else:
                    datgood = True
            for T in Tm:
                if T > Tcrit or datgood == False:
                    datgood = False
                    break
                else:
                    datgood = True
                    
            if datgood:
                good = good + 1
            else:
                bad = bad + 1    
                    
        samp[0] = samp[0] + [np.floor(i/1.0e+4)]
        samp[1] = samp[1] + [good]
        samp[2] = samp[2] + [bad]
     
    np.savetxt(datdir+'crit_data_'+str(day)+'.txt',np.c_[ \
            samp[0], samp[1], samp[2]])
    


