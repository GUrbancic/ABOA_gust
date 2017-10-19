import matplotlib.pylab as plt
import numpy as np
import os.path
from netCDF4 import Dataset
import glob
from timeit import default_timer as timer  
from scipy import stats as st

#----------------------------------------------------------------------------------------
# This script provides the spike detection method used in
# Suomi et al. (2017, QJRMS): http://onlinelibrary.wiley.com/doi/10.1002/qj.3059/pdf
#----------------------------------------------------------------------------------------

# READ DATA
datdir = '/Users/gabin/Documents/FMI/ABOA_data/ORIG_10m/'

# ascii file names
#fnames = np.sort(np.array(glob.glob(datdir+'10m-*')))
fnames = [datdir + '10m-1101' +i+ '.dat' for i in ['09','10','11','12',\
                                                        '13','14','15','16','17','18']]

print("number of data files:" + str(len(fnames)))
print("                     ")

size = 40                                                # window size for statistics
 

# Loop over ascii files:
for fname in fnames:
    
    spike_count=[[],[],[],[],[],[],[],[],[],[],[],[],[]]
    #
    # Read data
    [time,x,y,z,T] = np.transpose(np.genfromtxt(fname,usecols=(0,3,5,7,9)))
    day = np.floor(time[10]/1.0e+6)
    n = 1
    j = 0
    for D in [x,y,z,T]:
        
        print(fname + '---'+ ['x','y','z','T'][j])
        # Interpolate over possible nan values:
        msk_nans = np.isnan(D)
        if len(msk_nans[msk_nans])>0:
            xdata = np.arange(len(D))
            D[msk_nans] = np.interp(xdata[msk_nans], xdata[~msk_nans],D[~msk_nans])
            
        dat = np.array([ D[x:x+size] for x in range( len(D) - size + 1 ) ])    
            
        #calculate statistical parameters
        mn = np.mean(dat,axis = -1)
        std = np.std(dat,axis = -1)
        sN = len(dat[:,0])
        
        cov = np.array([np.corrcoef(dat[i,1:],dat[i,:-1])[1,0] for i in range(sN)])
      
        #calculate the forecasted value
        obs_old = D[(size-2):-1]
        fcst = cov*obs_old + (1.-cov)*mn 
                
        #select the observed values 
        obs_new = D[size-1:]  
        
        time1 = time[size-1:]
        
        #split into 1h sequences 
        
        u1hr = np.unique(np.floor(time1[:]/1.0e+4)*1.0e+4)            # YYMMHH0000
        
        for h in u1hr:
            if n == 1:
                spike_count[0] = spike_count[0] + [np.floor(h/1.0e+4)]
    
            msk1 = np.floor(time1[:]/1.0e+4)*1.0e+4 == h
            m = 0
            for C_spike in [2.0,3.0,4.0]:
            
                #(forcast-observed)/standard deviation Test 
                msk_spike = (np.absolute(fcst[msk1]-obs_new[msk1])/std[msk1] > C_spike)
                spk = len(msk_spike[msk_spike])
                tot = len(msk1[msk1])
                per = spk/tot *100 
                spike_count[m+n]=spike_count[m+n] + [per]
                m = m+1
        n = n+3
        j = j + 1       
  
    # Save the Data 

    np.savetxt(datdir+'low_spik_cnt_'+str(day)+'.txt',np.c_[ \
            spike_count[0], spike_count[1], spike_count[2],\
            spike_count[3], spike_count[4], spike_count[5],\
            spike_count[6], spike_count[7], spike_count[8],\
            spike_count[9], spike_count[10],spike_count[11],\
            spike_count[12]])



