##---------------------------------------------------------------------------
# R1_realistic.py 
# G.Urbancic 
#
# Takes raw .dat files: 
#    (1) removes NaN values and replace by interpolation
#    (2) determines weather or not there is bad data during the day  
#        [x,y,z,T]< [35,35,35,25]
##---------------------------------------------------------------------------
 
import numpy as np
import glob

# path to ascii raw data 
datdir = '/Users/gabin/Documents/ABOA/DATA/CSAT_2m/'
 
# ascii file names
fnames = np.sort(np.array(glob.glob(datdir+'2m*.dat')))

print("number of data files:" + str(len(fnames)))

size = 20
maxm = [35,35,35,25]


for fname in fnames:

    Ddate = fname.replace(datdir,"")
    Ddate = Ddate.replace("_0000.dat","")
    
    print("*************************************")
    print(Ddate)
    print("*************************************")
     
    # Read data
    [time,x,y,z,T] = np.transpose(np.genfromtxt(fname,usecols=(0,3,5,7,9)))
    
    b_std = [0,0,0,0]
    b_mean = [0,0,0,0]
    j = 0
    for D in [x,y,z,T]:
        # Calculate mean  
        
        # correct for NaN Values 
        msk_nans = np.isnan(D)
        if len(msk_nans[msk_nans])>0:
            xdata = np.arange(len(D))
            D[msk_nans] = np.interp(xdata[msk_nans], xdata[~msk_nans],D[~msk_nans])
        
        dat = np.array([ D[x:x+size] for x in range( len(D) - size + 1 ) ])
        std = np.std(dat,axis = -1)
        mean = np.mean(dat,axis = -1)
        for i in std:
            if i < 0.0001:
                b_std[j] = b_std[j]+ 1
        b_std[j] = b_std[j]/len(std) * 100
        for i in mean:
            if i > maxm[j]: 
                b_mean[j] = b_mean[j]+1 
        b_mean[j] = b_mean[j]/len(mean) * 100
        j = j + 1
    
    print('std max error: '  + str(max([b_std[0],b_std[1],b_std[2],b_std[3]])))
    print('mean max error: '  + str(max([b_mean[0],b_mean[1],b_mean[2],b_mean[3]])))

    
      
     