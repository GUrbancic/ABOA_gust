##---------------------------------------------------------------------------
# R1_filter_spike.py 
# G.Urbancic 
#
# Takes NetCDF files produced by R1_conv_ncdf: 
#    (1) Filters spikes 
#    (2) Builds groups: spk_min  
#    (3) Saves clean data in group CLEAN
##---------------------------------------------------------------------------

import numpy as np
import matplotlib.pyplot as plt
import os.path
from netCDF4 import Dataset
import glob

# path to NetCDF data 
datdir = '/Users/gabin/Documents/ABOA/DATA/CSAT_2m/'

# NetCDF file names
fnames = np.sort(np.array(glob.glob(datdir+'*.nc')))

print("number of NetCDF files:" + str(len(fnames)))

# Spike correction parameters:
size = 80                 
C_spike=[4.0,4.0,4.0,5.0]

for fname in fnames:

    Ddate = fname.replace(datdir,"")
    Ddate = Ddate.replace("_0000.dat","")
    
    print("*************************************")
    print(Ddate)
    print("*************************************")  
    
    # open NetCDF file 
    nc = Dataset(fname,'r',format='NETCDF4')
    gkeys = nc.groups.keys()
   
    # Read in RAW data 
    gg = nc.groups['RAW']
    tr = gg.variables['time']
    Xr = gg.variables['X']
    Yr = gg.variables['Y']
    Zr = gg.variables['Z']
    Tr = gg.variables['T']
    Fr = gg.variables['F']
    
    t = tr[:]
    X = Xr[:]
    Y = Yr[:]
    Z = Zr[:]
    T = Tr[:]
    F = Fr[:]
    
    # -------------------------------------------------------------------------
    # Despiking routine + spike count: 
    
    u1min = np.unique(np.floor(tr[:]/1.0e+2))
    spike_count=[u1min,[],[],[],[]]
        
    j = 0
    for D in [X,Y,Z,T]:
        
        # Calculate statistical parameters
        # index 1 represents main averaging (forward) until len(D) - size + 1
        # index 2 represents end averaging (backward) after len(D) - size + 1
        
        dat1 = np.array([ D[x:x+size] for x in range( len(D) - size + 1 ) ])
        dat2 = np.array([ D[x-size:x] for x in range( len(D) )[len(D)-size +1:]])
        
        # Calculate statistics on 'size' intervals
        
        mean1 = np.mean(dat1,axis = -1)
        mean2 = np.mean(dat2,axis = -1)
        
        std1 = np.std(dat1,axis = -1)
        std2 = np.std(dat2,axis = -1)
        
        sN1 = len(dat1[:,0])
        sN2 = len(dat2[:,0])

        cov1 = np.array([np.corrcoef(dat1[i,1:],dat1[i,:-1])[1,0] for i in range(sN1)])
        cov2 = np.array([np.corrcoef(dat2[i,1:],dat2[i,:-1])[1,0] for i in range(sN2)])        
      
        #calculate the forecasted value 
        prev1 = D[:len(D)-size]
        prev2 = D[len(D)-size:-1]
        mean1 = mean1[1:]
        mean2 = mean2
        fcst1 = cov1[1:]*prev1 + (1.-cov1[1:])*mean1
        fcst2 = cov2*prev2 + (1.-cov2)*mean2
        obs1 = D[1:len(D)-size+1]
        obs2 = D[len(D)-size+1:]
        
        # Spike test
        msk_spikes1 = np.absolute(fcst1[:]-obs1[:])/std1[1:] > C_spike[j]
        msk_spikes2 = np.absolute(fcst2-obs2)/std2 > C_spike[j] 
        
        # Count spikes per minute
        spike_tot = np.append(msk_spikes1,msk_spikes2)
        for ti in u1min:
            msktime = np.floor(tr[:]/1.0e+2) == ti 
            spike_min = spike_tot[msktime[1:]]
            spike_count[j+1]=spike_count[j+1]+ [len(spike_min[spike_min])/len(msktime[msktime]) * 100]
            
        
        # Calculate the number of new observed spikes
        N_spikes1 = len(msk_spikes1[msk_spikes1]) 
        N_spikes2 = len(msk_spikes2[msk_spikes2])         
       
        # Replace spikes by interpolation
        xdata1 = np.arange(len(D[1:len(D)-size+1]))
        xdata2 = np.arange(len(D[len(D)-size+1:]))
        if (N_spikes1>0):
            D[1:len(D)-size+1][msk_spikes1] = \
            np.interp(xdata1[msk_spikes1], xdata1[~msk_spikes1], D[1:len(D)-size+1][~msk_spikes1])
            F[1:len(D)-size+1][msk_spikes1] = [-2]*int(len(msk_spikes1[msk_spikes1]))
        if (N_spikes2>0):
            D[len(D)-size+1:][msk_spikes2] = \
            np.interp(xdata2[msk_spikes2], xdata2[~msk_spikes2], D[len(D)-size+1:][~msk_spikes2]) 
            F[len(D)-size+1:][msk_spikes2] = [-2]*int(len(msk_spikes2[msk_spikes2]))
        j = j + 1
       
    maxx = max(spike_count[1])  
    maxy = max(spike_count[2])
    maxz = max(spike_count[3])
    maxT = max(spike_count[4])
    
    print("                   ")
    print("WORST MINUTE STATS:")
    print("max % spike" + str(max([maxx,maxy,maxz,maxT])))
        
        
    #--------------------------------------------------------------------------
    # Build or edit CLEAN, spk_min 
    #--------------------------------------------
    
    # open NetCDF file 
    #nc2 = Dataset(fname,'r+',format='NETCDF4')
    
    if 'CLEAN' in gkeys:

    else:
        cc = nc.createGroup('CLEAN')
        cc.createDimension('dim',len(tr))
        tc = cc.createVariable('time', 'i8', 'dim')
        Xc = cc.createVariable('X', 'f4', 'dim')
        Yc = cc.createVariable('Y', 'f4', 'dim')
        Zc = cc.createVariable('Z', 'f4', 'dim')
        Tc = cc.createVariable('T', 'f4', 'dim')
        Fc = cc.createVariable('F', 'i8', 'dim')

    # save data into nc-file:
    tc[:] =  t
    Xc[:] =  X
    Yc[:] =  Y
    Zc[:] =  Z
    Tc[:] =  T
    Fc[:] =  F

    #--------------------------------------------
    if 'spk_min' in gkeys: 
        gspk = nc.groups['spk_min']
        ts = gspk.variables['t']
        Xs = gspk.variables['X']
        Ys = gspk.variables['Y']
        Zs = gspk.variables['Z']
        Ts = gspk.variables['T'] 
        
    else:
        gspk = nc.createGroup('spk_min')
        gspk.createDimension('dim2',len(spike_count[0]))
        ts = gspk.createVariable('t', 'i8', 'dim2')
        Xs = gspk.createVariable('X', 'f4', 'dim2')
        Ys = gspk.createVariable('Y', 'f4', 'dim2')
        Zs = gspk.createVariable('Z', 'f4', 'dim2')
        Ts = gspk.createVariable('T', 'f4', 'dim2')  
        
    ts[:] = spike_count[0]
    Xs[:] = spike_count[1]
    Ys[:] = spike_count[2]
    Zs[:] = spike_count[3]
    Ts[:] = spike_count[4]
        
    nc.close()  

        
        