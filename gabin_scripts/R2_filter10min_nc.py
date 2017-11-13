import numpy as np
import os.path
from netCDF4 import Dataset
import glob

# path to ascii raw data 
datdir = '/Users/gabin/Documents/FMI/ABOA_data/ORIG_10m/'
# path to netcdf data
outdir1 = '/Users/gabin/Documents/FMI/ABOA_data/G_NCDF_10m/'
outdir2 = '/Users/gabin/Documents/FMI/ABOA_data/B_NCDF_10m/'

# ascii file names
fnames = np.sort(np.array(glob.glob(datdir+'10m*.dat')))
print("number of data files:" + str(len(fnames)))

# file count 
good_file = 0
bad_file = 0

# window size - 2 seconds
size = 40                 

# correlation-based (SUOMI)  forecasting parameters. C_sp.. = [x,y,z,T]
C_spike=[4.0,4.0,4.0,5.0]
xcrit = 30.0
ycrit = 30.0
zcrit = 30.0
Tcrit = 10.0

# Loop over ascii files:

for fname in fnames:
    print(fname)
    
    # Read data
    [time,x,y,z,T] = np.transpose(np.genfromtxt(fname,usecols=(0,3,5,7,9)))
        
    # initialize mean value array mn= [x,y,z,T] 
    mn1=[[],[],[],[]]
    mn2=[[],[],[],[]]
    
    # -------------------------------------------------------------------------
    # Despiking routine: 
    
    j = 0
    for D in [x,y,z,T]:

        # Interpolate nan values 
        msk_nans = np.isnan(D)
        if len(msk_nans[msk_nans])>0:
            xdata = np.arange(len(D))
            D[msk_nans] = np.interp(xdata[msk_nans], xdata[~msk_nans],D[~msk_nans])
        
        # Calculate statistical parameters
        # index 1 represents main averaging (forward) until len(D) - size + 1
        # index 2 represents end averaging (backward) after len(D) - size + 1
        
        
        dat1 = np.array([ D[x:x+size] for x in range( len(D) - size + 1 ) ])
        dat2 = np.array([ D[x-size:x] for x in range( len(D) )[len(D)-size +1:]])
        
        # Calculate statistics on 'size' intervals
        
        mn1[j] = np.mean(dat1,axis = -1)
        mn2[j] = np.mean(dat2,axis = -1)
        
        std1 = np.std(dat1,axis = -1)
        std2 = np.std(dat2,axis = -1)
        
        sN1 = len(dat1[:,0])
        sN2 = len(dat2[:,0])

        cov1 = np.array([np.corrcoef(dat1[i,1:],dat1[i,:-1])[1,0] for i in range(sN1)])
        cov2 = np.array([np.corrcoef(dat2[i,1:],dat2[i,:-1])[1,0] for i in range(sN2)])        
      
        #calculate the forecasted value 
        prev1 = D[:len(D)-size]
        prev2 = D[len(D)-size:-1]
        mean1 = mn1[j][1:]
        mean2 = mn2[j]
        fcst1 = cov1[1:]*prev1 + (1.-cov1[1:])*mean1
        fcst2 = cov2*prev2 + (1.-cov2)*mean2
        obs1 = D[1:len(D)-size+1]
        obs2 = D[len(D)-size+1:]
        
        # Spike test
        msk_spikes1 = np.absolute(fcst1-obs1)/std1[1:] > C_spike[j]
        msk_spikes2 = np.absolute(fcst2-obs2)/std2 > C_spike[j] 
        
        # Calculate the number of new observed spikes
        N_spikes1 = len(msk_spikes1[msk_spikes1]) 
        N_spikes2 = len(msk_spikes2[msk_spikes2])         
       
        # Replace spikes by interpolation
        xdata1 = np.arange(len(D[1:len(D)-size+1]))
        xdata2 = np.arange(len(D[len(D)-size+1:]))
        if (N_spikes1>0):
            D[1:len(D)-size+1][msk_spikes1] = \
            np.interp(xdata1[msk_spikes1], xdata1[~msk_spikes1], D[1:len(D)-size+1][~msk_spikes1])    
        if (N_spikes2>0):
            D[len(D)-size+1:][msk_spikes2] = \
            np.interp(xdata2[msk_spikes2], xdata2[~msk_spikes2], D[len(D)-size+1:][~msk_spikes2])       
        j = j + 1
  
    #--------------------------------------------------------------------------
    #--------------------------------------------------------------------------
    # Build 10 minute bins 
    
    # Split data into 10 min sequencies: Good and Bad
    # u10min = TTMMDDhhm
    
    u10min = np.unique(np.floor(time[:]/1.0e+3)) 

    tim1 = np.array([ time[i] for i in range( len(D) - size + 1 ) ])
    tim2 = np.array([ time[i] for i in range( len(D) )[len(D)-size +1:] ])
    
    # Check correct sampling for 10 min, interval 
    for m in u10min:
        
        msk1 = np.floor(time[:]/1.0e+3) == m
        sec = np.unique(time[msk1])
            
        for ti in sec:
            msk2 = time[msk1] == ti
            if len(msk2[msk2]) == 20 and len(sec) == 600:
                datgood = True
            else: 
                datgood = False
                break
        
        msk3 = np.floor(tim1[:]/1.0e+3) == m
        msk3b = np.floor(tim2[:]/1.0e+3) == m 
        
        # Means: 
        xm1 = mn1[0][msk3]
        ym1 = mn1[1][msk3]
        zm1 = mn1[2][msk3]
        Tm1 = mn1[3][msk3]
        
        xm2 = mn2[0][msk3b]
        ym2 = mn2[1][msk3b]
        zm2 = mn2[2][msk3b]
        Tm2 = mn2[3][msk3b]      
        
        # Check 2s means are physically reasonable 
        for xi in np.append(xm1,xm2) :
            if xi > xcrit or datgood == False:
                datgood = False 
                break
            else: 
                datgood = True
        for yi in np.append(ym1,ym2):
            if yi > ycrit or datgood == False:
                datgood = False
                break
            else: 
                datgood = True 
        for zi in np.append(zm1,zm2):
            if zi > zcrit or datgood == False:
                datgood = False
                break
            else:
                datgood = True
        for Ti in np.append(Tm1,Tm2):
            if Ti > Tcrit or datgood == False:
                datgood = False
                break
            else:
                datgood = True
        
        # Save or update netCDF files for good and bad 10min intervals                     
        if datgood:
            
            good_file = good_file + 1
            
            t10 = time[msk1]
            x10 = x[msk1]
            y10 = y[msk1]
            z10 = z[msk1]
            T10 = T[msk1]
    
            # create/open nc-file
            
            if os.path.isfile(outdir1+'10min_AWS5_20hz_'+ str(m)+'.nc'):
                nc = Dataset(outdir1+'10min_AWS5_20hz_'+ str(m)+'.nc','r+',format='NETCDF4')
            else:
                nc = Dataset(outdir1+'10min_AWS5_20hz_'+str(m)+'.nc','w',format='NETCDF4')
                
            # create/open group & variables
            gkeys = nc.groups.keys()
            if ('10m') in gkeys:
                gg = nc.groups['10m']
                st = gg.variables['time']
                sX = gg.variables['X']
                sY = gg.variables['Y']
                sZ = gg.variables['Z']
                sT = gg.variables['T']

            else:
                gg = nc.createGroup(str('10m'))
                gg.createDimension('dim',12000)
                st = gg.createVariable('time', 'i8', ('dim'))
                sX = gg.createVariable('X', 'f4', ('dim'))
                sY = gg.createVariable('Y', 'f4', ('dim'))
                sZ = gg.createVariable('Z', 'f4', ('dim'))
                sT = gg.createVariable('T', 'f4', ('dim'))
                
            # save data into nc-file:
            st[:] =  t10.copy()
            sX[:] =  x10.copy()
            sY[:] =  y10.copy()
            sZ[:] =  z10.copy()
            sT[:] =  T10.copy()
            nc.close()
            
        else: 
            
            bad_file = bad_file + 1
            
            t10 = time[msk1]
            x10 = x[msk1]
            y10 = y[msk1]
            z10 = z[msk1]
            T10 = T[msk1]

            # create/open nc-file
            
            if os.path.isfile(outdir2+'10m_AWS5_20hz_'+ str(m) + '.nc'):
                nc = Dataset(outdir2+'10m_AWS5_20hz_'+ str(m) +'.nc','r+',format='NETCDF4')
            else:
                nc = Dataset(outdir2+'10m_AWS5_20hz_'+ str(m) +'.nc','w',format='NETCDF4')
                
            # create/open group & variables
            gkeys = nc.groups.keys()
            if ('10m_BAD') in gkeys:
                gg = nc.groups['10m_BAD']
                st = gg.variables['time']
                sX = gg.variables['X']
                sY = gg.variables['Y']
                sZ = gg.variables['Z']
                sT = gg.variables['T']

            else:
                gg = nc.createGroup('10m_BAD')
                gg.createDimension('dim',len(t10))
                st = gg.createVariable('time', 'i8', ('dim'))
                sX = gg.createVariable('X', 'f4', ('dim'))
                sY = gg.createVariable('Y', 'f4', ('dim'))
                sZ = gg.createVariable('Z', 'f4', ('dim'))
                sT = gg.createVariable('T', 'f4', ('dim'))
            #
            # save data into nc-file:
            st[:] =  t10.copy()
            sX[:] =  x10.copy()
            sY[:] =  y10.copy()
            sZ[:] =  z10.copy()
            sT[:] =  T10.copy()
            nc.close()

            
print( str(np.floor(good_file/(bad_file + good_file)*100)) +'% Success'  )
