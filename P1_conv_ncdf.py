##---------------------------------------------------------------------------
# P1_conv_ncdf.py 
# G.Urbancic 
#
# Takes NetCDF files produced by R1_conv_ncdf: 
#    (1) Determines- "great days", (largest % NaN/ min) < 1%
#                  - "good days" , 1% < (largest % NaN/ min) < 5%
#                  - "prob days",   5% < (largest % NaN/ min) 
#    (2) Plots figures of non-perfect days 
#    (3) Salvage or Toss?  
##---------------------------------------------------------------------------

import numpy as np
import matplotlib.pyplot as plt
import os.path
from netCDF4 import Dataset
import glob

# path to ascii raw data 
datdir = '/Users/gabin/Documents/ABOA/DATA/CSAT_2m/'

# NetCDF file names
fnames = np.sort(np.array(glob.glob(datdir+'*.nc')))

print("number of NetCDF files:" + str(len(fnames)))

count = [0,0,0]                                         # [great,good,prob]
prob_list = []

for fname in fnames:
    
    Ddate = fname.replace(datdir,"")
    Ddate = Ddate.replace("_2m","")

    # Read in NetCDF files for 2m - NaN_min group 
    
    nc2 = Dataset(fname,'r',format='NETCDF4')
    gnan2 = nc2.groups['NaN_min']
    t1_2 = gnan2.variables['t']
    Xn_2 = gnan2.variables['X']
    Yn_2 = gnan2.variables['Y']
    Zn_2 = gnan2.variables['Z']
    Tn_2 = gnan2.variables['T'] 
    
    max_X = max(Xn_2)
    max_Y = max(Yn_2)
    max_Z = max(Zn_2)
    max_T = max(Tn_2)
    
    maxtot2 = max([max_X,max_Y,max_Z,max_T])
    
    # Read in NetCDF files for 10m - NaN_min group 
    
    # 10m directory 
    fname2 = fname.replace("_2m","_10m")
    
    nc10 = Dataset(fname2,'r',format='NETCDF4')
    gnan10 = nc10.groups['NaN_min']
    t1_10 = gnan10.variables['t']
    Xn_10 = gnan10.variables['X']
    Yn_10 = gnan10.variables['Y']
    Zn_10 = gnan10.variables['Z']
    Tn_10 = gnan10.variables['T'] 
    
    max_X10 = max(Xn_10)
    max_Y10 = max(Yn_10)
    max_Z10 = max(Zn_10)
    max_T10 = max(Tn_10)
    
    maxtot10 = max([max_X10,max_Y10,max_Z10,max_T10])
    
    if (maxtot2 >= 5) or (maxtot10 >= 5):
        print(Ddate + " : prob")
        count[2] = count[2] + 1
        prob_list = prob_list + [fname]
        
    elif (maxtot2 >= 1) or (maxtot10 >= 1):
        print(Ddate + " : good")
        count[1] = count[1] + 1
    else: 
        print(Ddate + " : great")
        count[0] = count[0] + 1
        
print("Great: " + str(count[0]) + " Good: " + str(count[1]) + " Prob: " + str(count[2]))     
    

    
    
    