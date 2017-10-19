import numpy as np
from matplotlib import pyplot as plt
from netCDF4 import Dataset
import glob

###########################################################################################

# data path
datdir = '

# Read data file names
fnames = np.sort(np.array(glob.glob(datdir+'AWS5_sonic_20hz_*.nc')))
print(len(fnames))

hgt = 10

# loop over data files
for fname in fnames:
  #
  nc = Dataset(fname,'r+',format='NETCDF4')
  #
  keys = nc.groups.keys() # loop over heights
  for k in keys:
    #
    # open group
    gg = nc.groups[k]
    #
    # read variables
    x = gg.variables['X'][:]
    y = gg.variables['Y'][:]
    z = gg.variables['Z'][:]
    #
    # COORDINATE TRANSFORMATION: 2D method
    #
    # yaw-correction: <v> = 0
    Xmn = x.mean()
    Ymn = y.mean()
    Zmn = z.mean()
    # components of the untilt tensor:
    A11 = Xmn/np.sqrt(Xmn**2.0+Ymn**2.0)
    A21 = Ymn/np.sqrt(Xmn**2.0+Ymn**2.0)
    A12 = -Ymn/np.sqrt(Xmn**2.0+Ymn**2.0)
    A22 = Xmn/np.sqrt(Xmn**2.0+Ymn**2.0)
    Ayaw = np.array([[A11,A12,0.0],[A21,A22,0.0],[0.0,0.0,1.0]])
    xx,yy,zz = np.dot( Ayaw.T, np.array([x, y, z]))
    #
    # pitch-correction: <w> = 0
    Xmn = np.mean(xx)
    Ymn = np.mean(yy)
    Zmn = np.mean(zz)
    # components of the untilt tensor:
    A11 = Xmn/np.sqrt(Xmn**2.0+Zmn**2.0)
    A31 = Zmn/np.sqrt(Xmn**2.0+Zmn**2.0)
    A13 = -Zmn/np.sqrt(Xmn**2.0+Zmn**2.0)
    A33 = Xmn/np.sqrt(Xmn**2.0+Zmn**2.0)
    Apitch = np.array([[A11,0.0,A13],[0.0,1.0,0.0],[A31,0.0,A33]])
    u,v,w = np.dot( Apitch.T, np.array([xx, yy, zz]))
    #
    print fname[-26:-3], k, u.mean(), v.mean(), w.mean()
    #
    # save results to nc-file:
    vnams = gg.variables.keys()
    if 'u' in vnams:
      su = gg.variables['u']
    else:
      su = gg.createVariable('u', 'f4', ('dim'))
    if 'v' in vnams:
      sv = gg.variables['v']
    else:
      sv = gg.createVariable('v', 'f4', ('dim'))
    if 'w' in vnams:
      sw = gg.variables['w']
    else:
      sw = gg.createVariable('w', 'f4', ('dim'))
    su[:] = u
    sv[:] = v
    sw[:] = w
    #
  nc.close()
