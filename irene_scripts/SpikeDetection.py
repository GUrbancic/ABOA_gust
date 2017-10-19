import numpy as np
from matplotlib import pyplot as plt

#----------------------------------------------------------------------------------------
# This script provides the spike detection method used in
# Suomi et al. (2017, QJRMS): http://onlinelibrary.wiley.com/doi/10.1002/qj.3059/pdf
#----------------------------------------------------------------------------------------

# READ DATA
datdir = '/Users/Gabin/Documents/FMI/Data_for_Gabin/AWS5_3D_sonic'
fname = '/10m-101229_0000.dat'

[time,x,y,z,T] = np.transpose(np.genfromtxt(datdir+fname,usecols=(0,3,5,7,9)))

#----------------------------------------------------------------------------------------
# SPIKE DETECTION

# Threshold for spike detection
C_spike = 4.5 # THIS MUST BE TESTED FOR EACH VARIABLE

# window length for spike detection
size = 50 # THIS CAN ALSO BE TUNED

# Loop over parameters
for D in [x]:#[x,y,z,T]:
  D = D[:24000]
  D_orig = D.copy()
  N = len(D)
  #
  # interpolate over possible nan values:
  msk_nans = np.isnan(D)

  if len(msk_nans[msk_nans])>0:
    xdata = np.arange(len(D))
    D[msk_nans] = np.interp(xdata[msk_nans], xdata[~msk_nans], D[~msk_nans])
  #        
  # Initialize spike mask
  msk_spikes = np.zeros(len(D),dtype=bool)
  #
  # initialize number of spikes
  N_spikes = 1
  #
  # count iterations
  count = 0
  while (N_spikes>0):
    #
    print( "***",count+1,"***")
    #
    # reshape data array to enable moving average calculation
    dat = np.array([ D[x:x+size] for x in range( len(D) - size + 1 ) ])
    #
    # ----------
    # calculate statistical parameters
    mn = np.mean(dat,axis=-1) # moving average
    std = np.std(dat,axis=-1) # moving std
    sN = len(dat[:,0]) # length of the dat array
    # moving covariance from each row in dat array:
    cov = np.array([np.corrcoef(dat[i,1:],dat[i,:-1])[1,0] for i in range(sN)])
    #
    # ----------
    # Calculate the forecasted value
    obs_old = D[(size-2):-1]
    fcst = cov*obs_old + (1.-cov)*mn
    #
    # ----------
    # Define the observed value
    obs_new = D[size-1:]
    #
    # ----------
    # Create a temporary array for the spike mask
    msk_spikes_new = (np.absolute(fcst-obs_new)/std > (C_spike + count*0.1)) # HERE is the increase (0.1) in the threshold
    #
    # Calculate the number of new observed spikes
    N_spikes = len(msk_spikes_new[msk_spikes_new])
    #
    # Include the new spikes in the spike array
    msk_spikes[size-1:] = msk_spikes[size-1:] | msk_spikes_new
    #
    # Replace possible new spike(s) by linear interpolation using the neighbouring points
    xdata = np.arange(len(D))
    if (N_spikes>0):
      D[msk_spikes] = np.interp(xdata[msk_spikes], xdata[~msk_spikes], D[~msk_spikes])
      
    count += 1
    #
  # Check detected spikes:
  DoPlot=True
  if DoPlot:
    xdata = np.arange(N-size+1)
    plt.figure(figsize=(20,6))
    plt.plot(xdata,np.absolute((fcst-obs_new)/std),'b.')
    plt.plot(plt.gca().get_xlim(),[C_spike,C_spike],'m-',lw=2)
    plt.plot(xdata,D_orig[size-1:],'k-o')
    plt.plot(xdata,D[size-1:],'r-s')
    plt.tight_layout()

plt.show()
