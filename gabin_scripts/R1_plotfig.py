from datetime import datetime 
import matplotlib.pylab as plt
import glob
import matplotlib.dates as mdates
import numpy as np

# R1 Plot Figure 
# Script to build figures from data saved in R1_*.py routines.
# See R1_dataqual.py for information on data saves. 

# READ DATA
datdir = '/Users/gabin/Documents/FMI/ABOA_data/ORIG_10m/'

# ascii file names
fil_cnt = np.sort(np.array(glob.glob(datdir+'fil_cnt*')))
spk_cnt = np.sort(np.array(glob.glob(datdir+'spik_cnt*')))
crit_data = np.sort(np.array(glob.glob(datdir+'crit_data*')))
spk_cnt_g = [datdir + 'low_spik_cnt_1101' +i+ '.0.txt' for i in ['09','10','11','12',\
                                                        '13','14','15','16','17','18']]



##------------------------------------------------------------------------------
##------------------------------------------------------------------------------
# Plot # of valid 10 minute sequences per hour for the campaign.

sampsize = 0

if sampsize:
    samp= [[],[],[]]
    for fil in fil_cnt:
        [s0,s1,s2] = np.transpose(np.genfromtxt(fil,usecols=(0,1,2)))
        samp[0] = np.append(samp[0], s0)
        samp[1] = np.append(samp[1], s1)
        samp[2] = np.append(samp[2], s2)
    
    # samp[0] is the time in format YYMMDDhh
    # convert to datetime format: 
        
    dt = [0]*len(samp[0])
    tot = [0]*len(samp[0])
    idx = 0

    for i in samp[0]:
        year = int(2000 + np.floor(i/1.0e+6))
        month = int(np.floor(i/1.0e+4) - np.floor(i/1.0e+6)*100)
        day =  int(np.floor(i/1.0e+2) - np.floor(i/1.0e+4)*100)
        hour = int(i - np.floor(i/1.0e+2)*100)
        dt[idx] = datetime(year,month,day,hour)
        tot[idx] = samp[2][idx]+samp[1][idx]
        idx = idx + 1

    
    plt.figure(1,figsize=(20,6))

    ax1 = plt.subplot(211)
    plt.title('10m - 20Hz valid: 10 minutes sequences per hour',fontweight='bold')
    plt.plot(dt[0:int(np.floor(len(dt)/2.0))],tot[0:int(np.floor(len(dt)/2.0))],'bo' \
                ,dt[0:int(np.floor(len(dt)/2.0))],samp[1][0:int(np.floor(len(dt)/2.0))],'go')
    plt.grid(True)
    ax1.xaxis.set_major_locator(mdates.HourLocator(interval=24))
    ax1.xaxis.set_major_formatter(mdates.DateFormatter('%m/%d-%H'))
    ax1.xaxis.grid(True, which="minor")
    ax1.axis([min(dt),dt[int(np.floor(len(dt)/2.0))],-0.1,6.1])

    ax2 = plt.subplot(212)
    plt.plot(dt[int(np.floor(len(dt)/2.0)):],tot[int(np.floor(len(dt)/2.0)):],'bo',
            dt[int(np.floor(len(dt)/2.0)):],samp[1][int(np.floor(len(dt)/2.0)):],'go')
    plt.grid(True)
    ax2.xaxis.set_major_locator(mdates.HourLocator(interval=24))
    ax2.xaxis.set_major_formatter(mdates.DateFormatter('%m/%d-%H'))
    ax2.xaxis.grid(True, which="minor")
    ax2.axis([dt[int(np.floor(len(dt)/2.0))],max(dt),-0.1,6.1])

    plt.savefig('/Users/gabin/Documents/FMI/20Hzvalid_10m')






##---------------------------------------------------------------------------------
##---------------------------------------------------------------------------------
# Plot distribution of valid 10 minute sequence for time of day.
  
diurl_distro = 0

if diurl_distro:
    samp= [[],[],[]]
    for fil in fil_cnt:
        [s0,s1,s2] = np.transpose(np.genfromtxt(fil,usecols=(0,1,2)))
        samp[0] = np.append(samp[0], s0)
        samp[1] = np.append(samp[1], s1)
        samp[2] = np.append(samp[2], s2) 
    hr =[]
    for ti in samp[0]:
        hr = np.append(hr, np.array([ti-np.floor(ti/100)*100]))
        
    diurl_distro = np.zeros(24)
    
    m = 0    
    for h in hr:  
        diurl_distro[int(h)]=diurl_distro[int(h)] + samp[1][m]
        m = m + 1
    plt.figure(2)
    plt.plot(np.arange(24),diurl_distro,'bo')
    plt.title('Valid 10-minute sample per time of day')
    plt.xticks(np.arange(24))
    plt.grid('on')
    plt.xlabel('time of day')
    plt.ylabel('# of valid 10 minute samples')
    plt.savefig('/Users/gabin/Documents/FMI/20Hzvalid_10m_diurnal')





##------------------------------------------------------------------------------
##------------------------------------------------------------------------------
# Plot spikes per hour for duration of the campaign

spike_cnt = 0

if spike_cnt:
    spk= [[],[],[],[],[],[],[],[],[],[],[],[],[]]
    for fil in spk_cnt:
        [s0,s1,s2,s3,s4,s5,s6,s7,s8,s9,s10,s11,s12] = np.transpose(np.genfromtxt(fil, \
        usecols=(0,1,2,3,4,5,6,7,8,9,10,11,12)))
        spk[0] = np.append(spk[0], s0)
        spk[1] = np.append(spk[1], s1)
        spk[2] = np.append(spk[2], s2)
        spk[3] = np.append(spk[3], s3)
        spk[4] = np.append(spk[4], s4)
        spk[5] = np.append(spk[5], s5)
        spk[6] = np.append(spk[6], s6)
        spk[7] = np.append(spk[7], s7)
        spk[8] = np.append(spk[8], s8)
        spk[9] = np.append(spk[9], s9)
        spk[10] = np.append(spk[10], s10)
        spk[11] = np.append(spk[11], s11)
        spk[12] = np.append(spk[12], s12)
    
    # samp[0] is the time in format YYMMDDhh
    # convert to datetime format: 
        
    dt = [0]*len(spk[0])
  
    idx = 0

    for i in spk[0]:
        year = int(2000 + np.floor(i/1.0e+6))
        month = int(np.floor(i/1.0e+4) - np.floor(i/1.0e+6)*100)
        day =  int(np.floor(i/1.0e+2) - np.floor(i/1.0e+4)*100)
        hour = int(i - np.floor(i/1.0e+2)*100)
        dt[idx] = datetime(year,month,day,hour)
        idx = idx + 1
        
 #--------------------------------    
    plt.figure(3,figsize=(20,6))

    ax1 = plt.subplot(211)
    plt.title('10m - spk/h for x velocity',fontweight='bold')
    plt.plot(dt[0:int(np.floor(len(dt)/2.0))],spk[1][0:int(np.floor(len(dt)/2.0))],'bo' \
                ,dt[0:int(np.floor(len(dt)/2.0))],spk[2][0:int(np.floor(len(dt)/2.0))],'go' \
                ,dt[0:int(np.floor(len(dt)/2.0))],spk[3][0:int(np.floor(len(dt)/2.0))],'ro')
    plt.grid(True)
    ax1.xaxis.set_major_locator(mdates.HourLocator(interval=24))
    ax1.xaxis.set_major_formatter(mdates.DateFormatter('%m/%d-%H'))
    ax1.xaxis.grid(True, which="minor")
    ax1.axis([min(dt),dt[int(np.floor(len(dt)/2.0))],-0.1,0.4])

    ax2 = plt.subplot(212)
    plt.plot(dt[int(np.floor(len(dt)/2.0)):],spk[1][int(np.floor(len(dt)/2.0)):],'bo' \
            ,dt[int(np.floor(len(dt)/2.0)):],spk[2][int(np.floor(len(dt)/2.0)):],'go' \
            ,dt[int(np.floor(len(dt)/2.0)):],spk[3][int(np.floor(len(dt)/2.0)):],'ro')
    plt.grid(True)
    ax2.xaxis.set_major_locator(mdates.HourLocator(interval=24))
    ax2.xaxis.set_major_formatter(mdates.DateFormatter('%m/%d-%H'))
    ax2.xaxis.grid(True, which="minor")
    ax2.axis([dt[int(np.floor(len(dt)/2.0))],max(dt),-0.1,0.4])

    plt.savefig('/Users/gabin/Documents/FMI/spk_hr_x_10m')
    
#---------------------------------
    plt.figure(4,figsize=(20,6))

    ax1 = plt.subplot(211)
    plt.title('10m - spk/h for y velocity',fontweight='bold')
    plt.plot(dt[0:int(np.floor(len(dt)/2.0))],spk[4][0:int(np.floor(len(dt)/2.0))],'bo' \
                ,dt[0:int(np.floor(len(dt)/2.0))],spk[5][0:int(np.floor(len(dt)/2.0))],'go' \
                ,dt[0:int(np.floor(len(dt)/2.0))],spk[6][0:int(np.floor(len(dt)/2.0))],'ro')
    plt.grid(True)
    ax1.xaxis.set_major_locator(mdates.HourLocator(interval=24))
    ax1.xaxis.set_major_formatter(mdates.DateFormatter('%m/%d-%H'))
    ax1.xaxis.grid(True, which="minor")
    ax1.axis([min(dt),dt[int(np.floor(len(dt)/2.0))],-0.1,0.4])

    ax2 = plt.subplot(212)
    plt.plot(dt[int(np.floor(len(dt)/2.0)):],spk[4][int(np.floor(len(dt)/2.0)):],'bo' \
            ,dt[int(np.floor(len(dt)/2.0)):],spk[5][int(np.floor(len(dt)/2.0)):],'go' \
            ,dt[int(np.floor(len(dt)/2.0)):],spk[6][int(np.floor(len(dt)/2.0)):],'ro')
    plt.grid(True)
    ax2.xaxis.set_major_locator(mdates.HourLocator(interval=24))
    ax2.xaxis.set_major_formatter(mdates.DateFormatter('%m/%d-%H'))
    ax2.xaxis.grid(True, which="minor")
    ax2.axis([dt[int(np.floor(len(dt)/2.0))],max(dt),-0.1,0.4])

    plt.savefig('/Users/gabin/Documents/FMI/spk_hr_y_10m')

#---------------------------------
    plt.figure(5,figsize=(20,6))

    ax1 = plt.subplot(211)
    plt.title('10m - spk/h for z velocity',fontweight='bold')
    plt.plot(dt[0:int(np.floor(len(dt)/2.0))],spk[7][0:int(np.floor(len(dt)/2.0))],'bo' \
                ,dt[0:int(np.floor(len(dt)/2.0))],spk[8][0:int(np.floor(len(dt)/2.0))],'go' \
                ,dt[0:int(np.floor(len(dt)/2.0))],spk[9][0:int(np.floor(len(dt)/2.0))],'ro')
    plt.grid(True)
    ax1.xaxis.set_major_locator(mdates.HourLocator(interval=24))
    ax1.xaxis.set_major_formatter(mdates.DateFormatter('%m/%d-%H'))
    ax1.xaxis.grid(True, which="minor")
    ax1.axis([min(dt),dt[int(np.floor(len(dt)/2.0))],-0.1,0.4])

    ax2 = plt.subplot(212)
    plt.plot(dt[int(np.floor(len(dt)/2.0)):],spk[7][int(np.floor(len(dt)/2.0)):],'bo' \
            ,dt[int(np.floor(len(dt)/2.0)):],spk[8][int(np.floor(len(dt)/2.0)):],'go' \
            ,dt[int(np.floor(len(dt)/2.0)):],spk[9][int(np.floor(len(dt)/2.0)):],'ro')
    plt.grid(True)
    ax2.xaxis.set_major_locator(mdates.HourLocator(interval=24))
    ax2.xaxis.set_major_formatter(mdates.DateFormatter('%m/%d-%H'))
    ax2.xaxis.grid(True, which="minor")
    ax2.axis([dt[int(np.floor(len(dt)/2.0))],max(dt),-0.1,0.4])

    plt.savefig('/Users/gabin/Documents/FMI/spk_hr_z_10m')

#---------------------------------
    plt.figure(6,figsize=(20,6))

    ax1 = plt.subplot(211)
    plt.title('10m - spk/h for Temperature',fontweight='bold')
    plt.plot(dt[0:int(np.floor(len(dt)/2.0))],spk[10][0:int(np.floor(len(dt)/2.0))],'bo' \
                ,dt[0:int(np.floor(len(dt)/2.0))],spk[11][0:int(np.floor(len(dt)/2.0))],'go' \
                ,dt[0:int(np.floor(len(dt)/2.0))],spk[12][0:int(np.floor(len(dt)/2.0))],'ro')
    plt.grid(True)
    ax1.xaxis.set_major_locator(mdates.HourLocator(interval=24))
    ax1.xaxis.set_major_formatter(mdates.DateFormatter('%m/%d-%H'))
    ax1.xaxis.grid(True, which="minor")
    ax1.axis([min(dt),dt[int(np.floor(len(dt)/2.0))],-0.1,0.4])

    ax2 = plt.subplot(212)
    plt.plot(dt[int(np.floor(len(dt)/2.0)):],spk[10][int(np.floor(len(dt)/2.0)):],'bo' \
            ,dt[int(np.floor(len(dt)/2.0)):],spk[11][int(np.floor(len(dt)/2.0)):],'go' \
            ,dt[int(np.floor(len(dt)/2.0)):],spk[12][int(np.floor(len(dt)/2.0)):],'ro')
    plt.grid(True)
    ax2.xaxis.set_major_locator(mdates.HourLocator(interval=24))
    ax2.xaxis.set_major_formatter(mdates.DateFormatter('%m/%d-%H'))
    ax2.xaxis.grid(True, which="minor")
    ax2.axis([dt[int(np.floor(len(dt)/2.0))],max(dt),-0.1,0.4])

    plt.savefig('/Users/gabin/Documents/FMI/spk_hr_T_10m')




##------------------------------------------------------------------------------
##------------------------------------------------------------------------------
# Plot # of valid 10 minute sequences per hour for the campaign.

crt_data = 0

if crt_data:
    samp2= [[],[],[]]
    for fil in crit_data:
        [s0,s1,s2] = np.transpose(np.genfromtxt(fil,usecols=(0,1,2)))
        samp2[0] = np.append(samp2[0], s0)
        samp2[1] = np.append(samp2[1], s1)
        samp2[2] = np.append(samp2[2], s2)
    
    # samp[0] is the time in format YYMMDDhh
    # convert to datetime format: 
        
    dt2 = [0]*len(samp2[0])
    tot2 = [0]*len(samp2[0])
    idx = 0

    for i in samp2[0]:
        year = int(2000 + np.floor(i/1.0e+6))
        month = int(np.floor(i/1.0e+4) - np.floor(i/1.0e+6)*100)
        day =  int(np.floor(i/1.0e+2) - np.floor(i/1.0e+4)*100)
        hour = int(i - np.floor(i/1.0e+2)*100)
        dt2[idx] = datetime(year,month,day,hour)
        tot2[idx] = samp2[2][idx]+samp2[1][idx]
        idx = idx + 1

    
    plt.figure(7,figsize=(20,6))

    ax1 = plt.subplot(211)
    plt.title('10m - reasonable: 10 minutes sequences per hour',fontweight='bold')
    plt.plot(dt2[0:int(np.floor(len(dt)/2.0))],tot2[0:int(np.floor(len(dt)/2.0))],'bo' \
                ,dt2[0:int(np.floor(len(dt)/2.0))],samp2[1][0:int(np.floor(len(dt)/2.0))],'go')
    plt.grid(True)
    ax1.xaxis.set_major_locator(mdates.HourLocator(interval=24))
    ax1.xaxis.set_major_formatter(mdates.DateFormatter('%m/%d-%H'))
    ax1.xaxis.grid(True, which="minor")
    ax1.axis([min(dt),dt[int(np.floor(len(dt)/2.0))],-0.1,6.1])

    ax2 = plt.subplot(212)
    plt.plot(dt2[int(np.floor(len(dt)/2.0)):],tot2[int(np.floor(len(dt)/2.0)):],'bo',
            dt2[int(np.floor(len(dt)/2.0)):],samp2[1][int(np.floor(len(dt)/2.0)):],'go')
    plt.grid(True)
    ax2.xaxis.set_major_locator(mdates.HourLocator(interval=24))
    ax2.xaxis.set_major_formatter(mdates.DateFormatter('%m/%d-%H'))
    ax2.xaxis.grid(True, which="minor")
    ax2.axis([dt[int(np.floor(len(dt)/2.0))],max(dt),-0.1,6.1])

    plt.savefig('/Users/gabin/Documents/FMI/20Hzreason_10mB')
    
    
    

##---------------------------------------------------------------------------------
##---------------------------------------------------------------------------------
# Plot distribution of valid spikes per time of day.
  
diurl_distro_spk = 1

if diurl_distro_spk:
    samp= [[],[],[],[],[],[],[],[],[],[],[],[],[]]
    for fil in spk_cnt_g:
        [t0,x1,x2,x3,y1,y2,y3,z1,z2,z3,T1,T2,T3] = \
        np.transpose(np.genfromtxt(fil,usecols=(0,1,2,3,4,5,6,7,8,9,10,11,12)))
        samp[0] = np.append(samp[0], t0)
        samp[1] = np.append(samp[1], x1)
        samp[2] = np.append(samp[2], x2)
        samp[3] = np.append(samp[3], x3)
        samp[4] = np.append(samp[4], y1)
        samp[5] = np.append(samp[5], y2)
        samp[6] = np.append(samp[6], y3)        
        samp[7] = np.append(samp[7], z1)
        samp[8] = np.append(samp[8], z2)
        samp[9] = np.append(samp[9], z3)
        samp[10] = np.append(samp[10], T1)    
        samp[11] = np.append(samp[11], T2)
        samp[12] = np.append(samp[12], T3)
        
    hr = []    
    for ti in samp[0]:
        hr = np.append(hr, np.array([ti-np.floor(ti/100)*100]))
        
    diurl_distro = [np.zeros(24)]*12  
    m = 0    
    for h in hr:  
        diurl_distro[0][int(h)]=diurl_distro[0][int(h)] + samp[1][m]
        diurl_distro[1][int(h)]=diurl_distro[1][int(h)] + samp[2][m]
        diurl_distro[2][int(h)]=diurl_distro[2][int(h)] + samp[3][m]
        diurl_distro[3][int(h)]=diurl_distro[3][int(h)] + samp[4][m]
        diurl_distro[4][int(h)]=diurl_distro[4][int(h)] + samp[5][m]
        diurl_distro[5][int(h)]=diurl_distro[5][int(h)] + samp[6][m]
        diurl_distro[6][int(h)]=diurl_distro[6][int(h)] + samp[7][m]
        diurl_distro[7][int(h)]=diurl_distro[7][int(h)] + samp[8][m]
        diurl_distro[8][int(h)]=diurl_distro[8][int(h)] + samp[9][m]
        diurl_distro[9][int(h)]=diurl_distro[9][int(h)] + samp[10][m]
        diurl_distro[10][int(h)]=diurl_distro[10][int(h)] + samp[11][m]
        diurl_distro[11][int(h)]=diurl_distro[11][int(h)] + samp[12][m]
        m = m + 1
        
    diurl_distro[0]= [i/(len(spk_cnt_g)) for i in diurl_distro[0]] 
    diurl_distro[1]= [i/(len(spk_cnt_g)) for i in diurl_distro[1]]
    diurl_distro[2]= [i/(len(spk_cnt_g)) for i in diurl_distro[2]]
    diurl_distro[3]= [i/(len(spk_cnt_g)) for i in diurl_distro[3]] 
    diurl_distro[4]= [i/(len(spk_cnt_g)) for i in diurl_distro[4]]
    diurl_distro[5]= [i/(len(spk_cnt_g)) for i in diurl_distro[5]]
    diurl_distro[6]= [i/(len(spk_cnt_g)) for i in diurl_distro[6]] 
    diurl_distro[7]= [i/(len(spk_cnt_g)) for i in diurl_distro[7]]
    diurl_distro[8]= [i/(len(spk_cnt_g)) for i in diurl_distro[8]]
    diurl_distro[9]= [i/(len(spk_cnt_g)) for i in diurl_distro[9]] 
    diurl_distro[10]= [i/(len(spk_cnt_g)) for i in diurl_distro[10]]
    diurl_distro[11]= [i/(len(spk_cnt_g)) for i in diurl_distro[11]]
    
    hr = np.arange(24)
    
    plt.figure(8,figsize=(10,15))

    ax1 = plt.subplot(411)
    plt.title('10m - Spike Daily Distribution',fontweight='bold')
    plt.plot(hr,diurl_distro[0],'bo' \
                ,hr,diurl_distro[1],'go'\
                ,hr,diurl_distro[2],'ro')
    plt.grid(True)
    #ax1.axis([min(dt),dt[int(np.floor(len(dt)/2.0))],-0.1,6.1])

    ax2 = plt.subplot(412)
    plt.plot(hr,diurl_distro[3],'bo' \
                ,hr,diurl_distro[4],'go'\
                ,hr,diurl_distro[5],'ro')
    plt.grid(True)
    #ax1.axis([min(dt),dt[int(np.floor(len(dt)/2.0))],-0.1,6.1])

    ax3 = plt.subplot(413)
    plt.plot(hr,diurl_distro[6],'bo' \
                ,hr,diurl_distro[7],'go'\
                ,hr,diurl_distro[8],'ro')
    plt.grid(True)
    #ax1.axis([min(dt),dt[int(np.floor(len(dt)/2.0))],-0.1,6.1])
    
     ax4 = plt.subplot(414)
    plt.plot(hr,diurl_distro[9],'bo' \
                ,hr,diurl_distro[10],'go'\
                ,hr,diurl_distro[11],'ro')
    plt.grid(True)
    #ax1.axis([min(dt),dt[int(np.floor(len(dt)/2.0))],-0.1,6.1])
    
    
    
