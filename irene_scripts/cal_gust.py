import numpy as np

#------------------------------------------------------------------------------------
# moving average
def movingaverage(interval, window_size):
    window= np.ones(int(window_size))/float(window_size)
    return np.convolve(interval, window, 'same')

#------------------------------------------------------------------------------------
# wind gust from the horizontal wind speed, tg = gust duration in [s]
def cal_gust(U,tg):
    u_move = movingaverage(U, tg*20.0)[int(tg*10):-int(tg*10)]
    Umax = np.ma.max(u_move)
    Umin = np.ma.min(u_move)
    return Umax, Umin

#------------------------------------------------------------------------------------
# wind gust and wind minimum from the horizontal wind speed
# + wind components during the gust and the minimum,
# tg = gust duration in [s] (= moving average window width)
def cal_gust_components(u,v,w,tg):
    #
    # horizontal wind speed array
    U = np.sqrt(u**2.0+v**2.0)
    #
    # moving average of the horizontal wind speed, with a window length = gust length tg
    U_move = movingaverage(U, tg*20.0)[int(tg*10):-int(tg*10)]
    #
    # maximum and minimum gusts
    Umax = np.ma.max(U_move)
    #
    # minimum gusts
    Umin = np.ma.min(U_move)
    #
    # the location of the first occurrence of the gusts
    loc_max = np.where(U_move==Umax)[0][0]
    loc_min = np.where(U_move==Umin)[0][0]
    #
    # moving averages of the wind components:
    u_move = movingaverage(u, tg*20.0)[int(tg*10):-int(tg*10)]
    v_move = movingaverage(v, tg*20.0)[int(tg*10):-int(tg*10)]
    w_move = movingaverage(w, tg*20.0)[int(tg*10):-int(tg*10)]
    #
    # return the maximum and the minimum gusts + the wind components 
    return Umax, u_move[loc_max], v_move[loc_max], w_move[loc_max], Umin, u_move[loc_min], v_move[loc_min], w_move[loc_min]

