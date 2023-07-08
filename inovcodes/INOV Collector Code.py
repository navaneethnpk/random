##Date=Mar 29, 2022
import numpy as np
import matplotlib.pyplot as plt
import numpy.ma as ma
from astropy.stats import sigma_clip
import glob
import math
import datetime as dt
from scipy.stats import f
from astropy.table import Table
from astropy.io import ascii
def mjd_to_date(mjd):
    mjd = mjd + 2400000.5
    jd = mjd + 0.5
    F, I = math.modf(jd)
    I = int(I)
    A = math.trunc((I - 1867216.25)/36524.25)
    if I > 2299160:
        B = I + 1 + A - math.trunc(A / 4.)
    else:
        B = I  
    C = B + 1524
    D = math.trunc((C - 122.1) / 365.25)
    E = math.trunc(365.25 * D)
    G = math.trunc((C - E) / 30.6001)
    day = C - E + F - math.trunc(30.6001 * G)
    day = int(day)
    if G < 13.5:
        month = G - 1
    else:
        month = G - 13
    if month > 2.5:
        year = D - 4716
    else:
        year = D - 4715
    date = f"{day}-{month}-{year}"
    return date
def dur_to_time(val):
    mi,hr = math.modf(val)
    hr = int(hr)
    mi = mi*60
    mi = round(mi,2)
    time = f"{hr} Hrs, {mi} Mins"
    return time

cal_files = glob.glob('/home/ubuntu/INOV/calibrator/*r.dat', recursive = True)

eta_vals = []
for c in range(len(cal_files)):
    cal_obsid,cal_mjd,cal_mag,cal_magerr,cal_catflags,cal_ra,cal_dec = np.loadtxt(cal_files[c],unpack=True,usecols=[0,3,4,5,6,8,9],skiprows=54)
    #catflag
    cal_cat_ind = np.where(cal_catflags == 0)[0]
    cal_obsid,cal_mjd,cal_mag,cal_magerr,cal_catflags,cal_ra,cal_dec = cal_obsid[cal_cat_ind],cal_mjd[cal_cat_ind],cal_mag[cal_cat_ind],cal_magerr[cal_cat_ind],cal_catflags[cal_cat_ind],cal_ra[cal_cat_ind],cal_dec[cal_cat_ind]
    #sorted
    cal_mjd_ind = np.argsort(cal_mjd, axis = 0)
    cal_obsid,cal_mjd,cal_mag,cal_magerr,cal_catflags,cal_ra,cal_dec = cal_obsid[cal_mjd_ind],cal_mjd[cal_mjd_ind],cal_mag[cal_mjd_ind],cal_magerr[cal_mjd_ind],cal_catflags[cal_mjd_ind],cal_ra[cal_mjd_ind],cal_dec[cal_mjd_ind]
    #sigma clipped
    cal_clipped_mag = sigma_clip(cal_mag,sigma=3,maxiters=3,masked=True)
    cal_unmasked_id = ma.nonzero(cal_clipped_mag)
    cal_obsid,cal_mjd,cal_mag,cal_magerr,cal_catflags,cal_ra,cal_dec = np.array(cal_obsid)[cal_unmasked_id],np.array(cal_mjd)[cal_unmasked_id],np.array(cal_mag)[cal_unmasked_id],np.array(cal_magerr)[cal_unmasked_id],np.array(cal_catflags)[cal_unmasked_id],np.array(cal_ra)[cal_unmasked_id],np.array(cal_dec)[cal_unmasked_id],
    
    diff_arr = []
    #print(len(cal_mjd))
    #plt.errorbar(cal_mjd,cal_mag,yerr=cal_magerr,color='red', ls='none',linewidth=1, marker='.',mfc='white', capsize=2, capthick=1, markersize=5, ecolor='black')
    #plt.gca().invert_yaxis()
    #plt.title(f"calibrator{c} ; total datapoints={len(cal_mjd)}")        
    #plt.show()
    for i in range(len(cal_mjd)):
        if(i==0):
            diff_arr.append(float(0))
        else:
            cdiff = abs(cal_mjd[i]-cal_mjd[i-1])
            diff_arr.append(cdiff)
    gap_ind = np.where(np.array(diff_arr)>= (15/(60*24)))[0] #15 Minutes
    all_session_mjds,all_session_mags,all_session_errs = [],[],[]
    for j in range(len(gap_ind)):
        if(j == len(gap_ind)-1):
            continue
        session_mjd = cal_mjd[gap_ind[j]:gap_ind[j+1]]
        session_mjd_int = np.array(session_mjd,dtype=int)
        session_mjd_UT = session_mjd  - np.min(session_mjd_int)
        session_mag = cal_mag[gap_ind[j]:gap_ind[j+1]]
        session_err = cal_magerr[gap_ind[j]:gap_ind[j+1]]
        session_dur = (np.max(session_mjd)-np.min(session_mjd))*24
        session_counts = gap_ind[j+1]-gap_ind[j]
        if(session_counts>=10 and session_dur>=0): 
            date = mjd_to_date(np.min(session_mjd_int))
            time = dur_to_time(session_dur)
            #print(f"Calibrator{c} - session date={date}; duration={time}; datapoints={session_counts}")
            
            var = np.sum((session_mag - np.mean(session_mag))**2)/(len(session_mag)-1)
            sig = np.mean(session_err**2)
            eta = var/sig
            #print(eta)
            eta_vals.append(eta)
#print(eta_vals)
eeta = np.mean(eta_vals)
print("Avg eta = ",eeta)           