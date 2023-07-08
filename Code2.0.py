##Date=Mar 31, 2022
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
import wget
from tqdm import tqdm
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
qsr_RA, qsr_DEC = np.loadtxt("/home/ubuntu/INOV/mywork14/qsrdata.dat",unpack=True,usecols=[1,2])
for q in tqdm(range(len(qsr_RA))):
    #qlink = f"https://irsa.ipac.caltech.edu/cgi-bin/ZTF/nph_light_curves?POS=CIRCLE+{qsr_RA[q]}+{qsr_DEC[q]}+0.000416667&BANDNAME=r&BAD_CATFLAGS_MASK=32768&FORMAT=ipac_table"
    #qfilename = f"{round(qsr_RA[q],5)}{round(qsr_DEC[q],5)}_qsr.dat"
    #wget.download(qlink, qfilename)
    qsr_files = glob.glob('/home/ubuntu/INOV/mywork18/*r.dat', recursive = True) #loading ZTF data of all quasars
    qsrdata = []
    for doc1 in qsr_files:
        file1 = open(doc1, 'r')
        lines1 = file1.readlines()
        #print(len(lines1))
        if(len(lines1)<=54):
            continue
        else:
            qsrdata.append(doc1) #Quasars having data is saved in this array
col1,col2,col3,col4,col5,col6 = [],[],[],[],[],[]
variability = []
for m in range(len(qsrdata)):
    qsr_obsid,qsr_mjd,qsr_mag,qsr_magerr,qsr_catflags,qsr_ra,qsr_dec = np.loadtxt(qsrdata[m],unpack=True,usecols=[0,3,4,5,6,8,9],skiprows=54)
    qsr_cat_ind = np.where(qsr_catflags == 0)[0] #catflag criteria
    qsr_obsid,qsr_mjd,qsr_mag,qsr_magerr,qsr_catflags,qsr_ra,qsr_dec = qsr_obsid[qsr_cat_ind],qsr_mjd[qsr_cat_ind],qsr_mag[qsr_cat_ind],qsr_magerr[qsr_cat_ind],qsr_catflags[qsr_cat_ind],qsr_ra[qsr_cat_ind],qsr_dec[qsr_cat_ind]
    qsr_mjd_ind = np.argsort(qsr_mjd, axis = 0) #mjd sorting
    qsr_clipped_mag = sigma_clip(qsr_mag,sigma=3,maxiters=3,masked=True) #sigma clipping
    qsr_unmasked_id = ma.nonzero(qsr_clipped_mag)
    qsr_obsid,qsr_mjd,qsr_mag,qsr_magerr,qsr_catflags,qsr_ra,qsr_dec = np.array(qsr_obsid)[qsr_unmasked_id],np.array(qsr_mjd)[qsr_unmasked_id],np.array(qsr_mag)[qsr_unmasked_id],np.array(qsr_magerr)[qsr_unmasked_id],np.array(qsr_catflags)[qsr_unmasked_id],np.array(qsr_ra)[qsr_unmasked_id],np.array(qsr_dec)[qsr_unmasked_id]
    diff_arr = []
    for i in range(len(qsr_mjd)):
        if(i==0):
            diff_arr.append(float(0))
        else:
            cdiff = abs(qsr_mjd[i]-qsr_mjd[i-1])
            diff_arr.append(cdiff)
    gap_ind = np.where(np.array(diff_arr)>= (15/(60*24)))[0] #15 Minutes
    all_session_mjds,all_session_mags,all_session_errs = [],[],[]
    for j in range(len(gap_ind)):
        if(j == len(gap_ind)-1):
            continue
        session_mjd = qsr_mjd[gap_ind[j]:gap_ind[j+1]]
        session_mjd_int = np.array(session_mjd,dtype=int)
        session_mjd_UT = session_mjd  - np.min(session_mjd_int)
        session_mag = qsr_mag[gap_ind[j]:gap_ind[j+1]]
        session_err = qsr_magerr[gap_ind[j]:gap_ind[j+1]]
        session_dur = (np.max(session_mjd)-np.min(session_mjd))*24
        session_counts = gap_ind[j+1]-gap_ind[j]
        if(session_counts>=10 and session_dur>=0):
            date = mjd_to_date(np.min(session_mjd_int))
            time = dur_to_time(session_dur)
            #print(f"Source{m} - session date={date}; duration={time}; datapoints={session_counts}")
            var = np.sum((session_mag - np.mean(session_mag))**2)/(len(session_mag)-1)
            sig = np.mean(session_err**2)
            eta = 1.0
            f_eta = var/(eta*sig)
            psi = (np.max(session_mag) - np.min(session_mag) )**2 - 2*( eta**2 * np.mean(session_err**2) ) 
            #critical values
            F_c_95 = (f.ppf(0.95, len(session_mag)-1, len(session_mag)-1))
            F_c_99 = (f.ppf(0.999, len(session_mag)-1, len(session_mag)-1))
            if f_eta>=F_c_99:
                variability.append("Variable")
            elif f_eta > F_c_95 and f_eta <=F_c_99:
                variability.append("Probable Variable")
            else:
                variability.append("Non Variable")
            col1.append(f"QSR_{m}")
            col2.append(date)
            col3.append(time)
            col4.append(len(session_mag))
            col5.append(round(f_eta,6))
            col6.append(round(psi,6))
result = Table()
result['Quasar Name'] = col1
result['Session Date'] = col2
result['Session Duration'] = col3
result['Session Datapoints'] = col4
result['f_eta'] = col5
result['Psi'] = col6
result['Variability Status'] = variability
print(result)
#ascii.write(result, 'result.dat', overwrite=True)