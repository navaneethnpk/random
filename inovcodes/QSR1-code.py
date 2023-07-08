##Date=Mar 17, 2022
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

qsr_obsid,qsr_mjd,qsr_mag,qsr_magerr,qsr_catflags,qsr_ra,qsr_dec = np.loadtxt("/home/ubuntu/INOV/mywork4/110.7486667-75263278.dat",unpack=True,usecols=[0,3,4,5,6,8,9],skiprows=54)
qsr_cat_ind = np.where(qsr_catflags == 0)[0]
qsr_obsid,qsr_mjd,qsr_mag,qsr_magerr,qsr_catflags,qsr_ra,qsr_dec = qsr_obsid[qsr_cat_ind],qsr_mjd[qsr_cat_ind],qsr_mag[qsr_cat_ind],qsr_magerr[qsr_cat_ind],qsr_catflags[qsr_cat_ind],qsr_ra[qsr_cat_ind],qsr_dec[qsr_cat_ind]
qsr_mjd_ind = np.argsort(qsr_mjd, axis = 0)
qsr_obsid,qsr_mjd,qsr_mag,qsr_magerr,qsr_catflags,qsr_ra,qsr_dec = qsr_obsid[qsr_mjd_ind],qsr_mjd[qsr_mjd_ind],qsr_mag[qsr_mjd_ind],qsr_magerr[qsr_mjd_ind],qsr_catflags[qsr_mjd_ind],qsr_ra[qsr_mjd_ind],qsr_dec[qsr_mjd_ind]
diff_arr = []

for i in range(len(qsr_mjd)):
    if(i==0):
        diff_arr.append(float(0))
    else:
        cdiff = abs(qsr_mjd[i]-qsr_mjd[i-1])
        diff_arr.append(cdiff)
gap_ind = np.where(np.array(diff_arr)>= (15/(60*24)))[0]
all_session_mjds, all_session_mjds_UT, all_session_mags,all_session_errs = [],[],[],[]
all_session_durs,all_session_counts = [],[]
for i in range(len(gap_ind)):
    if(i == len(gap_ind)-1):
        continue
    session_mjd = qsr_mjd[gap_ind[i]:gap_ind[i+1]]
    session_mjd_int = np.array(session_mjd,dtype=int)
    session_mjd_UT = session_mjd  - np.min(session_mjd_int)
    session_mag = qsr_mag[gap_ind[i]:gap_ind[i+1]]
    session_err = qsr_magerr[gap_ind[i]:gap_ind[i+1]]
    session_dur = (np.max(session_mjd)-np.min(session_mjd))*24
    session_counts = gap_ind[i+1]-gap_ind[i]
    if(session_counts>=10 and session_dur>=0):
        date = mjd_to_date(np.min(session_mjd_int))
        time = dur_to_time(session_dur)
        #print(session_counts, 'datapoints in current session', time)
        #print('session max mjd=',np.max(session_mjd),'and session min mjd=',np.min(session_mjd),"; Date =",date)
        all_session_mjds.append(session_mjd)
        all_session_mjds_UT.append(session_mjd_UT)
        all_session_durs.append(session_dur)
        all_session_counts.append(session_counts)
        all_session_mags.append(session_mag)
        all_session_errs.append(session_err)       
#quasar data
ses_mjd,ses_mag,ses_magerr,ses_mjd_UT = np.array(all_session_mjds[0]),np.array(all_session_mags[0]),np.array(all_session_errs[0]),np.array(all_session_mjds_UT[0])
ses_clipped_mag = sigma_clip(ses_mag,sigma=5,maxiters=3,masked=True)
ses_unmasked_id = ma.nonzero(ses_clipped_mag)
ses_mjd,ses_mag,ses_magerr = np.array(ses_mjd)[ses_unmasked_id],np.array(ses_mag)[ses_unmasked_id],np.array(ses_magerr)[ses_unmasked_id]
ses_mjd_int = np.array(ses_mjd,dtype=int)
ses_mjd_UT = ses_mjd  - np.min(ses_mjd_int)
#star data
str_files = glob.glob('/home/ubuntu/INOV/mywork10/*r.dat', recursive = True)
data = []
for doc in str_files:
    file = open(doc, 'r')
    lines = file.readlines()
    #print(len(lines))
    if(len(lines)<=54):
        continue
    else:
        data.append(doc)
chi_val = []
test_mjd = []
test_mag = []
test_err = []
all_star_mjd,all_star_mjd_UT,all_star_mag,all_star_err = [],[],[],[]
for j in range(len(data)):
    str_obsid,str_mjd,str_mag,str_magerr,str_catflags,str_ra,str_dec = np.loadtxt(data[j],unpack=True,usecols=[0,3,4,5,6,8,9],skiprows=54)
    str_cat_ind = np.where(str_catflags == 0)[0]
    str_obsid,str_mjd,str_mag,str_magerr,str_catflags,str_ra,str_dec = str_obsid[str_cat_ind],str_mjd[str_cat_ind],str_mag[str_cat_ind],str_magerr[str_cat_ind],str_catflags[str_cat_ind],str_ra[str_cat_ind],str_dec[str_cat_ind]
    str_mjd_ind = np.argsort(str_mjd, axis = 0)
    str_obsid,str_mjd,str_mag,str_magerr,str_catflags,str_ra,str_dec = str_obsid[str_mjd_ind],str_mjd[str_mjd_ind],str_mag[str_mjd_ind],str_magerr[str_mjd_ind],str_catflags[str_mjd_ind],str_ra[str_mjd_ind],str_dec[str_mjd_ind]
    str_clipped_mag = sigma_clip(str_mag,sigma=5,maxiters=3,masked=True)
    str_unmasked_id = ma.nonzero(str_clipped_mag)
    str_obsid,str_mjd,str_mag,str_magerr,str_catflags,str_ra,str_dec = np.array(str_obsid)[str_unmasked_id],np.array(str_mjd)[str_unmasked_id],np.array(str_mag)[str_unmasked_id],np.array(str_magerr)[str_unmasked_id],np.array(str_catflags)[str_unmasked_id],np.array(str_ra)[str_unmasked_id],np.array(str_dec)[str_unmasked_id] 
    str_mjd_int = np.array(str_mjd,dtype=int)
    str_mjd_UT = str_mjd  - np.min(str_mjd_int)
    
    nodpts = int(len(ses_mjd) * 0.9) #90% 
    str_mean_mag = np.mean(str_mag)
    qsr_gaia_mag = 17.79
    
    chi_squared = (sum((str_mag-np.mean(str_mag))**2/str_magerr**2))/len(str_mag-1)
    chi_minlimit = 0.8
    chi_maxlimit = 1.5
    
    if(len(str_mjd)<=nodpts):
        continue
    if(str_mean_mag>qsr_gaia_mag):
        continue
    if(chi_squared < chi_maxlimit and chi_squared > chi_minlimit):
        chi_val.append(chi_squared)
        test_mjd.append(str_mjd)
        test_mag.append(str_mag)
        test_err.append(str_magerr)
chi_sort = np.argsort(chi_val)
i = chi_sort[0]
j = chi_sort[1]

C1 = np.where(test_mjd[i][np.in1d(test_mjd[i], test_mjd[j])])
C2 = np.where(test_mjd[j][np.in1d(test_mjd[j], test_mjd[i])])
D1 = np.where(ses_mjd[np.in1d(ses_mjd,test_mjd[i])])
D2 = np.where(test_mjd[i][np.in1d(test_mjd[i],ses_mjd)])
E1 = np.where(ses_mjd[np.in1d(ses_mjd,test_mjd[j])])
E2 = np.where(test_mjd[j][np.in1d(test_mjd[j],ses_mjd)])

QS1 = ses_mag[D1] - test_mag[i][D2]
QS2 = ses_mag[E1] - test_mag[j][E2]
S1S2 = test_mag[i][C1] - test_mag[j][C2]
QS1_e = np.sqrt(ses_magerr[D1]**2 + test_err[i][D2]**2)
QS2_e = np.sqrt(ses_magerr[E1]**2 + test_err[j][E2]**2)
S1S2_e = np.sqrt(test_err[i][C1]**2 + test_err[j][C2]**2)

#critical values
F_c_95 = (f.ppf(0.95, len(QS1)-1, len(QS1)-1))
F_c_99 = (f.ppf(0.999, len(QS1)-1, len(QS1)-1))
#variance of LCs
var_QS1 = ( (QS1 - np.mean(QS1))**2 ).sum()/(len(QS1)-1) 
var_QS2 = ( (QS2 - np.mean(QS2))**2 ).sum()/(len(QS2)-1) 
var_S1S2 =( (S1S2 - np.mean(S1S2))**2 ).sum()/(len(S1S2)-1)
#standard deviation of LCs 
sig_QS1 = np.mean(QS1_e**2)
sig_QS2 = np.mean(QS2_e**2)
sig_S1S2 = np.mean(S1S2_e**2)
eta = 1.5
eta_s = 1.5*1.5
si1 = (max(QS1) - min(QS1) )**2 - 2*( eta**2 * np.mean(QS1_e**2) ) 
si2 = (max(QS2) - min(QS2) )**2 - 2*( eta**2 * np.mean(QS2_e**2) )
var_am = (si1 + si2)/2.0 
precision = 0.5*( np.sqrt( eta**2 * np.mean(QS1_e**2) ) + np.sqrt( eta**2 * np.mean( QS2_e**2) ) ) 
#F-eta test    
f1_eta = var_QS1 / (eta_s * sig_QS1)
f2_eta = var_QS2 / (eta_s * sig_QS2)
variability = []
if f1_eta>=F_c_99:
    V = "Variable"
    variability.append(V)
elif f1_eta > F_c_95 and f1_eta <=F_c_99: 
    V = "Probable Variable"
    variability.append(V)
else:
    V = "Non Variable"
    variability.append(V)
if f2_eta>=F_c_99: 
    V = "Variable"
    variability.append(V)
elif f2_eta > F_c_95 and f2_eta <=F_c_99:
    V = "Probable Variable"
    variability.append(V)
else:
    V = "Non Variable"
    variability.append(V)
result = Table()
result['Blazars SDSS NAME'] = [f"QSR{round(qsr_ra[0],5)}{round(qsr_dec[0],5)}"]
result['Session Date'] = [date]
result['Session Duration'] = [time]
result['Session Datapoints'] = [len(ses_mjd)]
result['f1_eta'] = [round(f1_eta,5)]
result['f2_eta'] = [round(f2_eta,5)]
result['Psi1'] = [round(si1,5)]
result['Psi2'] = [round(si2,5)]
VAR = f"{variability[0]},{variability[1]}"
result['Variability Status']  = [VAR]
#print(result)
ascii.write(result, 'result.dat', overwrite=True)

