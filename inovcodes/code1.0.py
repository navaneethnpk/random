##Date=Mar 09, 2022
import numpy as np
import matplotlib.pyplot as plt
import numpy.ma as ma
from astropy.stats import sigma_clip
import glob
import math
import datetime as dt
import wget
from astroquery.gaia import Gaia
from astropy import units as u
from astropy.coordinates import Angle
from scipy.stats import f

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

for i in range(len(qsr_RA)):
    link = f"https://irsa.ipac.caltech.edu/cgi-bin/ZTF/nph_light_curves?POS=CIRCLE+{qsr_RA[i]}+{qsr_DEC[i]}+0.000416667&BANDNAME=r&BAD_CATFLAGS_MASK=32768&FORMAT=ipac_table"
    RA_round = round(qsr_RA[i],5)
    DEC_round = round(qsr_DEC[i],5)
    filename = f"{RA_round}{DEC_round}r.dat"
    wget.download(link, filename)

qsr_files = glob.glob('/home/ubuntu/INOV/mywork14/*r.dat', recursive = True)
qsrdata = []
for doc in qsr_files:
    file = open(doc, 'r')
    lines = file.readlines()
    #print(len(lines))
    if(len(lines)<=54):
        continue
    else:
        qsrdata.append(doc)

qsr_obsid,qsr_mjd,qsr_mag,qsr_magerr,qsr_catflags,qsr_ra,qsr_dec = np.loadtxt(qsrdata[0],unpack=True,usecols=[0,3,4,5,6,8,9],skiprows=54)
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
        #fig = plt.figure()
        #plt.errorbar(session_mjd_UT,session_mag,yerr=session_err,color='red', ls='none',linewidth=1, marker='.',mfc='white', capsize=2, capthick=1, markersize=5, ecolor='black')
        #plt.gca().invert_yaxis()
        #plt.minorticks_on()
        #plt.tick_params(axis ='both', which ='major', labelsize = 10, pad = 5, colors ='k')
        #plt.tick_params(axis ='both', which ='minor', labelsize = 8, colors ='b')
        #plt.xlabel("MJD")
        #plt.ylabel("MAG")
        #plt.title(f"Session {len(all_session_mjds)}\nDate= {date}, Duration= {time}")
        #plt.savefig(f"session{len(all_session_mjds)}.png",bbox_inches ="tight",pad_inches = 0.3,facecolor ="white",edgecolor ='white',orientation ='landscape',dpi=100)
        #plt.show()

#quasar data
ses_mjd,ses_mag,ses_magerr,ses_mjd_UT = np.array(all_session_mjds[0]),np.array(all_session_mags[0]),np.array(all_session_errs[0]),np.array(all_session_mjds_UT[0])

qsr_radius = Angle(2/60, u.arcminute)
qsr_rlimit = qsr_radius.degree 
qsr_gaia_ra = qsr_RA[0]
qsr_gaia_dec = qsr_DEC[0]

qsr_job = Gaia.launch_job("SELECT TOP 2000 "
                      "gaia_source.source_id,gaia_source.ra,gaia_source.dec,gaia_source.parallax,gaia_source.parallax_error,gaia_source.pm,gaia_source.pmra,gaia_source.pmra_error,gaia_source.pmdec,gaia_source.pmdec_error,gaia_source.phot_g_mean_mag "
                      "FROM gaiaedr3.gaia_source "
                      f"WHERE CONTAINS(POINT('ICRS',gaiaedr3.gaia_source.ra,gaiaedr3.gaia_source.dec),CIRCLE('ICRS',{qsr_gaia_ra},{qsr_gaia_dec},{qsr_rlimit}))=1 ",
                      dump_to_file=True, output_format='votable')
#print(qsr_job.outputFile)
gaia_qsr = qsr_job.get_results()
qsr_mag_g = gaia_qsr["phot_g_mean_mag"]

str_radius = Angle(20, u.arcminute)
str_rlimit = str_radius.degree 
str_ra = 110.748665 #qsr_RA[0]
str_dec = -7.5263276 #qsr_DEC[0]
pmlimit = 20
magupperlimit = 17.99 #qsr_mag_g+0.5
maglowerlimit = 15.49 #qsr_mag_g-1.5

str_job = Gaia.launch_job("SELECT TOP 2000 "
                      "gaia_source.source_id,gaia_source.ra,gaia_source.dec,gaia_source.parallax,gaia_source.parallax_error,gaia_source.pm,gaia_source.pmra,gaia_source.pmra_error,gaia_source.pmdec,gaia_source.pmdec_error,gaia_source.phot_g_mean_mag,gaia_source.dr2_radial_velocity,gaia_source.dr2_radial_velocity_error "
                      "FROM gaiaedr3.gaia_source "
                      f"WHERE CONTAINS(POINT('ICRS',gaiaedr3.gaia_source.ra,gaiaedr3.gaia_source.dec),CIRCLE('ICRS',{str_ra},{str_dec},{str_rlimit}))=1 "
                      f"AND (gaiaedr3.gaia_source.pm>={pmlimit} "
                      f"AND gaiaedr3.gaia_source.phot_g_mean_mag>={maglowerlimit} "
                      f"AND gaiaedr3.gaia_source.phot_g_mean_mag<={magupperlimit}) ", 
                      dump_to_file=True, output_format='votable')

#print(str_job.outputFile)
gaia_str = str_job.get_results()

ztf_ra = np.array(gaia_str['ra'])
ztf_dec = np.array(gaia_str['dec'])

mjdmin = np.min(ses_mjd)
mjdmax = np.max(ses_mjd)

for i in range(len(ztf_ra)):
    link = f"https://irsa.ipac.caltech.edu/cgi-bin/ZTF/nph_light_curves?POS=CIRCLE+{ztf_ra[i]}+{ztf_dec[i]}+0.000416667&TIME={mjdmin}+{mjdmax}&BANDNAME=r&BAD_CATFLAGS_MASK=32768&FORMAT=ipac_table"
    RA_round = round(ztf_ra[i],5)
    DEC_round = round(ztf_dec[i],5)
    filename = f"{RA_round}{DEC_round}r.dat"
    wget.download(link, filename)

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
#print("Number of sources with data =",len(data))

all_star_mjd,all_star_mjd_UT,all_star_mag,all_star_err = [],[],[],[]
for j in range(len(data)):
    str_obsid,str_mjd,str_mag,str_magerr,str_catflags,str_ra,str_dec = np.loadtxt(data[j],unpack=True,usecols=[0,3,4,5,6,8,9],skiprows=54)
    str_cat_ind = np.where(str_catflags == 0)[0]
    str_obsid,str_mjd,str_mag,str_magerr,str_catflags,str_ra,str_dec = str_obsid[str_cat_ind],str_mjd[str_cat_ind],str_mag[str_cat_ind],str_magerr[str_cat_ind],str_catflags[str_cat_ind],str_ra[str_cat_ind],str_dec[str_cat_ind]
    str_mjd_ind = np.argsort(str_mjd, axis = 0)
    str_obsid,str_mjd,str_mag,str_magerr,str_catflags,str_ra,str_dec = str_obsid[str_mjd_ind],str_mjd[str_mjd_ind],str_mag[str_mjd_ind],str_magerr[str_mjd_ind],str_catflags[str_mjd_ind],str_ra[str_mjd_ind],str_dec[str_mjd_ind]
    str_clipped_mag = sigma_clip(str_mag,sigma=3,maxiters=3,masked=True)
    str_unmasked_id = ma.nonzero(str_clipped_mag)
    str_obsid,str_mjd,str_mag,str_magerr,str_catflags,str_ra,str_dec = np.array(str_obsid)[str_unmasked_id],np.array(str_mjd)[str_unmasked_id],np.array(str_mag)[str_unmasked_id],np.array(str_magerr)[str_unmasked_id],np.array(str_catflags)[str_unmasked_id],np.array(str_ra)[str_unmasked_id],np.array(str_dec)[str_unmasked_id] 
    str_mjd_int = np.array(str_mjd,dtype=int)
    str_mjd_UT = str_mjd  - np.min(str_mjd_int)
    
    chi_squared = (sum((str_mag-np.mean(str_mag))**2/str_magerr**2))/len(str_mag-1)
    #print(f'chi_square of star {j} = {chi_squared}')
    chi_minlimit = 0.5
    chi_maxlimit = 1.5
    if(chi_squared < chi_maxlimit and chi_squared > chi_minlimit):
        #print(f'chi_square of star{j} = {chi_squared}')
        c1 = np.where(np.intersect1d(ses_mjd,str_mjd))
        str_mjd,str_mjd_UT,str_mag,str_magerr = str_mjd[c1],str_mjd_UT[c1],str_mag[c1],str_magerr[c1]
        
        all_star_mjd.append(str_mjd)
        all_star_mjd_UT.append(str_mjd_UT)
        all_star_mag.append(str_mag)
        all_star_err.append(str_magerr)
c = []
for k in range(len(all_star_mjd)):
    if(k == len(all_star_mjd)-1):
        continue
    c2 = np.where(np.intersect1d(all_star_mjd[k],all_star_mjd[k+1]))  
    c[k].append(all_star_mjd[c2])
    

#Common MJD values 
c1 = np.where(np.intersect1d(ses_mjd,sr1_mjd))
c2 = np.where(np.intersect1d(ses_mjd,sr2_mjd))
c3 = np.where(np.intersect1d(sr1_mjd,sr2_mjd))

QS1 = ses_mag[c1] - sr1_mag[c1]
QS2 = ses_mag[c2] - sr2_mag[c2]
S1S2 = sr1_mag[c3] - sr2_mag[c3]

QS1_e = np.sqrt(ses_magerr[c1]**2 + sr1_magerr[c1]**2)
QS2_e = np.sqrt(ses_magerr[c2]**2 + sr2_magerr[c2]**2)
S1S2_e = np.sqrt(sr1_magerr[c3]**2 + sr2_magerr[c3]**2)


###### critical values
F_c_95 = (f.ppf(0.95, len(QS1)-1, len(QS1)-1))
F_c_99 = (f.ppf(0.999, len(QS1)-1, len(QS1)-1))
###### variance of LCs
var_QS1 = ( (QS1 - np.mean(QS1))**2 ).sum()/(len(QS1)-1) 
var_QS2 = ( (QS2 - np.mean(QS2))**2 ).sum()/(len(QS2)-1) 
var_S1S2 =( (S1S2 - np.mean(S1S2))**2 ).sum()/(len(S1S2)-1)
###### standard deviation of LCs 
sig_QS1 = np.mean(QS1_e**2)
sig_QS2 = np.mean(QS2_e**2)
sig_S1S2 = np.mean(S1S2_e**2)
eta = 1.5
eta_s = 1.5*1.5
# Take eta_s = 1 for our case
si1, si2 = (max(QS1) - min(QS1) )**2 - 2*( eta**2 * np.mean(QS2_e**2) ), (max(QS2) - min(QS2) )**2 - 2*( eta**2 * np.mean(QS2_e**2) )
var_am = (si1 + si2)/2.0 
precision = 0.5*( np.sqrt( eta**2 * np.mean(QS1_e**2) ) + np.sqrt( eta**2 * np.mean( QS2_e**2) ) ) 
########### F-eta test ######################################    
f1_eta = var_QS1 / (eta_s * sig_QS1)
f2_eta = var_QS2 / (eta_s * sig_QS2)
 
if f1_eta>=F_c_99:    
        print("Variable")
elif f1_eta > F_c_95 and f1_eta <=F_c_99: 
        print("Probable Variable")   
else:
        print("Non Variable")
if f2_eta>=F_c_99:    
        print("Variable")  
elif f2_eta > F_c_95 and f2_eta <=F_c_99:
        print("Probable Variable")
else:
        print("Non Variable") 

fig = plt.figure()
#plt.plot(ses_mjd[c1], QS1, marker = '.', markersize = 2, linewidth=0.5, color='k')
plt.errorbar(ses_mjd[c1],QS1,yerr=QS1_e,color='red', ls='none',linewidth=1, marker='.',mfc='white', capsize=2, capthick=1, markersize=5, ecolor='black')
plt.gca().invert_yaxis()
plt.minorticks_on()
plt.tick_params(axis ='both', which ='major', labelsize = 10, pad = 5, colors ='k')
plt.tick_params(axis ='both', which ='minor', labelsize = 8, colors ='b')
plt.xlabel("MJD")
plt.ylabel("MAG")
plt.title(f"Q-S1")
#plt.savefig(f"session-Q-S2.png",bbox_inches ="tight",pad_inches = 0.3,facecolor ="white",edgecolor ='white',orientation ='landscape',dpi=100)
plt.show()









