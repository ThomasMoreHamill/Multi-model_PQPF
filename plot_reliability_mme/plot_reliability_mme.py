""" 
"""

from mpl_toolkits.basemap import Basemap
import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import rcParams
import numpy as np
from numpy import ma
import os, sys
from netCDF4 import Dataset
from verify_relia_bss_mme import verify_relia_bss_mme
import pygrib
from dateutils import daterange, dateshift
rcParams['legend.fontsize']='small'
rcParams['legend.fancybox']=True
rcParams['xtick.labelsize']='medium'

cleade = sys.argv[1]
cthresh = sys.argv[2]  

cyyyymmddhh_start = '2016040100'
cyyyymmddhh_end = '20160701000'
date_list = daterange(cyyyymmddhh_start, cyyyymmddhh_end, 24)

ileade = int(cleade)
ileadb = ileade - 12
cleadb = str(ileadb)
#if ileadb < 100: cleadb='0'+cleadb
if cthresh == 'POP':
    rthresh = 0.254
    cthresh_plot = 'POP'
else:
    rthresh = float(cthresh)
    cthresh_plot = r'$\geq$ '+cthresh+' mm'

nclasses = 21 # 0 to 100% probability by 5%
nxa = 464 # CCPA 1/8 degree grid over CONUS
nya = 224 #

# ---- read in precipitation analyses

date_list_anal = daterange(cyyyymmddhh_start, cyyyymmddhh_end, 24)
print 'date_list_anal = ',date_list_anal
ndates = len(date_list_anal)
apcp_anal_t = np.zeros((nxa,nya,ndates), dtype=np.float32)
mninetynine = -99.99*np.ones((nya,nxa), dtype=np.float32)

for idate, date in zip(range(ndates), date_list_anal):

    print '------------------------- getting precipitation for idate = ',idate,date
    date_fearly = dateshift(date, ileade-6)
    date_flate = dateshift(date, ileade)
    print 'date_fearly, date_flate = ', date_fearly, date_flate

    # --- read in the first of the two precip analysis files, 6 hourly
    
    infile = '/data/thamill/Rf2_tests/ccpa_v1/0.125d/ccpa.'+date_fearly[0:8]+\
        '/18/ccpa.t18z.06h.0p125.conus.gb2'
    #print infile
    fexist1 = os.path.exists(infile)
    if fexist1:
        afile1 = pygrib.open(infile)
        grb1 = afile1.select()[0]    # --- read in the first of the two precip analysis files, 6 hourly
    
    infile2 = '/data/thamill/Rf2_tests/ccpa_v1/0.125d/ccpa.'+date_flate[0:8]+\
        '/00/ccpa.t00z.06h.0p125.conus.gb2'
    #print infile2
    fexist2 = os.path.exists(infile2)
    if fexist2:
        afile2 = pygrib.open(infile2)
        grb2 = afile2.select()[0]
    
    if fexist1 and fexist2:
        apcp_anal = grb1.values + grb2.values
        apcp_anal = np.where(apcp_anal > 500., mninetynine, apcp_anal)
    else:
        print 'no analysis data for this date'
        apcp_anal = -99.99*np.ones((nya,nxa),dtype=np.float32)
    print np.shape(apcp_anal)

    apcp_anal_t[:,:,idate] = np.transpose(apcp_anal[:,:])
    #print 'min, max apcp_anal = ', np.min(apcp_anal), np.max(apcp_anal)

    afile1.close()
    afile2.close()

# ---- call fortran routine to do the reading in of forecasts and the generation
#      of reliability, frequency of use, and bss information.

relia = np.zeros(21, dtype=np.float32)
relia_NCEP = np.zeros(21, dtype=np.float32)
relia_CMC = np.zeros(21, dtype=np.float32)
relia_ECMWF = np.zeros(21, dtype=np.float32)
relia_MME = np.zeros(21, dtype=np.float32)

relia_05 = np.zeros(21, dtype=np.float32)
relia_NCEP_05 = np.zeros(21, dtype=np.float32)
relia_CMC_05 = np.zeros(21, dtype=np.float32)
relia_ECMWF_05 = np.zeros(21, dtype=np.float32)
relia_MME_05 = np.zeros(21, dtype=np.float32)

relia_95 = np.zeros(21, dtype=np.float32)
relia_NCEP_95 = np.zeros(21, dtype=np.float32)
relia_CMC_95 = np.zeros(21, dtype=np.float32)
relia_ECMWF_95 = np.zeros(21, dtype=np.float32)
relia_MME_95 = np.zeros(21, dtype=np.float32)

frequse_95 = np.zeros(21, dtype=np.float32)
frequse_NCEP_95 = np.zeros(21, dtype=np.float32)
frequse_CMC_95 = np.zeros(21, dtype=np.float32)
frequse_ECMWF_95 = np.zeros(21, dtype=np.float32)
frequse_MME_95 = np.zeros(21, dtype=np.float32)

bss_NCEP_daily = np.zeros(ndates,dtype=np.float32)
bss_CMC_daily = np.zeros(ndates,dtype=np.float32)
bss_ECMWF_daily = np.zeros(ndates,dtype=np.float32)
bss_MME_daily = np.zeros(ndates,dtype=np.float32)

bss = 0.
bss_NCEP = 0.
bss_CMC = 0.
bss_ECMWF = 0.
bss_MME = 0.

nxa = 464 # CCPA 1/8 degree grid over CONUS
nya = 224 #

relia_NCEP, relia_NCEP_05, relia_NCEP_95, frequse_NCEP, bss_NCEP, \
    relia_CMC, relia_CMC_05, relia_CMC_95, frequse_CMC, bss_CMC, \
    relia_ECMWF, relia_ECMWF_05, relia_ECMWF_95, frequse_ECMWF, bss_ECMWF, \
    relia_MME, relia_MME_05, relia_MME_95, frequse_MME, bss_MME = \
    verify_relia_bss_mme(cleade, nclasses, rthresh, \
    date_list_anal, apcp_anal_t, nxa, nya, ndates)
    
# ---- now make multi-panel reliability diagrams for each of the forecasts

fig = plt.figure(figsize=(6.,6.15))
plt.suptitle(cthresh_plot+' reliability diagrams for\n+'+cleadb+' to +'+\
    cleade+' hour post-processed forecasts',fontsize=15)
rcParams['xtick.labelsize']='xx-small'
rcParams['ytick.labelsize']='xx-small'

for itype in range(4):
    if itype == 0:
        ov_ax = [.08,.52,.36,.34]
        a1 = fig.add_axes(ov_ax)
        relia_out = relia_NCEP
        relia05 = relia_NCEP_05
        relia95 = relia_NCEP_95
        frequse_out = frequse_NCEP
        bss_out = bss_NCEP
        extra = False
        ctitle = '(a) NCEP'
    elif itype == 1:
        ov_ax = [.6,.52,.36,.34]
        a1 = fig.add_axes(ov_ax)
        relia_out = relia_CMC
        relia05 = relia_CMC_05
        relia95 = relia_CMC_95
        frequse_out = frequse_CMC
        bss_out = bss_CMC
        extra = False
        ctitle = '(b) CMC'
    elif itype == 2:
        ov_ax = [.08,.06,.36,.34]
        a1 = fig.add_axes(ov_ax)
        relia_out = relia_ECMWF
        relia05 = relia_ECMWF_05
        relia95 = relia_ECMWF_95
        frequse_out = frequse_ECMWF
        bss_out = bss_ECMWF
        extra = False
        ctitle = '(c) ECMWF'
    elif itype == 3:
        ov_ax = [.6,.06,.36,.34]
        a1 = fig.add_axes(ov_ax)
        relia_out = relia_MME
        relia05 = relia_MME_05
        relia95 = relia_MME_95
        frequse_out = frequse_MME
        bss_out = bss_MME
        extra = False
        ctitle = '(d) MME weighted average'

    # -- make reliability diagram for statistically downscaled and RAW

    a1.set_title(ctitle,fontsize=11)

    # --- add basic reliability diagram 

    yerrs = np.squeeze(np.array([[100.*(relia_out-relia05)],[100.*(relia95-relia_out)]]))
    
    #rcParams['xtick.labelsize']='medium'
    #rcParams['ytick.labelsize']='medium'
    strbss = 'BSS = %0.3f' %(bss_out)
    relia_m = ma.array(relia_out)
    relia_m = ma.masked_where(relia_m < 0.0, relia_m)
    probs = np.arange(nclasses) * 100./np.real(nclasses-1)
    a1.plot(probs,100.*relia_m,'o-',color='r',markersize=3)
    a1.errorbar(probs,100.*relia_m,yerr=yerrs,fmt='-',color='r')
    if extra == True: a1.plot(probs,100.*relia_mo,'o-',color='b')
    a1.plot([0,100],[0,100],'--',color='k')
    a1.set_ylabel('Observed Relative Frequency (%)',fontsize=9)
    a1.set_xlabel('Forecast Probability (%)',fontsize=10)
    a1.set_ylim(-1,101)
    a1.set_xlim(-1,101)

    # -- BSS inserted here

    a1.text(51,7,strbss,fontsize=10)
    
    # --- Frequency of usage inset diagram

    h0 = ov_ax[0] + ov_ax[2]*0.2
    h1 = ov_ax[1] + ov_ax[3]*0.6
    h2 = ov_ax[2]*0.3
    h3 = ov_ax[3]*0.3
    a2 = fig.add_axes([h0,h1,h2,h3])

    a2.bar(probs,frequse_out,width=5,bottom=0.001,log=True,color='red',\
        edgecolor='black',align='center')
    if extra == True:
        for i in range(len(frequse_out_o)):
            a2.plot([probs[i]-2.5, probs[i]+2.5], \
                [frequse_out_o[i],frequse_out_o[i]], color='Blue', linewidth=2)
    a2.set_xticks([0,25,50,75,100])
    a2.set_xlim(-3,103)
    a2.set_ylim(0.0,1.)
    a2.set_title('Frequency of usage',fontsize=7)
    a2.set_xlabel('Probability',fontsize=7)
    a2.set_ylabel('Frequency',fontsize=7)
    a2.hlines([.01,.1],0,100,linestyles='dashed',colors='black')


plot_title = 'relia_MME_'+cthresh+'_hour'+cleade+'.pdf'
print 'saving plot to file = ',plot_title
plt.savefig(plot_title)
print 'Plot done'


print 'relia_MME = ',relia_MME
print 'relia_MME_05 = ',relia_MME_05
print 'relia_MME_95 = ',relia_MME_95

    
