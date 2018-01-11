""" this routine will: (1) generate fitted gamma distribution parameters 
    for best-member dressing for individual forecast systems.  These
    dressing parameters are saved to a netcdf file, and diagnostic plots 
    are generated: (2) generate closest-member histograms, make diagnostic
    plots, and save to a netCDF file.
"""
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from numpy import ma
import os, sys
from netCDF4 import Dataset
from dateutils import daterange, dateshift
from matplotlib import rcParams

from scipy.optimize import curve_fit
from matplotlib import rcParams
import scipy.signal as signal
import scipy.stats as stats

rcParams['legend.fontsize']='medium'
rcParams['legend.fancybox']=True
rcParams['xtick.labelsize']='medium'

# --- COMMAND LINE INPUT

cleade = sys.argv[1] # ending of periods lead time in hours, e.g., '24', '108'
ileade = int(cleade)
date_forecast = sys.argv[2]
cempirical = sys.argv[3]
cleadb = str(int(cleade)-12)

# ======================================================================
# Part 1:  reading in data
# ======================================================================
# 

data_directory = '/Users/thamill/precip/ecmwf_data/'   
if cempirical == '1':
    infile = data_directory+'NCEP/closest_histogram_NCEP_date='+\
        date_forecast+'_lead='+cleade+'.nc'
else:
    infile = data_directory+'NCEP/closest_histogram_NCEP_date='+\
        date_forecast+'_lead='+cleade+'_gammaqmap.nc'
print infile
nc = Dataset(infile)    
precip_thresholds = nc.variables['precip_thresholds'][:]  
closest_histogram_NCEP = nc.variables['closest_histogram'][:,:]
nc.close()

data_directory = '/Users/thamill/precip/ecmwf_data/'   
if cempirical == '1':
    infile = data_directory+'CMC/closest_histogram_CMC_date='+\
        date_forecast+'_lead='+cleade+'.nc'
else:
    infile = data_directory+'CMC/closest_histogram_CMC_date='+\
        date_forecast+'_lead='+cleade+'_gammaqmap.nc'
print infile
nc = Dataset(infile)    
closest_histogram_CMC = nc.variables['closest_histogram'][:,:]
nc.close()

data_directory = '/Users/thamill/precip/ecmwf_data/'   
if cempirical == '1':
    infile = data_directory+'ECMWF/closest_histogram_ECMWF_date='+\
        date_forecast+'_lead='+cleade+'.nc'
else:
    infile = data_directory+'ECMWF/closest_histogram_ECMWF_date='+\
        date_forecast+'_lead='+cleade+'_gammaqmap.nc'
print infile
nc = Dataset(infile)    
closest_histogram_ECMWF = nc.variables['closest_histogram'][:,:]
nc.close()

print 'shape histogram = ',np.shape(closest_histogram_ECMWF)

for i in range(3):
    csavgol = np.copy(closest_histogram_NCEP[i,:])
    work = np.copy(csavgol[4:-4])
    csavgol[4:-4] = signal.savgol_filter(work, 9, 2, mode='interp')
    closest_histogram_NCEP[i,:] = csavgol[:] / np.sum(csavgol)
    
    csavgol = np.copy(closest_histogram_CMC[i,:])
    work = np.copy(csavgol[4:-4])
    csavgol[4:-4] = signal.savgol_filter(work, 9, 2, mode='interp')
    closest_histogram_CMC[i,:] = csavgol[:] / np.sum(csavgol)
    
    csavgol = np.copy(closest_histogram_ECMWF[i,:])
    work = np.copy(csavgol[4:-4])
    csavgol[4:-4] = signal.savgol_filter(work, 9, 2, mode='interp')
    closest_histogram_ECMWF[i,:] = csavgol[:] / np.sum(csavgol)
        
# ---- now make closest histogram plots

fig = plt.figure(figsize=(6.,9.))
plt.suptitle('Closest-member histograms for\n'+date_forecast+\
    ', lead = +'+cleadb+'-'+cleade+' h',fontsize=18)
for i in range(3):
    if i == 0:
        ctitle = r'(a) 0.01 mm $\leq$ mean precip. < 2 mm '
        axlocn = [0.15,0.66,0.8,0.21]
    elif i == 1:
        ctitle = r'(b) 2 mm $\leq$ mean precip. < 10 mm '
        axlocn = [0.15,0.36,0.8,0.21]
    else:
        ctitle = r'(c) mean precip. $\geq$ 10 mm'
        axlocn = [0.15,0.06,0.8,0.21]

    a1 = fig.add_axes(axlocn)
    a1.set_title(ctitle,fontsize=14)
    a1.set_xlabel('Fraction between lowest and highest rank',fontsize=12)
    a1.set_ylabel('Fraction',fontsize=12)
    a1.set_xlim(0,1)
    a1.set_ylim(0.0003,0.1)
    a1.set_yscale('log')
    a1.grid(color='Gray',lw=0.2,linestyle='--')

    ctext_low_N = "{:.3f}".format(closest_histogram_NCEP[i,0])
    ctext_high_N = "{:.3f}".format(closest_histogram_NCEP[i,-1])    
    a1.plot(np.real(range(1,501))/500., closest_histogram_NCEP[i,:],'-',\
        color='Red',label=r'$\leftarrow \ $'+ctext_low_N+' NCEP '+ctext_high_N+r'$\ \rightarrow$',linewidth=2)

    ctext_low_C = "{:.3f}".format(closest_histogram_CMC[i,0])
    ctext_high_C = "{:.3f}".format(closest_histogram_CMC[i,-1])
    a1.plot(np.real(range(1,501))/500., closest_histogram_CMC[i,:],'-',\
        color='LimeGreen',label=r'$\leftarrow \ $'+ctext_low_C+' CMC '+ctext_high_C+r'$\ \rightarrow$',linewidth=2)
    #a1.text(0.01,closest_histogram_CMC[i,0],ctext_low,color='LightGreen')
    #a1.text(0.88,closest_histogram_CMC[i,-1],ctext_high,color='LightGreen')
    
    ctext_low_E = "{:.3f}".format(closest_histogram_ECMWF[i,0])
    ctext_high_E = "{:.3f}".format(closest_histogram_ECMWF[i,-1])
    a1.plot(np.real(range(1,1251))/1250., closest_histogram_ECMWF[i,:],'-',\
        color='RoyalBlue',label=r'$\leftarrow \ $'+ctext_low_E+' ECMWF '+ctext_high_E+r'$\ \rightarrow$',linewidth=2)
    #a1.text(0.01,closest_histogram_ECMWF[i,0],ctext_low,color='RoyalBlue')
    #a1.text(0.88,closest_histogram_ECMWF[i,-1],ctext_high,color='RoyalBlue')
 
    a1.legend(loc=9)

outfile = 'closest_member_histogram_'+date_forecast+'_'+cleade+'h.pdf'
fig.savefig(outfile)
print 'saving plot to file = ',outfile










