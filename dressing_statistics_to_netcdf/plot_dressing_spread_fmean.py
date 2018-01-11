""" this routine will: (1) generate fitted gamma distribution parameters 
    for best-member dressing for individual forecast systems.  These
    dressing parameters are saved to a netcdf file, and diagnostic plots 
    are generated: (2) generate closest-member histograms, make diagnostic
    plots, and save to a netCDF file.
"""

from mpl_toolkits.basemap import Basemap
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from numpy import ma
import os, sys
from netCDF4 import Dataset
import pygrib
from dateutils import daterange, dateshift
from scipy.optimize import curve_fit
from matplotlib import rcParams
import scipy.signal as signal
import scipy.stats as stats
from scipy.stats import gamma
from matplotlib.backends.backend_pdf import PdfPages
rcParams['legend.fontsize']='small'
rcParams['legend.fancybox']=True
rcParams['xtick.labelsize']='medium'


# --- COMMAND LINE INPUT

cmodel = 'NCEP' # sys.argv[1] # NCEP, ECMWF, or CMC
cleade = '48' # sys.argv[2] # ending of periods lead time in hours, e.g., '24', '108'
cempirical = '1'
ileade = int(cleade)
date_forecast = '2016050100' # sys.argv[3]
date_end = dateshift(date_forecast,-ileade)
date_begin = dateshift(date_end, -61*24)


date_list = daterange(date_begin, date_end, 24)
ndates = len(date_list)

# ======================================================================
# Part 1:  reading in the daily data and summing up the information over
#          all case days
# ======================================================================

# ---- read in the array dimensions, precip values where gamma 
#      distribution parameters and fraction zero are calculated, and 
#      precipitation thresholds that are boundaries for light, 
#      medium, and heavy precipitation from netcdf file.  These
#      were previously generated in tally_gamma_statsfull_weighted.f90,
#      called by generate_gammadressing_stats_anymodel.f90


data_directory = '/Users/thamill/precip/ecmwf_data/'   
#data_directory = '/Projects/Reforecast2/netcdf/NationalBlend/'   

if cempirical == '1': 
    infile = data_directory+'gamma_and_closest_hist_stats_' + \
        cmodel + '_2016040100_fhour'+cleade+'.nc'
else:
    infile = data_directory+'gamma_and_closest_hist_stats_' + \
        cmodel + '_2016040100_fhour'+cleade+'_gammaqmap.nc'
print infile
nc = Dataset(infile)    
gamma_threshes = nc.variables['gamma_threshes'][:]  
closest_histogram = nc.variables['closest_histogram'][:]
npv, nmembersx25 = np.shape(closest_histogram)
print 'nmembersx25 = ', nmembersx25
nthreshes = len(gamma_threshes)
climo_pop_thresholds = nc.variables['climo_pop_thresholds'][:] 
thresh_light = nc.variables['thresh_light'][0]
thresh_mod = nc.variables['thresh_mod'][0]
thresh_high = nc.variables['thresh_high'][0]
output_thresh = nc.variables['output_threshes'][:]
nout_thresh = len(output_thresh)
nhist = len(closest_histogram)
climo_pop_thresholds_save = np.copy(climo_pop_thresholds)
n_climocats = len(climo_pop_thresholds) + 1
nc.close()

# ---- declare the arrays that will hold the summary information
#      tallied over many days.  

sum_of_pertval_total = np.zeros((3,3,nthreshes),dtype=np.float64)
sum_of_pertval2_total = np.zeros((3,3,nthreshes),dtype=np.float64)
npertval_count_total = np.zeros((3,3,nthreshes),dtype=np.float64)

# ---- loop over dates, read in data, and add the latest day's data
#      to the overall summation arrays

for idate, date in zip(range(ndates), date_list):
    infile = data_directory+'gamma_and_closest_hist_stats_' + cmodel + '_' + \
            date + '_fhour'+cleade+'.nc'
    print infile
    fexist  = os.path.exists(infile)
    if fexist:
        try: 
            
            #  --- read in 
            
            nc = Dataset(infile)               
            sum_of_pertval = nc.variables['sum_of_pertval'][:,:,:]
            sum_of_pertval2 = nc.variables['sum_of_pertval2'][:,:,:]
            npertval_count = nc.variables['npertval_count'][:,:,:]
            
            
            # ---- add to sum
            
            sum_of_pertval_total = sum_of_pertval_total + sum_of_pertval
            sum_of_pertval2_total = sum_of_pertval2_total + sum_of_pertval2
            npertval_count_total = npertval_count_total + npertval_count
            
            print 'sum today overall = ', np.sum(npertval_count), np.sum(npertval_count_total)
            
            nc.close()
                
        except (IOError, ValueError, RuntimeError):
            print 'Error reading ', infile
            

# ======================================================================
# Part 2:  process the data to spread of dressing distribution as a 
#          function of mean precipitation amount
# ======================================================================

# ---- calculate the empirical gamma distribution shape, scale, and fraction zero
#      for selected values of precipitation amount, defined by gamma_threshold
#      previously read in.  These parameters also vary with the precipitation 
#      amount category (4 categories, with boundaries defined by thresh_light, 
#      thresh_mod, thresh_high) and 3 closest-histogram ranks (lowest, interior,
#      highest).
 
spread_fgamma_threshes = -99.99 * np.ones((3,3,nthreshes),dtype=np.float32) # empirical
# dimid_3d =  (/ nthreshes_dimid, nprecipcats_dimid, nlomidhirank_dimid /)
#! ---- Here are the variables we intend to populate with the forecast, analysis information.
#!      dimension 2 below is for light, mod, heavy precip, with thresholds between
#!      them set by inputs thresh_light, thresh_mod, thresh_high.
#!      dimension 3 below is for lowest sorted member(1), intermediate(2), highest(3) 

for ipamount in range(3):
    for ilomidhi in range(0,3):
        print 'ilomidhi, ipamount, npertval_count_total[ilomidhi,ipamount,:] = ', \
            ilomidhi, ipamount, npertval_count_total[ilomidhi,ipamount,:]
        for ithresh in range(nthreshes):
            if npertval_count_total[ilomidhi,ipamount,ithresh] > 25:
                xbar = sum_of_pertval_total[ilomidhi,ipamount,ithresh] / \
                    npertval_count_total[ilomidhi,ipamount,ithresh]
                spread_fgamma_threshes[ilomidhi,ipamount,ithresh] = \
                    (sum_of_pertval2_total[ilomidhi,ipamount,ithresh] -
                    npertval_count_total[ilomidhi,ipamount,ithresh]*xbar**2)/ \
                    (npertval_count_total[ilomidhi,ipamount,ithresh]-1.)
                spread_fgamma_threshes[ilomidhi,ipamount,ithresh] = \
                    np.sqrt(spread_fgamma_threshes[ilomidhi,ipamount,ithresh])
            

# ---- create masked array with missing data.

spread_fgamma_threshes_m = ma.masked_less(spread_fgamma_threshes, 0.0 )

# ---- plot the dressing spread as a function of the precip amount

cpamount = ['light mean Precipitation (0.01 - 2 mm)', \
    'moderate mean precipitation (2-10 mm)', 'heavy mean precipitation (> 10 mm)']
cpamount_short = ['lightmean','modmean','heavymean']

clomidhi = ['lowest member','intermediate members','highest member']   
clomidhi_short = ['lowmbr','intmbrs','hibmr'] 
    
for ipamount in range(3): # for light, mod, heavy precip, with thresholds 
    for ilomidhi in range(0,3): # lowest sorted member(1), intermediate(2), highest(3)

        fig = plt.figure(figsize=(6.5,7.))
        ctitle = 'Dressing spread for '+clomidhi[ilomidhi]+'\n'+\
            cpamount[ipamount]
        ytitle = 'Average dressing uncertainty (mm)'
        xtitle = 'Mean precipitation amount (mm)'
        axlocn = [0.15,0.12,0.8,0.75]
        a1 = fig.add_axes(axlocn) 
        a1.set_title(ctitle,fontsize=15)
        a1.set_ylabel(ytitle,fontsize=12)
        a1.set_xlabel(xtitle,fontsize=12)
        xmax = np.max(spread_fgamma_threshes_m[ilomidhi,ipamount,:])
        a1.plot(gamma_threshes, spread_fgamma_threshes_m[ilomidhi,ipamount,:], '-', lw=2)
        #a1.plot([0,0],[xmax,xmax],linestyle='--',linewidth=0.5)
        a1.grid(color='Gray',lw=0.2,linestyle='--')
    
        outfile = 'dressing_spread_fmean_'+cmodel+'_'+\
            cpamount_short[ipamount]+'_'+clomidhi_short[ilomidhi]+'_'+cleade+'h.pdf'        
        fig.savefig(outfile)
        print 'saving plot to file = ',outfile



fig = plt.figure(figsize=(9,3.5))
plt.suptitle('NCEP +36 to +48 h dressing distribution spread', fontsize=17)
for ipamount in range(3):
    if ipamount == 0:
        ctitle = '(a) Light mean precip'
        axlocn = [0.06,0.14,0.24,0.65]
        xrange = [-1,26]
        xticks = [0,5,10,15,20,25]
    elif ipamount == 1:
        ctitle = '(b) Moderate mean precip'
        axlocn = [0.39,0.14,0.24,0.65]
        xrange = [-1,51]
        xticks = [0,10,20,30,40,50]
    else:
        ctitle = '(c) Heavy mean precip'
        axlocn = [0.72,0.14,0.24,0.65]
        xrange = [-1,121]
        xticks = [0,20,40,60,80,100,120]
        
        
    ytitle = 'Dressing spread (mm)'
    xtitle = 'Member precip amount (mm)'
    a1 = fig.add_axes(axlocn) 
    a1.set_title(ctitle,fontsize=13)
    a1.set_ylabel(ytitle,fontsize=12)
    a1.set_xlabel(xtitle,fontsize=11)
    a1.set_xlim(xrange)
    a1.set_xticks(xticks)
    a1.grid(color='Gray',lw=0.2,linestyle='--')
    a1.plot(gamma_threshes, spread_fgamma_threshes_m[0,ipamount,:], '-', color='Red',lw=2,label='Lowest member')
    a1.plot(gamma_threshes, spread_fgamma_threshes_m[1,ipamount,:], '-', color='RoyalBlue',lw=2,label='Intermediate')
    a1.plot(gamma_threshes, spread_fgamma_threshes_m[2,ipamount,:], '-', color='LimeGreen',lw=2,label='Highest member')
    if ipamount <= 1: 
        a1.legend(loc=2)
    else:
        a1.legend(loc=1)

outfile = 'NCEP_48h_dressing_spread_fmean.pdf'        
fig.savefig(outfile)
print 'saving plot to file = ',outfile    