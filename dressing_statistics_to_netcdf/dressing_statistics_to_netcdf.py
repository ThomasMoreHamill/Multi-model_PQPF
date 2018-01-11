""" this routine will: (1) generate fitted gamma distribution parameters 
    for best-member dressing for individual forecast systems.  These
    dressing parameters are saved to a netcdf file, and diagnostic plots 
    are generated: (2) generate closest-member histograms, make diagnostic
    plots, and save to a netCDF file.  

    There is now the option of using either empirical CDFs or fitted
    Gamma CDFs as input information.   It is possible the dressing statistics
    may be different depending on which type of CDF is used for the
    quantile mapping.
"""

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

rcParams['legend.fontsize']='medium'
rcParams['legend.fancybox']=True
rcParams['xtick.labelsize']='medium'

# ---- DEFINITION OF FUNCTIONS

def compute_gamma_parameters (nzeros, npositive, \
    gamma_sum, gamma_ln_sum, weight):
    
    # --- gamma parameters estmation technique follows Wilks text.
    #     Statistical Methods in the Atmospheric Sciences (3rd ed.)
    #     Here adaptation using eqs. 4.40 and 4.41 and 4.42, 
    #     the Thom (1958) estimator for the shape parameter.
    #     As opposed to the routine compute_gamma_parameters_fclim
    #     below, this routine is used for the parameter estimation
    #     in the situation when the ensemble-mean forecast is 
    #     nonzero.

    if weight >= 100.:
        rntot = nzeros + npositive
        fraction_zero = nzeros / rntot
        xmean = (gamma_sum/weight)
        D = np.log(xmean) - (gamma_ln_sum/weight)
        rnumer2 = 1. + 4.*D/3.
        gamma_shape = (1. + np.sqrt(rnumer2)) / (4.*D)
        gamma_scale = xmean / gamma_shape
    else:
        xmean = -99.99
        fraction_zero = -99.99 # not enough data to estimate well
        gamma_shape = -99.99
        gamma_scale = -99.99
    return fraction_zero, gamma_shape, gamma_scale, xmean

def compute_gamma_parameters_fclim (nzeros, npositive, \
    gamma_sum, gamma_ln_sum):
    
    # --- gamma parameters estmation technique follows Wilks text.
    #     Statistical Methods in the Atmospheric Sciences (3rd ed.)
    #     Here using eqs. 4.40 and 4.41 and 4.42, the Thom (1958)
    #     estimatotr for the shape parameter.  Used in the situation
    #     when the ensemble-mean forecast is zero.  

    rntot = nzeros + npositive
    fraction_zero = nzeros / rntot
    xmean = gamma_sum / float(npositive)  
    D = np.log(xmean) - (gamma_ln_sum / float(npositive))
    if D < 0.01 : D = 0.01
    rnumer2 = 1. + 4.*D/3.
    gamma_shape = (1. + np.sqrt(rnumer2)) / (4.*D)
    gamma_scale = xmean / gamma_shape
    
    return fraction_zero, gamma_shape, gamma_scale

# --- COMMAND LINE INPUT

cmodel = sys.argv[1] # NCEP, ECMWF, or CMC
cleade = sys.argv[2] # ending of periods lead time in hours, e.g., '24', '108'
ileade = int(cleade)
date_forecast = sys.argv[3] # yyyymmddhh format, date of initial condition
cempirical = sys.argv[4]  # 1 for empirical, 0 for gamma distributions in estimating precip CDFs
date_end = dateshift(date_forecast,-ileade)
date_begin = dateshift(date_end, -61*24)

# I've hardcoded here to generate dates of interest to me.  You'd change
# to flexibly input the dates spanning the previous 60 days.

#date_begin = sys.argv[3]
#date_end = sys.argv[4]
#date_begin = '2016040100'  # replace with command line input when
#date_end = '2016070100'    # used at MDL

date_list = daterange(date_begin, date_end, 24) # makes a character list of dates 
ndates = len(date_list)

# ======================================================================
# Part 1:  reading in the daily data and summing up the information over
#          all case days
# ======================================================================

# ---- from a sample file, read in the array dimensions, precip values where 
#      gamma distribution parameters and fraction zero are calculated, and 
#      precipitation thresholds that are boundaries for light, 
#      medium, and heavy precipitation from netcdf file.  These
#      were previously generated in tally_gamma_statsfull_weighted.f90,
#      called by generate_gammadressing_stats_anymodel.f90

data_directory = '/Users/thamill/precip/ecmwf_data/'   
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

weights_gamma_total = np.zeros((3,3,nthreshes),dtype=np.float64)
gamma_sum_total = np.zeros((3,3,nthreshes),dtype=np.float64)
gamma_ln_sum_total = np.zeros((3,3,nthreshes),dtype=np.float64)
nzeros_total = np.zeros((3,3,nthreshes),dtype=np.float64)
npositive_total = np.zeros((3,3,nthreshes),dtype=np.float64)
nsamps_total = np.zeros((3,3,nthreshes),dtype=np.float64)
gamma_sum_ensmeanzero_fclim_total = \
    np.zeros((n_climocats),dtype=np.float64)
gamma_ln_sum_ensmeanzero_fclim_total = \
    np.zeros((n_climocats),dtype=np.float64)
closest_histogram_total = np.zeros((3,nmembersx25),dtype=np.float32)
nzeros_fclim_total = np.zeros((n_climocats),dtype=np.float64)
npositive_fclim_total = np.zeros((n_climocats),dtype=np.float64)
nsamps_fclim_total = np.zeros((n_climocats),dtype=np.float64)
exceed_yes_fclimPOP_total = np.zeros((nout_thresh, 
    n_climocats),dtype=np.float64)
exceed_no_fclimPOP_total = np.zeros((nout_thresh, \
    n_climocats),dtype=np.float64)

# ---- loop over dates, read in data, and add the latest day's data
#      to the overall summation arrays

for idate, date in zip(range(ndates), date_list):
    if cempirical == '1':
        infile = data_directory+'gamma_and_closest_hist_stats_' + \
            cmodel + '_' + date + '_fhour'+cleade+'.nc'
    else:
        infile = data_directory+'gamma_and_closest_hist_stats_' + \
            cmodel + '_' + date + '_fhour'+cleade+'_gammaqmap.nc'        
    print infile
    fexist  = os.path.exists(infile)
    if fexist:
        try: 
            
            #  --- read in 
            
            nc = Dataset(infile)               
            gamma_sum = nc.variables['gamma_sum'][:,:,:]
            gamma_ln_sum = nc.variables['gamma_ln_sum'][:,:,:]
            weights_gamma = nc.variables['weights'][:,:,:]
            nzeros = nc.variables['nzeros'][:,:,:]
            npositive = nc.variables['npositive'][:,:,:]
            gamma_sum_ensmeanzero_fclim = \
                nc.variables['gamma_sum_ensmeanzero_fclim'][:]
            gamma_ln_sum_ensmeanzero_fclim = \
                nc.variables['gamma_ln_sum_ensmeanzero_fclim'][:]
            nzeros_fclim = nc.variables['nzeros_fclim'][:]
            npositive_fclim = nc.variables['npositive_fclim'][:]
            closest_histogram = nc.variables['closest_histogram'][:,:]
            exceed_yes_fclimPOP = nc.variables['exceed_yes_fclimPOP'][:,:]
            exceed_no_fclimPOP = nc.variables['exceed_no_fclimPOP'][:,:]
            
            # ---- add to sum over many days
            
            gamma_sum_total = gamma_sum_total + gamma_sum
            gamma_ln_sum_total = gamma_ln_sum_total + gamma_ln_sum
            nzeros_total = nzeros_total + nzeros
            npositive_total = npositive_total + npositive
            nsamps_total = nsamps_total + npositive + nzeros
            weights_gamma_total = weights_gamma_total + weights_gamma
            
            closest_histogram_total = closest_histogram_total + \
                closest_histogram
            
            gamma_sum_ensmeanzero_fclim_total = \
                gamma_sum_ensmeanzero_fclim_total + \
                gamma_sum_ensmeanzero_fclim
            gamma_ln_sum_ensmeanzero_fclim_total = \
                gamma_ln_sum_ensmeanzero_fclim_total + \
                gamma_ln_sum_ensmeanzero_fclim
            nzeros_fclim_total = nzeros_fclim_total + nzeros_fclim
            npositive_fclim_total = npositive_fclim_total + npositive_fclim 
            nsamps_fclim_total = nsamps_fclim_total + \
                nzeros_fclim + npositive_fclim
            exceed_yes_fclimPOP_total = exceed_yes_fclimPOP_total + \
                exceed_yes_fclimPOP 
            exceed_no_fclimPOP_total = exceed_no_fclimPOP_total + \
                exceed_no_fclimPOP           
                
        except (IOError, ValueError, RuntimeError):
            print 'Error reading ', infile
            

# ======================================================================
# Part 2:  process the data to generate the dressing statistics
# ======================================================================

# ---- calculate the empirical gamma distribution shape, scale, and fraction zero
#      for selected values of precipitation amount, defined by gamma_threshold
#      previously read in.  These parameters also vary with the precipitation 
#      amount category (4 categories, with boundaries defined by thresh_light, 
#      thresh_mod, thresh_high) and 3 closest-histogram ranks (lowest, interior,
#      highest).
 
fraction_zero = -99.99 * np.ones((3,3,nthreshes),dtype=np.float32) # empirical
gamma_shape = -99.99 * np.ones((3,3,nthreshes),dtype=np.float32) 
gamma_scale = -99.99 * np.ones((3,3,nthreshes),dtype=np.float32)

ilowest  = np.zeros((3,3),dtype=np.int32)
ihighest = np.zeros((3,3),dtype=np.int32)

iptitles = ['Light','Moderate','Heavy']
ilomidhi_titles = ['LowestMbr','MiddleMbr','HighestMbr']
for ipamount in range(3):
    for ilomidhi in range(0,3):
                
        print 'processing ',iptitles[ipamount], ilomidhi_titles[ilomidhi]
        
        # ---- step 1: find the threshold amount that has the largest 
        #      number of samples
        
        imax = 0
        kmax = 0
        for ithresh in range(nthreshes):
            nsamps = nzeros_total[ilomidhi,ipamount,ithresh] + \
                npositive_total[ilomidhi,ipamount,ithresh] 
            if nsamps > 100:
                imax = ithresh
                kmax = int(nsamps)

        # ---- step 2: working downward from max, estimate the
        #      parameters if we have enough samples
        
        ienough = True # enough samples data to be trustworthy?
        for ithresh in range(nthreshes):
            
            # --- see if these indices are ones where there
            #     are not enough samples to reliably estimate parameter 
            #     values any more.
    
            fraction_zero[ilomidhi,ipamount,ithresh], \
                gamma_shape[ilomidhi,ipamount,ithresh], \
                gamma_scale[ilomidhi,ipamount,ithresh], xmean = \
                compute_gamma_parameters(\
                nzeros_total[ilomidhi,ipamount,ithresh], 
                npositive_total[ilomidhi,ipamount,ithresh], \
                gamma_sum_total[ilomidhi,ipamount,ithresh], \
                gamma_ln_sum_total[ilomidhi,ipamount,ithresh], \
                weights_gamma_total[ilomidhi,ipamount,ithresh])  
                
        # ---- for the intermediate members, the fitted Gamma distribution
        #      should be centered on the ensemble-mean amount, with the
        #      spread relatively tiny, though larger when we have fewer
        #      ensemble members.  We will re-set these in an ad-hoc 
        #      fashion by assuming we know the mean and spread of the
        #      sample data and then setting the gamma distribution parameters
        #      by the method of moments.
        
        if ilomidhi == 1: 
            for ithresh in range(nthreshes):
                xmean = gamma_threshes[ithresh]
                stddev = (float(nmembersx25)/180.) * (xmean/20.)
                gamma_shape[ilomidhi,ipamount,ithresh] = xmean**2/stddev**2
                gamma_scale[ilomidhi,ipamount,ithresh] = stddev**2 / xmean
    
        # ---- Determine plotting bounds with good data
        
        #print 'fraction_zero[ilomidhi,ipamount,:] = ', \
        #    fraction_zero[ilomidhi,ipamount,:]
        ithresh = 0
        while fraction_zero[ilomidhi,ipamount,ithresh] < 0. and ithresh < nthreshes-1:
            ithresh = ithresh+1
        if ithresh == nthreshes-1: 
            ilowest[ilomidhi,ipamount] = 0
        else:
            ilowest[ilomidhi,ipamount] = ithresh
            
        ithresh = nthreshes-1
        while fraction_zero[ilomidhi,ipamount,ithresh] < 0. and ithresh > 0:
            ithresh = ithresh-1
        if ithresh == 0: 
            ihighest[ilomidhi,ipamount] = 1
        else:
            ihighest[ilomidhi,ipamount] = ithresh+1   #???? 
        
        #print 'ilowest, ihighest = ', ilowest[ilomidhi,ipamount],ihighest[ilomidhi,ipamount]

# --- in situations where the ensemble mean was near zero, practical experience
#     has shown that the distribution of nonzero analyzed values depends on the climatology.
#     here we will estimate the gamma shape and scale parameters in this situation
#     as well as the fraction of samples with zero analyzed precip

gamma_shape_meanzero_fclim = np.zeros((n_climocats), dtype=np.float32)
gamma_scale_meanzero_fclim = np.zeros((n_climocats), dtype=np.float32)
fraction_zero_meanzero_fclim = np.zeros((n_climocats), dtype=np.float32)
fz_meanzero_fclim_empirical = np.zeros((nout_thresh,n_climocats), dtype=np.float32)

for i in range(n_climocats):
    fraction_zero_meanzero_fclim[i], gamma_shape_meanzero_fclim[i], \
        gamma_scale_meanzero_fclim[i] = \
        compute_gamma_parameters_fclim(nzeros_fclim_total[i], \
        npositive_fclim_total[i], gamma_sum_ensmeanzero_fclim_total[i], \
        gamma_ln_sum_ensmeanzero_fclim_total[i])
    fz_meanzero_fclim_empirical[:,i] = exceed_no_fclimPOP_total[:,i] / \
        (exceed_yes_fclimPOP_total[:,i] + exceed_no_fclimPOP_total[:,i])

    
# ---- now we turn our attention to generating the closest histogram statistics
#      for the sorted ensemble.  Normalize the closest histograms so they sum 
#      to 1.0.   Then Savitzky-Golay smooth the interior values. 

for icat in range(3):
    closest_histogram_total[icat,:] = closest_histogram_total[icat,:] / \
        np.sum(closest_histogram_total[icat,:])

for i in range(180):
    print i,closest_histogram_total[0,i], closest_histogram_total[1,i],\
        closest_histogram_total[2,i]
        
closest_histogram_savgol = np.copy(closest_histogram_total)
for i in range(3):
    csavgol = np.copy(closest_histogram_savgol[i,:])
    work = np.copy(csavgol[4:-4])
    csavgol[4:-4] = signal.savgol_filter(work, 9, 2, mode='interp')
    closest_histogram_savgol[i,:] = csavgol[:] / np.sum(csavgol)

# ======================================================================
# ---- PART 3:  Save the gamma distribution 
#               information to netcdf files
# ======================================================================

# ---- (3a) define netCDF file name for gamma distribution parameters and open
#      the file
if cempirical == '1':
    outfile_nc = data_directory+cmodel+'/gamma_fraction_zero_dressing_'+\
        cmodel+'_date='+date_forecast+'_lead='+cleade+'.nc'
else:
    outfile_nc = data_directory+cmodel+'/gamma_fraction_zero_dressing_'+\
        cmodel+'_date='+date_forecast+'_lead='+cleade+'_gammaqmap.nc'
print 'writing netCDF dressing statistics to ',outfile_nc
rootgrp = Dataset(outfile_nc,'w',format='NETCDF4_CLASSIC')

# ---- (3b) declare array dimensions, indices for arrays

npv = len(gamma_threshes)
npvals = rootgrp.createDimension('npvals',npv)
npvals_indices = rootgrp.createVariable('npvals_indices','i4',('npvals',))
npvals_indices.long_name = "index for gamma_threshes array, "+\
    "the precip amounts where gamma distribution parameter information is stored"
npvals_indices.units = "n/a"

nclim = len(climo_pop_thresholds)
nclim_vals = rootgrp.createDimension('nclim_vals',nclim)
nclim_indices = rootgrp.createVariable('nclim_indices','i4',('nclim_vals',))
nclim_indices.long_name = \
    'index for climo_pop_thresholds, the climatological POP values that '+\
    'define boundaries for arrays of fraction zero and gamma parameters '+\
    'when forecast mean ~ = 0'
nclim_indices.units = "n/a"

nclim_p1 = nclim+1
nclim_p1_vals = rootgrp.createDimension('nclim_p1_vals',nclim_p1)
nclim_p1_indices = \
    rootgrp.createVariable('nclim_p1_indices','i4',('nclim_p1_vals',))
nclim_p1_indices.long_name = 'index into arrays gamma_shape_fclimpop, '+\
    'gamma_scale_fclimpop, fraction_zeros_fclimpop'
nclim_p1_indices.units = "n/a"

nlo_int_hi = 3
nlo_int_hi_vals = rootgrp.createDimension('nlo_int_hi_vals',nlo_int_hi)
nlo_int_hi_indices = \
    rootgrp.createVariable('nlo_int_hi_indices','i4',('nlo_int_hi_vals',))
nlo_int_hi_indices.long_name = \
    "indices of for lower/intermediate/highest sorted members"
nlo_int_hi_indices.units = "n/a"

npcats = 3 
npcatvals = rootgrp.createDimension('npcatvals',npcats)
npcat_indices = rootgrp.createVariable('npcat_indices','i4',('npcatvals',))
npcat_indices.long_name = "indices for discretization of arrays by "+\
    "ensemble-mean precipitation forecast amounts"
npcat_indices.units = "n/a"

# ---- (3c) initialize the arrays that will hold the data 

precip_values = rootgrp.createVariable('precip_values','f4',('npvals',), \
    zlib=True,least_significant_digit=3)
precip_values.units = "mm"
precip_values.long_name = \
    "precipitation values where fraction zero, "+\
    "Gamma shape, Gamma scale are estimated"
precip_values.valid_range = [0.,500.]
precip_values.missing_value = -99.99

fraction_zeros = rootgrp.createVariable('fraction_zeros','f4',\
    ('nlo_int_hi_vals','npcatvals','npvals',), \
    zlib=True,least_significant_digit=4)
fraction_zeros.units = "n/a"
fraction_zeros.long_name = \
    "Fraction of samples with analyzed precipitation = 0"
fraction_zeros.valid_range = [0.,1.]
fraction_zeros.missing_value = -99.99

fraction_zeros_fclimpop = rootgrp.createVariable('fraction_zeros_fclimpop',\
    'f4',('nclim_p1_vals',), zlib=True,least_significant_digit=4)
fraction_zeros_fclimpop.units = "n/a"
fraction_zeros_fclimpop.long_name = \
    "Fraction of samples where analyzed precipitation = 0 "+\
    "when mean forecast ~=0. (as function of climatological POP)"
fraction_zeros_fclimpop.valid_range = [0.,1.]
fraction_zeros_fclimpop.missing_value = -99.99

climo_pop_thresholds = rootgrp.createVariable('climo_pop_thresholds',\
    'f4',('nclim_vals',), zlib=True,least_significant_digit=4)
climo_pop_thresholds.units = "n/a"
climo_pop_thresholds.long_name = "Threshold of POP climatology for "+\
    "defining fraction zero and gamma distribution "+\
    "parameters when ens-mean forecast precip ~= 0"
climo_pop_thresholds.valid_range = [0.,1.]
climo_pop_thresholds.missing_value = -99.99

gamma_shapes = rootgrp.createVariable('gamma_shapes','f4',\
    ('nlo_int_hi_vals','npcatvals','npvals',), \
    zlib=True,least_significant_digit=3)
gamma_shapes.units = "n/a"
gamma_shapes.long_name = "Gamma distribution shape parameter (alpha)"
gamma_shapes.valid_range = [0.,1200.]
gamma_shapes.missing_value = -99.99

gamma_scales = rootgrp.createVariable('gamma_scales','f4',\
    ('nlo_int_hi_vals','npcatvals','npvals',), \
    zlib=True,least_significant_digit=3)
gamma_scales.units = "n/a"
gamma_scales.long_name =  "Gamma distribution scale parameter (beta)"
gamma_scales.valid_range = [0.,1200.]
gamma_scales.missing_value = -99.99

nsamps = rootgrp.createVariable('nsamps','f4',('nlo_int_hi_vals',\
    'npcatvals','npvals',), zlib=True,least_significant_digit=1)
nsamps.units = "n/a"
nsamps.long_name =  "Number of samples"
nsamps.valid_range = [0.,10000000000.]
nsamps.missing_value = -99.99

gamma_shape_fclimpop = rootgrp.createVariable('gamma_shape_fclimpop',\
    'f4',('nclim_p1_vals',), zlib=True,least_significant_digit=4)
gamma_shape_fclimpop.units = "n/a"
gamma_shape_fclimpop.long_name = \
    "Gamma shape parameter for nonzero analyzed precipitation" + \
    " where mean precip ~=0 as function of climatological POP"
gamma_shape_fclimpop.valid_range = [0.,100.]
gamma_shape_fclimpop.missing_value = -99.99

gamma_scale_fclimpop = rootgrp.createVariable('gamma_scale_fclimpop',\
    'f4',('nclim_p1_vals',), zlib=True,least_significant_digit=4)
gamma_scale_fclimpop.units = "n/a"
gamma_scale_fclimpop.long_name = \
    "Gamma scale parameter for nonzero analyzed precipitation " + \
    " where mean precip ~=0 as function of climatological POP"
gamma_scale_fclimpop.valid_range = [0.,100.]
gamma_scale_fclimpop.missing_value = -99.99

nsamps_fclimpop = rootgrp.createVariable('nsamps_fclimpop','f4',\
    ('nclim_p1_vals',), zlib=True,least_significant_digit=1)
nsamps_fclimpop.units = "n/a"
nsamps_fclimpop.long_name = \
    "Number of samples as function of climatological POP"
nsamps_fclimpop.valid_range = [0.,10000000000000.]
nsamps_fclimpop.missing_value = -99.99

precip_thresholds =  rootgrp.createVariable('precip_thresholds',\
    'f4',('npcatvals',), zlib=True,least_significant_digit=4)
precip_thresholds.units = "n/a"
precip_thresholds.long_name = \
    "Precipitation amounts that define borders "+\
    "between no, light, mod, heavy precip"
precip_thresholds.valid_range = [0.,100.]
precip_thresholds.missing_value = -99.99

# --- (3d) write calculated values for netcdf file indices and variables

nclim_indices[:] = range(nclim)
npvals_indices = range(npv)
nclim_p1_indices = range(nclim_p1)
nlo_int_hi_indices = range(nlo_int_hi) 
npcat_indices = range(npcats)

precip_values[:] = gamma_threshes[:]
fraction_zeros[:,:,:] = fraction_zero[:,:,:]
gamma_shapes[:,:,:] = gamma_shape[:,:,:]
gamma_scales[:,:,:] = gamma_scale[:,:,:]
nsamps[:,:,:] = nsamps_total[:,:,:]
fraction_zeros_fclimpop[:] = fraction_zero_meanzero_fclim[:]
gamma_shape_fclimpop[:] = gamma_shape_meanzero_fclim[:]
gamma_scale_fclimpop[:] = gamma_scale_meanzero_fclim[:]
nsamps_fclimpop[:] = nsamps_fclim_total[:]

climo_pop_thresholds[:] = climo_pop_thresholds_save[:]
precip_thresholds[:] = np.array([np.float(thresh_light), \
    np.float(thresh_mod), np.float(thresh_high)])

# ---- (3e) more metadata

rootgrp.stream = "s4" # ????
rootgrp.title = "Best-member dressing statistics for National Blend "+\
    cmodel+" multi-model ensemble"
rootgrp.Conventions = "CF-1.0"  # ????
rootgrp.history = "Created ~ Sep 2017 by Tom Hamill, tom.hamill@noaa.gov"
rootgrp.institution = "NOAA/ESRL/PSD"
rootgrp.platform = "Model"
rootgrp.references = "None"

# ---- (3f) close the file

rootgrp.close()


# ======================================================================
# ---- PART 4:  Save the closest-member histogram info to netcdf files
# ======================================================================

# ---- (4a) define netCDF file name for gamma distribution parameters and open
#      the file

if cempirical == '1':
    outfile_nc = data_directory+cmodel+'/closest_histogram_'+\
        cmodel+'_date='+date_forecast+'_lead='+cleade+'.nc'
else:
    outfile_nc = data_directory+cmodel+'/closest_histogram_'+\
        cmodel+'_date='+date_forecast+'_lead='+cleade+'_gammaqmap.nc'
        
print 'writing netCDF closest histogram data to ',outfile_nc
rootgrp = Dataset(outfile_nc,'w',format='NETCDF4_CLASSIC')

# ---- (4b) declare array dimensions, indices for arrays

npcats = 3 
npcatvals = rootgrp.createDimension('npcatvals',npcats)
npcat_indices = rootgrp.createVariable('npcat_indices','i4',('npcatvals',))
npcat_indices.long_name = "indices for discretization of arrays by '+\
    'ensemble-mean precipitation forecast amounts"
npcat_indices.units = "n/a"

nmembervals = rootgrp.createDimension('nmembervals',nmembersx25)
nmember_indices = rootgrp.createVariable('nmember_indices','i4',('nmembervals',))
nmember_indices.long_name = "sorted ensemble member number"
nmember_indices.units = "n/a"

# ---- (4c) declare arrays that will hold the dressing data

precip_thresholds =  rootgrp.createVariable('precip_thresholds','f4',('npcatvals',), \
    zlib=True,least_significant_digit=4)
precip_thresholds.units = "n/a"
precip_thresholds.long_name = \
    "Precipitation amounts that define borders between no, light, mod, heavy precip"
precip_thresholds.valid_range = [0.,100.]
precip_thresholds.missing_value = -99.99

closest_histogram = rootgrp.createVariable('closest_histogram','f4',('npcatvals','nmembervals',), \
    zlib=True,least_significant_digit=6)
closest_histogram.units = "n/a"
closest_histogram.long_name = \
    "histogram of how often this sorted member is the closet to analyzed value"
closest_histogram.valid_range = [0.,1.]
closest_histogram.missing_value = -99.99

# --- (4d) write calculated values for netcdf file indices and variables

npcat_indices = range(npcats)
precip_thresholds[:] = np.array([np.float(thresh_light), np.float(thresh_mod), np.float(thresh_high)])
nmember_indices[:] = range(nmembersx25)
closest_histogram[:,:] = closest_histogram_savgol[:,:]

# ---- (4e) more metadata

rootgrp.stream = "s4" # ????
rootgrp.title = "histograms of how often this sorted member is the closet to analyzed value for "+cmodel
rootgrp.Conventions = "CF-1.0"  # ????
rootgrp.history = "Created ~ Sep 2017 by Tom Hamill, tom.hamill@noaa.gov"
rootgrp.institution = "NOAA/ESRL/PSD"
rootgrp.platform = "Model"
rootgrp.references = "None"

# ---- (5f) close the file, and done

rootgrp.close()

sys.exit()


# ======================================================================
# ---- PART 5:  now make some diagnostic plots, save these to pdf file
# ======================================================================

#  first plot for the fraction zero, the gamma shape and scale for the empirical
#  for the various categories of ensemble-mean precip, and then as a function
#  of the best-member amount.

cmember = ['lowest sorted member', 'intermediate sorted members','highest sorted member']
camount = ['light mean precipitation','moderate mean precipitation',\
    'heavy mean precipitation']
cmember_short = ['lowest', 'intermed','highest']
camount_short = ['light_precip','mod_precip','heavy_precip']
with PdfPages(cmodel+'_'+cleade+'h_dressing_diagnostics.pdf') as pdf:
    for ilomidhi in range(0,3,2):   # lowest mbr, highest sorted
        for iamount in range(3):   # light, moderate, heavy
        
            print 'processing ilomidhi, amount', ilomidhi, camount_short[iamount]
            print 'ilowest, ihighest = ', ilowest[ilomidhi,iamount], \
                ihighest[ilomidhi,iamount]
            xlimits = [ gamma_threshes[ilowest[ilomidhi,iamount]], \
                gamma_threshes[ihighest[ilomidhi,iamount]]]
            if xlimits[0] == xlimits[1]: xlimits = [0,5]   
    
            # ---- make plots of the fraction zero and gamma shape and scale parameters.
    
            fig = plt.figure(figsize=(6.5,9.))
            fig.suptitle(cmodel+' +'+cleade+' h dressing statistics for\n'+\
            cmember[ilomidhi]+' and '+camount[iamount], fontsize=15)
                
            for ipanel in range(3):
                if ipanel == 0:
                    ctitle = '(a) Fraction of analyzed samples with zero precipitation'
                    ytitle = 'Fraction zero'
                    plotdata = fraction_zero[ilomidhi,iamount,\
                        ilowest[ilomidhi,iamount]:ihighest[ilomidhi,iamount]]
                    ylimits = [0.,1.]
                    axlocn = [0.15,0.68,0.8,0.21]
                elif ipanel == 1:
                    ctitle = r'(b) Gamma distribution shape parameter $\alpha$'
                    ytitle = r'Shape parameter $\alpha$'
                    plotdata = gamma_shape[ilomidhi,iamount,\
                        ilowest[ilomidhi,iamount]:ihighest[ilomidhi,iamount]]
                    print 'gamma_shape',plotdata
                    ylimits = [np.min(plotdata), np.max(plotdata)]
                    axlocn = [0.15,0.36,0.8,0.21]
                else:
                    ctitle = r'(c) Gamma distribution scale parameter $\alpha$'
                    ytitle = r'Scale parameter $\beta$'
                    plotdata = gamma_scale[ilomidhi,iamount,\
                        ilowest[ilomidhi,iamount]:ihighest[ilomidhi,iamount]]
                    print 'gamma_scale ', plotdata
                    ylimits = [np.min(plotdata), np.max(plotdata)]
                    axlocn = [0.15,0.06,0.8,0.21]

                a1 = fig.add_axes(axlocn)
                a1.set_title(ctitle,fontsize=13)
                a1.plot(gamma_threshes[ilowest[ilomidhi,iamount]:\
                    ihighest[ilomidhi,iamount]],plotdata,'r-',lw=3)
                a1.set_xlabel('Best-member precipitation amount (mm)',fontsize=11)      
                a1.set_ylabel(ytitle,fontsize=11)
                a1.set_xlim(xlimits)
                a1.set_ylim(ylimits)
                a1.grid(color='Gray',lw=0.2,linestyle='--')
    
            pdf.savefig()
            plt.close()

    # --- now make plots of sample dressing distributions
 
    for ilomidhi in range(0,3,2):   # lowest mbr, highest 
        clomidhi = ['lowest','intermediate','highest']
        fig = plt.figure(figsize=(6.5,9.))
        fig.suptitle('Dressing distribution for '+clomidhi[ilomidhi]+\
            ' sorted member,\nmodel = '+cmodel+', lead = '+cleade+' h', fontsize=16)
        for iamount in range(3):   # light, moderate, heavy
        
            ytitle = 'Probability density'
            xtitle = 'Precipitation amount (mm)'
            if ilomidhi == 0 and iamount == 0:  # 0.1 to 2.0
                plot_at_these_gamma_threshes = [0.0, 0.2, 0.4, 0.6]
                index_gamma = [0,3,5,7]
                xvalues = np.arange(0.0,2.0,0.01)
                ctitle = r'(a) 0.01 $\leq$ mean precip. < 2 mm'
                axlocn = [0.15,0.67,0.8,0.22]
            elif ilomidhi == 0 and iamount == 1: # 2.0 to 10
                plot_at_these_gamma_threshes = [1.0, 2.0, 4.0, 6.0]
                index_gamma =[11, 16, 22, 26]
                xvalues = np.arange(0.0,10.0,0.02)
                ctitle = r'(b) 2 $\leq$ mean precip. < 10 mm'
                axlocn = [0.15,0.36,0.8,0.22]
            elif ilomidhi == 0 and iamount == 2:
                plot_at_these_gamma_threshes = [4.0, 8.0, 12.0, 16.0]
                index_gamma = [22,30,35,39]
                xvalues = np.arange(0.0,30.0,0.1)
                ctitle = r'(c) Mean precip. $\geq$ 10 mm  '
                axlocn = [0.15,0.05,0.8,0.22]
                
            if ilomidhi == 1 and iamount == 0:
                plot_at_these_gamma_threshes = [0.1, 0.5, 1.0, 2.0]
                index_gamma = [2,6,11,16]
                xvalues = np.arange(0.0,4.0,0.02)
                ctitle = r'(a) 0.01 $\leq$ mean precip. < 2 mm'
                axlocn = [0.15,0.67,0.8,0.22]
            elif ilomidhi == 1 and iamount == 1: # 2.0 to 10
                plot_at_these_gamma_threshes = [2.0, 4.0, 7.0, 10.0]
                index_gamma = [16, 22, 28, 33]
                xvalues = np.arange(0.0,25.0,0.05)
                ctitle = r'(b) 2 $\leq$ mean precip. < 10 mm'
                axlocn = [0.15,0.36,0.8,0.22]
            elif ilomidhi == 1 and iamount == 2:
                plot_at_these_gamma_threshes = [10.0, 20.0, 30.0, 40.0]
                index_gamma = [33,41,45,48]
                xvalues = np.arange(0.0,100.0,0.1)
                ctitle = r'(c) Mean precip. $\geq$ 10 mm  '
                axlocn = [0.15,0.05,0.8,0.22]               
                
            if ilomidhi == 2 and iamount == 0:
                plot_at_these_gamma_threshes = [1.0, 2.0, 4.0, 7.0]
                index_gamma = [11,16,22,28]
                xvalues = np.arange(0.0,25.0,0.2)
                ctitle = r'(a) 0.01 $\leq$ mean precip. < 2 mm'
                axlocn = [0.15,0.67,0.8,0.22]
            elif ilomidhi == 2 and iamount == 1: # 2.0 to 10
                plot_at_these_gamma_threshes = [5.0, 10.0, 20.0, 30.0]
                index_gamma = [24,33,41,45]
                xvalues = np.arange(0.0,75.0,0.1)
                ctitle = r'(b) 2 $\leq$ mean precip. < 10 mm'
                axlocn = [0.15,0.36,0.8,0.22]
            elif ilomidhi == 2 and iamount == 2: # > 10
                plot_at_these_gamma_threshes = [15.0, 30.0, 50.0, 75.0]
                index_gamma = [38,45,50,55]
                xvalues = np.arange(0.0,150.0,0.1)
                ctitle = r'(c) Mean precip. $\geq$ 10 mm  '
                axlocn = [0.15,0.05,0.8,0.22]                
            
            xlimits = [0,xvalues[-1]]                
            a1 = fig.add_axes(axlocn)
            a1.set_title(ctitle,fontsize=13)
            a1.set_ylabel(ytitle,fontsize=11)
            a1.set_xlabel(xtitle,fontsize=11)
            a1.set_xlim(xlimits)
            
            nplots = len(plot_at_these_gamma_threshes)
            colors = ['Red','RoyalBlue','LimeGreen','Purple']
            for i in range(0,nplots):
                plotdata = gamma.pdf(xvalues, gamma_shape[ilomidhi,iamount,index_gamma[i]], \
                    0., gamma_scale[ilomidhi,iamount,index_gamma[i]])
                ctext1 = "{:.2f}".format(gamma_threshes[index_gamma[i]])
                a1.plot(xvalues,plotdata,'-',color=colors[i],lw=2,\
                    label='Forecast = '+ctext1)
            a1.grid(color='Gray',lw=0.2,linestyle='--')
            a1.legend(loc=0)
            
        pdf.savefig()
        plt.close()
        
    # ---- now make closest histogram plots

    fig = plt.figure(figsize=(6.,9.))
    plt.suptitle('Closest-member histograms for lead = +'+cleade+' h,\nmodel = '+cmodel,fontsize=16)
    for i in range(3):
        if i == 0:
            ctitle = r'(a) 0.01 mm $\leq$ mean precip. < 2 mm '
            plotdata1 = closest_histogram_total[0,:]/2
            plotdata2 = closest_histogram_savgol[0,:]
            axlocn = [0.15,0.68,0.8,0.21]
        elif i == 2:
            ctitle = r'(b) 2 mm $\leq$ mean precip. < 10 mm '
            plotdata1 = closest_histogram_total[1,:]/2
            plotdata2 = closest_histogram_savgol[1,:]
            axlocn = [0.15,0.36,0.8,0.21]
        else:
            ctitle = r'(c) mean precip. $\geq$ 10 mm'
            plotdata1 = closest_histogram_total[2,:]/2
            plotdata2 = closest_histogram_savgol[2,:]
            axlocn = [0.15,0.06,0.8,0.21]

        a1 = fig.add_axes(axlocn)
        a1.set_title(ctitle,fontsize=14)
        a1.set_xlabel('Rank of sorted closest member to analyzed',fontsize=12)
        a1.set_ylabel('Fraction',fontsize=12)
        a1.set_xlim(0,nmembersx25+2)
        a1.set_ylim(0.0001,0.5)
        a1.set_yscale('log')
        a1.grid(color='Gray',lw=0.2,linestyle='--')
        print 'plotdata1 = ',plotdata1
        print 'plotdata2 = ',plotdata2
        a1.plot(range(1,nmembersx25+1), plotdata1,'-',\
            color='Red',label='Raw / 2',linewidth=2)
        a1.plot(range(1,nmembersx25+1), plotdata2,'-',\
            color='RoyalBlue',label='Smoothed',linewidth=2)
        ctext_low = "{:.3f}".format(plotdata2[0])
        ctext_high = "{:.3f}".format(plotdata2[-1])
        a1.text(nmembersx25*0.01,plotdata2[0],ctext_low,color='RoyalBlue')
        a1.text(nmembersx25*0.88,plotdata2[-1],ctext_high,color='RoyalBlue')
        if i == 0: a1.legend(loc=9)

    pdf.savefig()
    plt.close()
        
    # ---- make plots of the fraction zero, gamma shape, and gamma 
    #      scale as f(POP climatology) when mean precip is near zero.    
        
    fig = plt.figure(figsize=(6.5,9.))
    fig.suptitle(cmodel+' +'+cleade+' h dressing statistics as f(climatological POP)\n'+\
        ' when mean forecast ~= 0', fontsize=15)
      
    probclim_to_plot = np.zeros((n_climocats),dtype=np.float32)
    for iclim in range(n_climocats):
        if iclim == 0:
            plow = 0.
            phigh = climo_pop_thresholds[0]
        elif iclim == n_climocats-1:
            plow = climo_pop_thresholds[iclim-1]
            phigh = 2.*climo_pop_thresholds[iclim-1] - climo_pop_thresholds[iclim-2] 
        else:
            plow = climo_pop_thresholds[iclim-1]
            phigh = climo_pop_thresholds[iclim]
        probclim_to_plot[iclim] = (plow + phigh) / 2.
                
    for ipanel in range(3):
        if ipanel == 0:
            ctitle = '(a) Fraction of analyzed samples with zero precipitation'
            ytitle = 'Fraction zero'
            plotdata = fraction_zero_meanzero_fclim
            ylimits = [0.,1.]
            axlocn = [0.15,0.68,0.8,0.21]
        elif ipanel == 1:
            ctitle = r'(b) Gamma distribution shape parameter $\alpha$'
            ytitle = r'Shape parameter $\alpha$'
            plotdata = gamma_shape_meanzero_fclim
            print 'gamma_shape',plotdata
            ylimits = [np.min(plotdata), np.max(plotdata)]
            axlocn = [0.15,0.36,0.8,0.21]
        else:
            ctitle = r'(c) Gamma distribution scale parameter $\alpha$'
            ytitle = r'Scale parameter $\beta$'
            plotdata = gamma_scale_meanzero_fclim
            print 'gamma_scale ', plotdata
            ylimits = [np.min(plotdata), np.max(plotdata)]
            axlocn = [0.15,0.06,0.8,0.21]

        a1 = fig.add_axes(axlocn)
        a1.set_title(ctitle,fontsize=13)
        #a1.semilogy(gamma_threshes[ilowest[ilomidhi,ipamount]:\
        #    ihighest[ilomidhi,ipamount]],plotdata,'r-',lw=3)
        a1.plot(probclim_to_plot,plotdata,'r-',lw=3)
        a1.set_xlabel('Climatological probability of non-zero precipitation',fontsize=11)      
        a1.set_ylabel(ytitle,fontsize=11)
        a1.set_xlim([0,probclim_to_plot[iclim]+0.01])
        a1.set_ylim(ylimits)
        a1.grid(color='Gray',lw=0.2,linestyle='--')
    
    pdf.savefig()
    plt.close()

    # --- now make plots of the probability distributions when mean forecast is zero
    #     as function of the climatological POP
 
    fig = plt.figure(figsize=(6.5,4.5))
    ctitle = 'Dressing distributions for model = '+cmodel+\
        ', lead = '+cleade+' h\nwhen mean forecast is ~= 0.0'
    ytitle = 'Non-exceedance probability'
    xtitle = 'Precipitation amount (mm)'
    axlocn = [0.15,0.15,0.8,0.73]
    a1 = fig.add_axes(axlocn)
    xlimits = [0.,2.]
    xvalues = np.arange(0.0,10.0,0.01)                
    a1.set_title(ctitle,fontsize=13)
    a1.set_ylabel(ytitle,fontsize=11)
    a1.set_xlabel(xtitle,fontsize=11)
    a1.set_xlim(xlimits)
    colors = ['Red','RoyalBlue','LimeGreen','Purple','Black','Orange','Gray','DarkRed']
    for i in range(n_climocats):  
        plotdata = gamma.cdf(xvalues, gamma_shape_meanzero_fclim[i], \
            0., gamma_scale_meanzero_fclim[i])
        ctext1 = "Climatological POP ~ {:.2f}".format(probclim_to_plot[i])
        a1.plot(xvalues,plotdata,'-',color=colors[i],lw=2,\
            label='Forecast = '+ctext1) 
    a1.grid(color='Gray',lw=0.2,linestyle='--')
    a1.legend(loc=0)
            
    pdf.savefig()
    plt.close() 
    
    # --- now make plots of the probability distributions when mean forecast is zero
    #     as function of the climatological POP, log scaling.
 
    fig = plt.figure(figsize=(6.5,4.5))
    ctitle = 'Dressing distributions for model = '+cmodel+\
        ', lead = '+cleade+' h\nwhen mean forecast is ~= 0.0'
    ytitle = 'Probability density'
    xtitle = 'Precipitation amount (mm)'
    axlocn = [0.15,0.15,0.8,0.73]
    a1 = fig.add_axes(axlocn)
    xlimits = [0.0001,10.]
    xvalues = np.arange(0.0,10.0,0.01)                
    a1.set_title(ctitle,fontsize=13)
    a1.set_ylabel(ytitle,fontsize=11)
    a1.set_xlabel(xtitle,fontsize=11)
    a1.set_xlim(xlimits)
    colors = ['Red','RoyalBlue','LimeGreen','Purple','Black','Orange','Gray','DarkRed']
    for i in range(n_climocats):  
        plotdata = gamma.pdf(xvalues, gamma_shape_meanzero_fclim[i], \
            0., gamma_scale_meanzero_fclim[i])
        ctext1 = "Climatological POP ~ {:.2f}".format(probclim_to_plot[i])
        a1.semilogy(xvalues,plotdata,'-',color=colors[i],lw=2,\
            label='Forecast = '+ctext1)
    a1.grid(color='Gray',lw=0.2,linestyle='--')
    a1.legend(loc=0)
            
    pdf.savefig()
    plt.close()       
    
    # --- now plot the fraction zero when mean forecast is zero
    #     as function of the climatological POP
 
    fig = plt.figure(figsize=(6.5,4.5))
    ctitle = 'Fraction zero for model = '+cmodel+\
        ', lead = '+cleade+' h\nwhen mean forecast is ~= 0.0'
    ytitle = 'Fraction zero'
    xtitle = 'Climatological probability of precipitation > 0.254 mm / 12 h'
    axlocn = [0.15,0.15,0.8,0.73]
    a1 = fig.add_axes(axlocn)
    xlimits = [0.0, 0.25]
    xvalues = np.arange(0.0,10.0,0.01)                
    a1.set_title(ctitle,fontsize=13)
    a1.set_ylabel(ytitle,fontsize=11)
    a1.set_xlabel(xtitle,fontsize=11)
    a1.set_xlim(xlimits)
    plotdata = fraction_zero_meanzero_fclim
    xvalues = [0.015, 0.045, 0.075, 0.105, 0.135, 0.165, 0.195, 0.225]
    a1.plot(xvalues,plotdata,'-',color='Red',lw=2)

    
    a1.grid(color='Gray',lw=0.2,linestyle='--')   
    pdf.savefig()
    plt.close()        
            
print 'Plotting to '+cmodel+'_'+cleade+'h_dressing_diagnostics.pdf done.'









