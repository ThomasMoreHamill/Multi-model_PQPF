from netCDF4 import Dataset
import matplotlib
matplotlib.use('TKAgg')
import matplotlib.pyplot as plt
from matplotlib import rcParams
from mpl_toolkits.basemap import Basemap, addcyclic
#import scipy.ndimage
import numpy as np
from numpy import ma
import math
import os, sys
#import cm as cm2
#from read_ppn_analog_locns_ccpa8 import read_ppn_analog_locns_ccpa8
from read_ppn_analog_locns_ccpa9 import read_ppn_analog_locns_ccpa9
from matplotlib import rcParams


cmonth = sys.argv[1] # Jan, Feb, etc
cthresh = sys.argv[2] # threshold for penalty function
rthresh = float(cthresh)

# ---- read in analog locations

nxa = 515
nya = 262
nxa_small = 464
nya_small = 224
nsupp = 100

xlocationa = np.zeros((nxa,nya,nsupp),dtype=np.int32)
ylocationa = np.zeros((nxa,nya,nsupp),dtype=np.int32)
conusmaska = np.zeros((nxa,nya),dtype=np.int16)
practical_mask = np.zeros((nxa,nya),dtype=np.int16)
penalty = np.zeros((nxa,nya,nsupp),dtype=np.float32)
penalty_pa = np.zeros((nxa,nya,nsupp),dtype=np.float32)
penalty_ter = np.zeros((nxa,nya,nsupp),dtype=np.float32)
penalty_facet = np.zeros((nxa,nya,nsupp),dtype=np.float32)
penalty_dist = np.zeros((nxa,nya,nsupp),dtype=np.float32)
rlona = np.zeros((nxa,nya),dtype=np.float)
rlata = np.zeros((nxa,nya),dtype=np.float)
#xlocationa, ylocationa, practical_mask, conusmaska, rlona, rlata, \
#    penalty, penalty_pa, penalty_ter, penalty_facet, penalty_dist = \
#    read_ppn_analog_locns_ccpa8 (cmonth, nxa, nya, nsupp)
xlocationa, ylocationa, conusmaska, rlona, rlata, \
    penalty, penalty_pa, penalty_ter, penalty_facet, penalty_dist = \
    read_ppn_analog_locns_ccpa9 (cmonth, nxa, nya, nsupp)

zeros = np.zeros((nxa,nya), dtype=np.int16)
conusmaska = np.where(conusmaska < 0, zeros, conusmaska)

numsupp_t = np.zeros((nya,nxa), dtype=np.int16)
penalties = np.zeros((nsupp),dtype=np.float32)
for ixa in range(nxa):
    for jya in range(nya):
        if conusmaska[ixa,jya] == 1:
            penalties[:] = penalty[ixa,jya,:]
            k=1
            while penalties[k-1] < rthresh and k < nsupp:
                k=k+1
            numsupp_t[jya,ixa]=np.max([k,50])
        else:
            numsupp_t[jya,ixa]=0

#numsupp_m = ma.masked_less(numsupp, 0)

#  ---- prep arrays for dumping to netCDF by switching order.

xlocationa_t = np.transpose(xlocationa)
ylocationa_t = np.transpose(ylocationa)
rlatsa_t     = np.transpose(rlata)
rlonsa_t     = np.transpose(rlona)
conusmaska_t = np.transpose(conusmaska)
#practical_mask_t = np.transpose(practical_mask)

print 'xlocationa_t[0,nya/2,:] = ', xlocationa_t[0,nya/2,:]
print 'ylocationa_t[0,nya/2,:] = ', ylocationa_t[0,nya/2,:]
print 'rlatsa_t[nya/2,:] = ', rlatsa_t[nya/2,:] 
print 'rlonsa_t[nya/2,:] = ', rlonsa_t[nya/2,:] 
print 'conusmaska_t[nya/2,:] = ',conusmaska_t[nya/2,:] 
#print 'practical_mask_t[nya/2,:] = ',practical_mask_t[nya/2,:] 
print 'numsupp_t[nya/2-1,:] = ',numsupp_t[nya/2-1,:]

#sys.exit()

# ---- set up output netCDF file metadata particulars

outfile_nc = 'supplemental_locations_eighth_degree_'+cmonth+'_v9.nc'
print 'writing netCDF precipitation analysis data to ',outfile_nc
rootgrp = Dataset(outfile_nc,'w',format='NETCDF4_CLASSIC')

xsupp = rootgrp.createDimension('xsupp',nsupp)
xvsupp = rootgrp.createVariable('xsupp','i4',('xsupp',))
xvsupp.long_name = "supplemental location number"
xvsupp.units = "n/a"

xa = rootgrp.createDimension('xa',nxa)
xva = rootgrp.createVariable('xa','i4',('xa',))
xva.long_name = "eastward grid point number"
xva.units = "n/a"

ya = rootgrp.createDimension('ya',nya)
yva = rootgrp.createVariable('ya','i4',('ya',))
yva.long_name = "northward grid point number"
yva.units = "n/a"

lons = rootgrp.createVariable('longitudes','f4',('ya','xa'))
lons.long_name = "longitude"
lons.units = "degrees_east"

lats = rootgrp.createVariable('latitudes','f4',('ya','xa',))
lats.long_name = "latitude"
lats.units = "degrees_north"

xlocations = rootgrp.createVariable('xlocations','i4',('xsupp','ya','xa',), \
    zlib=True,least_significant_digit=3)
xlocations.units = "none"
xlocations.long_name = "fortran-array grid point locations"
xlocations.valid_range = [1,nxa]
xlocations.missing_value = -99

ylocations = rootgrp.createVariable('ylocations','i4',('xsupp','ya','xa',), \
    zlib=True,least_significant_digit=3)
ylocations.units = "none"
ylocations.long_name = "fortran-array grid point locations"
ylocations.valid_range = [1,nya]
ylocations.missing_value = -99

conusmask  = rootgrp.createVariable('conusmask','i2',('ya','xa',), \
    zlib=True,least_significant_digit=1)
conusmask.units = "none"
conusmask.long_name = "1 = inside CONUS (or Columbia Basin), 0 = outside"
conusmask.valid_range = [0,1]
conusmask.missing_value = -99

#practical_mask  = rootgrp.createVariable('practical_mask','i2',('ya','xa',), \
#    zlib=True,least_significant_digit=1)
#practical_mask.units = "none"
#practical_mask.long_name = "1 = inside CONUS (or Columbia Basin), 0 = outside"
#practical_mask.valid_range = [0,1]
#practical_mask.missing_value = -99

nsupplemental  = rootgrp.createVariable('nsupplemental','i2',('ya','xa',), \
    zlib=True,least_significant_digit=1)
nsupplemental.units = "none"
nsupplemental.long_name = "number of supplemental locations to use for this grid point"
nsupplemental.valid_range = [0,100]
nsupplemental.missing_value = 0

rootgrp.stream = "s4" 
rootgrp.title = "Supplemental 1/8 degree grid point locations to be used in POP12 calibration for "+\
    "NOAA's National Blend of Models POP12 project, version 1.0"
rootgrp.Conventions = "CF-1.0"  
rootgrp.history = "Created ~Dec 2015 by Tom Hamill, tom.hamill@noaa.gov"
rootgrp.institution = "NOAA/OAR/ESRL/PSD"
rootgrp.platform = "Model"
rootgrp.references = "None"

# ---- now initialize the netCDF variables.

xvsupp[:] = np.arange(nsupp,dtype=np.int16)
xva[:] = np.arange(nxa,dtype=np.int16)
yva[:] = np.arange(nya,dtype=np.int16)
lons[:,:] = rlonsa_t[:,:]
lats[:,:] = rlatsa_t[:,:]
xlocations[:,:,:] = xlocationa_t[:,:,:]
ylocations[:,:,:] = ylocationa_t[:,:,:]
conusmask[:,:] = conusmaska_t[:,:]
#practical_mask[:,:] = practical_mask_t[:,:]
nsupplemental[:,:] = numsupp_t[:,:]

print 'min(nsupp_t), max(numsupp_t) =  ', np.min(numsupp_t), np.max(numsupp_t)


# ---- close the file, and done

rootgrp.close()

