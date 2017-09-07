from netCDF4 import Dataset
import matplotlib
matplotlib.use('TKAgg')
import matplotlib.pyplot as plt
from matplotlib import rcParams
from mpl_toolkits.basemap import Basemap, addcyclic
import numpy as np
from numpy import ma
import math
import os, sys
from read_supp_locns_ndfd2p5 import read_supp_locns_ndfd2p5
from matplotlib import rcParams

cmonth = sys.argv[1] # Jan, Feb, etc

# ---- read in supplemental locations on 2.5-km grid

nxa = 2345
nya = 1597
nsupp = 50

xlocationa = np.zeros((nxa,nya,nsupp),dtype=np.int32)
ylocationa = np.zeros((nxa,nya,nsupp),dtype=np.int32)
conusmaska = np.zeros((nxa,nya),dtype=np.int16)
practical_mask = np.zeros((nxa,nya),dtype=np.int16)
rlona = np.zeros((nxa,nya),dtype=np.float)
rlata = np.zeros((nxa,nya),dtype=np.float)

xlocationa, ylocationa, conusmaska, rlona, rlata = \
    read_supp_locns_ndfd2p5 (cmonth, nxa, nya, nsupp)
    
print 'rlona[nya/2, 0:nxa:20] = ', rlona[nya/2, 0:nxa:20]
print 'rlata[nya/2, 0:nxa:20] = ', rlata[nya/2, 0:nxa:20]

zeros = np.zeros((nxa,nya), dtype=np.int16)
conusmaska = np.where(conusmaska < 0, zeros, conusmaska)

print 'conusmaska[0:-1:20,nya/2] = ',conusmaska[0:-1:20,nya/2]

numsupp_t = np.zeros((nya,nxa), dtype=np.int16)
for ixa in range(nxa):
    if ixa%100 == 0: print 'processing ixa = ',ixa,' of ',nxa
    for jya in range(nya):
        if conusmaska[ixa,jya] == 1:
            vec = xlocationa[ixa,jya,:]
            isuppmin = np.argmin(vec)
            if xlocationa[ixa,jya,isuppmin] > 0:
                numsupp_t[jya,ixa] = nsupp
            else:
                numsupp_t[jya,ixa] = isuppmin + 1
        else:
            numsupp_t[jya,ixa] = -99

#  ---- prep arrays for dumping to netCDF by switching order.

xlocationa_t = np.transpose(xlocationa)
ylocationa_t = np.transpose(ylocationa)
rlatsa_t     = np.transpose(rlata)
rlonsa_t     = np.transpose(rlona)
conusmaska_t = np.transpose(conusmaska)

for i in range(nxa):
    print i, numsupp_t[nya/2,i], xlocationa_t[0,nya/2,i], \
        ylocationa_t[0,nya/2,i], xlocationa_t[49,nya/2,i], ylocationa_t[49,nya/2,i] 


#print 'numsupp_t[0:-1,nya/2] = ',numsupp_t[0:-1,nya/2]
#print 'xlocationa_t[0:-1,nya/2,0] = ', xlocationa_t[0:-1,nya/2,0]

# ---- set up output netCDF file metadata particulars

outfile_nc = 'supplemental_locations_ndfd2p5_'+cmonth+'.nc'
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

nsupplemental  = rootgrp.createVariable('nsupplemental','i2',('ya','xa',), \
    zlib=True,least_significant_digit=1)
nsupplemental.units = "none"
nsupplemental.long_name = "number of supplemental locations to use for this grid point"
nsupplemental.valid_range = [0,100]
nsupplemental.missing_value = -99

rootgrp.stream = "s4" 
rootgrp.title = "Supplemental NDFD 2.5-km grid locations to be used in POP12 and QPF06 calibration for "+\
    "NOAA's National Blend of Models project, version 3.1"
rootgrp.Conventions = "CF-1.0"  
rootgrp.history = "Created ~ May 2017 by Tom Hamill, tom.hamill@noaa.gov"
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
nsupplemental[:,:] = numsupp_t[:,:]

# ---- close the file, and done

rootgrp.close()

