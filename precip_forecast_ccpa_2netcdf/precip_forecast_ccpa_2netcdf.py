""" this script reads in grib forecast data and writes back out an equivalent
    netCDF file. """

import numpy as np
import sys
#import /Users/jwhitaker/python/pygrib.git/pygrib
import pygrib
import os
import time as timey
from netCDF4 import Dataset
from dateutils import hrstodate, daterange, dayofyear, \
     splitdate, datetohrs, dateshift, dateto_hrs_since_day1CE
from mpl_toolkits.basemap import interp
import cPickle

def read_gribdata(gribfilename):
    grb=[]
    istat = -1
    for iter in range(10):
        fexist_grib = False
        fexist_grib = os.path.exists(gribfilename)
        if fexist_grib :
            try:
                fcstfile = pygrib.open(gribfilename)
                grb = fcstfile.select()[0]
                istat = 0
                fcstfile.close()
                return istat,grb
            except (IOError, ValueError, RuntimeError):
                print 'Error reading ', gribfilename
                istat = -1
    return istat,grb

# =====================================================================================

center = sys.argv[1] # ECMWF, NCEP, or CMC currently (all ensembles)
cleade = sys.argv[2] # use 24 for 24 h end time, not 024

# ---- read in sample forecast in order to define the forecast lat/lon array

infilename = '/Users/thamill/precip/ecmwf_data/ECMWF_2016010100_fhour12_pertno1.grb'
grib = pygrib.open(infilename)
grb2 = grib.select()[0]
latsf, lonsf = grb2.latlons()
lonsf = lonsf - 360.
lonsf_1d = lonsf[0,:]
latsf_1d = latsf[:,0]
latsf_1d = np.flip(latsf_1d,0)
nyf, nxf = np.shape(latsf)
grib.close()

# ---- define a list of dates to read in, set number of members to read in

date_begin = '2016010100'
date_end = '2016063000'
date_list = daterange(date_begin, date_end, 24)
ndates = len(date_list)
ileade = int(cleade)

if center == 'ECMWF':
    nmembers = 50
else:
    nmembers = 20 # NCEP and CMC

# ---- read in sample conus mask, lat, lon

infile = '/Users/thamill/precip/conusmask_ccpa.nc'
nc = Dataset(infile)
conusmask_in = nc.variables['conusmask'][:,:]
lons = nc.variables['lons'][:,:]
lats = nc.variables['lats'][:,:]
nya, nxa = np.shape(lons)
nc.close()

# --- loop over dates

for idate, date in zip(range(ndates), date_list):

    # --- set up the output netCDF file information
    
    outfile = '/Users/thamill/precip/ecmwf_data/'+center+'_'+date+'_leadtime'+cleade+'h.nc'
    print outfile, timey.asctime()
    rootgrp = Dataset(outfile,'w',format='NETCDF4_CLASSIC')
    
    xa = rootgrp.createDimension('xa',nxa)
    xva = rootgrp.createVariable('xa','f4',('xa',))
    xva.long_name = "eastward grid point number, precip forecasts" 
    xva.units = "n/a" 
    
    ya = rootgrp.createDimension('ya',nya)
    yva = rootgrp.createVariable('ya','f4',('ya',))
    yva.long_name = "northward grid point number, precip forecasts" 
    yva.units = "n/a" 

    ens = rootgrp.createDimension('ens',nmembers)
    ensv = rootgrp.createVariable('ensv','i4',('ens',))
    ensv.long_name = "Ensemble member number" 
    ensv.units = " " 
    
    time = rootgrp.createDimension('time',None)
    timev = rootgrp.createVariable('time','i4',('time',))
    timev.units = "index into the file for time dimension"

    lonsa = rootgrp.createVariable('lons_anal','f4',('ya','xa',))
    lonsa.long_name = "longitude" 
    lonsa.units = "degrees_east" 

    latsa = rootgrp.createVariable('lats_anal','f4',('ya','xa',))
    latsa.long_name = "latitude" 
    latsa.units = "degrees_north" 
    
    conusmask = rootgrp.createVariable('conusmask','i2',('ya','xa',))
    conusmask.long_name = "CONUS mask (1=inside, 0=outside)" 
    conusmask.units = ""
    conusmask.valid_range = [0,1]
    conusmask.missing_value = -99

    apcp_fcst_ens = rootgrp.createVariable('apcp_fcst_ens','f4',('ens','ya','xa',),
        zlib=True,least_significant_digit=2)  
    apcp_fcst_ens.units = "mm" 
    apcp_fcst_ens.long_name = "Ensemble member forecast accumulated precipitation" 
    apcp_fcst_ens.valid_range = [0.,1000.]
    apcp_fcst_ens.missing_value = -99.99

    rootgrp.stream = "s4" # ????
    rootgrp.title = \
        "real-time ensemble forecast data interpolated to 1/8-degree CCPA grid over CONUS"
    rootgrp.Conventions = "CF-1.0"  # ????
    rootgrp.history = "Created 8 May 2017 by Tom Hamill" 
    rootgrp.institution = \
        "CCPA from NCEP/EMC, forecast data via TIGGE ECMWF portal"
    rootgrp.platform = "Model" 
    rootgrp.references = "" 

    # ---- set values for constant fields in netcdf file
    
    xva[:] = np.arange(nxa)
    yva[:] = np.arange(nya)
    ensv[:] = range(nmembers)

    lonsa[:] = lons
    latsa[:] = lats
    conusmask[:] = conusmask_in

    # --- process, reading in the grib data and writing netCDF records

    apcpf_out = np.zeros((nya,nxa),dtype=np.float)
    for imem in range(nmembers):
    
        # ---- build the name of the file to read in, then read in forecast.
        #      input data has n to s data, interpolation requires s to n,
        #      so also flip the array upside down
        
        cmem = str(imem+1)
        infile = '/Users/thamill/precip/ecmwf_data/'+\
            center+'_'+date+'_fhour'+cleade+'_pertno'+cmem+'.grb'
        istat, grb = read_gribdata(infile)
        if istat != -1:
            fmax = np.max(grb.values)
            print 'reading from ',infile,' max = ',fmax
            if fmax > 0.0:
                fcst = np.flipud(grb.values)
            else:
                fcst = -99.99*np.ones((nyf,nxf), dtype=np.float32)
        else:
            fcst = -99.99*np.ones((nyf,nxf), dtype=np.float32)
        
        # ---- interpolate from coarser 1/2-degree forecast grid to 1/8-degree CCPA grid
 
        apcpf_out[:,:] = interp(fcst, lonsf_1d, latsf_1d, \
                lons, lats, checkbounds=False, masked=False, order=1)
        #print 'apcpf_out[nya/2,0:nxa:12] = ',apcpf_out[nya/2,0:nxa:12]
        
        # ---- write out this netcdf record.
        
        apcp_fcst_ens[imem] = apcpf_out
        
    print 'writing to ',outfile
    rootgrp.close()

