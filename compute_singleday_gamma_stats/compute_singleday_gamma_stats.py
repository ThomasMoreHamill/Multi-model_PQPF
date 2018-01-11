""" this script synthesizes the important information from a particular
day's forecast data and analyzed data, writing the synthesized information
to a netCDF file.   At a later point (not in this script) a program will
read 60 days' worth of this synthesized data and from it will generate
Gamma-distributed CDFs of the forecast and analyzed, used for quantile
mapping.

written by Tom Hamill, tom.hamill@noaa.gov, (303) 497-3060 """


import numpy as np
import numpy.ma as ma
from netCDF4 import Dataset
import sys
import pygrib
from dateutils import hrstodate, daterange, dayofyear, \
     splitdate, datetohrs, dateshift, dateto_hrs_since_day1CE
from mpl_toolkits.basemap import interp
import os
import time as timey

# --- queries from command line

cyyyymmddhh = sys.argv[1]  # the initial date/time of the forecast
center = sys.argv[2] # ECMWF, NCEP, or CMC currently
cleade = sys.argv[3] # the ending lead time of the forecast in hours, e.g. 24 = 24h

if center == 'CMCE' or center == 'CMC':
    exchangeable = False  # need to tally the fcst stats member by member
else:
    exchangeable = True   # ok to tally fcst stats over all members

ileade = int(cleade)
ileadb = ileade - 12

print 'processing day = ',cyyyymmddhh,' center = ',center,' lead = ', cleade

# ---- read in the associated analysis
    
iyyyymmddhh_anal_end = int(dateshift(cyyyymmddhh,ileade))
iyyyymmddhh_anal_begin = int(dateshift(cyyyymmddhh,ileade-6))
    
# ---- read in ccpa precip analysis data for this date 
#      (12-h accum, generated from two 6-hourly files)

infilename = '/Projects/Reforecast2/netcdf/NationalBlend/precip_analyses_ccpa_v1_'+\
            '2002010100_to_2016123100.nc'
print infilename
nc = Dataset(infilename)
yyyymmddhh_anal_list = nc.variables['yyyymmddhh_anal_end'][:]
lons_anal_in = nc.variables['lons_anal'][:,:]
lats_anal_in = nc.variables['lats_anal'][:,:]
conusmask_in = nc.variables['conusmask'][:,:]
nya, nxa = lons_anal_in.shape
            
idx_today_begin = int(np.where(yyyymmddhh_anal_list == iyyyymmddhh_anal_begin)[0])
idx_today_end = int(np.where(yyyymmddhh_anal_list == iyyyymmddhh_anal_end)[0])

if idx_today_begin >= 0 and idx_today_end > 0:
    apcp_anal = nc.variables['apcp_anal'][idx_today_begin,:,:] + \
        nc.variables['apcp_anal'][idx_today_end,:,:]
else:
    apcp_anal = -99.99*np.ones((nya, nxa), dtype=np.float32)
pmax_anal = np.max(apcp_anal)
    
# ---- read in the forecast data for this date, and for the date 12 h previous
#      then subtract to get accumulated precip during this period
    
infile = '/Users/thamill/precip/ecmwf_data/'+center+'_'+\
    cyyyymmddhh+'_leadtime'+cleade+'h.nc'
print infile
nc2 = Dataset(infile)
apcp_fcst_ens_late = nc2.variables['apcp_fcst_ens'][:,:,:] # mbr,x,y
nc2.close()
    
cleadb = str(int(cleade)-12)
infile = '/Users/thamill/precip/ecmwf_data/'+center+'_'+\
    cyyyymmddhh+'_leadtime'+cleadb+'h.nc'
print infile
nc2 = Dataset(infile)
apcp_fcst_ens_early = nc2.variables['apcp_fcst_ens'][:,:,:] # mbr,x,y
nc2.close()
    
#   --- here there is a little trick so that we can use arrays with
#       the same numbers of dimensions for exchangeable ensembles (NCEP, ECMWF)
#       vs. not (CMC, and in the future SREF).   We'll have a third dimension
#       for the forecast array, dimensioned nmembers_qmap.  This is a dummy
#       dimension set to 1 for exchangeable systems.

nmembers = np.shape(apcp_fcst_ens_early)[0]
if exchangeable == False:
    nmembers_qmap = nmembers
else:
    nmembers_qmap = 1
    
apcp_fcst_ens = apcp_fcst_ens_late - apcp_fcst_ens_early
pmax_fcst = np.max(apcp_fcst_ens)
#print 'maximum forecast = ', pmax_fcst
#print 'maximum early forecast = ', np.max(apcp_fcst_ens_early)
#print 'maximum late forecast = ', np.max(apcp_fcst_ens_late)
    
# --- define arrays necessary to calculate fraction zero and Gamma distribution
#     parameters.   This is the compressed information that we write out to 
#     a netCDF file, later to read in over a 60-day period to calculate
#     fraction zero, gamma shape and scale parameters.

sumx_forecast_out = np.zeros((nmembers_qmap,nya,nxa), dtype = np.float32)
sumlnx_forecast_out = np.zeros((nmembers_qmap,nya,nxa), dtype = np.float32)
npositive_forecast_out = np.zeros((nmembers_qmap,nya,nxa), dtype = np.int16)
nzeros_forecast_out = np.zeros((nmembers_qmap,nya,nxa), dtype = np.int16)
        
ones = np.ones((nya,nxa), dtype = np.int16)
zeros= np.zeros((nya,nxa), dtype = np.int16)
    
sumx_analysis_out = np.zeros((nya,nxa), dtype = np.float32)
sumlnx_analysis_out = np.zeros((nya,nxa), dtype = np.float32)
npositive_analysis_out = np.zeros((nya,nxa), dtype = np.int16)
nzeros_analysis_out = np.zeros((nya,nxa), dtype = np.int16)
    
# ---- now let's increment the arrays assuming data is available

if pmax_anal >= 0 and pmax_fcst >=0 :
            
    sumx_analysis_out = np.where(apcp_anal*conusmask_in.astype(float) > 0.0001, apcp_anal, zeros)
    sumlnx_analysis_out = np.where(apcp_anal*conusmask_in.astype(float) > 0.0001, \
        np.log(apcp_anal), zeros)
    npositive_analysis_out = np.where(apcp_anal*conusmask_in.astype(float) > 0, ones, zeros)
    nzeros_analysis_out = np.where(apcp_anal*conusmask_in.astype(float) == 0, ones, zeros)
    for imem in range(nmembers):
        if exchangeable == False:
            imem_qmap = imem
        else:
            imem_qmap = 0
        apcp_fcst = apcp_fcst_ens[imem,:,:]
        apcp_fcst = np.where(apcp_fcst >= 0.0, apcp_fcst, zeros)
        sumx_forecast_out[imem_qmap,:,:] = sumx_forecast_out[imem_qmap,:,:] + \
            np.where(apcp_fcst*conusmask_in.astype(float) > 0.0001, apcp_fcst, zeros)
        sumlnx_forecast_out[imem_qmap,:,:] = sumlnx_forecast_out[imem_qmap,:,:] + \
            np.where(apcp_fcst*conusmask_in.astype(float) > 0.0001, np.log(apcp_fcst), zeros)
        npositive_forecast_out[imem_qmap,:,:] = npositive_forecast_out[imem_qmap,:,:] + \
            np.where(apcp_fcst*conusmask_in.astype(float) > 0, ones, zeros)
        nzeros_forecast_out[imem_qmap,:,:] = nzeros_forecast_out[imem_qmap,:,:] + \
            np.where(apcp_fcst*conusmask_in.astype(float) == 0, ones, zeros)

    # ---- open and initialize the netCDF file we'll be writing to.

    outfilename = '/Users/thamill/precip/ecmwf_data/'+center+'/'+center+\
        '_D_statistics_data_IC'+cyyyymmddhh+'_'+cleade+'h.nc'
        
    print 'writing to ',outfilename
    rootgrp = Dataset(outfilename,'w',format='NETCDF4_CLASSIC')
    
    xa = rootgrp.createDimension('xa',nxa)
    xav = rootgrp.createVariable('xa','f4',('xa',))
    xav.long_name = "analysis grid eastward distance from southwest corner of domain in grid points" 
    xav.units = "grid index (dimensionless)" 
    
    ya = rootgrp.createDimension('ya',nya)
    yav = rootgrp.createVariable('ya','f4',('ya',))
    yav.long_name = "analysis grid northward distance from southwest corner of domain in grid points" 
    yav.units = "grid index (dimensionless)"
    
    ens = rootgrp.createDimension('ens',nmembers_qmap)
    ensv = rootgrp.createVariable('ens','i4',('ens',))
    ensv.long_name = "effective member number for quantile mapping (0=exchangeable, 0 to nmembers-1 if not)" 
    ensv.units = "" 
    
    lons_anal = rootgrp.createVariable('lons_anal','f4',('ya','xa',))
    lons_anal.long_name = "longitude" 
    lons_anal.units = "degrees_east" 
    
    lats_anal = rootgrp.createVariable('lats_anal','f4',('ya','xa',))
    lats_anal.long_name = "latitude" 
    lats_anal.units = "degrees_north" 

    conusmask = rootgrp.createVariable('conusmask','i2',('ya','xa',))
    conusmask.long_name = "mask for grid points inside CONUS (1=yes,0=no)"
    conusmask.units=""

    sumx_analysis = rootgrp.createVariable('sumx_analysis','f4',('ya','xa',),
        zlib=True,least_significant_digit=5)  
    sumx_analysis.units = "" 
    sumx_analysis.long_name = "Precipitation amount, if nonzero (mm)" 
    sumx_analysis.valid_range = [0.0,1000.0]
    sumx_analysis.missing_value = -99.99
    
    sumlnx_analysis = rootgrp.createVariable('sumlnx_analysis','f4',('ya','xa',),
        zlib=True,least_significant_digit=5)  
    sumlnx_analysis.units = "" 
    sumlnx_analysis.long_name = "Natural log of precipitation amount, if nonzero (mm)" 
    sumlnx_analysis.valid_range = [-50,10.0]
    sumlnx_analysis.missing_value = -99.99

    npositive_analysis = rootgrp.createVariable('npositive_analysis','i2',('ya','xa',),
        zlib=True,least_significant_digit=5)  
    npositive_analysis.units = "" 
    npositive_analysis.long_name = "Flagged to 1 if sample had positive precip" 
    npositive_analysis.valid_range = [0,1]
    npositive_analysis.missing_value = -99
    
    nzeros_analysis = rootgrp.createVariable('nzeros_analysis','i2',('ya','xa',),
        zlib=True,least_significant_digit=5)  
    nzeros_analysis.units = "" 
    nzeros_analysis.long_name = "Flagged to 1 if sample had zero precip" 
    nzeros_analysis.valid_range = [0,1]
    nzeros_analysis.missing_value = -99

    sumx_forecast = rootgrp.createVariable('sumx_forecast','f4',('ens','ya','xa',),
        zlib=True,least_significant_digit=5)  
    sumx_forecast.units = "" 
    sumx_forecast.long_name = "Forecast precipitation amount for this member, if nonzero (mm)" 
    sumx_forecast.valid_range = [0.0,1000.0]
    sumx_forecast.missing_value = -99.99
    
    sumlnx_forecast = rootgrp.createVariable('sumlnx_forecast','f4',('ens','ya','xa',),
        zlib=True,least_significant_digit=5)  
    sumlnx_forecast.units = "" 
    sumlnx_forecast.long_name = "Natural log of forecast precipitation amount for this member, if nonzero (mm)" 
    sumlnx_forecast.valid_range = [-500.,100.0]
    sumlnx_forecast.missing_value = -99.99

    npositive_forecast = rootgrp.createVariable('npositive_forecast','i2',('ens','ya','xa',),
        zlib=True,least_significant_digit=5)  
    npositive_forecast.units = "" 
    npositive_forecast.long_name = "1 if positive precip" 
    npositive_forecast.valid_range = [0,nmembers]
    npositive_forecast.missing_value = -99
    
    nzeros_forecast = rootgrp.createVariable('nzeros_forecast','i2',('ens','ya','xa',),
        zlib=True,least_significant_digit=5)  
    nzeros_forecast.units = "" 
    nzeros_forecast.long_name = "1 if zero precip" 
    nzeros_forecast.valid_range = [0,nmembers]
    nzeros_forecast.missing_value = -99
 
    rootgrp.latcorners = [lats_anal_in[0,0], lats_anal_in[0,-1], lats_anal_in[-1,0], lats_anal_in[-1,-1]]
    rootgrp.loncorners = [lons_anal_in[0,0], lons_anal_in[0,-1], lons_anal_in[-1,0], lons_anal_in[-1,-1]]
    rootgrp.stream = "s4" # ????
    rootgrp.title = "analysis and forecast statistics underlying Gamma distribution parameters"
    rootgrp.Conventions = "CF-1.0"  # ????
    rootgrp.history = "Created Aug 2017 by Tom Hamill" 
    rootgrp.institution = "ESRL/PSD using CCPA data from NCEP/EMC"
    rootgrp.platform = "Model" 
    rootgrp.references = "" 

    if exchangeable == False:
        ensv[:] = range(nmembers)
    else:
        ensv[0] = 0

    # ---- define the corner lat/lons.

    xav[:]   = np.arange(nxa)
    yav[:]   = np.arange(nya)
    lons_anal[:,:]  = lons_anal_in
    lats_anal[:,:]  = lats_anal_in
    conusmask[:,:] = conusmask_in

    # --- write out data
        
    sumx_forecast[:,:,:] = sumx_forecast_out[:,:,:]
    sumlnx_forecast[:,:,:] = sumlnx_forecast_out[:,:,:]   
    npositive_forecast[:,:,:] = npositive_forecast_out[:,:,:]
    nzeros_forecast[:,:,:] = nzeros_forecast_out[:,:,:]

    npositive_analysis[:,:] = npositive_analysis_out[:,:]
    nzeros_analysis[:,:] = nzeros_analysis_out[:,:]
    sumx_analysis[:,:] = sumx_analysis_out[:,:]
    sumlnx_analysis[:,:] = sumlnx_analysis_out[:,:]

    rootgrp.close()

nc.close()