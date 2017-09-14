"""
This python script is intended to compute the parameters of the Gamma distributions
for CCPA (Climatology-Calibrated Precipitation Analyses; Hou et al 2014; 
http://journals.ametsoc.org/doi/abs/10.1175/JHM-D-11-0140.1 ).  Specifically, we
compute the parameters FZ, alpha, and beta, where FZ is the fraction of samples 
with zero precipitation, and alpha and beta are the shape and scale parameters
of the Gamma distribution.  To calculate these parameters, we follow the methodology
of Thom (1958) that is described in the Wilks Statistical Methods in the Atmospheric
Sciences (3rd ed.).  See eqs. 4.40 - 4.42 therein.
"""

import numpy as np
from netCDF4 import Dataset
import sys
from dateutils import daterange, dateshift
import os

cmonths = ['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']
cleadbs = ['00','12']
cleades = ['12','00']
data_directory = '/Projects/Reforecast2/netcdf/NationalBlend/'
cmonthsin = ['01','02','03','04','05','06','07','08','09','10','11','12']

# ---- get the shape of the input arrays, and infer conus mask

infilename = data_directory + 'ccpa2p5.201607.nc'
print infilename
nc = Dataset(infilename)
apcp_anal = nc.variables['apcp_anal_12h'][0,:,:] # 1st day in 2002 may be corrupt
lats = nc.variables['latitude'][:,:] 
lons = nc.variables['longitude'][:,:] 
nc.close()
nya, nxa = np.shape(apcp_anal)
print 'nya, nxa = ', nya, nxa
ones = np.ones((nya,nxa), dtype = np.int16)
zeros = np.zeros((nya,nxa), dtype = np.int16) 
conusmask_in = np.where(apcp_anal < -99., zeros, ones)

# ---- loop over months of year.

for imonth in range(12):
    cmonthname = cmonths[imonth]
    print '************** processing month = ',cmonthname
    
    # ---- perhaps there is a dependence of the climatology from night to day.
    #      to account for this possibility, we will process 00-12 UTC accumulations
    #      and the process 12-00 UTC accumulations
    
    for cleadb, cleade in zip(cleadbs, cleades):
    #for cleadb, cleade in zip(cleadbs[1], cleades[1]):
        
        print '  processing ',cleadb,' to ',cleade, ' UTC'

        # --- define arrays necessary to calculate fraction zero and Gamma distribution
        #     parameters

        sumx = np.zeros((nya,nxa), dtype = np.float64)
        sumlnx = np.zeros((nya,nxa), dtype = np.float64)
        npositive = np.zeros((nya,nxa), dtype = np.int64)
        nzeros = np.zeros((nya,nxa), dtype = np.int64)
        ones = np.ones((nya,nxa), dtype = np.int64)
        zeros = np.zeros((nya,nxa), dtype = np.int64)

        # ---- loop over years and days.  Read in the relevant precipitation analyses.
        #      Increment arrays that collect information relevant to calculating 
        #      fraction zero, alpha, beta.

        ndaysomo = [31,28,31,30,31,30, 31,31,30,31,30,31]
        for iyear in range(2002, 2017):
            print 'processing year = ', iyear
            cyearmo = str(iyear) + cmonthsin[imonth]            
            infilename = data_directory + 'ccpa2p5.'+cyearmo+'.nc'
            nc = Dataset(infilename) 
            yyyymmddhh_anal_end = nc.variables['anal_date_12h'][:]
            ndays = ndaysomo[imonth]
            if iyear%4 == 0 and imonth == 1: ndays = 29
            for iday in range(1,ndays+1):
                if iday < 10:
                    cday = '0'+str(iday)
                else:
                    cday = str(iday)
                cyyyymmddhh = str(iyear) + cmonthsin[imonth] + cday + '00'
                if cleadb == '00':
                    ishift = 12
                else:
                    ishift = 24
                iyyyymmddhh = int(dateshift(cyyyymmddhh,ishift))
                itemindex = np.where(yyyymmddhh_anal_end == iyyyymmddhh)
                print 'reading precip data, index = ',itemindex[0]
                if itemindex[0] >= 0:
            
                    # ---- 12-hourly precip analysis is sum of the two 6-hourly analyses.
                    
                    apcp_anal = np.squeeze(nc.variables['apcp_anal_12h'][itemindex[0],:,:])
                    if np.max(apcp_anal) >= 0:
                        sumx = np.where(apcp_anal*conusmask_in > 0, sumx + apcp_anal, sumx)
                        sumlnx = np.where(apcp_anal*conusmask_in > 0, sumlnx + np.log(apcp_anal), sumlnx)
                        npositive = np.where(apcp_anal*conusmask_in > 0, npositive+ones, npositive)
                        nzeros = np.where(apcp_anal*conusmask_in == 0, nzeros+ones, nzeros)

        # ---- now calculate the D statistic for estimation of Gamma distribution 
        #      parameters following Wilks text.

        D = np.where(nzeros*conusmask_in > 0, np.log(sumx/npositive)-sumlnx/npositive, zeros)
        alphahat = np.where(D>0, (1.0 + np.sqrt(1.0 + 4.0*D/3))/(4.0*D), ones)
        betahat = np.where(D>0, (sumx/npositive)/alphahat, ones)
        fraction_zeros = np.where(conusmask_in > 0, \
            nzeros.astype(np.float32)/(npositive+nzeros).astype(np.float32), np.real(zeros))

        # ---- store distribution parameters to netCDF file.

        outfilename = data_directory + \
            'climatology_gamma_parameters_ndfd2p5_'+cleadb+'_to_'+cleade+'_'+cmonthname+'.nc'
        print '   writing gamma climatology information to :'
        print '   ', outfilename

        rootgrp = Dataset(outfilename,'w',format='NETCDF4_CLASSIC')

        xa = rootgrp.createDimension('xa',nxa)
        xav = rootgrp.createVariable('xa','f4',('xa',))
        xav.long_name = "analysis grid eastward distance from southwest corner of domain in grid points" 
        xav.units = "grid index (dimensionless)" 

        ya = rootgrp.createDimension('ya',nya)
        yav = rootgrp.createVariable('ya','f4',('ya',))
        yav.long_name = "analysis grid northward distance from southwest corner of domain in grid points" 
        yav.units = "grid index (dimensionless)"

        lonsa = rootgrp.createVariable('lonsa','f4',('ya','xa',))
        lonsa.long_name = "longitude" 
        lonsa.units = "degrees_east" 

        latsa = rootgrp.createVariable('latsa','f4',('ya','xa',))
        latsa.long_name = "latitude" 
        latsa.units = "degrees_north" 

        conusmask = rootgrp.createVariable('conusmask','i2',('ya','xa',))
        conusmask.long_name = "mask for grid points inside CONUS (1=yes,0=no)"
        conusmask.units=""

        fraction_zero = rootgrp.createVariable('fraction_zero','f4',('ya','xa',),
            zlib=True,least_significant_digit=5)  
        fraction_zero.units = "" 
        fraction_zero.long_name = "Fraction of samples with zero precipitation" 
        fraction_zero.valid_range = [0.0,1.0]
        fraction_zero.missing_value = -99.99

        alpha = rootgrp.createVariable('alpha','f4',('ya','xa',),
            zlib=True,least_significant_digit=5)  
        alpha.units = "" 
        alpha.long_name = "alpha, fitted shape parameter for Gamma distribution of + values" 
        alpha.valid_range = [0.0,1000.]
        alpha.missing_value = -99.99

        beta = rootgrp.createVariable('beta','f4',('ya','xa',),
            zlib=True,least_significant_digit=5)  
        beta.units = "" 
        beta.long_name = "beta, fitted scale parameter for Gamma distribution of + values" 
        beta.valid_range = [0.0,1000.]
        beta.missing_value = -99.99

        rootgrp.latcorners = [lats[0,0], lats[0,-1], lats[-1,0], lats[-1,-1]]
        rootgrp.loncorners = [lons[0,0], lons[0,-1], lons[-1,0], lons[-1,-1]]

        rootgrp.stream = "s4" # ????
        rootgrp.title = "climatological precipitation fraction zero and Gamma distribution parameters for this month"
        rootgrp.Conventions = "CF-1.0"  # ????
        rootgrp.history = "Created Aug 2017 by Hamill" 
        rootgrp.institution = "Processing at NOAA ESRL/PSD"
        rootgrp.platform = "Model" 
        rootgrp.references = "" 
    
        # --- set values to output

        xav = range(nxa)
        yav = range(nya)
        lonsa[:,:]  = lons[:,:]
        latsa[:,:]  = lats[:,:]
        conusmask[:,:] = conusmask_in[:,:]
        fraction_zero[:,:] = fraction_zeros[:,:]
        alpha[:,:] = alphahat[:,:]
        beta[:,:] = betahat[:,:]

        rootgrp.close()
        
nc.close()
print 'Finished.'