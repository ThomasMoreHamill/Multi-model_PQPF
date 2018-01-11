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
from calculate_quantile_for_precip_values import calculate_quantile_for_precip_values
import os

cmonths = ['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']
cleadbs = ['00','12']
cleades = ['12','00']
data_directory = '/Projects/Reforecast2/netcdf/NationalBlend/'
cmonthsin = ['01','02','03','04','05','06','07','08','09','10','11','12']

# ---- get the shape of the input arrays, and infer conus mask

infile = data_directory + 'blend.precip_const.2p5.nc'
nc = Dataset(infile)
conusmask_in = nc.variables['conuslandmask'][:,:] 
validmask_in = nc.variables['validmask'][:,:] 
lats_terrain = nc.variables['latitude'][:,:] 
lons_terrain = nc.variables['longitude'][:,:] 
nyat, nxat = np.shape(lats_terrain)
print 'nyat, nxat = ', nyat, nxat
nc.close()
conusmask_in[61,1335] = 0

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

print 'lats_terrain[0,0], lats[0,0] = ', lats_terrain[0,0], lats[0,0] 
print 'lons_terrain[0,0], lons[0,0] = ', lons_terrain[0,0], lons[0,0]

alphahat = -99.99*np.ones((nya,nxa), dtype = np.float64) 
betahat = -99.99*np.ones((nya,nxa), dtype = np.float64) 
fraction_zeros = -99.99*np.ones((nya,nxa), dtype = np.float64) 

# ---- loop over months of year.

for imonth in range(12):
#for imonth in range(3,4):
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
        npositive = np.zeros((nya,nxa), dtype = np.float64)
        nzeros = np.zeros((nya,nxa), dtype = np.float64)
        ones = np.ones((nya,nxa), dtype = np.float64)
        zeros = np.zeros((nya,nxa), dtype = np.float64)

        # ---- loop over years and days.  Read in the relevant precipitation analyses.
        #      Increment arrays that collect information relevant to calculating 
        #      fraction zero, alpha, beta.

        ndaysomo = [31,28,31,30,31,30, 31,31,30,31,30,31]
        for iyear in range(2002, 2017):
        #for iyear in range(2002, 2003):
            print 'processing year = ', iyear
            cyearmo = str(iyear) + cmonthsin[imonth]            
            infilename = data_directory + 'ccpa2p5.'+cyearmo+'.nc' # 1st day in 2002 may be corrupt
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
                #print 'reading precip data, index = ',itemindex[0]
                if itemindex[0] >= 0:
            
                    # ---- 12-hourly precip analysis is sum of the two 6-hourly analyses.

                    apcp_anal = np.squeeze(nc.variables['apcp_anal_12h'][itemindex[0],:,:])
                    if np.max(apcp_anal) >= 0:
                        #print 'apcp_anal[61,1335] = ',apcp_anal[61,1335]
                        sumx = np.where(apcp_anal*conusmask_in > 0, sumx + apcp_anal, sumx)
                        sumlnx = np.where(apcp_anal*conusmask_in > 0, sumlnx + np.log(apcp_anal), sumlnx)
                        npositive = np.where(apcp_anal*conusmask_in > 0, npositive+ones, npositive)
                        nzeros = np.where(apcp_anal*conusmask_in == 0, nzeros+ones, nzeros)


        # ---- in locations with extremely dry climatologies, we may have only a few
        #      precipitation samples, not enough to reliably estimate distribution 
        #      characteristics.  In such a case, bogus in quasi-realistic values
        #      that amount to an exponential distribution with a mean precipitation 
        #      amount of 0.5 mm.   NOTE:   this is effectively tailored to the CONUS
        #      and 2002-2016 precip data.  If this is being applied to another training
        #      data set, the statistician may wish to reconsider these assumptions.
        #
        #      if not extremely dry, then calculate parameters following Wilks text
        #      and Thom (1958) estimator using D statistic.

        #print 'npositive[61,1335], nzeros[61,1335] = ', npositive[61,1335], nzeros[61,1335]
        #print 'conusmask_in[61,1335] = ', conusmask_in[61,1335]
        for j in range(nya):
            for i in range(nxa):
                if conusmask_in[j,i] == 1 :
                    if npositive[j,i] < 4: # not enough data
                        if nzeros[j,i] + npositive[j,i] == 0:
                            print 'conusmask:   no data at python index j,i,lat,lon = ',j,i,lats[j,i], lons[j,i]
                            alphahat[j,i] = -99.99
                            betahat[j,i] = -99.99
                            fraction_zeros[j,i]= -99.99
                            conusmask_in[j,i] = 0
                        else:
                            alphahat[j,i] = 1.0
                            betahat[j,i] = 0.5                          
                            fraction_zeros[j,i] = nzeros[j,i] / (npositive[j,i]+nzeros[j,i])
                    else:
                        D = np.log(sumx[j,i]/npositive[j,i]) - sumlnx[j,i]/npositive[j,i]
                        alphahat[j,i] = (1.0 + np.sqrt(1.0 + 4.0*D/3)) / (4.0*D)
                        betahat[j,i] = (sumx[j,i]/npositive[j,i]) / alphahat[j,i]
                        fraction_zeros[j,i] = nzeros[j,i] / (npositive[j,i]+nzeros[j,i])
                        if fraction_zeros[j,i] > 0.97 and alphahat[j,i] > 1.0:
                            alphahat[j,i] = 1.0 # not enough data
                            betahat[j,i] = 0.5                                  
                elif validmask_in[j,i] == 1 :
                    if npositive[j,i] < 4: # not enough data
                        if nzeros[j,i] + npositive[j,i] == 0:
                            #print 'validmask:   no data at python index j,i,lat,lon = ',j,i,lats[j,i], lons[j,i]
                            alphahat[j,i] = -99.99
                            betahat[j,i] = -99.99
                            fraction_zeros[j,i]= -99.99
                            validmask_in[j,i] = 0
                        else:
                            alphahat[j,i] = 1.0
                            betahat[j,i] = 0.5                          
                            fraction_zeros[j,i] = nzeros[j,i] / (npositive[j,i] + nzeros[j,i])
                    else:
                        D = np.log(sumx[j,i]/npositive[j,i]) - sumlnx[j,i]/npositive[j,i]
                        alphahat[j,i] = (1.0 + np.sqrt(1.0 + 4.0*D/3)) / (4.0*D)
                        betahat[j,i] = (sumx[j,i]/npositive[j,i]) / alphahat[j,i]
                        fraction_zeros[j,i] = nzeros[j,i] / (npositive[j,i]+nzeros[j,i])
                        if fraction_zeros[j,i] > 0.97 and alphahat[j,i] > 1.0:
                            alphahat[j,i] = 1.0 # not enough data
                            betahat[j,i] = 0.5
                else:
                    alphahat[j,i] = -99.99
                    betahat[j,i] = -99.99
                    fraction_zeros[j,i] = -99.99   
        
        #print 'fraction_zeros[61,1335] = ', fraction_zeros[61,1335]         
        #sys.exit()      
        
        rthreshes = np.array([0.25, 0.5, 1.0, 3.0, 5.0, 7.0, 10.0, 15.0], dtype=np.float32)
        npthreshes = len(rthreshes)
        quant_fpamt = np.zeros((npthreshes, nya, nxa), dtype=np.float32)
        print 'npthreshes, nya, nxa = ', npthreshes, nya, nxa 
        print calculate_quantile_for_precip_values.__doc__
        quant_fpamt = calculate_quantile_for_precip_values(npthreshes, nya, nxa,  \
            rthreshes, fraction_zeros, alphahat, betahat, conusmask_in)
            
        # ---- store distribution parameters to netCDF file.

        outfilename = data_directory + \
            'climatology_gamma_parameters_ndfd2p5_'+cleadb+'_to_'+cleade+'_'+cmonthname+'.nc'
        print '   writing gamma climatology information to :'
        print '   ', outfilename

        rootgrp = Dataset(outfilename,'w',format='NETCDF4_CLASSIC')

        xa = rootgrp.createDimension('xa',nxa)
        xav = rootgrp.createVariable('xa','i4',('xa',))
        xav.long_name = "analysis grid eastward distance from southwest corner of domain in grid points" 
        xav.units = "grid index (dimensionless)" 

        ya = rootgrp.createDimension('ya',nya)
        yav = rootgrp.createVariable('ya','i4',('ya',))
        yav.long_name = "analysis grid northward distance from southwest corner of domain in grid points" 
        yav.units = "grid index (dimensionless)"
        
        tha = rootgrp.createDimension('tha',npthreshes)
        thav = rootgrp.createVariable('thav','i4',('tha',))
        thav.long_name = "number index associated with this precip threshold used in CDF quantifications" 
        thav.units = "index (dimensionless)"

        lonsa = rootgrp.createVariable('lonsa','f4',('ya','xa',))
        lonsa.long_name = "longitude" 
        lonsa.units = "degrees_east" 

        latsa = rootgrp.createVariable('latsa','f4',('ya','xa',))
        latsa.long_name = "latitude" 
        latsa.units = "degrees_north" 

        conusmask = rootgrp.createVariable('conusmask','i2',('ya','xa',))
        conusmask.long_name = "mask for grid points inside CONUS (1=yes,0=no)"
        conusmask.units=""
        
        validmask = rootgrp.createVariable('validmask','i2',('ya','xa',))
        validmask.long_name = "mask for grid points inside domain w usable data (1=yes,0=no)"
        validmask.units=""

        pthreshes = rootgrp.createVariable('pthreshes','f4',('tha',),
            zlib=True,least_significant_digit=5)  
        pthreshes.units = "" 
        pthreshes.long_name = "precipitation thresholds CDF evaluated at" 
        pthreshes.valid_range = [0.0,50.]
        pthreshes.missing_value = -99.99

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

        quantile_fpamt = rootgrp.createVariable('quantile_fpamt','f4',('tha','ya','xa',),
            zlib=True,least_significant_digit=4)  
        quantile_fpamt.units = "" 
        quantile_fpamt.long_name = "quantile associated with precip amt, and grid pt location" 
        quantile_fpamt.valid_range = [0.0,1000.]
        quantile_fpamt.missing_value = -99.99

        rootgrp.latcorners = [lats[0,0], lats[0,-1], lats[-1,0], lats[-1,-1]]
        rootgrp.loncorners = [lons[0,0], lons[0,-1], lons[-1,0], lons[-1,-1]]

        rootgrp.stream = "s4" # ????
        rootgrp.title = "climatological precipitation fraction zero, Gamma dist parameters, '+\
             'and quantiles associated with precip amts for this month"
        rootgrp.Conventions = "CF-1.0"  # ????
        rootgrp.history = "Created Aug 2017 by Hamill" 
        rootgrp.institution = "Processing at NOAA ESRL/PSD"
        rootgrp.platform = "Model" 
        rootgrp.references = "" 
    
        # --- set values to output

        xav = range(nxa)
        yav = range(nya)
        thav = rthreshes[:]
        lonsa[:,:]  = lons[:,:]
        latsa[:,:]  = lats[:,:]
        conusmask[:,:] = conusmask_in[:,:]
        validmask[:,:] = validmask_in[:,:]
        fraction_zero[:,:] = fraction_zeros[:,:]
        alpha[:,:] = alphahat[:,:]
        beta[:,:] = betahat[:,:]
        quantile_fpamt[:,:,:] = quant_fpamt[:,:]
        pthreshes[:] = rthreshes[:]

        rootgrp.close()
        
nc.close()
print 'Finished.'