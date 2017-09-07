""" This python script will, for an input month, read in precipitation 
analyses for that month (data stored for every 6 h) and compute an 
exceedance probability for many common precipitation thresholds, 
storing that information to disk in a netCDF file. """

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

data_directory = '/Projects/Reforecast2/netcdf/NationalBlend/'

cmonths = ['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']
cleadbs = ['00','12']
cleades = ['12','00']

for cmonth in cmonths:
    print 'processing month = ', cmonth
    for cleadtime_begin, cleadtime_end in zip(cleadbs, cleades):
        print 'processing ',cleadtime_begin,' to ',cleadtime_end
        imonth = cmonths.index(cmonth) + 1
        if imonth < 10:
            cmono = '0' + str(imonth)
        else:
            cmono = str(imonth)

        # ---- determine lat/lon information and masks from a sample file
        #      conusmask is the strict mask for inside CONUS, practical_mask is the 
        #      set of points with CCPA data not missing (1=available, 0 = no)

        filename= data_directory + 'precip_analyses_ccpa_v1_2002010100_to_2016123100.nc'
        print 'reading ',filename
        nc = Dataset(filename)
        lats_anal = nc.variables['lats_anal'][:]
        lons_anal = nc.variables['lons_anal'][:]
        xavals = nc.variables['xa'][:]
        yavals = nc.variables['ya'][:]
        nia = len(xavals)
        nja = len(yavals)
        conusmask_in = nc.variables['conusmask'][:,:]
        iyyyymmddhh_list = nc.variables['yyyymmddhh_anal_end'][:]
        cyyyymmddhh_list = str(iyyyymmddhh_list)

        # ---- make a list of 1/0 for whether or not to read in this date in iyyyymmddhh_list
        #      for a given month we want data from that month's central date +/- 31 days

        ndates = len(iyyyymmddhh_list)
        usethisdate = np.zeros((ndates),dtype=np.int16)

        for iyear in range(2002, 2017):
            
            # ---- find yyyymmddhh that is approximately at the center of the month
            
            if cleadtime_begin == '12' and cleadtime_end == '00':
                cyyyymmddhh_center = str(iyear) + cmono + '1500'
            elif cleadtime_begin == '00' and cleadtime_end == '12':
                cyyyymmddhh_center = str(iyear) + cmono + '1512'
            else:
                print 'invalid begin, end times for accumulation period.  Stopping'
                sys.exit()
            cyyyymmddhh_begin = dateshift(cyyyymmddhh_center,-24*31)
            cyyyymmddhh_end = dateshift(cyyyymmddhh_center,24*31)
            date_list = daterange(cyyyymmddhh_begin, cyyyymmddhh_end, 24)
            
            # ---- go thru the list of dates that are in the data set, and 
            #      flag the elements that correspond to date_list.  These
            #      represent the indices for the first of two times we read in.
            
            for cdate in date_list:
                idate = int(cdate)
                idx = np.where(iyyyymmddhh_list == idate)
                if idx >= 0: usethisdate[idx] = 1

        # ---- estimate the climatology from relative frequency for selected thresholds

        ktr = 0    
        climo_POP = np.zeros((nja,nia), dtype=np.float32) # 0.254 mm
        climo_1mm = np.zeros((nja,nia), dtype=np.float32) 
        climo_2p5mm = np.zeros((nja,nia), dtype=np.float32) 
        climo_5mm = np.zeros((nja,nia), dtype=np.float32) 
        climo_10mm = np.zeros((nja,nia), dtype=np.float32) 
        climo_25mm = np.zeros((nja,nia), dtype=np.float32)
        climo_50mm = np.zeros((nja,nia), dtype=np.float32) 
        ones  = np.ones((nja,nia), dtype=np.float32)
        zeros = np.zeros((nja,nia), dtype=np.float32)       
        for idate in range(ndates):
            if usethisdate[idate] == 1 and idate < ndates-1:
                apcp_anal = nc.variables['apcp_anal'][idate,:,:] + \
                    nc.variables['apcp_anal'][idate+1,:,:]
                if apcp_anal[nja/2,nia/2] >= 0. :
                    work = np.where(apcp_anal > 0.254, ones, zeros)
                    climo_POP[:,:] = climo_POP[:,:] + work[:,:]
                    work = np.where(apcp_anal > 1.0, ones, zeros)
                    climo_1mm[:,:] = climo_1mm[:,:] + work[:,:]
                    work = np.where(apcp_anal > 2.5, ones, zeros)
                    climo_2p5mm[:,:] = climo_2p5mm[:,:] + work[:,:]
                    work = np.where(apcp_anal > 5.0, ones, zeros)
                    climo_5mm[:,:] = climo_5mm[:,:] + work[:,:]
                    work = np.where(apcp_anal > 10.0, ones, zeros)
                    climo_10mm[:,:] = climo_10mm[:,:] + work[:,:]
                    work = np.where(apcp_anal > 25.0, ones, zeros)
                    climo_25mm[:,:] = climo_25mm[:,:] + work[:,:]
                    work = np.where(apcp_anal > 50.0, ones, zeros)
                    climo_50mm[:,:] = climo_50mm[:,:] + work[:,:]
                    ktr = ktr+1

        nc.close()
        
        #print 'number of samples: ', ktr

        climo_POP = climo_POP / np.float(ktr)
        climo_1mm = climo_1mm / np.float(ktr)
        climo_2p5mm = climo_2p5mm / np.float(ktr)
        climo_5mm = climo_5mm / np.float(ktr)
        climo_10mm = climo_10mm / np.float(ktr)
        climo_25mm = climo_25mm / np.float(ktr)
        climo_50mm = climo_50mm / np.float(ktr)

        # ---- set points outside CONUS to -99.99

        mninetynine = -99.99*np.ones((nja,nia), dtype=np.float32)
        climo_POP = np.where(conusmask_in == 0, mninetynine, climo_POP)
        climo_1mm = np.where(conusmask_in == 0, mninetynine, climo_1mm)
        climo_2p5mm = np.where(conusmask_in == 0, mninetynine, climo_2p5mm)
        climo_5mm = np.where(conusmask_in == 0, mninetynine, climo_5mm)
        climo_10mm = np.where(conusmask_in == 0, mninetynine, climo_10mm)
        climo_25mm = np.where(conusmask_in == 0, mninetynine, climo_25mm)
        climo_50mm = np.where(conusmask_in == 0, mninetynine, climo_50mm)

        nthreshes = 7
        pthresh = np.array([0.254,1,2.5,5.0,10.0,25.0,50.0])
        climo_array = np.zeros((nthreshes, nja, nia))
        climo_array[0,:,:] = climo_POP[:,:]
        climo_array[1,:,:] = climo_1mm[:,:]
        climo_array[2,:,:] = climo_2p5mm[:,:]
        climo_array[3,:,:] = climo_5mm[:,:]
        climo_array[4,:,:] = climo_10mm[:,:]
        climo_array[5,:,:] = climo_25mm[:,:]
        climo_array[6,:,:] = climo_50mm[:,:]

        # ---- open and initialize the netCDF file we'll be writing to.

        if cleadtime_begin == '12' and cleadtime_end == '00':
            outfilename = data_directory + 'apcp_climatologies_12_to_00UTC_'+\
                cmonth+'_2002_to_2016.nc'      
        else: 
            outfilename = data_directory + 'apcp_climatologies_00_to_12UTC_'+\
                cmonth+'_2002_to_2016.nc' 
        print outfilename
        rootgrp = Dataset(outfilename,'w',format='NETCDF4_CLASSIC')
    
        xa = rootgrp.createDimension('xa',nia)
        xav = rootgrp.createVariable('xa','i4',('xa',))
        xav.long_name = "analysis grid eastward distance from southwest corner of domain in grid points" 
        xav.units = "grid index (dimensionless)" 

        ya = rootgrp.createDimension('ya',nja)
        yav = rootgrp.createVariable('ya','i4',('ya',))
        yav.long_name = "analysis grid northward distance from southwest corner of domain in grid points" 
        yav.units = "grid index (dimensionless)"

        thra = rootgrp.createDimension('thra',nthreshes)
        thrav = rootgrp.createVariable('thra','f4',('thra',))
        thrav.long_name = "index for thresholds of precipitation amount" 
        thrav.units = "mm / 12h" 
 
        pthreshes = rootgrp.createVariable('pthreshes','f4',('thra',))
        pthreshes.long_name = "precip amounts where climatological probabilities calc'd" 
        pthreshes.units = "mm / 12h" 

        lonsa = rootgrp.createVariable('lonsa','f4',('ya','xa',))
        lonsa.long_name = "longitude" 
        lonsa.units = "degrees_east" 

        latsa = rootgrp.createVariable('latsa','f4',('ya','xa',))
        latsa.long_name = "latitude" 
        latsa.units = "degrees_north" 

        conusmask = rootgrp.createVariable('conusmask','i2',('ya','xa',))
        conusmask.long_name = "mask for grid points inside CONUS (1=yes,0=no)"
        conusmask.units=""

        climo_probs = rootgrp.createVariable('climo_prob','f4',('thra','ya','xa',),
            zlib=True,least_significant_digit=4)  
        climo_probs.units = "" 
        climo_probs.long_name = "Climatological probability in accumulation period" 
        climo_probs.valid_range = [0.0,1.0]
        climo_probs.missing_value = -99.99

        rootgrp.latcorners = \
            [lats_anal[0,0], lats_anal[0,-1], lats_anal[-1,0], lats_anal[-1,-1]]
        rootgrp.loncorners = \
            [lons_anal[0,0], lons_anal[0,-1], lons_anal[-1,0], lons_anal[-1,-1]]

        rootgrp.stream = "s4" # ????
        rootgrp.title = \
            "12-h climatological probabilities of exceeding various precip thresholds"
        rootgrp.Conventions = "CF-1.0"  # ????
        rootgrp.history = "Created 5 May 2017 by Hamill" 
        rootgrp.institution = "ESRL/PSD using CCPA data from NCEP/EMC"
        rootgrp.platform = "Model" 
        rootgrp.references = "n/a"

        # --- now write data to records

        climo_probs[:,:,:] = climo_array[:,:,:]
        pthreshes[:] = pthresh[:]
        xav[:]   = np.arange(nia)
        yav[:]   = np.arange(nja)
        lonsa[:]  = lons_anal
        latsa[:]  = lats_anal
        conusmask[:] = conusmask_in

        rootgrp.close()
print 'done writing'

