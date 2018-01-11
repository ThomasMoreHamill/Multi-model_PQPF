"""
this python script controls the generation of empirical CDFs, whereby
the empirical CDF information for precipitation for quantile mapping 
is generated at a large number of precipitation amount thresholds.  Since
these files are large, rather than generating one new one with the 
latest 60 days of forecast & analyzed data, instead we generate one 
only every 10th day roughly.

"""
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

from get_cdf_precip_anal import get_cdf_precip_anal
from get_cdf_cmc import get_cdf_cmc
from get_cdf_ecmwf_ncep import get_cdf_ecmwf_ncep
from get_quantiles_linear import get_quantiles_linear
from get_quantiles_linear_cmc import get_quantiles_linear_cmc

# --- queries from command line

center = sys.argv[1]
cleade = sys.argv[2]
date_end = sys.argv[3]
idaysbefore = -1 -int(int(cleade)/24)
ihoursbefore = idaysbefore*24
date_end_shift = dateshift(date_end, ihoursbefore)
date_begin_shift = dateshift(date_end_shift, -24*61)

# --- initialize stuff

date_list = daterange(date_begin_shift, date_end_shift, 24)
print 'date_list = ', date_list
ndates = len(date_list)
date_middle = date_list[ndates/2]
cmonthno = date_middle[4:6]
imo = int(cmonthno)-1
cmonths = ['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']
cmonth = cmonths[imo]
if center == 'ECMWF':
    nmembers = 50
else:
    nmembers = 20
    
nmembers2process = 5  # for ecmwf, ncep only process 1st five mbrs

# ---- define the precipitation amount thresholds that we will calculate CDF at.

npct = 99 + 8 # quantiles, from 0.01 to 0.99 by 0.01, also (.0001, .005, .001, .0005), 
              # and 1 - (.0001, .005, .001, .0005)

xc = range(90)
thresh = [.001,.003,.005,.01,.03, .05,.07,.1,.2,.3,  .4,.5,.6,.7,.8,  \
      .9,1.0,1.2, 1.4, 1.6,    1.8, 2.0, 2.25, 2.5, 2.75,   3.0, 3.5, 4.0, 4.5, 5.0, \
      6.0,7.0,8.0,9.0,10.0,  11.0,12.0,13.0,14.0,15.0,  16.0,17.0,18.0,19.0,19.5,  \
      20.0,22.5,25.,27.5,30.0,   32.5,35.0,37.5,40.0,42.5,  45.0,50.0,55.0,60.0,65.0,  \
      70.0,75.0,80.0,85.0,90.0,   95.0,100.0,105.0,110.0,120.0,  130.0,140.0,150.0,160.0,170.0,  \
      180.0,190.0,200.0,220.0,240.0,   260.0,280.0,300.0,325.0,350.0,   400.0,500.0,600.0,700.0,1000.]
pctvalues = [.0001, .0005, .001, .005, \
          .01, .02, .03, .04, .05, .06, .07, .08, .09, .10, \
          .11, .12, .13, .14, .15, .16, .17, .18, .19, .20, \
          .21, .22, .23, .24, .25, .26, .27, .28, .29, .30, \
          .31, .32, .33, .34, .35, .36, .37, .38, .39, .40, \
          .41, .42, .43, .44, .45, .46, .47, .48, .49, .50, \
          .51, .52, .53, .54, .55, .56, .57, .58, .59, .60, \
          .61, .62, .63, .64, .65, .66, .67, .68, .69, .70, \
          .71, .72, .73, .74, .75, .76, .77, .78, .79, .80, \
          .81, .82, .83, .84, .85, .86, .87, .88, .89, .90, \
          .91, .92, .93, .94, .95, .96, .97, .98, .99, \
          .995, .999, .9995, .9999]
nthresh = len(xc)
print 'len xc, thresh = ',len(xc), len(thresh)

# ---- read in the supplemental locations information for the CCPA grid.

infile = '/Users/thamill/precip/supplemental_locations_eighth_degree_'+cmonth+'_v9.nc'
print infile
nc = Dataset(infile)
conusmask_lg_in = nc.variables['conusmask'][:,:]
lons_lg = nc.variables['longitudes'][:,:]
lats_lg = nc.variables['latitudes'][:,:]
xlocations_lg = nc.variables['xlocations'][:,:,:]
ylocations_lg = nc.variables['ylocations'][:,:,:]
nsupplemental_lg = nc.variables['nsupplemental'][:,:]
nsuppmax = 100
nsupp, nya_lg, nxa_lg = np.shape(xlocations_lg)
nc.close()

# ---- extract subset on smaller CCPA grid.

nxa = 464
nya = 224
xoffset = 21  
yoffset = 22  
xlocations = np.zeros((nsupp,nya,nxa), dtype=np.int32)
ylocations = np.zeros((nsupp,nya,nxa), dtype=np.int32)
nsupplemental = np.zeros((nya,nxa), dtype=np.int16)
conusmask_in = np.zeros((nya,nxa), dtype=np.int16)
lons_anal = np.zeros((nya,nxa), dtype=np.float32)
lats_anal = np.zeros((nya,nxa), dtype=np.float32)

xlocations[:,:,:] = xlocations_lg[:,yoffset:yoffset+nya,xoffset:xoffset+nxa]
ylocations[:,:,:] = ylocations_lg[:,yoffset:yoffset+nya,xoffset:xoffset+nxa]
nsupplemental[:,:] = nsupplemental_lg[yoffset:yoffset+nya,xoffset:xoffset+nxa]
conusmask_in[:,:] = conusmask_lg_in[yoffset:yoffset+nya,xoffset:xoffset+nxa]
lons_anal[:,:] = lons_lg[yoffset:yoffset+nya,xoffset:xoffset+nxa]
lats_anal[:,:] = lats_lg[yoffset:yoffset+nya,xoffset:xoffset+nxa]

# ---- set up the CDF arrays

CDFa = np.zeros((nthresh,nya,nxa), dtype=np.float64) # analyzed
CDFworka = np.zeros((nthresh,nya,nxa), dtype=np.float64)
if center == 'CMC':
    CDFf = np.zeros((nthresh,nmembers,nya,nxa), dtype=np.float64) # forecast
    CDFworkf = np.zeros((nthresh,nmembers,nya,nxa), dtype=np.float64)
else:
    CDFf = np.zeros((nthresh,nya,nxa), dtype=np.float64) # forecast
    CDFworkf = np.zeros((nthresh,nya,nxa), dtype=np.float64)
    
# --- loop thru forecasts and analyses for each date in date_list    
    
icdf_a = np.zeros((nya,nxa), dtype=np.float32)
icdf_f = np.zeros((nya,nxa), dtype=np.float32)
icounta = np.zeros((nya,nxa), dtype=np.int32)
icountf = np.zeros((nya,nxa), dtype=np.int32)

for idate, date in zip(range(ndates), date_list):
    
    # ---- read in ccpa precip analysis data for this date (12-h accum)

    print 'processing date = ', date
    if idate == 0: 
        infilename = '/data/thamill/Rf2_tests/ccpa_v1/precip_ccpav1_'+\
            '2002010200_to_2016123100.nc'            
        print infilename
        afile1 = Dataset(infilename,"r")
        yyyymmddhh = afile1.variables['yyyymmddhh'][:]
    
    fdate_today = int(dateshift(date, int(cleade)))
    idx_today = int(np.where(yyyymmddhh == fdate_today)[0])
        
    if idx_today >= 0:
        apcp_anal = afile1.variables['apcp_anal'][idx_today,:,:]
    else:
        apcp_anal = -99.99*np.ones((nya, nxa), dtype=np.float32)
            
    # ---- read in the forecast data for this date, and for the date 12 h previous
    #      then subtract to get accumulated precip during this period
    
    infile1 = '/Users/thamill/precip/ecmwf_data/'+center+'_'+date+'_leadtime'+cleade+'h.nc'
    cleadb = str(int(cleade)-12)
    infile2 = '/Users/thamill/precip/ecmwf_data/'+center+'_'+date+'_leadtime'+cleadb+'h.nc'
    fexist1  = os.path.exists(infile1)
    fexist2  = os.path.exists(infile1)
    
    if fexist1 and fexist2:
        try:
            print infile1
            nc = Dataset(infile1)
            apcp_fcst_ens_late = nc.variables['apcp_fcst_ens'][:,:,:] # mbr,x,y
            print infile2
            nc.close()
            nc = Dataset(infile2)
            apcp_fcst_ens_early = nc.variables['apcp_fcst_ens'][:,:,:] # mbr,x,y
            nc.close()
            apcp_fcst_ens = apcp_fcst_ens_late - apcp_fcst_ens_early
            pmax = np.max(apcp_fcst_ens)
        except (IOError, ValueError, RuntimeError):
            print 'Error reading either ', infile1
            print '  or ', infile2
            pmax = -99.99
    else:
        pmax = -99.99

    # ---- populate work arrays with CDF data for this date.

    if pmax > 0:
        
        print 'get_cdf_precip_anal '
        CDFworka, icounta, istata = get_cdf_precip_anal(nthresh, nsuppmax, nya, nxa, \
            yoffset, xoffset, thresh, apcp_anal, xlocations, ylocations, \
            nsupplemental, conusmask_in)
    
        if center == 'CMC':
            print 'get_cdf_precip_cmc '
            CDFworkf, icountf, istatf = get_cdf_cmc(nthresh, nmembers, nya, nxa, nsuppmax, \
                yoffset, xoffset, thresh, apcp_fcst_ens, xlocations, ylocations, \
                nsupplemental, conusmask_in)
        
        else:
            print 'get_cdf_precip_ecmwf_ncep '
            CDFworkf, icountf, istatf = get_cdf_ecmwf_ncep(nthresh, nmembers, nya, nxa, nsuppmax, \
                yoffset, xoffset, nmembers2process, thresh, apcp_fcst_ens, xlocations, ylocations, \
                nsupplemental, conusmask_in)     
    
        # ---- if both forecast and analysis data exist for this date, increment CDF arrays

        print 'adding to cdfs, istata, istatf ',istata,istatf
        if istata == 1 and istatf == 1 :
            CDFa = CDFa + CDFworka
            icdf_a = icdf_a + icounta
            CDFf = CDFf + CDFworkf
            icdf_f = icdf_f + icountf

# ---- divide thru by the number of samples.  keep track of cdfs for each member
#      in the Canadian system, since that doesn't have members with exchangeable
#      error statistics.

ones = np.ones((nya,nxa),dtype=np.float32)
icdf_a = np.where(icdf_a == 0, ones, icdf_a)
icdf_f = np.where(icdf_f == 0, ones, icdf_f) # make sure not divide by zero

for it in range(nthresh):
    CDFa[it,:,:] = CDFa[it,:,:] / icdf_a[:,:]
    if center == 'CMC':
        for imem in range(nmembers):
            CDFf[it,imem,:,:] = CDFf[it,imem,:,:] / icdf_f[:,:]
    else:
        CDFf[it,:,:] = CDFf[it,:,:] / icdf_f[:,:]
    
# --- get quantiles associated with the precipitation amounts, both for the forecast and analyzed

precip_qret_a = np.zeros((npct,nya,nxa),dtype=np.float)
precip_qret_a, istat = get_quantiles_linear(nthresh,npct,nya,nxa,pctvalues,thresh,CDFa)

if center != 'CMC':
    precip_qret_f = np.zeros((npct,nya,nxa),dtype=np.float)
    precip_qret_f, istat = get_quantiles_linear(nthresh,npct,nya,nxa,pctvalues,thresh,CDFf)
else:
    precip_qret_f = np.zeros((npct,nmembers,nya,nxa),dtype=np.float)
    precip_qret_f, istat = get_quantiles_linear_cmc(nthresh,npct,nmembers,nya,nxa,pctvalues,thresh,CDFf)

# ---- open and initialize the netCDF file we'll be writing to.

outfilename = '/Users/thamill/precip/ecmwf_data/'+center+'_CDF_flead'+cleade+\
    '_'+date_begin_shift+'_to_'+date_end_shift+'.nc'
print outfilename
rootgrp = Dataset(outfilename,'w',format='NETCDF4_CLASSIC')
    
xa = rootgrp.createDimension('xa',nxa)
xav = rootgrp.createVariable('xa','f4',('xa',))
xav.long_name = "analysis grid eastward distance from southwest corner of domain in grid points" 
xav.units = "grid index (dimensionless)" 
    
ya = rootgrp.createDimension('ya',nya)
yav = rootgrp.createVariable('ya','f4',('ya',))
yav.long_name = "analysis grid northward distance from southwest corner of domain in grid points" 
yav.units = "grid index (dimensionless)"
    
pct = rootgrp.createDimension('pct',npct)
pctv = rootgrp.createVariable('pct','f4',('pct'))
pctv.long_name = "quantiles of the distribution"
pctv.units = "fraction"

thrnum = rootgrp.createDimension('thrnum',nthresh)
thrnumv = rootgrp.createVariable('thrnum','i4',('thrnum',))
thrnumv.long_name = "Threshold iterator (0:nthresh)" 
thrnumv.units = " "  
    
thrval = rootgrp.createDimension('thrval',nthresh)
thrvalv = rootgrp.createVariable('thrval','f4',('thrval',))
thrvalv.long_name = "Precip thresholds (mm) that precip_CDFs evaluated at" 
thrvalv.units = "K" 
    
lonsa = rootgrp.createVariable('lonsa','f4',('ya','xa',))
lonsa.long_name = "longitude" 
lonsa.units = "degrees_east" 
    
latsa = rootgrp.createVariable('latsa','f4',('ya','xa',))
latsa.long_name = "latitude" 
latsa.units = "degrees_north" 

conusmask = rootgrp.createVariable('conusmask','i2',('ya','xa',))
conusmask.long_name = "mask for grid points inside CONUS (1=yes,0=no)"
conusmask.units=""

panal_CDF = rootgrp.createVariable('panal_CDF','f4',('thrnum','ya','xa',),
    zlib=True,least_significant_digit=3)  
panal_CDF.units = "" 
panal_CDF.long_name = "Cumulative distribution function of analyzed precip, defined at thrval" 
panal_CDF.valid_range = [0.0,1.0]
panal_CDF.missing_value = -99.99

panal_quantiles = rootgrp.createVariable('panal_quantiles','f4',('pct','ya','xa',),
    zlib=True,least_significant_digit=2)  
panal_quantiles.units = "mm" 
panal_quantiles.long_name = "Precip amount associated with quantiles of analyzed distribution" 
panal_quantiles.valid_range = [0.0,10000.]
panal_quantiles.missing_value = -99.99
  
# ----  set for CMC differently
  
if center == 'CMC': 
    
    ens = rootgrp.createDimension('ens',nmembers)
    ensv = rootgrp.createVariable('ensv','i4',('ens',))
    ensv.long_name = "Ensemble member number" 
    ensv.units = " " 

    pfcst_CDF = rootgrp.createVariable('pfcst_CDF','f4',('thrnum','ens','ya','xa',),
        zlib=True,least_significant_digit=3)  
    pfcst_CDF.units = "" 
    pfcst_CDF.long_name = "Cumulative distribution function of forecast precip, defined at thrval" 
    pfcst_CDF.valid_range = [0.0,1.0]
    pfcst_CDF.missing_value = -99.99
    
    pfcst_quantiles = rootgrp.createVariable('pfcst_quantiles','f4',('pct','ens','ya','xa',),
        zlib=True,least_significant_digit=2)  
    pfcst_quantiles.units = "mm" 
    pfcst_quantiles.long_name = "Precip amount associated with quantiles of forecast distribution" 
    pfcst_quantiles.valid_range = [0.0,10000.]
    pfcst_quantiles.missing_value = -99.99
    
else:
    
    pfcst_CDF = rootgrp.createVariable('pfcst_CDF','f4',('thrnum','ya','xa',),
        zlib=True,least_significant_digit=3)  
    pfcst_CDF.units = "" 
    pfcst_CDF.long_name = "Cumulative distribution function of forecast precip, defined at thrval" 
    pfcst_CDF.valid_range = [0.0,1.0]
    pfcst_CDF.missing_value = -99.99
    
    pfcst_quantiles = rootgrp.createVariable('pfcst_quantiles','f4',('pct','ya','xa',),
        zlib=True,least_significant_digit=2)  
    pfcst_quantiles.units = "mm" 
    pfcst_quantiles.long_name = "Precip amount associated with quantiles of forecast distribution" 
    pfcst_quantiles.valid_range = [0.0,10000.]
    pfcst_quantiles.missing_value = -99.99
    
rootgrp.latcorners = [lats_anal[0,0], lats_anal[0,-1], lats_anal[-1,0], lats_anal[-1,-1]]
rootgrp.loncorners = [lons_anal[0,0], lons_anal[0,-1], lons_anal[-1,0], lons_anal[-1,-1]]

rootgrp.stream = "s4" # ????
rootgrp.title = "anal CDF"
rootgrp.Conventions = "CF-1.0"  # ????
rootgrp.history = "Created 8 May 2017 by Tom Hamill" 
rootgrp.institution = "ESRL/PSD using CCPA data from NCEP/EMC"
rootgrp.platform = "Model" 
rootgrp.references = "" 

pctv[:] = pctvalues[:]
    
thrnumv[:] = xc[:]
thrvalv[:] = thresh[:]
if center == 'CMC':
    ensv[:] = range(nmembers)

# ---- define the corner lat/lons.

llcrnrlat = lats_anal[0,0]
llcrnrlon = lons_anal[0,0]
urcrnrlat = lats_anal[-1,-1]
urcrnrlon = lons_anal[-1,-1]

xav[:]   = np.arange(nxa)
yav[:]   = np.arange(nya)
lonsa[:]  = lons_anal
latsa[:]  = lats_anal

conusmask[:] = conusmask_in

# --- write out the CDF record

for ithresh in range(nthresh):
    panal_CDF[ithresh] = CDFa[ithresh,:,:]
        
if center == 'CMC':
    for ithresh in range(nthresh):
        pfcst_CDF[ithresh] = CDFf[ithresh,:,:,:]
else:
    for ithresh in range(nthresh):
        pfcst_CDF[ithresh] = CDFf[ithresh,:,:]
    
# --- write out the quantiles record

for ipct in range(npct):
    panal_quantiles[ipct] = precip_qret_a[ipct,:,:] 
     
if center == 'CMC':
    for ipct in range(npct):
        pfcst_quantiles[ipct] = precip_qret_f[ipct,:,:,:]
else:
    for ipct in range(npct):
        pfcst_quantiles[ipct] = precip_qret_f[ipct,:,:]

rootgrp.close()
