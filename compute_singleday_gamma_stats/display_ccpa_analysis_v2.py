import numpy as np
import numpy.ma as ma
import sys
import pygrib
import os
import time as timey
import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap, interp
from dateutils import hrstodate, daterange, dayofyear, \
     splitdate, datetohrs, dateshift, dateto_hrs_since_day1CE
from netCDF4 import Dataset
import cPickle

# --- Read input data from command line.  Get beginning and ending lead time 
#     in hours and the name of the file with the dates we want

# --- queries from command line

cyyyymmddhh = sys.argv[1]  # the initial date/time of the forecast

# ---- read in the associated analysis
    
iyyyymmddhh_anal_end = int(dateshift(cyyyymmddhh,0))
iyyyymmddhh_anal_begin = int(dateshift(cyyyymmddhh,-6))
print 'iyyyymmddh_anal_end, iyyyymmddh_anal_begin = ',\
    iyyyymmddhh_anal_end, iyyyymmddhh_anal_begin 
    
# ---- read in ccpa precip analysis data for this date 
#      (12-h accum, generated from two 6-hourly files)

infilename = '/Projects/Reforecast2/netcdf/NationalBlend/precip_analyses_ccpa_v1_'+\
            '2002010100_to_2016123100.nc'
print infilename
nc = Dataset(infilename)
yyyymmddhh_anal_list = nc.variables['yyyymmddhh_anal_end'][:]
lons_anal = nc.variables['lons_anal'][:,:]
lats_anal = nc.variables['lats_anal'][:,:]
conusmask = nc.variables['conusmask'][:,:]
nya, nxa = lons_anal.shape
            
idx_today_begin = int(np.where(yyyymmddhh_anal_list == iyyyymmddhh_anal_begin)[0])
idx_today_end = int(np.where(yyyymmddhh_anal_list == iyyyymmddhh_anal_end)[0])
print 'idx_today_begin, idx_today_end = ', idx_today_begin, idx_today_end
print 'yyyymmddhh_anal_list[idx_today_begin:idx_today_end+1] = ',\
    yyyymmddhh_anal_list[idx_today_begin:idx_today_end+1]

if idx_today_begin >= 0 and idx_today_end > 0:
    apcp_anal = nc.variables['apcp_anal'][idx_today_begin,:,:] + \
        nc.variables['apcp_anal'][idx_today_end,:,:]
else:
    apcp_anal = -99.99*np.ones((nya, nxa), dtype=np.float32)
pmax_anal = np.max(apcp_anal)
print 'pmax_anal = ',np.max(apcp_anal)
print 'apcp_anal[50,0:nxa] = ', apcp_anal[50,0:nxa]
print 'apcp_anal[nya/2,0:nxa] = ', apcp_anal[nya/2,0:nxa]
print 'apcp_anal[200,0:nxa] = ', apcp_anal[200,0:nxa]


iyyyymmddhh_end = int(cyyyymmddhh)
rthresh = [0.1,1.0,2.5,5.0,10.0,15.0,20.0,25.0,30.0,35.0,40.0] # list of event thresholds
nthresh = len(rthresh)

# ----- make plots of precipitation analyses and forecasts

colorst = ['White','#E4FFFF','#C4E8FF','#8FB3FF','#D8F9D8','#A6ECA6','#42F742','Yellow','Gold',\
    'Orange','#FCD5D9','#F6A3AE','#FA5257','Orchid','#AD8ADB','#A449FF','LightGray'] 
colorstblack=['White','Black','Black','Black','Black', 'Black','Black','Black',\
                      'Black','Black','Black','Black','Black','Black','Black','Black','Black']
colorstwhite=['White','Black','Black','White','White','White','White',\
                      'White','White','White','White','White','White','White','White','White']

fig1 = plt.figure(figsize=(9.,5.5))

# ---- plot 1: analysis data

clevs = [0.0, 0.1, 0.5, 1.0, 2.5, 5., 10., 25., 50., 100., 150., 200.]
ax = fig1.add_axes([0.01,.13,0.98,0.8])
ctitle = '(a) 12-h accumulated precipitation analysis ending '+cyyyymmddhh 
ax.set_title(ctitle,fontsize=16)
m = Basemap(projection='mill',llcrnrlon=-125.,llcrnrlat=25.,\
        urcrnrlon=-65,urcrnrlat=50.,resolution='l')  # resolution next step = i
x,y = m(lons_anal,lats_anal)
CS1 = m.contour(x,y,apcp_anal,clevs,colors=colorstblack,cmap=None,linewidths=0.3)
CS2 = m.contourf(x,y,apcp_anal,clevs,colors=colorst,cmap=None,extend='neither')
m.drawcoastlines(linewidth=.5)
m.drawstates(linewidth=.5)
m.drawcountries(linewidth=.5)

print clevs

cax = fig1.add_axes([0.01,0.10,0.98,0.03])
cbar = fig1.colorbar(CS2,extend='neither', \
   orientation='horizontal',cax=cax,drawedges=True,ticks=clevs,format='%g')
cax.set_xlabel('Analyzed precipitation amount (mm)',fontsize=14)

# ---- set plot title

plot_title = 'apcp_anal_'+cyyyymmddhh+'.pdf'
fig1.savefig(plot_title)
print 'saving plot to file = ',plot_title
#plt.show()
