from netCDF4 import Dataset
import matplotlib
matplotlib.use('TKAgg')
import matplotlib.pyplot as plt
from matplotlib import rcParams
from mpl_toolkits.basemap import Basemap, addcyclic
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
from numpy import ma
import math
import os, sys


rcParams['xtick.labelsize']='x-small'
rcParams['ytick.labelsize']='x-small'
rcParams['contour.negative_linestyle']='solid'
rcParams['path.simplify'] = False
rcParams['legend.fontsize']='small'
rcParams['legend.fancybox']=True

cmonths = ['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']

def add_analog_locns(ycoord,xcoord,ylocs,xlocs,lons,lats,map,symbolshape,color):
    istat=0
    nsize = len(ylocs)
    print 'ycoord = ',ycoord
    print 'xcoord = ',xcoord
    print 'ylocs = ',  ylocs
    print 'xlocs = ',  xlocs
    print 'nsize = ',nsize
    print 'np.shape(ylocs) = ',np.shape(ylocs)
    print 'np.shape(lons), lats', np.shape(lons), np.shape(lats)
    #print '**** xcoord, ycoord = ',xcoord-1,ycoord-1
    # --- first plot the original point in a larger font
    x,y = map(lons[ycoord,xcoord],lats[ycoord,xcoord])
    map.scatter(x,y,s=25,color=color,marker=symbolshape,zorder=12)
    for i in range(nsize):
        #print 'xlocs, ylocs, i = ',xlocs[i]-1,ylocs[i]-1,i
        x,y = map(lons[ylocs[i]-1,xlocs[i]-1],lats[ylocs[i]-1,xlocs[i]-1])
        if i < 20:
            map.scatter(x,y,s=10,color=color,marker=symbolshape,alpha=1)
        elif i >= 20 and i <=40:
            map.scatter(x,y,s=10,color='DimGray',marker=symbolshape,zorder=12,edgecolors='face',alpha=1)
        elif i >= 40 and i <=60:
            map.scatter(x,y,s=10,color='DarkGray',marker=symbolshape,zorder=12,edgecolors='face',alpha=1)
        else:
            map.scatter(x,y,s=10,color='Gainsboro',marker=symbolshape,zorder=12,edgecolors='face',alpha=1)
    return istat

def find_nearest_latlon(lonin, latin, rlon, rlat):
    istat = 0
    ilon = -99
    ilat = -99
    dmin = 999999.
    ny, nx = np.shape(rlon)
    for i in range(nx):
        for j in range(ny):
            distx = np.abs(lonin - rlon[j,i])
            disty = np.abs(latin - rlat[j,i])
            dist = np.sqrt(distx**2 + disty**2)
            if dist < dmin:
                dmin = dist
                ilon = i
                jlat = j
    if dmin >= 999999: 
        istat = -1
    else:
        print 'nearest lat, lon index = ',jlat, ilon
        print 'rlon, rlat = ',rlon[jlat,ilon], rlat[jlat,ilon]
    return istat, ilon, jlat

cmonths = ['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']
cmonth = sys.argv[1] # Jan
imonth = cmonths.index(cmonth)
imop1 = imonth+1
if imop1 < 10:
    cmop1 = '0'+str(imop1)
else:
    cmop1= str(imop1)
cleadb = sys.argv[2]
cleade = sys.argv[3]

# ---- make a dictionary of cities and their lon/lat

cities = {'POR':[-122.6, 45.54], 'BER':[-122.0, 37.95], 'BOZ':[-111.05, 45.67],\
    'BIS':[-101.7, 46.8], 'BOU':[-105.23, 40.02], 'LUB':[-101.88, 33.59],\
    'PHX':[-112.1, 33.6], 'LIT':[-92.23, 34.72], 'ATL': [-84.42, 33.76],\
    'CHI':[-87.73, 41.83],'CIN':[-84.54, 39.13], 'NYC': [-73.84, 40.7], \
    'OMA':[-96.00, 41.25],'ORL':[-81.37, 28.53], 'SAT':[-98.50, 29.41],\
    'BWI':[-76.67, 39.17]}

# --- read in the supplemental locations (huge file!)

data_directory = '/Projects/Reforecast2/netcdf/NationalBlend/'
infile = data_directory + 'supplemental_locations_CONUS_ndfd2p5_'+cleadb+'_to_'+cleade+'_'+cmonth+'.nc'
print 'reading supplemental locations file ',infile

ncfile = Dataset(infile,'r',format='NETCDF4_CLASSIC')
xlocation = ncfile.variables['xlocation_supp'][:,:,:]
ylocation = ncfile.variables['ylocation_supp'][:,:,:]
conusmask = ncfile.variables['conusmask'][:,:]
rlon = ncfile.variables['lons'][:,:]
rlat = ncfile.variables['lats'][:,:]
rlon = rlon - 360.
nsupp, ny_2p5, nx2p5 = np.shape(xlocation)
ncfile.close()
print 'min, max rlon = ', np.min(rlon), np.max(rlon)

#  ----  read in the fraction of samples with climatological zero precipitation

infile = data_directory + \
    'climatology_gamma_parameters_ndfd2p5_'+cleadb+'_to_'+cleade+'_'+cmonth+'.nc'
print 'reading from ', infile
nc = Dataset(infile)
fraction_zero = nc.variables['fraction_zero'][:]
nc.close()

# ---- get nearest NDFD 2.5-km gridpoint to input lon, lat.  Then same for CCPA

with PdfPages('supp_locations_'+cleadb+'_to_'+cleade+'_'+cmonth+'.pdf') as pdf:

    for city in cities:
        rlonin, rlatin = cities[city]
        print 'rlonin, rlatin = ',rlonin, rlatin
        istat, iloc, jloc = find_nearest_latlon(rlonin, rlatin, rlon, rlat)
        print city,' jloc, iloc = ',jloc, iloc
        print 'i,   iloc,   jloc,    lon(2.5),   lat(2.5) '
        for i in range(50):
            ix = xlocation[i,jloc,iloc]
            jy = ylocation[i,jloc,iloc]
            print i,ix,jy,rlon[jy,ix],rlat[jy,ix]

        colorst = ['White','#ECFFFF','#D9F7FF','#C4E8FF','#E8FBE8','#C7F4C7','#92F592','Yellow','Gold',\
            'Orange','#FFB2B2','#EC5B71','Red','Magenta','DarkOrchid','White']
        colorstblack=['White','Black','Black','Black','Black', 'Black','Black','Black',\
          'Black','Black','Black','Black','Black','Black','Black','Black']

        # --- Plot fraction nonzero and NDFD supplemental locations 

        clevs = [.03,.05,.07,.10,.15,.20,.25,.30,.35,.40,.50,.60,.70,.80,.90]
        symbolshapes = ['o']

        #  ---- make map plot of supplemental locations with fraction nonzero of climo superimposed

        fig1 = plt.figure(figsize=(9.,5.9))
        ax = fig1.add_axes([0.03,.12,0.94,.83])
        ctitle = cmonth+' NDFD supp. locations and fraction of climatology with nonzero precipitation'
        ax.set_title(ctitle,fontsize=12)

        m = Basemap(projection='mill',llcrnrlon=-125,llcrnrlat=22.,\
            urcrnrlon=-65.,urcrnrlat=53,resolution='l')
    
        x,y = m(rlon,rlat)
        CS1 = m.contourf(x,y,1.-fraction_zero,clevs,colors=colorst,cmap=None,extend='neither')
        m.drawcoastlines(linewidth=.5)
        m.drawstates(linewidth=.5)
        m.drawcountries(linewidth=.5)

        istat = add_analog_locns(ycoord=jloc,xcoord=iloc,\
            ylocs=ylocation[:,jloc,iloc],xlocs=xlocation[:,jloc,iloc],\
            lons=rlon,lats=rlat,map=m,symbolshape='o',color='Black') 

        cax = fig1.add_axes([0.1,0.07,0.8,0.03])
        cbar = fig1.colorbar(CS1,extend='neither', orientation='horizontal',\
            cax=cax,drawedges=True,ticks=clevs,format='%g')
        cax.set_xlabel('Fraction of climatological samples with nonzero'+\
            ' 12-h accumulated precipitation')

        pdf.savefig()
        plt.close()  
          
print 'saving plot to file = ','supp_locations_'+cleadb+'_to_'+cleade+'_'+cmonth+'.pdf'
print 'Plot done' 










