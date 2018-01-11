""" this routine will: (1) generate fitted gamma distribution parameters 
    for best-member dressing for individual forecast systems.  These
    dressing parameters are saved to a netcdf file, and diagnostic plots 
    are generated: (2) generate closest-member histograms, make diagnostic
    plots, and save to a netCDF file.
"""

from mpl_toolkits.basemap import Basemap
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from numpy import ma
import os, sys
from netCDF4 import Dataset
import pygrib
from dateutils import daterange, dateshift
from scipy.optimize import curve_fit
from matplotlib import rcParams
import scipy.signal as signal
import scipy.stats as stats
from scipy.stats import gamma
from matplotlib.backends.backend_pdf import PdfPages
rcParams['legend.fontsize']='x-small'
rcParams['legend.fancybox']=True
rcParams['xtick.labelsize']='medium'
rcParams['ytick.labelsize']='small'


# --- COMMAND LINE INPUT

cyyyymmddhh = sys.argv[1]
cleade = sys.argv[2] # ending of periods lead time in hours, e.g., '24', '108'
cempirical = sys.argv[3]
cleadb = str(int(cleade) - 12)
ileade = int(cleade)


# ======================================================================
# Part 1:  reading in dressing parameter data from netCDF files.
# ======================================================================

data_directory = '/Users/thamill/precip/ecmwf_data/NCEP/'   
if cempirical == '1':
    infile = data_directory+'gamma_fraction_zero_dressing_'+\
        'NCEP_date='+cyyyymmddhh+'_lead='+cleade+'.nc'
else:
    infile = data_directory+'gamma_fraction_zero_dressing_'+\
        'NCEP_date='+cyyyymmddhh+'_lead='+cleade+'_gammaqmap.nc'
print infile
nc = Dataset(infile)    
precip_values = nc.variables['precip_values'][:]  
fraction_zeros_NCEP = nc.variables['fraction_zeros'][:,:,:]
fraction_zeros_fclimpop_NCEP = nc.variables['fraction_zeros_fclimpop'][:]
climo_pop_thresholds = nc.variables['climo_pop_thresholds'][:] 
gamma_shapes_NCEP = nc.variables['gamma_shapes'][:,:,:]
gamma_scales_NCEP = nc.variables['gamma_scales'][:,:,:]
gamma_shape_fclimpop_NCEP = nc.variables['gamma_shape_fclimpop'][:] 
gamma_scale_fclimpop_NCEP = nc.variables['gamma_scale_fclimpop'][:] 
nsamps_NCEP = nc.variables['nsamps'][:,:,:]
nsamps_fclimpop_NCEP = nc.variables['nsamps_fclimpop'][:]
precip_thresholds_NCEP = nc.variables['precip_thresholds'][:] 
nc.close()

data_directory = '/Users/thamill/precip/ecmwf_data/CMC/'   
if cempirical == '1':
    infile = data_directory+'gamma_fraction_zero_dressing_'+\
        'CMC_date='+cyyyymmddhh+'_lead='+cleade+'.nc'
else:
    infile = data_directory+'gamma_fraction_zero_dressing_'+\
        'CMC_date='+cyyyymmddhh+'_lead='+cleade+'_gammaqmap.nc'
print infile
nc = Dataset(infile)    
precip_values = nc.variables['precip_values'][:]  
fraction_zeros_CMC = nc.variables['fraction_zeros'][:,:,:]
fraction_zeros_fclimpop_CMC = nc.variables['fraction_zeros_fclimpop'][:]
gamma_shapes_CMC = nc.variables['gamma_shapes'][:,:,:]
gamma_scales_CMC = nc.variables['gamma_scales'][:,:,:]
gamma_shape_fclimpop_CMC = nc.variables['gamma_shape_fclimpop'][:] 
gamma_scale_fclimpop_CMC = nc.variables['gamma_scale_fclimpop'][:] 
nsamps_CMC = nc.variables['nsamps'][:,:,:]
nsamps_fclimpop_CMC = nc.variables['nsamps_fclimpop'][:]
precip_thresholds_CMC = nc.variables['precip_thresholds'][:] 
nc.close()

data_directory = '/Users/thamill/precip/ecmwf_data/ECMWF/'   
if cempirical == '1':
    infile = data_directory+'gamma_fraction_zero_dressing_'+\
        'ECMWF_date='+cyyyymmddhh+'_lead='+cleade+'.nc'
else:
    infile = data_directory+'gamma_fraction_zero_dressing_'+\
        'ECMWF_date='+cyyyymmddhh+'_lead='+cleade+'_gammaqmap.nc'
print infile
nc = Dataset(infile)    
precip_values = nc.variables['precip_values'][:]  
fraction_zeros_ECMWF= nc.variables['fraction_zeros'][:,:,:]
fraction_zeros_fclimpop_ECMWF = nc.variables['fraction_zeros_fclimpop'][:]
gamma_shapes_ECMWF = nc.variables['gamma_shapes'][:,:,:]
gamma_scales_ECMWF = nc.variables['gamma_scales'][:,:,:]
gamma_shape_fclimpop_ECMWF = nc.variables['gamma_shape_fclimpop'][:] 
gamma_scale_fclimpop_ECMWF = nc.variables['gamma_scale_fclimpop'][:] 
nsamps_ECMWF = nc.variables['nsamps'][:,:,:]
nsamps_fclimpop_ECMWF = nc.variables['nsamps_fclimpop'][:]
precip_thresholds_ECMWF = nc.variables['precip_thresholds'][:] 
nc.close()

# ---- generate gamma pdfs based on shape and scale parameters.

# ======================================================================
# ---- PART 2: make some diagnostic plots, save these to pdf file
# ======================================================================

#  first plot for the fraction zero, the gamma shape and scale for the empirical
#  for the various categories of ensemble-mean precip, and then as a function
#  of the best-member amount.

# ---- plot the Fraction Zero and fitted Gamma distributions for nonzero amounts 
#      for the lowest ranked member

# precip_values = 0, 0.05, 0.1, 0.2, 0.3, \
#    0.4, 0.5, 0.6, 0.7, 0.8, \
#    0.9, 1, 1.2, 1.4,  1.6, 
#    1.8, 2, 2.3, 2.6, 3, 
#    3.3, 3.6, 4, 4.5, 5, 5.5, 
#    6, 6.5, 7, 7.5, 8, 
#    8.5, 9, 10, 11, 12, 
#    13, 14, 15, 16, 18, 
#    20, 22.5, 25, 27.5, 
#    30, 33, 36, 40, 45, 
#    50, 55, 60, 65, 70, 
#    75, 80, 85, 90, 95, 
#    100, 120, 140, 160, 180, 
#    200, 250, 300 ;


fig = plt.figure(figsize=(6.5,9.))

fig.suptitle('Dressing distributions, +'+cleadb+'-'+cleade+' h, lowest-ranked member', fontsize=15)
        
nplots = 4  
colors = ['Red','RoyalBlue','LimeGreen','Purple']     
for ipanel in range(12):
    if ipanel == 0:
        ctitle = r'(a) $M(\bar{\tilde{x}}^{f}) = 1$'+'\n[0.1 - 2 mm)'
        ytitle = 'Fraction zero'
        xtitle = 'Precipitation amount (mm)'
        yplotdata_N = fraction_zeros_NCEP[0,0,:] # 0 = lowest member, 0 = light precip, : = precip_thresholds
        yplotdata_C = fraction_zeros_CMC[0,0,:] # 0 = lowest member, 0 = light precip, : = precip_thresholds
        yplotdata_E = fraction_zeros_ECMWF[0,0,:] # 0 = lowest member, 0 = light precip, : = precip_thresholds
        print 'yplotdata_N = ',yplotdata_N[0:20]
        print 'yplotdata_C = ',yplotdata_C[0:20]
        print 'yplotdata_E = ',yplotdata_E[0:20]
        xplotdata = precip_values 
        print 'xplotdata = ', xplotdata
        #sys.exit()
        ylimits = [-0.01,1.01]
        xlimits = [-0.01,1.01]
        axlocn = [0.1,0.74,0.25,0.16]
        xticks = [0.0, 0.2, 0.4, 0.6, 0.8, 1.0]
    elif ipanel == 1:
        ctitle = r'(b) $M(\bar{\tilde{x}}^{f}) = 2$'+'\n[2 - 10 mm)'
        ytitle = ' '
        xtitle = 'Precipitation amount (mm)'
        yplotdata_N = fraction_zeros_NCEP[0,1,:] # 0 = lowest member, 1 = mod precip, : = precip_thresholds
        yplotdata_C = fraction_zeros_CMC[0,1,:] # 0 = lowest member, 1 = mod precip, : = precip_thresholds
        yplotdata_E = fraction_zeros_ECMWF[0,1,:] # 0 = lowest member, 1 = mod precip, : = precip_thresholds
        xplotdata = precip_values 
        ylimits = [-0.01,1.01]
        xlimits = [0.,10.]
        axlocn = [0.415,0.74,0.25,0.16]
        xticks = [0, 2, 4, 6, 8, 10]
    elif ipanel == 2:
        ctitle = r'(c) $M(\bar{\tilde{x}}^{f}) = 3$'+'\n[> 10 mm)'
        ytitle = ' '
        xtitle = 'Precipitation amount (mm)'
        yplotdata_N = fraction_zeros_NCEP[0,2,:] # 0 = lowest member, 2 = heavy precip, : = precip_thresholds
        yplotdata_C = fraction_zeros_CMC[0,2,:] # 0 = lowest member, 2 = heavy precip, : = precip_thresholds
        yplotdata_E = fraction_zeros_ECMWF[0,2,:] # 0 = lowest member, 2 = heavy precip, : = precip_thresholds
        xplotdata = precip_values 
        ylimits = [-0.01,1.01]
        xlimits = [0.,30.]
        axlocn = [0.72,0.74,0.25,0.16]
        xticks = [0,5,10,15,20,25,30]
    elif ipanel == 3:
        ctitle = '(d) NCEP [0.1 - 2 mm)'
        ytitle = 'Probability density'
        xtitle = 'Precipitation amount (mm)'
        gamma_shape = gamma_shapes_NCEP[0,0,:] # 0 = lowest member, 0 = light precip, : = precip_thresholds
        gamma_scale = gamma_scales_NCEP[0,0,:] # 0 = lowest member, 0 = light precip, : = precip_thresholds
        axlocn = [0.1,0.51,0.25,0.16]
        plot_at_these_gamma_threshes = [0.1, 0.2, 0.3, 0.4]
        index_gamma = [2,3,4,5]
        #plot_at_these_gamma_threshes = [0.2, 0.4, 0.6, 0.8]
        #index_gamma = [3,5,7,9]
        xvalues = np.arange(0.0,1.0,0.01)
        xlimits = [0.,1.]
        ylimits = [0.,25.]
        xticks = [0.0, 0.2, 0.4, 0.6, 0.8, 1.0]
    elif ipanel == 4:
        ctitle = '(e) NCEP [2 - 10 mm)'
        ytitle = ''
        xtitle = 'Precipitation amount (mm)'
        gamma_shape = gamma_shapes_NCEP[0,1,:] # 0 = lowest member, 1 = mod precip, : = precip_thresholds
        gamma_scale = gamma_scales_NCEP[0,1,:] # 0 = lowest member, 1 = mod precip, : = precip_thresholds
        axlocn = [0.415,0.51,0.25,0.16]          
        plot_at_these_gamma_threshes = [1,2,3,4] 
        index_gamma =[11,16,19,22]
        #plot_at_these_gamma_threshes = [1.0, 2.0, 4.0, 6.0]
        #index_gamma =[11, 16, 22, 26]
        xvalues = np.arange(0.0,10.0,0.02)
        xlimits = [0.,10.]
        ylimits = [0.,2.5]
        xticks = [0, 2, 4, 6, 8, 10]
    elif ipanel == 5:
        ctitle = '(f) NCEP [> 10 mm)'
        ytitle = ''
        xtitle = 'Precipitation amount (mm)'
        gamma_shape = gamma_shapes_NCEP[0,2,:] # 0 = lowest member, 2 = heavy precip, : = precip_thresholds
        gamma_scale = gamma_scales_NCEP[0,2,:] # 0 = lowest member, 2 = heavy precip, : = precip_thresholds
        axlocn = [0.72,0.51,0.25,0.16] 
        plot_at_these_gamma_threshes = [3,6,9,12]
        index_gamma = [19,25,31,34]
        #plot_at_these_gamma_threshes = [4.0, 8.0, 12.0, 16.0]
        #index_gamma = [22,30,35,39]
        xvalues = np.arange(0.0,30.0,0.1)
        xlimits = [0.,30.]  
        ylimits = [0.,0.5]
        xticks = [0,5,10,15,20,25,30]  
    elif ipanel == 6:
        ctitle = '(g) CMC [0.1 - 2 mm)'
        ytitle = 'Probability density'
        xtitle = 'Precipitation amount (mm)'
        gamma_shape = gamma_shapes_CMC[0,0,:] # 0 = lowest member, 0 = light precip, : = precip_thresholds
        gamma_scale = gamma_scales_CMC[0,0,:] # 0 = lowest member, 0 = light precip, : = precip_thresholds
        axlocn = [0.1,0.28,0.25,0.16]
        plot_at_these_gamma_threshes = [0.1, 0.2, 0.3, 0.4]
        index_gamma = [2,3,4,5]
        #plot_at_these_gamma_threshes = [0.2, 0.4, 0.6, 0.8]
        #index_gamma = [3,5,7,9]
        xvalues = np.arange(0.0,1.0,0.01)
        xlimits = [0.,1.]
        ylimits = [0.,25.]
        xticks = [0.0, 0.2, 0.4, 0.6, 0.8, 1.0]
    elif ipanel == 7:
        ctitle = '(h) CMC [2 - 10 mm)'
        ytitle = ''
        xtitle = 'Precipitation amount (mm)'
        gamma_shape = gamma_shapes_CMC[0,1,:] # 0 = lowest member, 1 = mod precip, : = precip_thresholds
        gamma_scale = gamma_scales_CMC[0,1,:] # 0 = lowest member, 1 = mod precip, : = precip_thresholds
        axlocn = [0.415,0.28,0.25,0.16]  
        plot_at_these_gamma_threshes = [1,2,3,4] 
        index_gamma =[11,16,19,22]        
        #plot_at_these_gamma_threshes = [1.0, 2.0, 4.0, 6.0]
        #index_gamma =[11, 16, 22, 26]
        xvalues = np.arange(0.0,10.0,0.02)
        xlimits = [0.,10.]
        ylimits = [0.,2.5]
        xticks = [0, 2, 4, 6, 8, 10]
    elif ipanel == 8:
        ctitle = '(i) CMC [> 10 mm)'
        ytitle = ''
        xtitle = 'Precipitation amount (mm)'
        gamma_shape = gamma_shapes_CMC[0,2,:] # 0 = lowest member, 2 = heavy precip, : = precip_thresholds
        gamma_scale = gamma_scales_CMC[0,2,:] # 0 = lowest member, 2 = heavy precip, : = precip_thresholds
        axlocn = [0.72,0.28,0.25,0.16] 
        plot_at_these_gamma_threshes = [3,6,9,12]
        index_gamma = [19,25,31,34] 
        #plot_at_these_gamma_threshes = [4.0, 8.0, 12.0, 16.0]
        #index_gamma = [22,30,35,39]
        xvalues = np.arange(0.0,30.0,0.1)
        xlimits = [0.,30.]  
        ylimits = [0.,0.5]
        xticks = [0,5,10,15,20,25,30]  
    elif ipanel == 9:
        ctitle = '(j) ECMWF [0.1 - 2 mm)'
        ytitle = 'Probability density'
        xtitle = 'Precipitation amount (mm)'
        gamma_shape = gamma_shapes_ECMWF[0,0,:] # 0 = lowest member, 0 = light precip, : = precip_thresholds
        gamma_scale = gamma_scales_ECMWF[0,0,:] # 0 = lowest member, 0 = light precip, : = precip_thresholds
        axlocn = [0.1,0.05,0.25,0.16]
        plot_at_these_gamma_threshes = [0.1, 0.2, 0.3, 0.4]
        index_gamma = [2,3,4,5]
        #plot_at_these_gamma_threshes = [0.2, 0.4, 0.6, 0.8]
        #index_gamma = [3,5,7,9]
        xvalues = np.arange(0.0,1.0,0.01)
        xlimits = [0.,1.]
        ylimits = [0.,25.]
        xticks = [0.0, 0.2, 0.4, 0.6, 0.8, 1.0]
    elif ipanel == 10:
        ctitle = '(k) ECMWF [2 - 10 mm)'
        ytitle = ''
        xtitle = 'Precipitation amount (mm)'
        gamma_shape = gamma_shapes_ECMWF[0,1,:] # 0 = lowest member, 1 = mod precip, : = precip_thresholds
        gamma_scale = gamma_scales_ECMWF[0,1,:] # 0 = lowest member, 1 = mod precip, : = precip_thresholds
        axlocn = [0.415,0.05,0.25,0.16]  
        plot_at_these_gamma_threshes = [1,2,3,4] 
        index_gamma =[11,16,19,22]        
        #plot_at_these_gamma_threshes = [1.0, 2.0, 4.0, 6.0]
        #index_gamma =[11, 16, 22, 26]
        xvalues = np.arange(0.0,10.0,0.02)
        xlimits = [0.,10.]
        ylimits = [0.,2.5]
        xticks = [0, 2, 4, 6, 8, 10]
    elif ipanel == 11:
        ctitle = '(l) ECMWF [> 10 mm)'
        ytitle = ''
        xtitle = 'Precipitation amount (mm)'
        gamma_shape = gamma_shapes_ECMWF[0,2,:] # 0 = lowest member, 2 = heavy precip, : = precip_thresholds
        gamma_scale = gamma_scales_ECMWF[0,2,:] # 0 = lowest member, 2 = heavy precip, : = precip_thresholds
        axlocn = [0.72,0.05,0.25,0.16]  
        plot_at_these_gamma_threshes = [3,6,9,12]
        index_gamma = [19,25,31,34]
        #plot_at_these_gamma_threshes = [4.0, 8.0, 12.0, 16.0]
        #index_gamma = [22,30,35,39]
        xvalues = np.arange(0.0,30.0,0.1)
        xlimits = [0.,30.] 
        ylimits = [0.,0.5]
        xticks = [0,5,10,15,20,25,30]  
        
    if ipanel <= 2:
            
        a1 = fig.add_axes(axlocn) 
        a1.set_title(ctitle,fontsize=11)
        a1.set_ylabel(ytitle,fontsize=10)
        #a1.set_xlabel(xtitle,fontsize=9)
        a1.set_xlim(xlimits)
        a1.set_ylim(ylimits)
        a1.plot(xplotdata, yplotdata_N,'-',color='Red',lw=2,label='NCEP')
        if ipanel >= 1:
            a1.plot(xplotdata, yplotdata_C,'-',color='LimeGreen',lw=2,label='CMC')
        else:
            a1.plot(xplotdata[0], yplotdata_C[0],'.',color='LimeGreen',lw=2,label='CMC',markersize=10)
        a1.plot(xplotdata, yplotdata_E,'-',color='RoyalBlue',lw=2,label='ECMWF')
        a1.grid(color='Gray',lw=0.2,linestyle='--')
        a1.set_xticks(xticks)
        a1.legend(loc=0)

    else:
        
        a1 = fig.add_axes(axlocn) 
        a1.set_title(ctitle,fontsize=11)
        a1.set_ylabel(ytitle,fontsize=10)
        if ipanel >= 9: a1.set_xlabel(xtitle,fontsize=9)
        a1.set_xlim(xlimits)
        a1.set_ylim(ylimits)
        for i in range(0,nplots):
            plotdata = gamma.pdf(xvalues, gamma_shape[index_gamma[i]], 0., gamma_scale[index_gamma[i]])
            ctext1 = "{:.2f}".format(precip_values[index_gamma[i]])
            a1.plot(xvalues, plotdata, '-', color=colors[i], lw=2, label='Forecast = '+ctext1)
        a1.grid(color='Gray',lw=0.2,linestyle='--')
        a1.legend(loc=0)   
        a1.set_xticks(xticks)
            
outfile = 'dressing_distributions_lowmbr_'+cyyyymmddhh+'_'+cleade+'h.pdf'        
plot_title = outfile
fig.savefig(plot_title)
print 'saving plot to file = ',plot_title


# ---- plot the fitted Gamma distributions for nonzero amounts 
#      for the highest ranked member



fig = plt.figure(figsize=(6.5,7.))
fig.suptitle('Dressing distributions, +'+cleadb+'-'+cleade+' h, highest-ranked member', fontsize=15)
        
nplots = 4
colors = ['Red','RoyalBlue','LimeGreen','Purple']     
for ipanel in range(9):
    if ipanel == 0:
        ctitle = '(a) NCEP [0.1 - 2 mm)'
        ytitle = 'Probability density'
        xtitle = 'Precipitation amount (mm)'
        gamma_shape = gamma_shapes_NCEP[2,0,:] # 2 = highest member, 0 = light precip, : = precip_thresholds
        gamma_scale = gamma_scales_NCEP[2,0,:] # 2 = highest member, 0 = light precip, : = precip_thresholds
        axlocn = [0.09,0.69,0.25,0.2]  
        plot_at_these_gamma_threshes = [1,2,3,4]
        index_gamma = [11,16,19,22]
        xvalues = np.arange(0.0,10.0,0.01)
        xlimits = [0.,10.]
        ylimits = [0,0.7]
        xticks = [0,2,4,6,8,10]
    elif ipanel == 1:
        ctitle = '(b) NCEP [2 - 10 mm)'
        ytitle = ''
        xtitle = 'Precipitation amount (mm)'
        gamma_shape = gamma_shapes_NCEP[2,1,:] # 2 = highest member, 1 = mod precip, : = precip_thresholds
        gamma_scale = gamma_scales_NCEP[2,1,:] # 2 = highest member, 1 = mod precip, : = precip_thresholds
        axlocn = [0.405,0.69,0.25,0.2]          
        plot_at_these_gamma_threshes = [5, 10, 15, 20]
        index_gamma =[24, 33, 38, 41]
        xvalues = np.arange(0.0,50.0,0.02)
        xlimits = [0.,50.]
        ylimits = [0., 0.4]
        xticks = [0,10,20,30,40,50]
    elif ipanel == 2:
        ctitle = '(c) NCEP [> 10 mm)'
        ytitle = ''
        xtitle = 'Precipitation amount (mm)'
        gamma_shape = gamma_shapes_NCEP[2,2,:] # 2 = highest member, 2 = heavy precip, : = precip_thresholds
        gamma_scale = gamma_scales_NCEP[2,2,:] # 2 = highest member, 2 = heavy precip, : = precip_thresholds
        axlocn = [0.73,0.69,0.25,0.2]  
        plot_at_these_gamma_threshes = [15, 20, 30, 40]
        index_gamma = [38,43,45,48]
        xvalues = np.arange(0.0,100.0,0.1)
        xlimits = [0.,100.]
        ylimits = [0.,0.15]
        xticks = [0,20,40,60,80,100]
    elif ipanel == 3:
        ctitle = '(d) CMC [0.1 - 2 mm)'
        ytitle = 'Probability density'
        xtitle = 'Precipitation amount (mm)'
        gamma_shape = gamma_shapes_CMC[2,0,:] # 2 = highest member, 0 = light precip, : = precip_thresholds
        gamma_scale = gamma_scales_CMC[2,0,:] # 2 = highest member, 0 = light precip, : = precip_thresholds
        axlocn = [0.09,0.38,0.25,0.2]  
        plot_at_these_gamma_threshes = [1,2,3,4]
        index_gamma = [11,16,19,22]
        xvalues = np.arange(0.0,10.0,0.01)
        xlimits = [0.,10.]
        ylimits = [0,0.7]
        xticks = [0,2,4,6,8,10]
    elif ipanel == 4:
        ctitle = '(e) CMC [2 - 10 mm)'
        ytitle = ''
        xtitle = 'Precipitation amount (mm)'
        gamma_shape = gamma_shapes_CMC[2,1,:] # 2 = highest member, 1 = mod precip, : = precip_thresholds
        gamma_scale = gamma_scales_CMC[2,1,:] # 2 = highest member, 1 = mod precip, : = precip_thresholds
        axlocn = [0.405,0.38,0.25,0.2]          
        plot_at_these_gamma_threshes = [5, 10, 15, 20]
        index_gamma =[24, 33, 38, 41]
        xvalues = np.arange(0.0,50.0,0.02)
        xlimits = [0.,50.]
        ylimits = [0., 0.4]
        xticks = [0,10,20,30,40,50]
    elif ipanel == 5:
        ctitle = '(f) CMC [> 10 mm)'
        ytitle = ''
        xtitle = 'Precipitation amount (mm)'
        gamma_shape = gamma_shapes_CMC[2,2,:] # 2 = highest member, 2 = heavy precip, : = precip_thresholds
        gamma_scale = gamma_scales_CMC[2,2,:] # 2 = highest member, 2 = heavy precip, : = precip_thresholds
        axlocn = [0.73,0.38,0.25,0.2]  
        plot_at_these_gamma_threshes = [15, 20, 30, 40]
        index_gamma = [38,43,45,48]
        xvalues = np.arange(0.0,100.0,0.1)
        xlimits = [0.,100.]
        ylimits = [0.,0.15]
        xticks = [0,20,40,60,80,100]                
    elif ipanel == 6:
        ctitle = '(g) ECMWF [0.1 - 2 mm)'
        ytitle = 'Probability density'
        xtitle = 'Precipitation amount (mm)'
        gamma_shape = gamma_shapes_ECMWF[2,0,:] # 2 = highest member, 0 = light precip, : = precip_thresholds
        gamma_scale = gamma_scales_ECMWF[2,0,:] # 2 = highest member, 0 = light precip, : = precip_thresholds
        axlocn = [0.09,0.07,0.25,0.2]  
        plot_at_these_gamma_threshes = [1,2,3,4]
        index_gamma = [11,16,19,22]
        xvalues = np.arange(0.0,10.0,0.01)
        xlimits = [0.,10.]
        ylimits = [0,0.7]
        xticks = [0,2,4,6,8,10]    
    elif ipanel == 7:
        ctitle = '(h) ECMWF [2 - 10 mm)'
        ytitle = ''
        xtitle = 'Precipitation amount (mm)'
        gamma_shape = gamma_shapes_ECMWF[2,1,:] # 2 = highest member, 1 = mod precip, : = precip_thresholds
        gamma_scale = gamma_scales_ECMWF[2,1,:] # 2 = highest member, 1 = mod precip, : = precip_thresholds
        axlocn = [0.405,0.07,0.25,0.2]          
        plot_at_these_gamma_threshes = [5, 10, 15, 20]
        index_gamma =[24, 33, 38, 41]
        xvalues = np.arange(0.0,50.0,0.02)
        xlimits = [0.,50.]
        ylimits = [0., 0.4]
        xticks = [0,10,20,30,40,50]
    elif ipanel == 8:
        ctitle = '(i) ECMWF [> 10 mm)'
        ytitle = ''
        xtitle = 'Precipitation amount (mm)'
        gamma_shape = gamma_shapes_ECMWF[2,2,:] # 2 = highest member, 2 = heavy precip, : = precip_thresholds
        gamma_scale = gamma_scales_ECMWF[2,2,:] # 2 = highest member, 2 = heavy precip, : = precip_thresholds
        axlocn = [0.73,0.07,0.25,0.2]  
        plot_at_these_gamma_threshes = [15, 20, 30, 40]
        index_gamma = [38,43,45,48]
        xvalues = np.arange(0.0,100.0,0.1)
        xlimits = [0.,100.]
        ylimits = [0.,0.15]
        xticks = [0,20,40,60,80,100]         

        
    a1 = fig.add_axes(axlocn) 
    a1.set_title(ctitle,fontsize=11)
    a1.set_ylabel(ytitle,fontsize=10)
    a1.set_xlabel(xtitle,fontsize=9)
    a1.set_xlim(xlimits)
    a1.set_ylim(ylimits)
    a1.set_xticks(xticks)
    for i in range(0,nplots):
        plotdata = gamma.pdf(xvalues, gamma_shape[index_gamma[i]], 0., gamma_scale[index_gamma[i]])
        ctext1 = "{:.2f}".format(precip_values[index_gamma[i]])
        a1.plot(xvalues, plotdata, '-', color=colors[i], lw=2, label='Forecast = '+ctext1)
    a1.grid(color='Gray',lw=0.2,linestyle='--')
    a1.legend(loc=0)
            
outfile = 'dressing_distributions_highmbr_'+cyyyymmddhh+'_'+cleade+'h.pdf'        
plot_title = outfile
fig.savefig(plot_title)
print 'saving plot to file = ',plot_title

# ---- make a plot of the fraction zero and gamma cdfs when ensemble mean is zero

rcParams['legend.fontsize']='medium'
fig = plt.figure(figsize=(8.0,6.0))
fig.suptitle('Dressing distributions, +'+cleadb+'-'+cleade+' h, '+r'$\bar{\tilde{x}}^{f} \leq 0.01\ mm$', fontsize=17)
        
xpamts = np.array([climo_pop_thresholds[0], climo_pop_thresholds[1], climo_pop_thresholds[2], \
    climo_pop_thresholds[3], climo_pop_thresholds[4], climo_pop_thresholds[5], \
    climo_pop_thresholds[6],2.*climo_pop_thresholds[6]-climo_pop_thresholds[5]])        
        
for ipanel in range(6):
    rcParams['xtick.labelsize']='medium'
    if ipanel == 0:
        axlocn = [0.1,0.55,0.25,0.32]
        ctitle = r'(a) Fraction zero, NCEP'
        ytitle = 'Fraction zero'
        xtitle = 'Climatological POP'
        xplotdata = xpamts 
        yplotdata = fraction_zeros_fclimpop_NCEP 
        ylimits = [0.94, 1.0]
        xlimits = [-0.01,0.26]
        xticks = [0.0, 0.05, 0.10, 0.15, 0.20, 0.25]
    elif ipanel == 1:
        axlocn = [0.415,0.55,0.25,0.32]
        ctitle = r'(b) Fraction zero, CMC'
        ytitle = 'Fraction zero'
        xtitle = 'Climatological POP'
        xplotdata = xpamts
        yplotdata = fraction_zeros_fclimpop_CMC 
        ylimits = [0.94, 1.0]
        xlimits = [-0.01,0.26]
        xticks = [0.0, 0.05, 0.10, 0.15, 0.20, 0.25]
    elif ipanel == 2:
        axlocn = [0.72,0.55,0.25,0.32]
        ctitle = r'(c) Fraction zero, ECMWF'
        ytitle = 'Fraction zero'
        xtitle = 'Climatological POP'
        xplotdata = xpamts 
        yplotdata = fraction_zeros_fclimpop_ECMWF
        ylimits = [0.94, 1.0]
        xlimits = [-0.01,0.26]
        xticks = [0.0, 0.05, 0.10, 0.15, 0.20, 0.25]
    elif ipanel == 3:
        axlocn = [0.1,0.08,0.25,0.32]
        ctitle = '(d) NCEP CDFs'
        ytitle = 'Non-exceedance probability'
        xtitle = 'Precipitation amount (mm)'
        gamma_shape = gamma_shape_fclimpop_NCEP 
        gamma_scale = gamma_scale_fclimpop_NCEP
        xvalues = np.arange(0.0,1.0,0.01)
        xlimits = [0.,3.]
        ylimits = [0.,1.]
        xticks = [0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0]
    elif ipanel == 4:
        axlocn = [0.415,0.08,0.25,0.32] 
        ctitle = '(e) CMC CDFs'
        ytitle = 'Non-exceedance probability'
        xtitle = 'Precipitation amount (mm)'
        gamma_shape = gamma_shape_fclimpop_CMC 
        gamma_scale = gamma_scale_fclimpop_CMC 
        xlimits = [0.,3.]
        ylimits = [0.,1.]
        xticks = [0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0]
    elif ipanel == 5:
        axlocn = [0.72,0.08,0.25,0.32] 
        ctitle = '(f) ECMWF CDFs'
        ytitle = 'Non-exceedance probability'
        xtitle = 'Precipitation amount (mm)'
        gamma_shape = gamma_shape_fclimpop_ECMWF
        gamma_scale = gamma_scale_fclimpop_ECMWF
        plot_at_these_gamma_threshes = [3,6,9,12]
        xlimits = [0.,3.]
        ylimits = [0.,1.]
        xticks = [0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0]  
    
    if ipanel <=2 :

        a1 = fig.add_axes(axlocn) 
        a1.set_title(ctitle,fontsize=12)
        if ipanel == 0: a1.set_ylabel(ytitle,fontsize=11)
        a1.set_xlabel(xtitle,fontsize=11)
        a1.set_xlim(xlimits)
        a1.set_ylim(ylimits)
        a1.set_xticks(xticks)
        print 'xplotdata = ', xplotdata
        print 'yplotdata = ', yplotdata
        a1.plot(xplotdata, yplotdata, 'r-', lw=2)
        a1.grid(color='Gray',lw=0.2,linestyle='--')
        
    else:
        
        a1 = fig.add_axes(axlocn) 
        a1.set_title(ctitle,fontsize=12)
        if ipanel == 3: a1.set_ylabel(ytitle,fontsize=11)
        a1.set_xlabel(xtitle,fontsize=11)
        a1.set_xlim(xlimits)
        a1.set_xticks(xticks)
        a1.grid(color='Gray',lw=0.2,linestyle='--')
        nt = len(gamma_shape_fclimpop_NCEP)
        colors = ['Red','RoyalBlue','LimeGreen','Purple','Black','Orange','Gray','DarkRed']
        for i in range(0,nt,2):
            xvals = np.arange(0.01, 3.01, 0.01)
            plotdata = gamma.cdf(xvals, gamma_shape[i], 0., gamma_scale[i])
            ctext1 = "{:.2f}".format(xpamts[i])
            a1.plot(xvals, plotdata, '-', color=colors[i], lw=2, label='POP = '+ctext1)
        a1.legend(loc=0)
            
outfile = 'fractionzero_dressdistns_'+cyyyymmddhh+'_'+cleade+'h.pdf'        
plot_title = outfile
fig.savefig(plot_title)
print 'saving plot to file = ',plot_title

