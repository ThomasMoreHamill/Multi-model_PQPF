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
cmodel = sys.argv[2] # NCEP, ECMWF, or CMC
cleade = sys.argv[3] # ending of periods lead time in hours, e.g., '24', '108'
cleadb = str(int(cleade) - 12)
ileade = int(cleade)


# ======================================================================
# Part 1:  reading in dressing parameter data from netCDF files.
# ======================================================================

data_directory = '/Projects/Reforecast2/netcdf/NationalBlend/'   
infile = data_directory+cmodel+'/gamma_fraction_zero_dressing_'+\
    cmodel+'_date='+cyyyymmddhh+'_lead='+cleade+'.nc'
print infile
nc = Dataset(infile)    
precip_values = nc.variables['precip_values'][:]  
fraction_zeros = nc.variables['fraction_zeros'][:,:,:]
fraction_zeros_fclimpop = nc.variables['fraction_zeros_fclimpop'][:]
climo_pop_thresholds = nc.variables['climo_pop_thresholds'][:] 
gamma_shapes = nc.variables['gamma_shapes'][:,:,:]
gamma_scales = nc.variables['gamma_scales'][:,:,:]
gamma_shape_fclimpop = nc.variables['gamma_shape_fclimpop'][:] 
gamma_scale_fclimpop = nc.variables['gamma_scale_fclimpop'][:] 
precip_thresholds = nc.variables['precip_thresholds'][:] 
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

fig = plt.figure(figsize=(6.5,5.))
fig.suptitle(cmodel+' +'+cleadb+'-'+cleade+' h, lowest-ranked member', fontsize=17)
        
nplots = 4  
colors = ['Red','RoyalBlue','LimeGreen','Purple']     
for ipanel in range(6):
    if ipanel == 0:
        ctitle = r'$M(\bar{\tilde{x}}^{f}) = 1$'+'\n[0.1 - 2 mm)'
        ytitle = 'Fraction zero'
        xtitle = 'Precipitation amount (mm)'
        yplotdata = fraction_zeros[0,0,:] # 0 = lowest member, 0 = light precip, : = precip_thresholds
        xplotdata = precip_values 
        ylimits = [0.,1.]
        xlimits = [0.,1.]
        axlocn = [0.1,0.52,0.25,0.31]
        xticks = [0.0, 0.2, 0.4, 0.6, 0.8, 1.0]
        cletter = '(a)'
        cposition = [0.1, 0.9]
    elif ipanel == 1:
        ctitle = r'$M(\bar{\tilde{x}}^{f}) = 2$'+'\n[2 - 10 mm)'
        ytitle = ' '
        xtitle = 'Precipitation amount (mm)'
        yplotdata = fraction_zeros[0,1,:] # 0 = lowest member, 1 = mod precip, : = precip_thresholds
        xplotdata = precip_values 
        ylimits = [0.,1.]
        xlimits = [0.,10.]
        axlocn = [0.415,0.52,0.25,0.31]
        xticks = [0, 2, 4, 6, 8, 10]
        cletter = '(b)'
        cposition = [1.0, 0.9]
    elif ipanel == 2:
        ctitle = r'$M(\bar{\tilde{x}}^{f}) = 3$'+'\n[> 10 mm]'
        ytitle = ' '
        xtitle = 'Precipitation amount (mm)'
        yplotdata = fraction_zeros[0,2,:] # 0 = lowest member, 2 = heavy precip, : = precip_thresholds
        xplotdata = precip_values 
        ylimits = [0.,1.]
        xlimits = [0.,30.]
        axlocn = [0.73,0.52,0.25,0.31]
        xticks = [0,5,10,15,20,25,30]
        cletter = '(c)'
        cposition = [3.0, 0.9]
    elif ipanel == 3:
        ctitle = ' '
        ytitle = 'Probability density'
        xtitle = 'Precipitation amount (mm)'
        gamma_shape = gamma_shapes[0,0,:] # 0 = lowest member, 0 = light precip, : = precip_thresholds
        gamma_scale = gamma_scales[0,0,:] # 0 = lowest member, 0 = light precip, : = precip_thresholds
        axlocn = [0.1,0.09,0.25,0.31]
        plot_at_these_gamma_threshes = [0.2, 0.4, 0.6, 0.8]
        index_gamma = [3,5,7,9]
        xvalues = np.arange(0.0,1.0,0.01)
        xlimits = [0.,1.]
        xticks = [0.0, 0.2, 0.4, 0.6, 0.8, 1.0]
        cletter = '(d)'
        cposition = [0.1,10.0]
    elif ipanel == 4:
        ctitle = ' '
        ytitle = ''
        xtitle = 'Precipitation amount (mm)'
        gamma_shape = gamma_shapes[0,1,:] # 0 = lowest member, 1 = mod precip, : = precip_thresholds
        gamma_scale = gamma_scales[0,1,:] # 0 = lowest member, 1 = mod precip, : = precip_thresholds
        axlocn = [0.415,0.09,0.25,0.31]          
        plot_at_these_gamma_threshes = [1.0, 2.0, 4.0, 6.0]
        index_gamma =[11, 16, 22, 26]
        xvalues = np.arange(0.0,10.0,0.02)
        xlimits = [0.,10.]
        xticks = [0, 2, 4, 6, 8, 10]
        cletter = '(e)'
        cposition = [1, 1.3]
    elif ipanel == 5:
        ctitle = ' '
        ytitle = ''
        xtitle = 'Precipitation amount (mm)'
        gamma_shape = gamma_shapes[0,2,:] # 0 = lowest member, 2 = heavy precip, : = precip_thresholds
        gamma_scale = gamma_scales[0,2,:] # 0 = lowest member, 2 = heavy precip, : = precip_thresholds
        axlocn = [0.73,0.09,0.25,0.31]  
        plot_at_these_gamma_threshes = [4.0, 8.0, 12.0, 16.0]
        index_gamma = [22,30,35,39]
        xvalues = np.arange(0.0,30.0,0.1)
        xlimits = [0.,30.]  
        xticks = [0,5,10,15,20,25,30]  
        cletter = '(f)'
        cposition = [5, 0.27]  

    if ipanel <= 2:
            
        a1 = fig.add_axes(axlocn) 
        a1.set_title(ctitle,fontsize=12)
        a1.set_ylabel(ytitle,fontsize=11)
        a1.set_xlabel(xtitle,fontsize=9)
        a1.set_xlim(xlimits)
        a1.set_ylim(ylimits)
        a1.plot(xplotdata, yplotdata,'-',color='Red',lw=2)
        a1.grid(color='Gray',lw=0.2,linestyle='--')
        a1.set_xticks(xticks)
        a1.annotate(cletter, cposition, xycoords = 'data',textcoords = 'data', horizontalalignment = 'center')

    else:
        
        a1 = fig.add_axes(axlocn) 
        a1.set_title(ctitle,fontsize=12)
        a1.set_ylabel(ytitle,fontsize=11)
        a1.set_xlabel(xtitle,fontsize=9)
        a1.set_xlim(xlimits)
        for i in range(0,nplots):
            plotdata = gamma.pdf(xvalues, gamma_shape[index_gamma[i]], 0., gamma_scale[index_gamma[i]])
            ctext1 = "{:.2f}".format(precip_values[index_gamma[i]])
            a1.plot(xvalues, plotdata, '-', color=colors[i], lw=2, label='Forecast = '+ctext1)
        a1.grid(color='Gray',lw=0.2,linestyle='--')
        a1.legend(loc=0)   
        a1.set_xticks(xticks)
        a1.annotate(cletter, cposition, xycoords = 'data',textcoords = 'data', horizontalalignment = 'center')
            
outfile = 'dressing_distributions_lowmbr_'+cyyyymmddhh+'_'+cmodel+'_'+cleade+'h.pdf'        
plot_title = outfile
fig.savefig(plot_title)
print 'saving plot to file = ',plot_title



# ---- plot the fitted Gamma distributions for nonzero amounts 
#      for the highest ranked member



# precip_values = 0, 0.04980469, 0.09960938, 0.2001953, 0.2998047, \
#    0.4003906, 0.5, 0.5996094, 0.7001953, 0.7998047, \
#    0.9003906, 1, 1.200195, 1.400391,  1.599609, 
#    1.799805, 2, 2.299805, 2.599609, 3, 
#    3.299805, 3.599609, 4, 4.5, 5, 5.5, 
#    6, 6.5, 7, 7.5, 8, 
#    8.5, 9, 10, 11, 12, 
#    13, 14, 15, 16, 18, 
#    20, 22.5, 25, 27.5, 30, 
#    33, 36, 40, 45, 50, 
#    55, 60, 65, 70, 75, 
#    80, 85, 90, 95, 100, 
#    120, 140, 160, 180, 200, 
#    250, 300 ;


fig = plt.figure(figsize=(6.5,3.))
fig.suptitle(cmodel+' +'+cleadb+'-'+cleade+' h, highest-ranked member', fontsize=17)
        
nplots = 4  
colors = ['Red','RoyalBlue','LimeGreen','Purple']     
for ipanel in range(3):
    if ipanel == 0:
        ctitle = r'$M(\bar{\tilde{x}}^{f}) = 1$'+'\n[0.1 - 2 mm)'
        ytitle = 'Probability density'
        xtitle = 'Precipitation amount (mm)'
        gamma_shape = gamma_shapes[2,0,:] # 2 = highest member, 0 = light precip, : = precip_thresholds
        gamma_scale = gamma_scales[2,0,:] # 2 = highest member, 0 = light precip, : = precip_thresholds
        axlocn = [0.08,0.15,0.25,0.58]
        plot_at_these_gamma_threshes = [1,2,3,4]
        index_gamma = [11,16,19,22]
        xvalues = np.arange(0.0,10.0,0.01)
        xlimits = [0.,10.]
        xticks = [0,2,4,6,8,10]
        cletter = '(a)'
        cposition = [0.7, 0.46]  
        
    elif ipanel == 1:
        ctitle = r'$M(\bar{\tilde{x}}^{f}) = 2$'+'\n[2 - 10 mm)'
        ytitle = ''
        xtitle = 'Precipitation amount (mm)'
        gamma_shape = gamma_shapes[2,1,:] # 2 = highest member, 1 = mod precip, : = precip_thresholds
        gamma_scale = gamma_scales[2,1,:] # 2 = highest member, 1 = mod precip, : = precip_thresholds
        axlocn = [0.405,0.15,0.25,0.58]          
        plot_at_these_gamma_threshes = [5, 10, 15, 20]
        index_gamma =[24, 33, 38, 41]
        xvalues = np.arange(0.0,50.0,0.02)
        xlimits = [0.,50.]
        xticks = [0,10,20,30,40,50]
        cletter = '(b)'
        cposition = [3, 0.24]
    elif ipanel == 2:
        ctitle = r'$M(\bar{\tilde{x}}^{f}) = 3$'+'\n[> 10 mm]'
        ytitle = ''
        xtitle = 'Precipitation amount (mm)'
        gamma_shape = gamma_shapes[2,2,:] # 2 = highest member, 2 = heavy precip, : = precip_thresholds
        gamma_scale = gamma_scales[2,2,:] # 2 = highest member, 2 = heavy precip, : = precip_thresholds
        axlocn = [0.73,0.15,0.25,0.58]  
        plot_at_these_gamma_threshes = [15, 20, 30, 40]
        index_gamma = [38,43,45,48]
        xvalues = np.arange(0.0,100.0,0.1)
        xlimits = [0.,100.]
        xticks = [0,20,40,60,80,100]
        cletter = '(c)'
        cposition = [7, 0.129]                 

        
    a1 = fig.add_axes(axlocn) 
    a1.set_title(ctitle,fontsize=12)
    a1.set_ylabel(ytitle,fontsize=10)
    a1.set_xlabel(xtitle,fontsize=9)
    a1.set_xlim(xlimits)
    a1.set_xticks(xticks)
    for i in range(0,nplots):
        plotdata = gamma.pdf(xvalues, gamma_shape[index_gamma[i]], 0., gamma_scale[index_gamma[i]])
        ctext1 = "{:.2f}".format(precip_values[index_gamma[i]])
        a1.plot(xvalues, plotdata, '-', color=colors[i], lw=2, label='Forecast = '+ctext1)
    a1.grid(color='Gray',lw=0.2,linestyle='--')
    a1.legend(loc=0)
    a1.annotate(cletter, cposition, xycoords = 'data',textcoords = 'data', horizontalalignment = 'center')
            
outfile = 'dressing_distributions_highmbr_'+cyyyymmddhh+'_'+cmodel+'_'+cleade+'h.pdf'        
plot_title = outfile
fig.savefig(plot_title)
print 'saving plot to file = ',plot_title



# ---- make a plot of the fraction zero and gamma pdfs when ensemble mean is zero

rcParams['legend.fontsize']='medium'
fig = plt.figure(figsize=(6.5,6.))
fig.suptitle(cmodel+' +'+cleadb+'-'+cleade+' h, '+r'$\bar{\tilde{x}}^{f} \leq 0.01\ mm$', fontsize=17)
        
nplots = 4  
colors = ['Red','RoyalBlue','LimeGreen','Purple']     

ctitle = '(a) Fraction zero'
ytitle = 'Fraction Zero'
xtitle = 'Climatological POP (> 0.254 mm / 12 h)'
axlocn = [0.12,0.55,0.82,0.33]
print climo_pop_thresholds
xvalues = np.array([0,climo_pop_thresholds[0], climo_pop_thresholds[1], climo_pop_thresholds[2], \
    climo_pop_thresholds[3], climo_pop_thresholds[4], climo_pop_thresholds[5], \
    climo_pop_thresholds[6],2.*climo_pop_thresholds[6]-climo_pop_thresholds[5]])

print 'xvalues = ', xvalues 
xvalues_avg = (xvalues[0:-1] + xvalues[1:])/2.
print 'xvalues_avg = ', xvalues_avg
yvalues = fraction_zeros_fclimpop 
a1 = fig.add_axes(axlocn) 
a1.set_title(ctitle,fontsize=14)
a1.set_ylabel(ytitle,fontsize=11)
a1.set_xlabel(xtitle,fontsize=11)
a1.set_xlim([0,0.25])
a1.set_xticks([0.0,0.05,0.10,0.15,0.20,0.25])
a1.plot(xvalues_avg, yvalues, 'r-', lw=2)
a1.grid(color='Gray',lw=0.2,linestyle='--')
a1.legend(loc=0)

ctitle = '(b) Sample Gamma dressing distributions'
ytitle = 'Probability Density'
xtitle = 'Precipitation amount (mm)'
axlocn = [0.12,0.08,0.82,0.33]
xvalues = np.arange(0.0, 0.2, 0.001)
a1 = fig.add_axes(axlocn) 
a1.set_title(ctitle,fontsize=14)
a1.set_ylabel(ytitle,fontsize=11)
a1.set_xlabel(xtitle,fontsize=11)
a1.set_xlim([0,0.2])
a1.set_xticks([0,0.05,0.1,0.15,0.2])
a1.grid(color='Gray',lw=0.2,linestyle='--')

plotdata1 = gamma.pdf(xvalues, gamma_shape_fclimpop[1], 0., gamma_scale_fclimpop[1])
ctext1 = "{:.2f}".format(climo_pop_thresholds[1])
a1.plot(xvalues, plotdata1, '-', color=colors[0], lw=2, label='Climo POP = '+ctext1)

plotdata2 = gamma.pdf(xvalues, gamma_shape_fclimpop[2], 0., gamma_scale_fclimpop[2])
ctext2 = "{:.2f}".format(climo_pop_thresholds[2])
a1.plot(xvalues, plotdata2, '-', color=colors[1], lw=2, label='Climo POP = '+ctext2)
 
plotdata3 = gamma.pdf(xvalues, gamma_shape_fclimpop[4], 0., gamma_scale_fclimpop[4])
ctext3 = "{:.2f}".format(climo_pop_thresholds[4])
a1.plot(xvalues, plotdata3, '-', color=colors[2], lw=2, label='Climo POP = '+ctext3)

plotdata4 = gamma.pdf(xvalues, gamma_shape_fclimpop[6], 0., gamma_scale_fclimpop[6])
ctext4 = "{:.2f}".format(climo_pop_thresholds[6])
a1.plot(xvalues, plotdata4, '-', color=colors[3], lw=2, label='Climo POP = '+ctext4)

a1.legend(loc=0)
            
outfile = 'fractionzero_dressdistns_'+cyyyymmddhh+'_'+cmodel+'_'+cleade+'h.pdf'        
plot_title = outfile
fig.savefig(plot_title)
print 'saving plot to file = ',plot_title

