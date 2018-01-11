""" this routine will generate a script file that can be used to 
    run python routine dressing_statistics_to_netcdf.py 
    for multiple leads, dates, and models.  The python routine
    generates the information that is later used to generate 
    fitted Gamma CDFs used for quantile mapping.   This is an
    alternative to the use of empirical CDFs in the 
    precip_forecast_CCPA_2netcdf folder.
"""

import numpy as np
from numpy import ma
import os, sys
from dateutils import daterange, dateshift

cmodels = ['NCEP','ECMWF','CMC'] # the models we seek to generate CDFs for
cleades = ['24','48','72','96','120','144','168'] # ending forecast lead times
date_forecasts = daterange('2015120100', '2016070100', 24) # list of dates

fname = 'compute_singleday_gamma_stats.deck'
outfile = open(fname,'w')
textstring = '#!/bin/tcsh \n'
outfile.write(textstring)
outfile.write(' \n')

for cmodel in cmodels:
    for cleade in cleades:
        outfile.write(' \n')
        for date_forecast in date_forecasts:
            textstring = 'python compute_singleday_gamma_stats.py '+date_forecast+' '+cmodel+' '+cleade+'\n'
            outfile.write(textstring)
            
outfile.close()