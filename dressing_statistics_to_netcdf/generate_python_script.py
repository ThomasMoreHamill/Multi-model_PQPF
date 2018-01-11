""" this routine will generate a text file that can be used to 
    run a python script dressing_statistics_to_netcdf.py 
    for multiple leads, dates, and models
"""

import numpy as np
from numpy import ma
import os, sys
from dateutils import daterange, dateshift

cempirical = sys.argv[1] # 1 for empirical, 0 for fitted CDFs
#cempirical = '1'  # if 1, empirical is true
cmodels = ['NCEP','ECMWF','CMC']
cleades = ['24','48','72','96','120','144','168']
date_forecasts = daterange('2016040100', '2016070100', 24)

fname = 'dressing_statistics_to_netcdf.deck'
outfile = open(fname,'w')
textstring = '#!/bin/tcsh \n'
outfile.write(textstring)
outfile.write(' \n')

for cmodel in cmodels:
    for cleade in cleades:
        outfile.write(' \n')
        for date_forecast in date_forecasts:
            textstring = 'python dressing_statistics_to_netcdf.py '+\
                cmodel+' '+cleade+' '+date_forecast+' '+cempirical+'\n'
            outfile.write(textstring)
            
outfile.close()