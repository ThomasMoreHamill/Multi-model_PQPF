""" this routine will generate a text file that can be used to 
    run a python script dressing_statistics_to_netcdf.py 
    for multiple leads, dates, and models
"""

import numpy as np
from numpy import ma
import os, sys
from dateutils import daterange, dateshift

cmodels = ['NCEP','ECMWF','CMC']
cleades = ['24','48','72','96','120','144','168']
date_forecasts = daterange('2016040100', '2016070100', 24)

for cmodel in cmodels:

    fname = 'blend_execute_'+cmodel+'.deck'
    outfile = open(fname,'w')
    textstring = '#!/bin/tcsh \n'
    outfile.write(textstring)
    outfile.write(' \n')
    
    for cleade in cleades:
        outfile.write(' \n')
        for date_forecast in date_forecasts:
            textstring = 'blend_precip_singlemodel_dressed_gammacdf.x '+date_forecast+' '+cleade+' '+cmodel+'\n'
            outfile.write(textstring)
            
    outfile.close()