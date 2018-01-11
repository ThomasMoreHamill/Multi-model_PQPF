""" this routine will generate a text file that can be used to 
    run a python script dressing_statistics_to_netcdf.py 
    for multiple leads, dates, and models
"""

import numpy as np
from numpy import ma
import os, sys
from dateutils import daterange, dateshift

cmodels = ['NCEP','ECMWF','CMC'] # list of models
cleades = ['12','24','36','48','60','72','84',\
    '96','108','120','132','144','156','168'] # list of lead times

fname = 'precip_forecast_ccpa_2netcdf.deck'
outfile = open(fname,'w')
textstring = '#!/bin/tcsh \n'
outfile.write(textstring)
outfile.write(' \n')

for cmodel in cmodels:
    for cleade in cleades:
        textstring = 'python precip_forecast_ccpa_2netcdf.py '+' '+cmodel+' '+cleade+'\n'
        outfile.write(textstring)
            
outfile.close()