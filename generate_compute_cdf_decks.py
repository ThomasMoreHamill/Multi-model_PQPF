import numpy as np
from dateutils import daterange, dateshift
import sys, os

cmodel = sys.argv[1]  # ECMWF, CMC, or NCEP

outfilename = 'control_cdf_generation_'+cmodel+'.deck'
cleads = ['24','48','72','96','120','144','168']
cdates = ['2016013000','0216021000','2016022000','2016022800','2016031000',\
    '2016032000','2016033000','2016041000','2016042000','2016043000',\
    '2016051000','2016052000','2016053000','2016061000','2016062000',\
    '2016063000'] # the list of dates for when we want to generate CDFs

outfile = open(outfilename,'w')  
print 'writing to ',outfilename
textstring ='#!/bin/tcsh'      
outfile.write(textstring)
outfile.write(' \n')    
for clead in cleads:
    outfile.write(' \n')
    for cdate in cdates:
        textstring = 'python compute_precip_CCPAgrid_cdfs.py '+cmodel+' '+clead+' '+cdate+' \n'
        outfile.write(textstring)       
outfile.close()

