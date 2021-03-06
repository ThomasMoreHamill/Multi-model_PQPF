import pygrib
from dateutils import daterange, dateshift, dayofyear, splitdate
import os, sys

# ---- I note (sadly) that for this project, to save space the grib inputs and outputs
#      were deleted.

# --- set up dictionary associating center_wmo_code with the more common name of 
#     the prediction center

dict_center = {"babj":"CMA", "cwao":"CMC", "ecmf":"ECMWF", "rjtd":"JMA", 
               "egrr":"UKMO", "kwbc":"NCEP" , "rksl":"KMA", "sbsj":"CPTEC"}


infile = sys.argv[1] # the name of the grib file with lots of records in it
fverif = pygrib.open(infile)
fverif.seek(0)
for grb in fverif:
   print grb.centre, grb.forecastTime, grb.dataDate, grb.endStep, grb.perturbationNumber
   cdataDate = str(grb.dataDate)+'00'
   cendStep = str(grb.endStep)
   cpertno = str(grb.perturbationNumber)
   if grb.shortName != 'lsm' and grb.shortName != 'orog':
      #if grb.forecastTime == 0:
      msg = grb.tostring()
      center = dict_center[grb.centre]
      outfilename = '/Users/thamill/precip/ecmwf_data/'+\
          center+'/'+center+'_'+cdataDate+'_fhour'+cendStep+'_pertno'+cpertno+'.grb'
      print outfilename
      grbout = open(outfilename,'ab')  # append, binary
      grbout.write(msg)
      grbout.close
         
