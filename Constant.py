###############################################################################################################################
#
#        Analysis constants: detector caracteristics and usefull constants
#
###############################################################################################################################

import future_builtins
SensorType=future_builtins.SensorType


if "Timepix" in SensorType :
# Timepix Specifications
    pitchX = 0.055
    pitchY = 0.055
    npix_X = 256
    npix_Y = 256
    
    print "Using Timepix detector"
    

elif SensorType=="CLICpix" : 
# CLICPix Specifications
    pitchX = 0.025
    pitchY = 0.025
    npix_X = 64
    npix_Y = 64
    
    print "Using CLICpix detector"


halfChip_X = npix_X*pitchX/2.
halfChip_Y = npix_Y*pitchY/2.

# sigma = 0.015
um = 1e-3
mm = 1
cm =10


someData_in_um = 1000*um

scaler =1



