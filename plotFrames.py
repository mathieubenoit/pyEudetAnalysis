# Will process a tbtrack file
# Will pop out a canvas showing the DUT pixel map
# For options, run:
# python plotFrames.py -h

from math import fsum
import time,os
from optparse import OptionParser
import future_builtins

parser = OptionParser()
parser.add_option("-r", "--run",
                  help="Run Number", dest="RUN", type="int")

parser.add_option("-n", "--nevent",
                  help="Number of events to process", dest="NEVENT")

parser.add_option("-f", "--fevent",
                  help="First event to process", dest="FEVENT")

parser.add_option("-d", "--data",
                  help="tbtrack Input Folder", dest="INPUT")

parser.add_option("-m", "--minpix",
                  help="minimum pixel in frame for display", dest="NPIX")

parser.add_option("-s", "--sensor",
                  help="Sensor type", dest="SENSOR", default="Timepix")

(options, args) = parser.parse_args()

if(options.RUN) :
    RunNumber = int(options.RUN)
else :
    print "Please provide a Run Number (-r [Run Number])"
    parser.print_help()
    exit()  
     
if(options.NPIX) :
    n_pix_min = int(options.NPIX)
else :
    n_pix_min = 0
       
if(options.INPUT):
    input_folder=options.INPUT
else :
    print "Please provide an input folder with tbtrack files (-d [PathToData] , put no / at the end )"
    parser.print_help()
    exit()

if(("Timepix" in options.SENSOR) or options.SENSOR=="CLICpix"):
    future_builtins.SensorType=options.SENSOR
else :
    print "Please provide known sensor name. Timepix/Timepix3 (default) or CLICpix"
    parser.print_help()
    exit()


from ROOT import *
import ROOT
from ROOT import gStyle
from ROOT import TMath
from ToolBox import *
import pyximport; pyximport.install(pyimport=True)
from EudetData import *
from array import array


gStyle.SetOptStat("nemruoi")
gStyle.SetOptFit(1111)


aDataSet = EudetData("%s/tbtrackrun%06i.root"%(input_folder,RunNumber),50000.0,0,1,RunNumber,"tbtrack")

if(options.FEVENT):
    first_event = int(options.FEVENT)
else:
    first_event = 0

if(options.NEVENT):
    n_proc= int(options.NEVENT)

    if n_proc >= aDataSet.t_nEntries - first_event:
        n_proc = -1
    if n_proc == -1:
        n_proc = aDataSet.t_nEntries - first_event

else :
    n_proc= aDataSet.t_nEntries - first_event

print "Running on run %i, will show maximum %i frames, starting at event %i" %(RunNumber,n_proc,first_event)

canvas = TCanvas("canvas","",0,0,800,800)

for i in range(first_event,first_event+n_proc) : 
    aDataSet.PlotFrame(i,canvas,n_pix_min)



