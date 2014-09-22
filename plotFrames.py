# Will process a tbtrack file
# Will produce an alignment file
# For options, run:
# python ComputeAlignment.py -h

from math import fsum
import time,os
from optparse import OptionParser

parser = OptionParser()
parser.add_option("-r", "--run",
                  help="Run Number", dest="RUN", type="int")

parser.add_option("-n", "--nevent",
                  help="Number of events to process", dest="NEVENT")

parser.add_option("-d", "--data",
                  help="tbtrack Input Folder", dest="INPUT")

parser.add_option("-m", "--minpix",
                  help="minimum pixel in frame for display", dest="NPIX")
(options, args) = parser.parse_args()

if(options.RUN) :
    RunNumber = int(options.RUN)
else :
    print "Please provide a Run Number (-r [Run Number])"
    parser.print_help()
    exit()  
     
if(options.RUN) :
    npix = int(options.NPIX)
else :
    npix=1
       
if(options.INPUT):
    input_folder=options.INPUT
else :
    print "Please provide an input folder with tbtrack files (-d [PathToData] , put no / at the end )"
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


scaler = 1
method_name="DigitalCentroid"

if(options.NEVENT):
    n_proc= int(options.NEVENT)

    if n_proc >= aDataSet.t_nEntries:
        n_proc = -1
    if n_proc == -1:
        n_proc = aDataSet.t_nEntries

else :
    n_proc= aDataSet.t_nEntries

print "Running on run %i, with Method %s, on %i Events"%(RunNumber,method_name,n_proc)


plot=0

for i in range(n_proc) : 
    aDataSet.PlotFrame(i,plot,npix)



