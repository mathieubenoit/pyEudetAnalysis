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

parser.add_option("-m", "--method",
                  help="Position Reconstruction Method, QWeighted,  DigitalCentroid, maxTOT, EtaCorrection", dest="METHOD", default="QWeighted")

parser.add_option("-d", "--data",
                  help="tbtrack Input Folder", dest="INPUT")

parser.add_option("-o", "--output",
                  help="Histograms and results output folder", dest="OUTPUT", default=".")

parser.add_option("-a", "--alignment",
                  help="alignement file", dest="ALIGNMENT", default="alignement.dat")

parser.add_option("-e", "--edge",
                  help="edge width", dest="EDGE", default=0.0, type="float")

(options, args) = parser.parse_args()

if(options.RUN) :
    RunNumber = int(options.RUN)
else :
    print "Please provide a Run Number (-r [Run Number])"
    parser.print_help()
    exit()  
     
if(options.EDGE) :
    edge_width = float(options.EDGE)
else : 
    edge_width = 0.0
    
if(options.METHOD) :
    if(options.METHOD=="QWeighted"):
        method_name=options.METHOD
    elif(options.METHOD=="maxTOT"):
        method_name=options.METHOD
    elif(options.METHOD=="DigitalCentroid"):
        method_name=options.METHOD
    elif(options.METHOD=="EtaCorrection"):
        method_name=options.METHOD
    else :
        print "Please provide a valid cluster position reconstruction method ( -m [method]  QWeighted,maxTOT,DigitalControid,EtaCorrection)"
        parser.print_help()
        exit()

else:
    print "Please provide a valid cluster position reconstruction method ( -m [method]  QWeighted,maxTOT,DigitalControid,EtaCorrection)"
    parser.print_help()
    exit()

if(options.INPUT):
    input_folder=options.INPUT
else :
    print "Please provide an input folder with tbtrack files (-d [PathToData] , put no / at the end )"
    parser.print_help()
    exit()

if(options.OUTPUT):
    PlotPath=options.OUTPUT
else :
    print "Please provide an output folder with tbtrack files (-o [PathToOutput] , put no / at the end )"
    parser.print_help()
    exit()

if(options.ALIGNMENT):
    AlignementPath = "%s"%(options.ALIGNMENT)
else :
    print "Please provide an Alignment File (-a [PathToFile]  0 0 0 0 0 if no alignement needed )"
    parser.print_help()
    exit()

os.system("mkdir %s/Run%i"%(PlotPath,RunNumber))
os.system("mkdir %s/Run%i/%s"%(PlotPath,RunNumber,method_name))

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


aDataSet = EudetData("%s/tbtrackrun%06i.root"%(input_folder,RunNumber),50000.0,edge_width,1,RunNumber,"tbtrack")


scaler = 1

if(options.NEVENT):
    n_proc= int(options.NEVENT)

    if n_proc >= aDataSet.t_nEntries:
        n_proc = -1
    if n_proc == -1:
        n_proc = aDataSet.t_nEntries

else :
    n_proc= aDataSet.t_nEntries

print "Running on run %i, with Method %s, on %i Events"%(RunNumber,method_name,n_proc)

histo_nhits,histo_hit,histo_hot,histo_freq = aDataSet.FindHotPixel(0.01,n_proc)

cannhits = TCanvas()
histo_nhits.Draw()

canhit = TCanvas()
histo_hit.Draw("colz")

canhot = TCanvas()
histo_hot.Draw("colz")

canfreq = TCanvas()
canfreq.SetLogx()
canfreq.SetLogy()
histo_freq.Draw("")


prev_pixel_xhits = []
last_time=time.time()

for i in range(0,n_proc) :
    aDataSet.getEvent(i)

    # is this a new pixel map?
    npixels_hit = len(aDataSet.p_col)
    pixel_x_hits = []
    for k in xrange(npixels_hit):
        pixel_x_hits.append(aDataSet.p_col[k])

    if (pixel_x_hits == prev_pixel_xhits):
        # print "same pixel map as before, will add clusters already computed"
        aDataSet.AllClusters.append(clusters_tmp)
    else:
        # print "this is a new event, will cluster"
        etacorr_sigma = 0.003
        aDataSet.ClusterEvent(i, method_name, etacorr_sigma, scaler)
        clusters_tmp = aDataSet.AllClusters[i]
    prev_pixel_xhits = pixel_x_hits

    aDataSet.GetTrack(i)
    if i%1000 ==0 :
        print "Event %d"%i
        print "Elapsed time/1000 Event. Clustering : %f s"%(time.time()-last_time)
        last_time = time.time()

last_time=time.time()

tccorx1,tccory1 = TrackClusterCorrelation(aDataSet,6,n_proc)
tccorx1.SetName("tccor1")
tccory1.SetName("tccory1")
cantccorx1 = TCanvas()
tccorx1.Draw("colz")
cantccory1 = TCanvas()
tccory1.Draw("colz")

print "Press any key to continue, ctrl-D to stop"
bla = raw_input()


print "Performing prealignment"
alignment_constants, prealix, prealiy = PerformPreAlignement(aDataSet,n_proc,1,AlignementPath,6,[0,0,0])
canprealix = TCanvas()
prealix.Draw()
canprealiy = TCanvas()
prealiy.Draw()

print "Press any key to continue, ctrl-D to stop"
b=raw_input() 

last_time = time.time()

for i in range(0,n_proc) :

    for alignement in alignment_constants :
        ApplyAlignment_at_event(i,aDataSet,[alignement[3],alignement[4],0],[alignement[0],alignement[1],alignement[2]])

    aDataSet.FindMatchedCluster(i,0.3,6)
    a,b=aDataSet.ComputeResiduals(i)
    if i%1000 ==0 :
        print "Event %d"%i
        print "Elapsed time/1000 Event. Apply Alignement and TrackMatching : %f s"%(time.time()-last_time)
        last_time = time.time()


tccorx2,tccory2 = TrackClusterCorrelation(aDataSet,6,n_proc)
tccorx2.SetName("tccorx2")
tccory2.SetName("tccory2")
cancorx2 = TCanvas()
tccorx2.Draw("colz")
cancory2 = TCanvas()
tccory2.Draw("colz")


niter = 2
for i in range(niter) :
    resr,rest = Perform3StepAlignment(aDataSet,[[0,360],[0,360],[0,360],[-0.5,0.5],[-0.5,0.5]],n_proc,1,0.05,AlignementPath,1e-3,[0,0,0])
    ApplyAlignment(aDataSet,rest,resr)


tccorx3,tccory3 = TrackClusterCorrelation(aDataSet,6,n_proc)
tccorx3.SetName("tccorx3")
tccory3.SetName("tccory3")
cancorx3 = TCanvas()
tccorx3.Draw("colz")
cancory3 = TCanvas()
tccory3.Draw("colz")


n_matched = 0

for i in range(0,n_proc) :

    aDataSet.FindMatchedCluster(i,0.3,6)
    a,b=aDataSet.ComputeResiduals(i)
    n_matched+=a
    if i%1000 ==0 :
        print "Event %d"%i
        print "Elapsed time/1000 Event. Apply Alignement and TrackMatching : %f s"%(time.time()-last_time)
        last_time = time.time()

      
# Write all histograms to output root file
out = TFile("%s/Run%i/%s/alignment_rootfile.root"%(PlotPath,RunNumber,method_name), "recreate")
out.cd()
histo_nhits.Write()
histo_hits.Write()
histo_hot.Write()
histo_freq.Write()
tccorx1.Write()
tccory1.Write()
tccorx2.Write()
tccory2.Write()
tccorx3.Write()
tccory3.Write()
prealix.Write()
prealiy.Write()

print "Found %i matched track-cluster binome"%n_matched
print "Produced Alignment file %s"%AlignementPath
