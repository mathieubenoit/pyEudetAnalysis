# Will process a tbtrack file
# Will produce an ntuple for further analysis with pyEudetAnalysisOnly.py
# For options, run:
# python pyEudetReconstructionOnly.py -h

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
                  help="alignment file", dest="ALIGNMENT", default="alignment.dat")

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
        print "Please provide a valid cluster position reconstruction method ( -m [method] QWeighted, maxTOT, DigitalControid, EtaCorrection)"
        parser.print_help()
        exit()
else:
    print "Please provide a valid cluster position reconstruction method ( -m [method] QWeighted, maxTOT, DigitalControid, EtaCorrection)"
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
    print "Please provide an output folder (-o [PathToOutput] , put no / at the end )"
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


import sys
sys.argv.append( '-b-' )
from ROOT import *
import ROOT
from ROOT import gStyle
from ROOT import TMath
from ToolBox import *
import pyximport; pyximport.install(pyimport=True)
from EudetData import *
from array import array

alignment_constants = ReadAlignment(AlignementPath)


gStyle.SetOptStat("nemruoi")
gStyle.SetOptFit(1111)


aDataSet = EudetData("%s/tbtrackrun%06i.root"%(input_folder,RunNumber),50000.0,edge_width,1,"tbtrack")


# Computing Chi2 cut and plotting Chi2 distribution
h_chi2,h_chi2ndof = aDataSet.GetChi2Cut()

can_chi2 = TCanvas()
h_chi2.Draw("")
can_chi2.SetLogx()
can_chi2.SetLogy()
can_chi2.SaveAs("%s/Run%i/Chi2_run%06i.root"%(PlotPath,RunNumber,RunNumber))

can_chi2ndof = TCanvas()
h_chi2ndof.Draw("")
can_chi2ndof.SetLogx()
can_chi2ndof.SetLogy()
can_chi2ndof.SaveAs("%s/Run%i/Chi2ndof_run%06i.root"%(PlotPath,RunNumber,RunNumber))

if(options.NEVENT):
    n_proc= int(options.NEVENT)

    if n_proc >= aDataSet.t_nEntries:
        n_proc = -1
    if n_proc == -1:
        n_proc = aDataSet.t_nEntries

else :
    n_proc= aDataSet.t_nEntries

print "Running on run %i, with Method %s, on %i Events"%(RunNumber,method_name,n_proc)


# Count pixel map repeats
histo_maprepeats = CountPixelMapRepeats(aDataSet,n_proc)
canreps = TCanvas()
histo_maprepeats.Draw()
histo_maprepeats.SaveAs("%s/Run%i/PixelMapRepeats_run%06i.root"%(PlotPath,RunNumber,RunNumber))


# Find Hot Pixels
histo_nhits,histo_hit,histo_hot,histo_freq = aDataSet.FindHotPixel(0.01,n_proc)

cannhits = TCanvas()
histo_nhits.Draw()
cannhits.SaveAs("%s/Run%i/NHit_pixels_run%06i.root"%(PlotPath,RunNumber,RunNumber))

canhit = TCanvas()
histo_hit.Draw("colz")
canhit.SaveAs("%s/Run%i/Hit_pixels_run%06i.root"%(PlotPath,RunNumber,RunNumber))

canhot = TCanvas()
histo_hot.Draw("colz")
canhot.SaveAs("%s/Run%i/Hot_pixels_run%06i.root"%(PlotPath,RunNumber,RunNumber))

canfreq = TCanvas()
canfreq.SetLogx()
canfreq.SetLogy()
histo_freq.Draw("")
canfreq.SaveAs("%s/Run%i/Firing_frequency_run%06i.root"%(PlotPath,RunNumber,RunNumber))


n_matched = 0
n_matched_edge = 0
last_time = time.time()
distances_histo = TH1F("distances_histo","",100,0.0,1.0)
prev_pixel_xhits = []

for i in range(0,n_proc) :
    aDataSet.getEvent(i)

    # is this a new pixel map?
    npixels_hit = len(aDataSet.p_col)
    pixel_x_hits = []
    for k in xrange(npixels_hit):
        pixel_x_hits.append(aDataSet.p_col[k])

    if (pixel_x_hits == prev_pixel_xhits):
        # same pixel map as before, will add clusters already computed
        aDataSet.AllClusters.append(clusters_tmp)
    else:
        # this is a new event, will cluster
        etacorr_sigma = 0.003
        aDataSet.ClusterEvent(i, method_name, etacorr_sigma)
        clusters_tmp = aDataSet.AllClusters[i]
    prev_pixel_xhits = pixel_x_hits

    aDataSet.GetTrack(i)

    for alignement in alignment_constants :
        ApplyAlignment_at_event(i,aDataSet,[alignement[3],alignement[4],0],[alignement[0],alignement[1],alignement[2]])

    aDataSet.FindMatchedCluster(i,0.1,6,distances_histo)
    m,me=aDataSet.ComputeResiduals(i)
    n_matched+=m
    n_matched_edge+=me
    if i%1000 ==0 :
        print "Event %d"%i
        print "Elapsed time/1000 Event : %f s"%(time.time()-last_time)
        last_time = time.time()


print "Found %i matched cluster"%(n_matched)

candist = TCanvas()
candist.SetLogy()
distances_histo.GetXaxis().SetTitle("Track-cluster distance (mm)")
distances_histo.Draw()
candist.SaveAs("%s/Run%i/%s/Cluster-track_dist.pdf"%(PlotPath,RunNumber,method_name))

root_file = "%s/Run%i/%s/pyEudetNtuple_run%i_%s.root"%(PlotPath,RunNumber,method_name,RunNumber,method_name)
os.system("rm %s"%root_file)

print "Writing reconstructed data to %s"%root_file
aDataSet.WriteReconstructedData(root_file,6)
