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

# method_name = "QWeighted"
# method_name = "DigitalCentroid"
# method_name = "maxTOT"
#method_name = "EtaCorrection"


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

#

os.system("mkdir %s/Run%i"%(PlotPath,RunNumber))
os.system("mkdir %s/Run%i/QWeighted"%(PlotPath,RunNumber))
os.system("mkdir %s/Run%i/maxTOT"%(PlotPath,RunNumber))
os.system("mkdir %s/Run%i/DigitalCentroid"%(PlotPath,RunNumber))
os.system("mkdir %s/Run%i/EtaCorrection"%(PlotPath,RunNumber))

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


###############################################################################################################################
# first compile the code in ROOT using ACLIC: (testLangauFit.C)
# .L File.C+

#from guppy import hpy
#h = hpy()
#gSystem.Load('testLangauFit_C.so')
# gSystem.Load('/afs/cern.ch/work/a/apequegn/public/DESY_TB_DATA_02_07-06-2013_results/pyEudetAnalysisPlots/testLangauFit_C.so')

# Now the function should be available in the ROOT namespace
#from ROOT import langaufun, langaufit, langaupro

###############################################################################################################################

#PlotPath = "/VertexScratch/TB_Data/DESY_TB_DATA_August2013_results/pyEudetAnalysis_plots"

# global n_sizeX2sizeY2


def h1_style(h, optstat=0) :
    h.SetStats(optstat)
    h.SetLabelFont(42,"X")
    h.SetLabelFont(42,"Y")
    h.SetLabelOffset(0.005,"X")
    h.SetLabelOffset(0.005,"Y")
    h.SetLabelSize(0.045,"X")
    h.SetLabelSize(0.045,"Y")
    h.SetTitleOffset(1.15,"X")
    h.SetTitleOffset(1.15,"Y")
    h.SetTitleSize(0.04,"X")
    h.SetTitleSize(0.04,"Y")
    #h.SetTitle(0)
    h.SetTitleFont(42, "XYZ")

def TGraph_style (h) :
    h.GetXaxis().SetLabelOffset(0.005)
    h.GetXaxis().SetLabelFont(42)
    h.GetXaxis().SetLabelSize(0.055)
    h.GetXaxis().SetTitleOffset(1.15)
    h.GetXaxis().SetTitleSize(0.04)
    h.GetXaxis().SetTitleFont(42)
#     h.GetYaxis().SetRangeUser(0.,1.)
    h.GetYaxis().SetLabelOffset(0.005)
    h.GetYaxis().SetLabelFont(42)
    h.GetYaxis().SetLabelSize(0.045)
    h.GetYaxis().SetTitleOffset(1.2)
    h.GetYaxis().SetTitleFont(42)
    h.GetYaxis().SetTitleSize(0.04)



gStyle.SetOptStat("nemruoi")
gStyle.SetOptFit(1111)



#aDataSet = EudetData("/VertexScratch/TB_Data/DESY_TB_DATA_02_07-06-2013_results/histo/tbtrackrun000062.root",500.0)
aDataSet = EudetData("%s/tbtrackrun%06i.root"%(input_folder,RunNumber),50000.0,edge_width,1,"tbtrack")




# Computing Chi2 cut and plotting Chi2 distribution
h_chi2,h_chi2ndof = aDataSet.GetChi2Cut()

can_chi2 = TCanvas()
h_chi2.Draw("")
can_chi2.SetLogx()
can_chi2.SetLogy()

can_chi2ndof = TCanvas()
h_chi2ndof.Draw("")
can_chi2ndof.SetLogx()
can_chi2ndof.SetLogy()


scaler = 1
#n_proc= 25000

if(options.NEVENT):
    n_proc= int(options.NEVENT)
else :
    n_proc= aDataSet.t_nEntries

print "Running on run %i, with Method %s, on %i Events"%(RunNumber,method_name,n_proc)
# aDataSet.PrintTBranchElement()

trackX_vs_trackY_plan3 = TH2D("trackX_vs_trackY_plan3","track_posX[3] wrt track_posY[3]",300,-20.,20.,300,-20.,20.)
# trackX_vs_trackY_plan3.GetXaxis().SetRangeUser(-0.,14.08)
trackX_vs_trackY_plan3.GetXaxis().SetTitle("Track X position within pixel [mm]")
# trackX_vs_trackY_plan3.GetYaxis().SetRangeUser(-0.,14.08)
trackX_vs_trackY_plan3.GetYaxis().SetTitle("Track Y position within pixel [mm]")

trackX_vs_trackY_plan0 = TH2D("trackX_vs_trackY_plan0","track_posX[0] wrt track_posY[0]",300,-20.,20.,300,-20.,20.)
# trackX_vs_trackY_plan0.GetXaxis().SetRangeUser(-0.,14.08)
trackX_vs_trackY_plan0.GetXaxis().SetTitle("Track X position within pixel [mm]")
# trackX_vs_trackY_plan0.GetYaxis().SetRangeUser(-0.,14.08)
trackX_vs_trackY_plan0.GetYaxis().SetTitle("Track Y position within pixel [mm]")

h1_style(trackX_vs_trackY_plan3)
h1_style(trackX_vs_trackY_plan0)

# Filter Hot Pixels
# histo_hot,histo_freq = aDataSet.FilterHotPixel(0.005,200)
histo_hot,histo_freq = aDataSet.FilterHotPixel(0.005,n_proc,15)

canhot = TCanvas()
histo_hot.Draw("colz")

canfreq = TCanvas()
canfreq.SetLogx()
canfreq.SetLogy()
histo_freq.Draw("")

canfreq.SaveAs("%s/Run%i/Firing_frequency_run%06i.root"%(PlotPath,RunNumber,RunNumber))

n_matched = 0
n_matched_edge = 0
#for i in range(aDataSet.p_nEntries) :

last_time=time.time()

print alignment_constants

for i in range(0,n_proc,scaler) :
    aDataSet.getEvent(i)
    aDataSet.ClusterEvent(i,method_name,0.003,scaler)
    #print "Event %i"%i
    for ind in range(i,i+scaler):
    #print "copying to event %i"%ind
        aDataSet.GetTrack(ind)
        #trackX_vs_trackY_plan3.Fill(aDataSet.t_posX[3],aDataSet.t_posY[3])
        #trackX_vs_trackY_plan0.Fill(aDataSet.t_posX[0],aDataSet.t_posY[0])

        for alignement in alignment_constants :
            ApplyAlignment_at_event(ind,aDataSet,[alignement[3],alignement[4],0],[alignement[0],alignement[1],alignement[2]])

        #ApplyAlignment_at_event(ind,aDataSet,[-0.062820313, 0.051976563 , 0],[ 0.000151986,0.000118889,0.266226043])
        #ApplyAlignment_at_event(ind,aDataSet,[-0.0622656250,0.0527500000 , 0],[0.0000048612,0.0000164430,0.2499066771 ])

        aDataSet.FindMatchedCluster(ind,0.5 ,6)
        m,me=aDataSet.ComputeResiduals(ind)
        n_matched+=m
        n_matched_edge+=me
        if ind%1000 ==0 :
            print "Event %d"%ind
            print "Elapsed time/1000 Event : %f s"%(time.time()-last_time)
            last_time = time.time()
            #print h.heap()

    #aDataSet.PrintClusters(i)
    #aDataSet.PrintEvent(i)
print "Found %i matched cluster"%(n_matched)
    
root_file = "%s/Run%i_pyEudetNtuple.root"%(PlotPath,RunNumber)
os.system("rm %s"%root_file)

print "Writing reconstructed data to %s"%root_file
aDataSet.WriteReconstructedData(root_file,6)
