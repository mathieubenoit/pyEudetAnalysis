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
###############################################################################################################################
# first compile the code in ROOT us
#ing ACLIC: (testLangauFit.C)
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




#aDataSet = EudetData("/VertexScratch/TB_Data/DESY_TB_DATA_02_07-06-2013_results/histo/tbtrackrun000062.root",500.0)
#aDataSet = EudetData("%s/Run%i_pyEudetNtuple.root"%(input_folder,RunNumber),50000.0,edge_width,1,RunNumber,"tbtrack")
#aDataSet.ReadReconstructedData()


aDataSet = EudetData("%s/tbtrackrun%06i.root"%(input_folder,RunNumber),50000.0,edge_width,1,RunNumber,"tbtrack")

#h_chi2,h_chi2ndof = aDataSet.GetChi2Cut(0.95)

#can_chi2 = TCanvas()
#h_chi2.Draw("")
#can_chi2.SetLogx()
#can_chi2.SetLogy()
#
#can_chi2ndof = TCanvas()
#h_chi2ndof.Draw("")
#can_chi2ndof.SetLogx()
#can_chi2ndof.SetLogy()


scaler = 1
#n_proc= 25000


histo_hot,histo_freq = aDataSet.FilterHotPixel(0.01,5000,15)

canhot = TCanvas()
histo_hot.Draw("colz")

canfreq = TCanvas()
canfreq.SetLogx()
canfreq.SetLogy()
histo_freq.Draw("")




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
# histo_hot,histo_freq = aDataSet.FilterHotPixel(0.005,25000)
#histo_hot,histo_freq = aDataSet.FilterHotPixel(0.1,5000,10)
#
#canhot = TCanvas()
#histo_hot.Draw("colz")
#
#canfreq = TCanvas()
#canfreq.SetLogx()
#canfreq.SetLogy()
#histo_freq.Draw("")

n_matched = 0
#for i in range(aDataSet.p_nEntries) :

last_time=time.time()


for i in range(0,n_proc) :
    aDataSet.getEvent(i)
    aDataSet.ClusterEvent(i,method_name,0.003,scaler)
    aDataSet.GetTrack(i)
    if i%1000 ==0 :
        print "Event %d"%i
        print "Elapsed time/1000 Event. Clustering : %f s"%(time.time()-last_time)
        #print h.heap()
        last_time = time.time()


#Rot = [0,0,0]
#Trans=[0,0,0]
#for i in range(0,n_proc-1) :
#    ApplyAlignment_at_event(i,aDataSet,Trans,Rot)


last_time=time.time()

hx,hy = TrackClusterCorrelation(aDataSet,6,50000)

cancorx = TCanvas()
hx.Draw("colz")

cancory = TCanvas()
hy.Draw("colz")

#bla = raw_input()

alignment_constants=PerformPreAlignement(aDataSet,n_proc,1,AlignementPath,6,[0,0,0])

for i in range(0,n_proc-1) :
#    aDataSet.getEvent(i)
#    aDataSet.ClusterEvent(i,method_name,0.003,scaler)
#    #print "Event %i"%i
#    for ind in range(i,i+scaler):
#    #print "copying to event %i"%ind
#        aDataSet.GetTrack(ind)
#        trackX_vs_trackY_plan3.Fill(aDataSet.t_posX[3],aDataSet.t_posY[3])
#        trackX_vs_trackY_plan0.Fill(aDataSet.t_posX[0],aDataSet.t_posY[0])

    for alignement in alignment_constants :
        ApplyAlignment_at_event(i,aDataSet,[alignement[3],alignement[4],0],[alignement[0],alignement[1],alignement[2]])

        #ApplyAlignment_at_event(ind,aDataSet,[0.97, 0, 0.],[0.0000000000, 0.0000000000,0.0])

    aDataSet.FindMatchedCluster(i,0.3 ,6,True)
    a,b=aDataSet.ComputeResiduals(i)
    n_matched+=a
    if i%1000 ==0 :
        print "Event %d"%i
        print "Elapsed time/1000 Event. Apply Alignement and TrackMatching : %f s"%(time.time()-last_time)
        last_time = time.time()
            #print h.heap()
#bleh=PerformPreAlignement(aDataSet,n_proc,1,AlignementPath,6,[0,0,0])


hx,hy = TrackClusterCorrelation(aDataSet,6,50000)

cancorx = TCanvas()
hx.Draw("colz")

cancory = TCanvas()
hy.Draw("colz")

#bla = raw_input()




niter = 2
for i in range(niter) :
    resr,rest = Perform3StepAlignment(aDataSet,[[0,360],[0,360],[0,360],[-0.5,0.5],[-0.5,0.5]],n_proc,1,0.05,AlignementPath,1e-3,[0,0,0])
    ApplyAlignment(aDataSet,rest,resr)


hx,hy = TrackClusterCorrelation(aDataSet,6,50000)

cancorx = TCanvas()
hx.Draw("colz")

cancory = TCanvas()
hy.Draw("colz")

#bla = raw_input()





n_matched = 0

for i in range(0,n_proc-1) :
#    aDataSet.getEvent(i)
#    aDataSet.ClusterEvent(i,method_name,0.003,scaler)
#    #print "Event %i"%i
#    for ind in range(i,i+scaler):
#    #print "copying to event %i"%ind
#        aDataSet.GetTrack(ind)
#        trackX_vs_trackY_plan3.Fill(aDataSet.t_posX[3],aDataSet.t_posY[3])
#        trackX_vs_trackY_plan0.Fill(aDataSet.t_posX[0],aDataSet.t_posY[0])

#    for alignement in alignment_constants :
#        ApplyAlignment_at_event(i,aDataSet,[alignement[3],alignement[4],0],[alignement[0],alignement[1],alignement[2]])

        #ApplyAlignment_at_event(ind,aDataSet,[0.97, 0, 0.],[0.0000000000, 0.0000000000,0.0])

    aDataSet.FindMatchedCluster(i,0.3 ,6,True)
    a,b=aDataSet.ComputeResiduals(i)
    n_matched+=a
    if i%1000 ==0 :
        print "Event %d"%i
        print "Elapsed time/1000 Event. Apply Alignement and TrackMatching : %f s"%(time.time()-last_time)
        last_time = time.time()


#bleh=PerformPreAlignement(aDataSet,n_proc,1,AlignementPath,6,[0,0,0])      
        
print "Found %i matched track-cluster binome"%n_matched
print "Produced Alignment file %s"%AlignementPath
