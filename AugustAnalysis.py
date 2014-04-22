from ROOT import *
import ROOT
from ROOT import gStyle
from ROOT import TMath
from math import fsum
from array import array
from EudetData import *
from ToolBox import *

PlotPath = "/afs/cern.ch/work/a/apequegn/public/DESY_TB_DATA_02_07-06-2013_results/pyEudetAnalysisPlots"

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

method_name = "QWeighted"
# method_name = "DigitalCentroid"
# method_name = "maxTOT"
#method_name = "EtaCorrection"

scaler = 1
n_proc= 15000


aDataSet = EudetData("/VertexScratch/TB_Data/DESY_TB_DATA_August2013_results/histo/tbtrackrun000048.root",50000.0,1)

for i in range(0,n_proc,scaler) :
    aDataSet.getEvent(i)
    aDataSet.ClusterEvent(i,method_name,0.003,scaler)
    for ind in range(i,i+scaler):
        aDataSet.GetTrack(ind)
        #aDataSet.FindMatchedCluster(ind, 0.350, 0.350,6)
        #n_matched+=aDataSet.ComputeResiduals(i)
        if ind%1000 ==0 :
            print "Event %d"%ind



#Align track to clusters
ApplyAlignment(aDataSet,[0.97, 0, 0.],[0.0000000000, 0.0000000000,0.0])
resr,rest = Perform3StepAlignment(aDataSet,[[0,360],[0,360],[0,360],[-0.5,0.5],[-0.5,0.5]],n_proc,1)
ApplyAlignment(aDataSet,rest,resr)

for i in range(n_proc) :
#for i in range(aDataSet.p_nEntries) :
    aDataSet.FindMatchedCluster(i,0.1 ,6)
    aDataSet.ComputeResiduals(i)

hx,hy = TrackClusterCorrelation(aDataSet)

#Correlation plot
h1_style(hx,1)
h1_style(hy,1)
cancorx = TCanvas()
hx.Draw("colz")
cancory = TCanvas()
hy.Draw("colz")



# Total Unbiased Residual
resX = TH1D("resX","Unbiased residual X",600,-0.3,0.300)
resY = TH1D("resY","Unbiased residual Y",600,-0.3,0.300)
resX.GetXaxis().SetTitle("X_{track} - X_{Timepix} (mm)")
resX.GetYaxis().SetTitle("Number of hits")
resY.GetXaxis().SetTitle("Y_{track} - Y_{Timepix} (mm)")
resY.GetYaxis().SetTitle("Number of hits")


for j,tracks in enumerate(aDataSet.AllTracks) :
    for track in tracks :
        if track.cluster!=-11 and len(aDataSet.AllClusters[j])!=0 :
            aCluster = aDataSet.AllClusters[j][track.cluster]
            resX.Fill(aCluster.resX)
            resY.Fill(aCluster.resY)

can3 = TCanvas()
resX.Draw("")
#resX.Fit("gaus","R","",-0.03,0.03)

#resX_calib.DrawNormalized("same")

can4 = TCanvas()
resY.Draw("")
#resY.Fit("gaus","R","",-0.03,0.03)
#resY_calib.DrawNormalized("same")
