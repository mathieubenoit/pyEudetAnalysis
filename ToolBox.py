import numpy as np
from scipy.optimize import minimize,basinhopping
from scipy.stats import skew
from ROOT import *
from math import *
from array import array
from EudetData import *
from Constant import *
import time

###############################################################################################################################
#
#                        Box with tools usefull for the analysis
#
###############################################################################################################################

def Chi2Distribution(x,par):
    res = []
    for X in x :
        A = (TMath.Power(TMath.Power(X,2),(-1+par[0]/2)))
        B = (1.0/((TMath.Power(2.,par[0]))*TMath.Gamma(par[0]/2.)))
        C = TMath.Exp(-(X*X)/2.)
        res.append(A*B*C)

    return res


def drange(start, stop, step):
    r = start
    while r < stop:
        yield r
        r += step


def CountClusterSize(dataSet):
    # Count the number of clusters for the different topologies and draw histograms of the results

    n_s1x1y1 = 0.
    n_s2x1y2 = 0.
    n_s2x2y1 = 0.
    n_s2x2y2 = 0.
    n_s3x2y2 = 0.
    n_s4x2y2 = 0.
    n_else = 0.

    hClusterSize = TH1D("hClusterSize","Cluster size distribution",4,1,5)
    hClusterSizeX = TH1D("hClusterSizeX","Cluster size X distribution",4,1,5)
    hClusterSizeY = TH1D("hClusterSizeY","Cluster size Y distribution",4,1,5)
    hClusterSizeXvsSizeY = TH2D("hClusterSizeXvsSizeY","Cluster size X vs Y distribution",4,1,5,4,1,5)

    last_time = time.time()

    for i,tracks in enumerate(dataSet.AllTracks) :
        for track in tracks :
            if track.cluster!=-11 and len(dataSet.AllClusters[i])!=0 :
                aCluster = dataSet.AllClusters[i][track.cluster]

                hClusterSize.Fill(aCluster.size)
                hClusterSizeX.Fill(aCluster.sizeX)
                hClusterSizeY.Fill(aCluster.sizeY)
                hClusterSizeXvsSizeY.Fill(aCluster.sizeX, aCluster.sizeY)

                if(aCluster.size==1) :
                    n_s1x1y1 = n_s1x1y1 + 1.
                elif(aCluster.size==2 and aCluster.sizeY==2 and aCluster.sizeX==1) :
                    n_s2x1y2 = n_s2x1y2 + 1.
                elif(aCluster.size==2 and aCluster.sizeY==1 and aCluster.sizeX==2) :
                    n_s2x2y1 = n_s2x2y1 + 1.
                elif(aCluster.size==2 and aCluster.sizeY==2 and aCluster.sizeX==2) :
                    n_s2x2y2 = n_s2x2y2 + 1.
                elif(aCluster.size==3 and aCluster.sizeY==2 and aCluster.sizeX==2) :
                    n_s3x2y2 = n_s3x2y2 + 1.
                elif(aCluster.size==4 and aCluster.sizeY==2 and aCluster.sizeX==2) :
                    n_s4x2y2 = n_s4x2y2 + 1.
                else :
                    n_else = n_else + 1.

                if(i%10000==0):
                    print "Elapsed time for Counting Cluster Size , event %i: %f s"%(i,(time.time()-last_time))


    n_tot = n_s1x1y1 + n_s2x1y2 + n_s2x2y1 + n_s2x2y2 + n_s3x2y2 + n_s4x2y2 + n_else
    allSizes = [n_s1x1y1 ,n_s2x1y2 , n_s2x2y1 , n_s2x2y2 , n_s3x2y2 , n_s4x2y2 , n_else]
    allSizesPercent = [n_s1x1y1/n_tot*100. ,n_s2x1y2/n_tot*100. , n_s2x2y1/n_tot*100. , n_s2x2y2/n_tot*100. , n_s3x2y2/n_tot*100. , n_s4x2y2/n_tot*100. , n_else/n_tot*100.]
    for size in allSizes :
        print"cluster sizes number..."
        print size
    for size in allSizesPercent :
        print"cluster sizes percentage..."
        print size

    hClusterSize.GetXaxis().SetTitle("Cluster size")
    hClusterSizeX.GetXaxis().SetTitle("Cluster size X")
    hClusterSizeY.GetXaxis().SetTitle("Cluster size Y")
    hClusterSizeXvsSizeY.GetXaxis().SetTitle("Cluster size X")
    hClusterSizeXvsSizeY.GetYaxis().SetTitle("Cluster size Y")

    hClusterSizeCounter = TH1D("hClusterSizeCounter","Number of clusters for different sizes",7,0.,7.)
    hClusterSizeCounter.GetXaxis().SetTitle("Cluster size")
    hClusterSizeCounter.GetYaxis().SetTitle("Number of events")

    hClusterSizeCounter.SetBinContent(1,n_s1x1y1)
    hClusterSizeCounter.SetBinContent(2,n_s2x1y2)
    hClusterSizeCounter.SetBinContent(3,n_s2x2y1)
    hClusterSizeCounter.SetBinContent(4,n_s2x2y2)
    hClusterSizeCounter.SetBinContent(5,n_s3x2y2)
    hClusterSizeCounter.SetBinContent(6,n_s4x2y2)
    hClusterSizeCounter.SetBinContent(7,n_else)
    hClusterSizeCounter.GetXaxis().SetNdivisions(7,kTRUE)
    hClusterSizeCounter.GetXaxis().SetBinLabel(1,"size 1 (1x1)")
    hClusterSizeCounter.GetXaxis().SetBinLabel(2,"size 2 (1x2)")
    hClusterSizeCounter.GetXaxis().SetBinLabel(3,"size 2 (2x1)")
    hClusterSizeCounter.GetXaxis().SetBinLabel(4,"size 2 (2x2)")
    hClusterSizeCounter.GetXaxis().SetBinLabel(5,"size 3 (2x2)")
    hClusterSizeCounter.GetXaxis().SetBinLabel(6,"size 4 (2x2)")
    hClusterSizeCounter.GetXaxis().SetBinLabel(7,"else")
    hClusterSizeCounter.SetStats(0)

    hClusterSizeCounter_percent = TH1D("hClusterSizeCounter_percent","Percentage of clusters for different sizes",7,0.,7.)
    hClusterSizeCounter_percent.GetXaxis().SetTitle("Cluster size")
    hClusterSizeCounter_percent.GetYaxis().SetTitle("Percent of events")

    hClusterSizeCounter_percent.SetBinContent(1,n_s1x1y1/n_tot*100.)
    hClusterSizeCounter_percent.SetBinContent(2,n_s2x1y2/n_tot*100.)
    hClusterSizeCounter_percent.SetBinContent(3,n_s2x2y1/n_tot*100.)
    hClusterSizeCounter_percent.SetBinContent(4,n_s2x2y2/n_tot*100.)
    hClusterSizeCounter_percent.SetBinContent(5,n_s3x2y2/n_tot*100.)
    hClusterSizeCounter_percent.SetBinContent(6,n_s4x2y2/n_tot*100.)
    hClusterSizeCounter_percent.SetBinContent(7,n_else/n_tot*100.)
    hClusterSizeCounter_percent.GetXaxis().SetNdivisions(7,kTRUE)
    hClusterSizeCounter_percent.GetXaxis().SetBinLabel(1,"size 1 (1x1)")
    hClusterSizeCounter_percent.GetXaxis().SetBinLabel(2,"size 2 (1x2)")
    hClusterSizeCounter_percent.GetXaxis().SetBinLabel(3,"size 2 (2x1)")
    hClusterSizeCounter_percent.GetXaxis().SetBinLabel(4,"size 2 (2x2)")
    hClusterSizeCounter_percent.GetXaxis().SetBinLabel(5,"size 3 (2x2)")
    hClusterSizeCounter_percent.GetXaxis().SetBinLabel(6,"size 4 (2x2)")
    hClusterSizeCounter_percent.GetXaxis().SetBinLabel(7,"else")
    hClusterSizeCounter_percent.SetStats(0)

    return hClusterSize, hClusterSizeX, hClusterSizeY, hClusterSizeXvsSizeY, hClusterSizeCounter,hClusterSizeCounter_percent





def RotationMatrix(theta):
    tx,ty,tz = theta
    tx = 2*pi*tx/360.
    ty = 2*pi*ty/360.
    tz = 2*pi*tz/360.
    Rx = np.array([[1,0,0], [0, cos(tx), -sin(tx)], [0, sin(tx), cos(tx)]])
    Ry = np.array([[cos(ty), 0, -sin(ty)], [0, 1, 0], [sin(ty), 0, cos(ty)]])
    Rz = np.array([[cos(tz), -sin(tz), 0], [sin(tz), cos(tz), 0], [0,0,1]])

    return np.dot(Rx, np.dot(Ry, Rz))


def rms(x):
    r"""Compute root mean square."""

    x = np.asanyarray(x)
    rms = np.sqrt(np.sum(x ** 2) / x.size)

    return rms

#
#compute the number of tracks falling in the detector acceptance
#first parameter: a data set (class EudetData)
#second parameter: position of the device under test in the list of planes (the timepix detector in our case)
#
def ComputeDetectorAcceptance(dataSet, dut=6, edges = 0):
    n_tracks_in = 0
    last_time = time.time()
    for i,tracks in enumerate(dataSet.AllTracks) :
        for track in tracks :
            #if (i%1000==0):
                #print "Elapsed time for track in acceptance, event %i: %f s"%(i,(time.time()-last_time))
            if (track.trackX[track.iden.index(dut)]>=(-(pitchX*npix_X)/2.-edges) and track.trackX[track.iden.index(dut)]<=((pitchX*npix_X)/2.+edges)) and (track.trackY[track.iden.index(dut)]>=(-(pitchY*npix_Y)/2.-edges) and track.trackY[track.iden.index(dut)]<=((pitchY*npix_Y)/2.+edges)) :
                n_tracks_in+=1
    return n_tracks_in



def ComputeChargeDistance(dataSet,d=0.005,dut=6):
    # Compute the distance between the track and the edge of the pixel with highest energy
    # d: the value is used to remove the tracks which are in the corner of the pixel

    QrelWrtMindistance = TH2D("QrelWrtMindistance","",110,0,0.0275,100,0.5,1)

    for i,tracks in enumerate(dataSet.AllTracks) :
        for track in tracks :
            if track.cluster!=-11 :
                if(dataSet.AllClusters[i][track.cluster].size==2 and dataSet.AllClusters[i][track.cluster].aspectRatio!=1) :

                    cluster = dataSet.AllClusters[i][track.cluster]

                    # computing the track positions X and Y in the pixel and the relative charge 
                    X = (track.trackX[track.iden.index(dut)])%pitchX
                    Y = (track.trackY[track.iden.index(dut)])%pitchY
                    Qrel = max(cluster.tot) / cluster.totalTOT

                    if ((X<=d and Y<=d) or (X>=(pitchX-d) and Y<=d) or (X>=(pitchX-d) and Y>=(pitchY-d)) or (X<=d and Y>=(pitchY-d))):
                        # removing tracks in the corners of the pixel
                        continue
                    else :
                        # computing the distance between the track and the closest pixel edge
                        if(Y<X and Y<(-X+pitchY)):
                            d = Y
                        elif(Y<X and Y>=(-X+pitchY)):
                            d = pitchX - X
                        elif(Y>=X and Y<(-X+pitchY)):
                            d = X
                        elif(Y>=X and Y>=(-X+pitchY)):
                            d = pitchY - Y

                        QrelWrtMindistance.Fill(d, Qrel)

    return QrelWrtMindistance


#
#compute the correlation in X i.e. X position of the tracks as a function of the X position of the clusters
#param 1: a data set (class EudetData)
#param 2: number of bins in the histograms
#param 3: position of the device under test in the list of planes
#
def HitProbCorrelationX(dataSet,nbin,dut=6):
    HitProb_1_correlationX = TH2D("HitProb_1_correlationX_nbin%i"%nbin,"Hit probability, cluster size 1",nbin,0.,0.055,nbin,0.,0.055)
    #HitProb_1_correlationX.GetXaxis().SetRangeUser(0.,0.055)
    HitProb_1_correlationX.GetXaxis().SetTitle("Track X position within pixel [mm]")
    #HitProb_1_correlationX.GetYaxis().SetRangeUser(0.,0.055)
    HitProb_1_correlationX.GetYaxis().SetTitle("Cluster X position within pixel [mm]")

    HitProb_2_correlationX = TH2D("HitProb_2_correlationX_nbin%i"%nbin,"Hit probability, cluster size 2",nbin,0.,0.055,nbin,0.,0.055)
    #HitProb_2_correlationX.GetXaxis().SetRangeUser(0.,0.055)
    HitProb_2_correlationX.GetXaxis().SetTitle("Track X position within pixel [mm]")
    #HitProb_2_correlationX.GetYaxis().SetRangeUser(0.,0.055)
    HitProb_2_correlationX.GetYaxis().SetTitle("Cluster X position within pixel [mm]")

    HitProb_3_correlationX = TH2D("HitProb_3_correlationX_nbin%i"%nbin,"Hit probability, cluster size 3",nbin,0.,0.055,nbin,0.,0.055)
    #HitProb_3_correlationX.GetXaxis().SetRangeUser(0.,0.055)
    HitProb_3_correlationX.GetXaxis().SetTitle("Track X position within pixel [mm]")
    #HitProb_3_correlationX.GetYaxis().SetRangeUser(0.,0.055)
    HitProb_3_correlationX.GetYaxis().SetTitle("Cluster X position within pixel [mm]")

    HitProb_4_correlationX = TH2D("HitProb_4_correlationX_nbin%i"%nbin,"Hit probability, cluster size 4",nbin,0.,0.055,nbin,0.,0.055)
    #HitProb_4_correlationX.GetXaxis().SetRangeUser(0.,0.055)
    HitProb_4_correlationX.GetXaxis().SetTitle("Track X position within pixel [mm]")
    #HitProb_4_correlationX.GetYaxis().SetRangeUser(0.,0.055)
    HitProb_4_correlationX.GetYaxis().SetTitle("Cluster X position within pixel [mm]")

    for i,tracks in enumerate(dataSet.AllTracks) :
        for track in tracks :
            if track.cluster!=-11 :
                if(dataSet.AllClusters[i][track.cluster].size==1) :
                    HitProb_1_correlationX.Fill((track.trackX[track.iden.index(dut)])%pitchX,(dataSet.AllClusters[i][track.cluster].absX)%pitchX)
                elif(dataSet.AllClusters[i][track.cluster].size==2) :
                    HitProb_2_correlationX.Fill((track.trackX[track.iden.index(dut)])%pitchX,(dataSet.AllClusters[i][track.cluster].absX)%pitchX)
                elif(dataSet.AllClusters[i][track.cluster].size==3 and (dataSet.AllClusters[i][track.cluster].sizeX==2 and dataSet.AllClusters[i][track.cluster].sizeY==2)) :
                    HitProb_3_correlationX.Fill((track.trackX[track.iden.index(dut)])%pitchX,(dataSet.AllClusters[i][track.cluster].absX)%pitchX)
                elif(dataSet.AllClusters[i][track.cluster].size==4 and (dataSet.AllClusters[i][track.cluster].sizeX==2 and dataSet.AllClusters[i][track.cluster].sizeY==2)) :
                    HitProb_4_correlationX.Fill((track.trackX[track.iden.index(dut)])%pitchX,(dataSet.AllClusters[i][track.cluster].absX)%pitchX)

    return HitProb_1_correlationX,HitProb_2_correlationX,HitProb_3_correlationX,HitProb_4_correlationX


#
#compute the correlation in Y i.e. Y position of the tracks as a function of the Y position of the clusters
#param 1: a data set (class EudetData)
#param 2: number of bins in the histograms
#param 3: position of the device under test in the list of planes
#
def HitProbCorrelationY(dataSet,nbin,dut=6):
    HitProb_1_correlationY = TH2D("HitProb_1_correlationY_nbin%i"%nbin,"Hit probability, cluster size 1",nbin,0.,0.055,nbin,0.,0.055)
    #HitProb_1_correlationY.GetXaxis().SetRangeUser(0.,0.055)
    HitProb_1_correlationY.GetXaxis().SetTitle("Track Y position within pixel [mm]")
    #HitProb_1_correlationY.GetYaxis().SetRangeUser(0.,0.055)
    HitProb_1_correlationY.GetYaxis().SetTitle("Cluster Y position within pixel [mm]")

    HitProb_2_correlationY = TH2D("HitProb_2_correlationY_nbin%i"%nbin,"Hit probability, cluster size 2",nbin,0.,0.055,nbin,0.,0.055)
    #HitProb_2_correlationY.GetXaxis().SetRangeUser(0.,0.055)
    HitProb_2_correlationY.GetXaxis().SetTitle("Track Y position within pixel [mm]")
    #HitProb_2_correlationY.GetYaxis().SetRangeUser(0.,0.055)
    HitProb_2_correlationY.GetYaxis().SetTitle("Cluster Y position within pixel [mm]")

    HitProb_3_correlationY = TH2D("HitProb_3_correlationY_nbin%i"%nbin,"Hit probability, cluster size 3",nbin,0.,0.055,nbin,0.,0.055)
    #HitProb_3_correlationY.GetXaxis().SetRangeUser(0.,0.055)
    HitProb_3_correlationY.GetXaxis().SetTitle("Track Y position within pixel [mm]")
    #HitProb_3_correlationY.GetYaxis().SetRangeUser(0.,0.055)
    HitProb_3_correlationY.GetYaxis().SetTitle("Cluster Y position within pixel [mm]")

    HitProb_4_correlationY = TH2D("HitProb_4_correlationY_nbin%i"%nbin,"Hit probability, cluster size 4",nbin,0.,0.055,nbin,0.,0.055)
    #HitProb_4_correlationY.GetXaxis().SetRangeUser(0.,0.055)
    HitProb_4_correlationY.GetXaxis().SetTitle("Track Y position within pixel [mm]")
    #HitProb_4_correlationY.GetYaxis().SetRangeUser(0.,0.055)
    HitProb_4_correlationY.GetYaxis().SetTitle("Cluster Y position within pixel [mm]")

    for i,tracks in enumerate(dataSet.AllTracks) :
        for track in tracks :
            if track.cluster!=-11 :
                if(dataSet.AllClusters[i][track.cluster].size==1) :
                    HitProb_1_correlationY.Fill((track.trackY[track.iden.index(dut)])%pitchY,(dataSet.AllClusters[i][track.cluster].absY)%pitchY)
                elif(dataSet.AllClusters[i][track.cluster].size==2) :
                    HitProb_2_correlationY.Fill((track.trackY[track.iden.index(dut)])%pitchY,(dataSet.AllClusters[i][track.cluster].absY)%pitchY)
                elif(dataSet.AllClusters[i][track.cluster].size==3 and (dataSet.AllClusters[i][track.cluster].sizeX==2 and dataSet.AllClusters[i][track.cluster].sizeY==2)) :
                    HitProb_3_correlationY.Fill((track.trackY[track.iden.index(dut)])%pitchY,(dataSet.AllClusters[i][track.cluster].absY)%pitchY)
                elif(dataSet.AllClusters[i][track.cluster].size==4 and (dataSet.AllClusters[i][track.cluster].sizeX==2 and dataSet.AllClusters[i][track.cluster].sizeY==2)) :
                    HitProb_4_correlationY.Fill((track.trackY[track.iden.index(dut)])%pitchY,(dataSet.AllClusters[i][track.cluster].absY)%pitchY)

    return HitProb_1_correlationY,HitProb_2_correlationY,HitProb_3_correlationY,HitProb_4_correlationY


#
#this function computes the position hit probability coming from the tracks in the pixel plan
#i.e. Track X and Track Y positions within pixel
#it returns 4 histograms for the different cluster sizes
#param 1: a data set (class EudetData)
#param 2: number of bins in the histograms
#param 3: position of the device under test in the list of planes
#
def TrackHitProb(dataSet,nbin,dut=6):
    HitProb_1_track = TH2D("HitProb_1_track_nbin%i"%nbin,"Hit probability, cluster size 1",nbin,0.,0.055,nbin,0.,0.055)
    #HitProb_1_track.GetXaxis().SetRangeUser(0.,0.055)
    HitProb_1_track.GetXaxis().SetTitle("Track X position within pixel [mm]")
    #HitProb_1_track.GetYaxis().SetRangeUser(0.,0.055)
    HitProb_1_track.GetYaxis().SetTitle("Track Y position within pixel [mm]")

    HitProb_2_track = TH2D("HitProb_2_track_nbin%i"%nbin,"Hit probability, cluster size 2",nbin,0.,0.055,nbin,0.,0.055)
    #HitProb_2_track.GetXaxis().SetRangeUser(0.,0.055)
    HitProb_2_track.GetXaxis().SetTitle("Track X position within pixel [mm]")
    #HitProb_2_track.GetYaxis().SetRangeUser(0.,0.055)
    HitProb_2_track.GetYaxis().SetTitle("Track Y position within pixel [mm]")

    HitProb_3_track = TH2D("HitProb_3_track_nbin%i"%nbin,"Hit probability, cluster size 3",nbin,0.,0.055,nbin,0.,0.055)
    #HitProb_3_track.GetXaxis().SetRangeUser(0.,0.055)
    HitProb_3_track.GetXaxis().SetTitle("Track X position within pixel [mm]")
    #HitProb_3_track.GetYaxis().SetRangeUser(0.,0.055)
    HitProb_3_track.GetYaxis().SetTitle("Track Y position within pixel [mm]")

    HitProb_4_track = TH2D("HitProb_4_track_nbin%i"%nbin,"Hit probability, cluster size 4",nbin,0.,0.055,nbin,0.,0.055)
    #HitProb_4_track.GetXaxis().SetRangeUser(0.,0.055)
    HitProb_4_track.GetXaxis().SetTitle("Track X position within pixel [mm]")
    #HitProb_4_track.GetYaxis().SetRangeUser(0.,0.055)
    HitProb_4_track.GetYaxis().SetTitle("Track Y position within pixel [mm]")

    last_time = time.time()

    for i,tracks in enumerate(dataSet.AllTracks) :
        if(i%1000==0):
            print "Elapsed time for Hitprob Calculation , event %i: %f s"%(i,(time.time()-last_time))
        for track in tracks :
            if track.cluster!=-11 :
                if(dataSet.AllClusters[i][track.cluster].size==1) :
                    HitProb_1_track.Fill((track.trackX[track.iden.index(dut)])%pitchX,(track.trackY[track.iden.index(dut)])%pitchY)
                elif(dataSet.AllClusters[i][track.cluster].size==2) :
                    HitProb_2_track.Fill((track.trackX[track.iden.index(dut)])%pitchX,(track.trackY[track.iden.index(dut)])%pitchY)
                elif(dataSet.AllClusters[i][track.cluster].size==3 and (dataSet.AllClusters[i][track.cluster].sizeX==2 and dataSet.AllClusters[i][track.cluster].sizeY==2)) :
                    HitProb_3_track.Fill((track.trackX[track.iden.index(dut)])%pitchX,(track.trackY[track.iden.index(dut)])%pitchY)
                elif(dataSet.AllClusters[i][track.cluster].size==4 and (dataSet.AllClusters[i][track.cluster].sizeX==2 and dataSet.AllClusters[i][track.cluster].sizeY==2)) :
                    HitProb_4_track.Fill((track.trackX[track.iden.index(dut)])%pitchX,(track.trackY[track.iden.index(dut)])%pitchY)

    return HitProb_1_track,HitProb_2_track,HitProb_3_track,HitProb_4_track

def TOTProfile(dataSet,nbin,dut=6):
    TOTProfileX_1 = TH2D("TOTProfileX_1_nbin%i"%nbin,"Hit probability, cluster size 1",nbin,0.,0.055,1000,0.,1000)
    #TOTProfileX_1.GetXaxis().SetRangeUser(0.,0.055)
    TOTProfileX_1.GetXaxis().SetTitle("Track X position within pixel [mm]")
    #TOTProfileX_1.GetYaxis().SetRangeUser(0.,0.055)
    TOTProfileX_1.GetYaxis().SetTitle("TOT (A.U.)]")

    TOTProfileX_2 = TH2D("TOTProfileX_2_nbin%i"%nbin,"Hit probability, cluster size 2",nbin,0.,0.055,1000,0.,1000)
    #TOTProfileX_2.GetXaxis().SetRangeUser(0.,0.055)
    TOTProfileX_2.GetXaxis().SetTitle("Track X position within pixel [mm]")
    #TOTProfileX_2.GetYaxis().SetRangeUser(0.,0.055)
    TOTProfileX_2.GetYaxis().SetTitle("TOT (A.U.)]")

    TOTProfileX_3 = TH2D("TOTProfileX_3_nbin%i"%nbin,"Hit probability, cluster size 3",nbin,0.,0.055,1000,0.,1000)
    #TOTProfileX_3.GetXaxis().SetRangeUser(0.,0.055)
    TOTProfileX_3.GetXaxis().SetTitle("Track X position within pixel [mm]")
    #TOTProfileX_3.GetYaxis().SetRangeUser(0.,0.055)
    TOTProfileX_3.GetYaxis().SetTitle("TOT (A.U.)]")

    TOTProfileX_4 = TH2D("TOTProfileX_4_nbin%i"%nbin,"Hit probability, cluster size 4",nbin,0.,0.055,1000,0.,1000)
    #TOTProfileX_4.GetXaxis().SetRangeUser(0.,0.055)
    TOTProfileX_4.GetXaxis().SetTitle("Track X position within pixel [mm]")
    #TOTProfileX_4.GetYaxis().SetRangeUser(0.,0.055)
    TOTProfileX_4.GetYaxis().SetTitle("TOT (A.U.)]")
  
 
    
    TOTProfileY_1 = TH2D("TOTProfileY_1_nbin%i"%nbin,"Hit probability, cluster size 1",nbin,0.,0.055,1000,0.,1000)
    #TOTProfileY_1.GetXaxis().SetRangeUser(0.,0.055)
    TOTProfileY_1.GetXaxis().SetTitle("Track X position within pixel [mm]")
    #TOTProfileY_1.GetYaxis().SetRangeUser(0.,0.055)
    TOTProfileY_1.GetYaxis().SetTitle("TOT (A.U.)]")

    TOTProfileY_2 = TH2D("TOTProfileY_2_nbin%i"%nbin,"Hit probability, cluster size 2",nbin,0.,0.055,1000,0.,1000)
    #TOTProfileY_2.GetXaxis().SetRangeUser(0.,0.055)
    TOTProfileY_2.GetXaxis().SetTitle("Track X position within pixel [mm]")
    #TOTProfileY_2.GetYaxis().SetRangeUser(0.,0.055)
    TOTProfileY_2.GetYaxis().SetTitle("TOT (A.U.)]")

    TOTProfileY_3 = TH2D("TOTProfileY_3_nbin%i"%nbin,"Hit probability, cluster size 3",nbin,0.,0.055,1000,0.,1000)
    #TOTProfileY_3.GetXaxis().SetRangeUser(0.,0.055)
    TOTProfileY_3.GetXaxis().SetTitle("Track X position within pixel [mm]")
    #TOTProfileY_3.GetYaxis().SetRangeUser(0.,0.055)
    TOTProfileY_3.GetYaxis().SetTitle("TOT (A.U.)]")

    TOTProfileY_4 = TH2D("TOTProfileY_4_nbin%i"%nbin,"Hit probability, cluster size 4",nbin,0.,0.055,1000,0.,1000)
    #TOTProfileY_4.GetXaxis().SetRangeUser(0.,0.055)
    TOTProfileY_4.GetXaxis().SetTitle("Track X position within pixel [mm]")
    #TOTProfileY_4.GetYaxis().SetRangeUser(0.,0.055)
    TOTProfileY_4.GetYaxis().SetTitle("TOT (A.U.)]")
    
    TOTProfileX = TH2D("TOTProfileX_nbin%i"%nbin,"Hit probability, cluster size 4",nbin,0.,0.055,1000,0.,2000)
    #TOTProfileY_4.GetXaxis().SetRangeUser(0.,0.055)
    TOTProfileX.GetXaxis().SetTitle("Track X position within pixel [mm]")
    #TOTProfileY_4.GetYaxis().SetRangeUser(0.,0.055)
    TOTProfileX.GetYaxis().SetTitle("TOT (A.U.)]")  
    
    TOTProfileY = TH2D("TOTProfileY_nbin%i"%nbin,"Hit probability, cluster size 4",nbin,0.,0.055,1000,0.,2000)
    #TOTProfileY_4.GetXaxis().SetRangeUser(0.,0.055)
    TOTProfileY.GetXaxis().SetTitle("Track X position within pixel [mm]")
    #TOTProfileY_4.GetYaxis().SetRangeUser(0.,0.055)
    TOTProfileY.GetYaxis().SetTitle("TOT (A.U.)]")     
    
    TOTProfile= TH2D("TOTProfile","2D TOT Profile",nbin,0.,0.055,nbin,0.,0.055)
    TOTProfileN= TH2D("TOTProfileN","2D TOT Profile N",nbin,0.,0.055,nbin,0.,0.055)  
    last_time = time.time()
    
    

    for i,tracks in enumerate(dataSet.AllTracks) :
        if(i%1000==0):
            print "Elapsed time for Hitprob Calculation , event %i: %f s"%(i,(time.time()-last_time))
        for track in tracks :
            if track.cluster!=-11 :
                TOTProfileX.Fill((track.trackX[track.iden.index(dut)])%pitchX,dataSet.AllClusters[i][track.cluster].totalTOT)                   
                TOTProfileY.Fill((track.trackY[track.iden.index(dut)])%pitchX,dataSet.AllClusters[i][track.cluster].totalTOT)
                TOTProfile.Fill((track.trackX[track.iden.index(dut)])%pitchX,(track.trackY[track.iden.index(dut)])%pitchY,dataSet.AllClusters[i][track.cluster].totalTOT) 
                TOTProfileN.Fill((track.trackX[track.iden.index(dut)])%pitchX,(track.trackY[track.iden.index(dut)])%pitchY)                              
                if(dataSet.AllClusters[i][track.cluster].size==1) :
                    TOTProfileX_1.Fill((track.trackX[track.iden.index(dut)])%pitchX,dataSet.AllClusters[i][track.cluster].totalTOT)
                    TOTProfileY_1.Fill((track.trackY[track.iden.index(dut)])%pitchX,dataSet.AllClusters[i][track.cluster].totalTOT)
                    
                                          
                elif(dataSet.AllClusters[i][track.cluster].size==2) :
                    TOTProfileX_2.Fill((track.trackX[track.iden.index(dut)])%pitchX,dataSet.AllClusters[i][track.cluster].totalTOT)
                    TOTProfileY_2.Fill((track.trackY[track.iden.index(dut)])%pitchX,dataSet.AllClusters[i][track.cluster].totalTOT)
                    
                    
                elif(dataSet.AllClusters[i][track.cluster].size==3 and (dataSet.AllClusters[i][track.cluster].sizeX==2 and dataSet.AllClusters[i][track.cluster].sizeY==2)) :
                    TOTProfileX_3.Fill((track.trackX[track.iden.index(dut)])%pitchX,dataSet.AllClusters[i][track.cluster].totalTOT)
                    TOTProfileY_3.Fill((track.trackY[track.iden.index(dut)])%pitchX,dataSet.AllClusters[i][track.cluster].totalTOT)
                    
                    
                elif(dataSet.AllClusters[i][track.cluster].size==4 and (dataSet.AllClusters[i][track.cluster].sizeX==2 and dataSet.AllClusters[i][track.cluster].sizeY==2)) :
                    TOTProfileX_4.Fill((track.trackX[track.iden.index(dut)])%pitchX,dataSet.AllClusters[i][track.cluster].totalTOT)
                    TOTProfileY_4.Fill((track.trackY[track.iden.index(dut)])%pitchX,dataSet.AllClusters[i][track.cluster].totalTOT)
                    
    TOTProfileAv = TOTProfile.Clone("efficiency")
    TOTProfileAv.Divide(TOTProfile,TOTProfileN,1.,1.,"B")             

    return TOTProfileX_1,TOTProfileX_2,TOTProfileX_3,TOTProfileX_4,TOTProfileY_1,TOTProfileY_2,TOTProfileY_3,TOTProfileY_4,TOTProfileX,TOTProfileY,TOTProfileAv
#
#this function computes the position hit probability reconstructed from the clusters in the pixel plan
#i.e. cluster X and cluster Y positions within pixel
#it returns 4 histograms for the different cluster sizes
#param 1: a data set (class EudetData)
#param 2: number of bins in the histograms
#param 3: position of the device under test in the list of planes
#
def ClusterHitProb(dataSet,nbin,dut=6):
    HitProb_1_cluster = TH2D("HitProb_1_cluster_nbin%i"%nbin,"Hit probability, cluster size 1",nbin,0.,0.055,nbin,0.,0.055)
    #HitProb_1_cluster.GetXaxis().SetRangeUser(0.,0.055)
    HitProb_1_cluster.GetXaxis().SetTitle("Cluster X position within pixel [mm]")
    #HitProb_1_cluster.GetYaxis().SetRangeUser(0.,0.055)
    HitProb_1_cluster.GetYaxis().SetTitle("Cluster Y position within pixel [mm]")

    HitProb_2_cluster = TH2D("HitProb_2_cluster_nbin%i"%nbin,"Hit probability, cluster size 2",nbin,0.,0.055,nbin,0.,0.055)
    #HitProb_2_cluster.GetXaxis().SetRangeUser(0.,0.055)
    HitProb_2_cluster.GetXaxis().SetTitle("Cluster X position within pixel [mm]")
    #HitProb_2_cluster.GetYaxis().SetRangeUser(0.,0.055)
    HitProb_2_cluster.GetYaxis().SetTitle("Cluster Y position within pixel [mm]")

    HitProb_3_cluster = TH2D("HitProb_3_cluster_nbin%i"%nbin,"Hit probability, cluster size 3",nbin,0.,0.055,nbin,0.,0.055)
    #HitProb_3_cluster.GetXaxis().SetRangeUser(0.,0.055)
    HitProb_3_cluster.GetXaxis().SetTitle("Cluster X position within pixel [mm]")
    #HitProb_3_cluster.GetYaxis().SetRangeUser(0.,0.055)
    HitProb_3_cluster.GetYaxis().SetTitle("Cluster Y position within pixel [mm]")

    HitProb_4_cluster = TH2D("HitProb_4_cluster_nbin%i"%nbin,"Hit probability, cluster size 4",nbin,0.,0.055,nbin,0.,0.055)
    #HitProb_4_cluster.GetXaxis().SetRangeUser(0.,0.055)
    HitProb_4_cluster.GetXaxis().SetTitle("Cluster X position within pixel [mm]")
    #HitProb_4_cluster.GetYaxis().SetRangeUser(0.,0.055)
    HitProb_4_cluster.GetYaxis().SetTitle("Cluster Y position within pixel [mm]")



    for i,tracks in enumerate(dataSet.AllTracks) :
        for track in tracks :
            if track.cluster!=-11 :
#                 print 'len(dataSet.AllClusters[i]) : '
#                 print len(dataSet.AllClusters[i])
#                 print 'track.cluster : '
#                 print track.cluster
                if(dataSet.AllClusters[i][track.cluster].size==1) :
                    HitProb_1_cluster.Fill((dataSet.AllClusters[i][track.cluster].absX)%pitchX,(dataSet.AllClusters[i][track.cluster].absY)%pitchY)
                elif(dataSet.AllClusters[i][track.cluster].size==2) :
                    HitProb_2_cluster.Fill((dataSet.AllClusters[i][track.cluster].absX)%pitchX,(dataSet.AllClusters[i][track.cluster].absY)%pitchY)
                elif((dataSet.AllClusters[i][track.cluster].size==3) and (dataSet.AllClusters[i][track.cluster].sizeX==2 and dataSet.AllClusters[i][track.cluster].sizeY==2)) :
                    HitProb_3_cluster.Fill((dataSet.AllClusters[i][track.cluster].absX)%pitchX,(dataSet.AllClusters[i][track.cluster].absY)%pitchY)
                elif(dataSet.AllClusters[i][track.cluster].size==4 and (dataSet.AllClusters[i][track.cluster].sizeX==2 and dataSet.AllClusters[i][track.cluster].sizeY==2)) :
                    HitProb_4_cluster.Fill((dataSet.AllClusters[i][track.cluster].absX)%pitchX,(dataSet.AllClusters[i][track.cluster].absY)%pitchY)

    return HitProb_1_cluster,HitProb_2_cluster,HitProb_3_cluster,HitProb_4_cluster


def TrackClusterCorrelation(dataSet,dut=6,imax=1000):

    histox = TH2D("corX","Track-cluster x correlation",npix_X,-(npix_X)*pitchX/2.,npix_X*pitchX/2.,npix_X,-(npix_X)*pitchX/2.,npix_X*pitchX/2.)
    histoy = TH2D("corY","Track-cluster y correlation",npix_Y,-(npix_Y)*pitchY/2.,npix_Y*pitchY/2.,npix_Y,-(npix_Y)*pitchY/2.,npix_Y*pitchY/2.)

    for h in [histox,histoy] :
        h.GetXaxis().SetTitle("Cluster position (mm)")
        h.GetYaxis().SetTitle("Track position (mm)")

    last_time = time.time()
    for i,tracks in enumerate(dataSet.AllTracks[0:imax]) :
        if(i%1000==0) :
            print "Correlation, event %i %f s elapsed"%(i,time.time()-last_time)

        for track in tracks :
            for cluster in dataSet.AllClusters[i]:
                histox.Fill(cluster.absX, track.trackX[track.iden.index(dut)])
                histoy.Fill(cluster.absY, track.trackY[track.iden.index(dut)])

    return histox,histoy



def TotalMeanFunctionX(Translations,Rotations,aDataDet,nevents,skip,cut = 0.1,dut=6):

    totaldist_evaluator = 0.
    n = 0
    dist_tmp_x = []
    dist_tmp_y = []
    rotationMatrix = RotationMatrix(Rotations)
    h_dist_x_1 = TH1D("h_dist_x_1","TotalMeanFunctionX: dist_x",8000,-4.,4.)
    h_dist_y_1 = TH1D("h_dist_y_1","TotalMeanFunctionX: dist_y",8000,-4.,4.)
    cut2 = cut**2
    for i,clusters in enumerate(aDataDet.AllClusters[0:nevents]) :
        for index,cluster in enumerate(clusters) :
            if i%skip==0 :
                for track in aDataDet.AllTracks[i] :
                    tmp=np.dot(rotationMatrix,[track.trackX[track.iden.index(dut)],track.trackY[track.iden.index(dut)],0])
                    tmp[0] = tmp[0] + Translations[0]
                    tmp[1] = tmp[1]
                    distx=cluster.absX -tmp[0]
                    disty=cluster.absY -tmp[1]
                    h_dist_x_1.Fill(distx)
                    h_dist_y_1.Fill(disty)
                    dist_tmp_x.append(distx)
                    dist_tmp_y.append(disty)

    maxx_bin = h_dist_x_1.GetMaximumBin()
    maxx = h_dist_x_1.GetXaxis().GetBinCenter(maxx_bin)
    maxy_bin = h_dist_y_1.GetMaximumBin()
    maxy = h_dist_y_1.GetXaxis().GetBinCenter(maxy_bin)

    for index,eventx in enumerate(dist_tmp_x) :
        if((eventx)**2 < cut2) :
            totaldist_evaluator+=eventx
            n+=1

    print "Evaluating for Trans : %.9f %.9f  [mm] metric = %.9f  n = %i"%(Translations[0],0,fabs(totaldist_evaluator/n),n)
    return fabs(totaldist_evaluator/n)


def TotalMeanFunctionY(Translations,Tx,Rotations,aDataDet,nevents,skip,cut = 0.1,dut=6):

    totaldist_evaluator = 0.
    n = 0
    dist_tmp_x = []
    dist_tmp_y = []
    rotationMatrix = RotationMatrix(Rotations)
    h_dist_x_2 = TH1D("h_dist_x_2","TotalMeanFunctionY: dist_x",8000,-4.,4.)
    h_dist_y_2 = TH1D("h_dist_y_2","TotalMeanFunctionY: dist_y",8000,-4.,4.)
    cut2 = cut**2
    for i,clusters in enumerate(aDataDet.AllClusters[0:nevents]) :
        for index,cluster in enumerate(clusters) :
            if i%skip==0 :
                for track in aDataDet.AllTracks[i] :
                    tmp=np.dot(rotationMatrix,[track.trackX[track.iden.index(dut)],track.trackY[track.iden.index(dut)],0])
                    tmp[0] = tmp[0] + Tx
                    tmp[1] = tmp[1] + Translations[0]
                    distx=cluster.absX -tmp[0]
                    disty=cluster.absY -tmp[1]
                    h_dist_x_2.Fill(distx)
                    h_dist_y_2.Fill(disty)
                    dist_tmp_x.append(distx)
                    dist_tmp_y.append(disty)

    maxx_bin = h_dist_x_2.GetMaximumBin()
    maxx = h_dist_x_2.GetXaxis().GetBinCenter(maxx_bin)
    maxy_bin = h_dist_y_2.GetMaximumBin()
    maxy = h_dist_y_2.GetXaxis().GetBinCenter(maxy_bin)

    for index,eventy in enumerate(dist_tmp_y) :
        if((eventy)**2 < cut2) :
            totaldist_evaluator+=eventy
            n+=1

    print "Evaluating for Trans : %.9f %.9f  [mm] metric = %.9f  n = %i"%(Tx,Translations[0],fabs(totaldist_evaluator/n),n)
    return fabs(totaldist_evaluator/n)


def TotalRotationFunction(Rotations,Translations,aDataDet,nevents,skip=1,cut = 0.1,dut=6):

    n = 0
    dist_tmp_x = []
    dist_tmp_y = []
    dist_good_x = []
    dist_good_y = []
    rotationMatrix = RotationMatrix(Rotations)
    h_dist_x_3 = TH1D("h_dist_x_3","TotalRotationFunction: dist_x",8000,-4.,4.)
    h_dist_y_3 = TH1D("h_dist_y_3","TotalRotationFunction: dist_y",8000,-4.,4.)

    for i,clusters in enumerate(aDataDet.AllClusters[0:nevents]) :
        for cluster in clusters :
            if i%skip==0 :
                for track in aDataDet.AllTracks[i] :

                    tmp=np.dot(rotationMatrix,[track.trackX[track.iden.index(dut)],track.trackY[track.iden.index(dut)],0])
                    tmp[0] = tmp[0] + Translations[0]
                    tmp[1] = tmp[1] + Translations[1]
                    distx=cluster.absX -tmp[0]
                    disty=cluster.absY -tmp[1]

                    if((distx*distx + disty*disty) < (cut*cut)):                        
                        dist_tmp_x.append(distx)
                        dist_tmp_y.append(disty)
                        n+=1
                    h_dist_x_3.Fill(distx)
                    h_dist_y_3.Fill(disty)


    result=sqrt(rms(dist_tmp_x)**2 + rms(dist_tmp_x)**2)/n
    print "Evaluating for Rotation : %.9f %.9f %.9f [deg] Trans : %f %f  [mm] metric = %.9f  n = %i"%(Rotations[0],Rotations[1],Rotations[2],Translations[0],Translations[1],result,n)
    return result


#
#return the X resolution (sigma of the X residuals distribution) for a given value of the sigma for the charge sharing (eta correction)
#param 1: value of the sigma for the charge sharing (eta correction)
#param 2: a data set (class EudetData)
#param 3: number of skiped events (compute residuals for 1 event over 'skip' events)
#param 4: position of the device under test in the list of planes
#
def TotalSigmaFunctionX(sigmaCharge_tmp_X,sigmaCharge_tmp_Y,dataSet,skip,dut=6):

    tmpx = TH1D("resX_2","Unbiased residual X, cluster size Y = 2",300,-0.150,0.150)
    
    
    for j,tracks in enumerate(dataSet.AllTracks) :
        if j%skip==0 :
            for track in tracks :
                if track.cluster!=-11 and len(dataSet.AllClusters[j])!=0 :
                    aCluster = dataSet.AllClusters[j][track.cluster]
                    if(aCluster.sizeX==2 and aCluster.sizeY==1) :
                        aCluster.GetEtaCorrectedQWeightedCentroid(fabs(sigmaCharge_tmp_X),fabs(sigmaCharge_tmp_Y))
                        dataSet.ComputeResiduals(j)
                        tmpx.Fill(aCluster.resX)

#     for i in range(dataSet.p_nEntries) :
#         if i%skip==0 :
#             dataSet.ClusterEvent(i,method_name,sigmaCharge_tmp_Y)
#             dataSet.ComputeResiduals(i)
#    resY_cs = []
#    n_cs = 3
#    for i in range(1,n_cs+2) : #n_cs+2 excluded
#        tmpy = TH1D("resY_%i"%i,"Unbiased residual Y, cluster size Y = %i"%i,300,-0.150,0.150)
#        tmpy.GetXaxis().SetTitle("Y_{track} - Y_{Timepix} (mm)")
#        tmpy.GetYaxis().SetTitle("Number of hits")
#        tmpy.SetLineColor(i)
#        resY_cs.append(tmpy)
#    for j,tracks in enumerate(dataSet.AllTracks) :
#        if j%skip==0 :
#            for track in tracks :
#                if track.cluster!=-11 and len(dataSet.AllClusters[j])!=0 :
#                    aCluster = dataSet.AllClusters[j][track.cluster]
#                    for i in range(1,n_cs+2) :
#                        if(aCluster.sizeY==i and aCluster.sizeX==1) :
#                            resY_cs[i-1].Fill(aCluster.resY)

    g2 = TF1("m1","gaus",-0.03,0.03)
    if tmpx.GetEffectiveEntries() > 10.:
        rX = tmpx.Fit(g2,"RS","")
        sigmaResX = rX.Parameter(2)
        sigmaResX_err = rX.ParError(2)
        print "resolution = %f for sigma=%f"%(sigmaResX,sigmaCharge_tmp_X)
    else:
        print "not enough events to fit."
        sigmaResX = 9.
    return sigmaResX




#
#return the Y resolution (sigma of the Y residuals distribution) for a given value of the sigma for the charge sharing (eta correction)
#param 1: value of the sigma for the charge sharing (eta correction)
#param 2: a data set (class EudetData)
#param 3: number of skiped events (compute residuals for 1 event over 'skip' events)
#param 4: position of the device under test in the list of planes
#
def TotalSigmaFunctionY(sigmaCharge_tmp_Y,sigmaCharge_tmp_X,dataSet,skip,dut=6):
    
    tmpy = TH1D("resY_1","Unbiased residual Y, cluster size Y = 2",300,-0.150,0.150)
    for j,tracks in enumerate(dataSet.AllTracks) :
        if j%skip==0 :
            for track in tracks :
                if track.cluster!=-11 and len(dataSet.AllClusters[j])!=0 :
                    aCluster = dataSet.AllClusters[j][track.cluster]
                    if(aCluster.sizeY==2 and aCluster.sizeX==1) :
                        aCluster.GetEtaCorrectedQWeightedCentroid(fabs(sigmaCharge_tmp_X),fabs(sigmaCharge_tmp_Y))
                        dataSet.ComputeResiduals(j)
                        tmpy.Fill(aCluster.resY)

#     for i in range(dataSet.p_nEntries) :
#         if i%skip==0 :
#             dataSet.ClusterEvent(i,method_name,sigmaCharge_tmp_Y)
#             dataSet.ComputeResiduals(i)
#    resY_cs = []
#    n_cs = 3
#    for i in range(1,n_cs+2) : #n_cs+2 excluded
#        tmpy = TH1D("resY_%i"%i,"Unbiased residual Y, cluster size Y = %i"%i,300,-0.150,0.150)
#        tmpy.GetXaxis().SetTitle("Y_{track} - Y_{Timepix} (mm)")
#        tmpy.GetYaxis().SetTitle("Number of hits")
#        tmpy.SetLineColor(i)
#        resY_cs.append(tmpy)
#    for j,tracks in enumerate(dataSet.AllTracks) :
#        if j%skip==0 :
#            for track in tracks :
#                if track.cluster!=-11 and len(dataSet.AllClusters[j])!=0 :
#                    aCluster = dataSet.AllClusters[j][track.cluster]
#                    for i in range(1,n_cs+2) :
#                        if(aCluster.sizeY==i and aCluster.sizeX==1) :
#                            resY_cs[i-1].Fill(aCluster.resY)

    g2 = TF1("m1","gaus",-0.03,0.03)

    rY = tmpy.Fit(g2,"RS","")
    sigmaResY = rY.Parameter(2)
    sigmaResY_err = rY.ParError(2)
    print "resolution = %f for sigma=%f"%(sigmaResY,sigmaCharge_tmp_Y)
    return sigmaResY



def TotalDistanceFunction(parameters,aDataDet,nevents,skip,cutx = 0.1, cuty = 0.1,dut=6):

    totaldist_evaluator = 0.
    n = 0
    dist_tmp_x = []
    dist_tmp_y = []
    for i,clusters in enumerate(aDataDet.AllClusters[0:nevents]) :
        for index,cluster in enumerate(clusters) :
            if i%skip==0 :

                for track in aDataDet.AllTracks[i] :
                    #print len(aDataDet.AllTracks[i]),track.cluster
                    #cluster = aDataDet.AllTracks[i][track.cluster]
                    tmp=np.dot(RotationMatrix(parameters[0:3]),[track.trackX[track.iden.index(dut)],track.trackX[track.iden.index(dut)],0])
                    tmp[0] = tmp[0] + parameters[3]
                    tmp[1] = tmp[1] + parameters[4]


                    #dist=sqrt(pow(cluster.absX-tmp[0],2)+pow(cluster.absY.-tmp[1],2))
                    distx=cluster.absX -tmp[0]
                    disty=cluster.absY -tmp[1]

                    if(fabs(distx)<cutx and fabs(disty)<cuty):
                        dist_tmp_x.append(distx)
                        dist_tmp_y.append(disty)
                        totaldist_evaluator+=distx
                        n+=1

    if(n!=0):
        result = fabs(totaldist_evaluator/n)
    else :
        result = 1000.
    print "Evaluating for Rotation : %f %f %f [deg] Trans : %f %f  [mm] metric = %f  n = %i"%(parameters[0],parameters[1],parameters[2],parameters[3],parameters[4],result,n)
    return result

def ReadAlignment(filename) :
    f = open(filename,'r')
    alignments = []

    for line in f.readlines() :
        result = line.split(" ")
        #print result
        alignments.append([float(result[2]),float(result[3]),float(result[4]),float(result[8]),float(result[9])])
    return alignments


def PerformPreAlignment(aDataSet,nevents,skip=1,filename='Alignment.txt',dut=6,Rotations=[0,0,0]):
    # will plot the x and y distance between every cluster-track pair in an event
    # returns the central values of the x, y bins with the most entries, and hists

    last_time=time.time()
    h_dist_x = TH1D("h_dist_x","Track-cluster distance x",800,-20.,20.)
    h_dist_y = TH1D("h_dist_y","Track-cluster distance y",800,-20.,20.)
                    
    for i in xrange(nevents):
        if i%skip == 0:
            clusters = aDataSet.AllClusters[i]
            for cluster in clusters:
                for track in aDataSet.AllTracks[i]:

                    tmp=[track.trackX[track.iden.index(dut)],track.trackY[track.iden.index(dut)],0]
                    tmp=np.dot(RotationMatrix(Rotations),tmp)
            
                    distx=cluster.absX -tmp[0]
                    disty=cluster.absY -tmp[1]

                    h_dist_x.Fill(distx)
                    h_dist_y.Fill(disty)

        if i%1000 == 0 :
            print "Event %d"%i
            print "Elapsed time/1000 Event. Prealignment : %f s"%(time.time()-last_time)
            last_time = time.time()

    maxx_bin = h_dist_x.GetMaximumBin()
    maxx = h_dist_x.GetXaxis().GetBinCenter(maxx_bin)
    maxy_bin = h_dist_y.GetMaximumBin()
    maxy = h_dist_y.GetXaxis().GetBinCenter(maxy_bin)
    
    f = open(filename,'w')
    f.write("Rotation : %f %f %f [deg] Trans : %f %f  [mm] \n"%(Rotations[0],Rotations[1],Rotations[2],maxx,maxy))
    f.close()

    print "Prealignment yield Translations : %.9f %.9f  [mm]  Rotation : %f %f %f [deg] "%(maxx,maxy,Rotations[0],Rotations[1],Rotations[2])
    print "Time for Prealignment : %f s"%(time.time()-last_time)

    return [[Rotations[0],Rotations[1],Rotations[2],maxx,maxy]], h_dist_x, h_dist_y



def PerformAlignement(aDataSet, boundary) :
    x0 = np.array([0.,0.,0.,0.,0])
    res = minimize(TotalDistanceFunction,x0,[aDataSet,3000],method='Nelder-Mead',options={'disp': True})
    return res.x[0:3],res.x[3:]


def Perform3StepAlignment(aDataSet,boundary,nevent,skip,cut = 0.1,filename='Alignment.txt',gtol=1e-5,Rotations=[0,0,0]) :
    x_tx = np.array([0.])
    x_ty = np.array([0.])
    xr= np.array(Rotations)

    sigmas = []
    rzs = drange(-1,1,0.01)
    theRs = [x for x in rzs]
    for rZ in theRs:
        aSigma=TotalRotationFunction([0,0,rZ],[0,0,0],aDataSet,nevent,skip,cut,dut=6)
        if (isnan(aSigma) or isinf(aSigma)) : 
            sigmas.append(1e7)
        else : 
            sigmas.append(aSigma)
  
    rZ = theRs[sigmas.index(min(sigmas))]
    
    print "Starting guess Z angle : %f"%rZ
      
    xr= np.array([0,0,rZ])
    
    argTuple = [x_tx,x_ty],aDataSet,nevent,skip,cut        
    resr = minimize(TotalRotationFunction,xr,argTuple,method='BFGS',options={'disp': True,'gtol': gtol , 'eps':0.5, 'maxiter' : 15 })
    print "resr", resr
    print "best guess for rotation matrix: resr.x", resr.x
    if resr.success == False:
        print "Minimisation didn't converge. Press any key to continue, ctrl D to exit"
        b=raw_input()


    argTuple = resr.x,aDataSet,nevent,skip,cut  
    rest = minimize(TotalMeanFunctionX,x_tx, argTuple,method='Nelder-Mead',options={'xtol': 1e-5,'disp': True})
    print "rest", rest
    print "best guess for x translation: rest.x", rest.x
    if rest.success == False:
        print "Minimisation didn't converge. Press any key to continue, ctrl D to exit"
        b=raw_input() 

    argTuple = rest.x[0],resr.x,aDataSet,nevent,skip,cut
    rest2= minimize(TotalMeanFunctionY,x_ty, argTuple,method='Nelder-Mead',options={'xtol': 1e-5,'disp': True})
    print "rest2", rest2
    print "best guess for y translation: rest2.x", rest2.x
    if rest2.success == False:
        print "Minimisation didn't converge. Press any key to continue, ctrl D to exit"
        b=raw_input() 

    f = open(filename,'a')
    f.write("Rotation : %f %f %f [deg] Trans : %f %f  [mm] \n"%(resr.x[0],resr.x[1],resr.x[2],rest.x[0],rest2.x[0]))
    f.close()

    return resr.x ,[rest.x[0],rest2.x[0],0]

def Perform2StepAlignment(aDataSet,boundary,nevent,skip,cut = 0.1,filename='Alignment.txt') :
    x_tx = np.array([0.])
    x_ty = np.array([0.])
    xr= np.array([0.,0.,0.])
    #resr = minimize(TotalRotationFunction,xr,[x_tx,aDataSet,nevent,skip],method='Nelder-Mead',options={'xtol': 1e-5,'disp': True})

    rest = minimize(TotalMeanFunctionX,x_tx,[xr,aDataSet,nevent,skip,cut],method='Nelder-Mead',options={'xtol': 1e-5,'disp': True})
    rest2 = minimize(TotalMeanFunctionY,x_ty,[rest.x[0],xr,aDataSet,nevent,skip,cut],method='Nelder-Mead',options={'xtol': 1e-5,'disp': True})
#     rest = minimize(TotalMeanFunctionX,x_tx,[xr,aDataSet,nevent,skip,cut],method='Nelder-Mead',options={'xtol': 1e-3,'disp': True})
#     rest2 = minimize(TotalMeanFunctionY,x_ty,[rest.x[0],xr,aDataSet,nevent,skip,cut],method='Nelder-Mead',options={'xtol': 1e-3,'disp': True})

    f = open(filename,'a')
    f.write("Rotation : %f %f %f [deg] Trans : %f %f  [mm] \n"%(xr[0],xr[1],xr[2],rest.x[0],rest.x[1]))
    f.close()

    return xr,[rest.x[0],rest2.x[0],0]


#
#optimize the sigma for the eta correction
#i.e. find the sigma for charge sharing which gives the best detector resolution (minimum value for the sigma of the residuals distributions)
#param 1: a data set (class EudetData)
#param 2: number of events we are running on
#param 3: number of skiped events
#
def FindSigmaMin(dataSet,nevent,skip=1) :
    
    #grX = TGraph(100)
    #grY = TGraph(100) 
    
    bestsigma=1000
    bestres=1000    
    for sigmaint in range(100) :
    	sigma=sigmaint*1e-4
	
	
	
    	resX=TotalSigmaFunctionX(sigma,sigma,dataSet,skip)
    	resY=TotalSigmaFunctionY(sigma,sigma,dataSet,skip)
	res=resX/2.0+resY/2.0
	#grX.SetPoint(sigmaint,sigma,resX)
	#grY.SetPoint(sigmaint,sigma,resY)
	
	if (res < bestres) : 
		bestres=res
		bestsigma=sigma
    #grX.Draw("AL")
    #grY.Draw("same")
    print "Best Sigma found at : %f um for resolution : %f um"%(bestsigma*1000,bestres*1000)
    #blah=raw_input()
     
    return bestres,bestres
     
     
	    
#    xsigmachargeX = np.array([0.005])
#    xsigmachargeY = np.array([0.005])
#    ressigmachargeX = minimize(TotalSigmaFunctionX,xsigmachargeX,[xsigmachargeY,dataSet,skip],method='BFGS',options={'gtol': 1e-5,'disp': True ,'eps': 0.001})
#    ressigmachargeY = minimize(TotalSigmaFunctionY,xsigmachargeY,[ressigmachargeX.x,dataSet,skip],method='BFGS',options={'gtol': 1e-5,'disp': True ,'eps': 0.001})

    return ressigmachargeX.x ,ressigmachargeY.x

def ApplyAlignment(dataSet,translations,rotations,dut=6,filename="Alignment.txt") :

    print "Applying Alignment with  Rotation : %0.10f %0.10f %0.10f [deg] Trans : %0.10f %0.10f  [mm]"%(rotations[0],rotations[1],rotations[2],translations[0],translations[1])
    RotMat = RotationMatrix(rotations)
    f = open(filename,'w')
    f.write("Rotation : %f %f %f [deg] Trans : %f %f  [mm] \n"%(rotations[0],rotations[1],rotations[2],translations[0],translations[1]))
    f.close()
    for Tracks in dataSet.AllTracks :
        for index,track in enumerate(Tracks) :
            tmp=np.dot(RotMat,[Tracks[index].trackX[track.iden.index(dut)],Tracks[index].trackY[track.iden.index(dut)],0])
            Tracks[index].trackX[track.iden.index(dut)] = tmp[0] + translations[0]
            Tracks[index].trackY[track.iden.index(dut)] = tmp[1] + translations[1]
#             track.trackZ[track.iden.index(dut)] = tmp[2] + translations[2]

def ApplyAlignment_at_event(i,dataSet,translations,rotations,dut=6) :

    #print "Applying Alignment with  Rotation : %0.10f %0.10f %0.10f [deg] Trans : %0.10f %0.10f  [mm]"%(rotations[0],rotations[1],rotations[2],translations[0],translations[1])
    RotMat = RotationMatrix(rotations)
    #f = open(filename,'w')
    #f.write("Rotation : %f %f %f [deg] Trans : %f %f  [mm] \n"%(rotations[0],rotations[1],rotations[2],translations[0],translations[1]))
    #f.close()
    for Tracks in dataSet.AllTracks[i:(i+1)] :
        for index,track in enumerate(Tracks) :
            tmp=np.dot(RotMat,[Tracks[index].trackX[track.iden.index(dut)],Tracks[index].trackY[track.iden.index(dut)],0])
            Tracks[index].trackX[track.iden.index(dut)] = tmp[0] + translations[0]
            Tracks[index].trackY[track.iden.index(dut)] = tmp[1] + translations[1]
#             track.trackZ[track.iden.index(dut)] = tmp[2] + translations[2]

#
#apply the eta correction i.e. compute the final residuals with the optimum charge sharing sigma
#param 1: a data set (class EudetData)
#param 2: value of the optimum charge sharing sigma coming from the X residuals study
#param 3: value of the optimum charge sharing sigma coming from the Y residuals study
#param 4: position of the device under test in the list of planes
#param 5: name of the file in which we store the results of the sigma optimisation
#
def ApplyEtaCorrection(dataSet,ressigmachargeX,ressigmachargeY,dut=6,filename="EtaCorrection.txt") :

    print "Applying Eta Correction with  chargeSigmaX : %f [mm] chargeSigmaY : %f [mm]"%(float(ressigmachargeX),float(ressigmachargeY))

#     f2 = open(filename,'w')
#     f2.write("ChargeSigmaX : %f [mm] ChargeSigmaY : %f [mm] \n"%(float(ressigmachargeX),float(ressigmachargeY))
#     f2.close()
    #ressigmachargeMean = (ressigmachargeX + ressigmachargeY)/2.
    for j,tracks in enumerate(dataSet.AllTracks) :
        for track in tracks :
            if track.cluster!=-11 and len(dataSet.AllClusters[j])!=0 :
                dataSet.AllClusters[j][track.cluster].GetEtaCorrectedQWeightedCentroid(ressigmachargeX,ressigmachargeY)
        dataSet.ComputeResiduals(j)


def EdgeEfficiency(aDataSet,dut) :

    TotalTrack = TH1D("TotalTrack","Track distribution in edge",50,0,aDataSet.edge)
    MatchedTrack = TH1D("MatchedTrack","Matched track distribution in edge",50,0,aDataSet.edge)
    TOT_vs_edge = TH2D("TOT_vs_edge","TOT vs track position in edge",50,0,aDataSet.edge,200,0,1000)

    edge_tracks = []
    edge_matched = []
    edge_efficiencies = []
    for i in range(4):
        TotalTrack_edge = TH1D("TotalTrack_edge_%i"%i,"Track distribution in edge %i" %i,50,0,aDataSet.edge)
        MatchedTrack_edge = TH1D("MatchedTrack_edge_%i"%i,"Matched track distribution in edge %i" %i,50,0,aDataSet.edge)
        edge_tracks.append(TotalTrack_edge)
        edge_matched.append(MatchedTrack_edge)


    for j,tracks in enumerate(aDataSet.AllTracks) :
        for track in tracks :
            if(aDataSet.IsInEdges(track,dut)) :
                if(fabs(track.trackX[track.iden.index(dut)])>halfChip_X and fabs(track.trackY[track.iden.index(dut)])<halfChip_Y):
                    TotalTrack.Fill(fabs(track.trackX[track.iden.index(dut)])-halfChip_X)
                    if(track.trackX[track.iden.index(dut)]>0) :
                        edge_tracks[0].Fill(fabs(track.trackX[track.iden.index(dut)])-halfChip_X)
                    else :
                        edge_tracks[2].Fill(fabs(track.trackX[track.iden.index(dut)])-halfChip_X)

                    if track.cluster!=-11 and len(aDataSet.AllClusters[j])!=0 :
                        MatchedTrack.Fill(fabs(track.trackX[track.iden.index(dut)])-halfChip_X)
                        TOT_vs_edge.Fill(fabs(track.trackX[track.iden.index(dut)])-halfChip_X,aDataSet.AllClusters[j][track.cluster].totalTOT)
                        if(track.trackX[track.iden.index(dut)]>0) :
                            edge_matched[0].Fill(fabs(track.trackX[track.iden.index(dut)])-halfChip_X)
                        else :
                            edge_matched[2].Fill(fabs(track.trackX[track.iden.index(dut)])-halfChip_X)

                if(fabs(track.trackX[track.iden.index(dut)])<halfChip_X and fabs(track.trackY[track.iden.index(dut)])>halfChip_Y):
                    TotalTrack.Fill(fabs(track.trackY[track.iden.index(dut)])-halfChip_Y)
                    if(track.trackY[track.iden.index(dut)]>0) :
                        edge_tracks[1].Fill(fabs(track.trackY[track.iden.index(dut)])-halfChip_Y)
                    else :
                        edge_tracks[3].Fill(fabs(track.trackY[track.iden.index(dut)])-halfChip_Y)

                    if track.cluster!=-11 and len(aDataSet.AllClusters[j])!=0 :
                        MatchedTrack.Fill(fabs(track.trackY[track.iden.index(dut)])-halfChip_Y)
                        TOT_vs_edge.Fill(fabs(track.trackY[track.iden.index(dut)])-halfChip_Y,aDataSet.AllClusters[j][track.cluster].totalTOT)
                        if(track.trackY[track.iden.index(dut)]>0) :
                            edge_matched[1].Fill(fabs(track.trackY[track.iden.index(dut)])-halfChip_Y)
                        else :
                            edge_matched[3].Fill(fabs(track.trackY[track.iden.index(dut)])-halfChip_Y)


    for i in range(4):
        edge_tracks[i].SetLineColor(i+1)
        edge_matched[i].SetLineColor(i+1)
        h = edge_matched[i].Clone("Efficiency")
        h.Divide(edge_matched[i],edge_tracks[i],1.,1.,"B")
        h.SetTitle("Efficiency in edge %i"%i)
        edge_efficiencies.append(h)

    Efficiency = MatchedTrack.Clone("Efficiency")
    Efficiency.Divide(MatchedTrack,TotalTrack,1.,1.,"B")
    Efficiency.SetTitle("Efficiency")

    return TotalTrack, MatchedTrack, Efficiency, edge_tracks, edge_matched, edge_efficiencies, TOT_vs_edge



def ComputeEfficiency(aDataSet,n_matched,n_matched_edge,edge,PlotPath):
    n_tracks_in_total = ComputeDetectorAcceptance(aDataSet,6,edge)
    n_tracks_in_main = ComputeDetectorAcceptance(aDataSet,6,0)
    n_tracks_in_edge = n_tracks_in_total - n_tracks_in_main

    print "n_tracks_in_total", n_tracks_in_total
    print "n_tracks_in_main", n_tracks_in_main
    print "n_tracks_in_edge", n_tracks_in_edge

    efficiency_in_main = 0.
    efficiency_in_edge = 0.
    efficiency_in_total = 0.

    if n_tracks_in_main != 0 :
        efficiency_in_main = float(n_matched)/n_tracks_in_main
    else:
        efficiency_in_main = 0.

    if n_tracks_in_edge != 0 :
        efficiency_in_edge = float(n_matched_edge)/n_tracks_in_edge
    else:
        efficiency_in_edge = 0.

    if n_tracks_in_total != 0. :
        efficiency_in_total = float(n_matched + n_matched_edge)/n_tracks_in_total
    else:
        efficiency_in_total = 0.
       
    print "Efficiency in main : %f %%"%(efficiency_in_main*100)
    print "Efficiency in edge : %f %%"%(efficiency_in_edge*100)
    print "Efficiency in total : %f %%"%(efficiency_in_total*100)

    f = open("%s/Efficiency.txt"%PlotPath,'w')
    f.write("Efficiency in main : %f %%"%(efficiency_in_main*100))
    f.write("Efficiency in edge : %f %%"%(efficiency_in_edge*100))
    f.write("Efficiency in total : %f %%"%(efficiency_in_total*100))
    f.close()


def ApplyGlobalEnergyCalibration(aDataSet,nevents,a,b,c,t):

    for i,clusters in enumerate(aDataSet.AllClusters[0:nevents]) :
        for j,cluster in enumerate(clusters) :
            for k,tot in enumerate(cluster.tot):
                
                cluster.tot[k] = ( t*a + tot - b + sqrt( ( b + t*a - tot )**2 + 4*a*c ) ) / (2*a) # energy in keV

def ReadCalibFile(calibFile):

    a = [[ 0. for x in xrange(npix_X)] for x in xrange(npix_Y)]
    b = [[ 0. for x in xrange(npix_X)] for x in xrange(npix_Y)]
    c = [[ 0. for x in xrange(npix_X)] for x in xrange(npix_Y)]
    t = [[ 0. for x in xrange(npix_X)] for x in xrange(npix_Y)]
    
    rootfile = TFile(calibFile)
    tree = rootfile.Get("fitPara")

    for entry in tree:
        col = tree.pixx
        row = tree.pixy
        a[col][row] = tree.a
        b[col][row] = tree.b
        c[col][row] = tree.c
        t[col][row] = tree.d

    return a,b,c,t

def ApplyPixelEnergyCalibration(aDataSet,nevents,calibFile):

    a,b,c,t = ReadCalibFile(calibFile)

    for i,clusters in enumerate(aDataSet.AllClusters[0:nevents]) :
        for j,cluster in enumerate(clusters) :
            for k in xrange(len(cluster.tot)):
                
                col = cluster.col[k]
                row = cluster.row[k]
                tot = cluster.tot[k]

                cluster.tot[k] = ( t[col][row]*a[col][row] + tot - b[col][row] + sqrt( ( b[col][row] + t[col][row]*a[col][row] - tot )**2 + 4*a[col][row]*c[col][row] ) ) / (2*a[col][row]) # energy in keV


def CountPixelMapRepeats(dataSet,n_proc):
    # counts the number of times a pixel map is repeated
    # returns a histogram of result

    histo_mapreps = TH1F("histo_mapreps","",20,0.5,20.5)
    histo_mapreps.GetXaxis().SetTitle("Number of repeated pixel maps")
    prev_pixel_xhits = []
    nrepeats = 1

    for i in range(0,n_proc) :
        dataSet.getEvent(i)

        npixels_hit = len(dataSet.p_col)
        pixel_x_hits = []
        for k in xrange(npixels_hit):
            pixel_x_hits.append(dataSet.p_col[k])

        if (pixel_x_hits == prev_pixel_xhits):
            # same pixel map as before
            nrepeats = nrepeats+1
        else:
            # new pixel map
            if i != 0:
                histo_mapreps.Fill(nrepeats)
            nrepeats = 1

        prev_pixel_xhits = pixel_x_hits

    return histo_mapreps


###############################################################################################################################
#
#                        landau * gauss fit tools
#
###############################################################################################################################



# def langaufun(x,par) :
#
#     #Fit parameters:
#     #par[0]=Width (scale) parameter of Landau density
#     #par[1]=Most Probable (MP, location) parameter of Landau density
#     #par[2]=Total area (integral -inf to inf, normalization constant)
#     #par[3]=Width (sigma) of convoluted Gaussian function
#     #
#     #In the Landau distribution (represented by the CERNLIB approximation),
#     #the maximum is located at x=-0.22278298 with the location parameter=0.
#     #This shift is corrected within this function, so that the actual
#     #maximum is identical to the MP parameter.
#
#     # Numeric constants
#     invsq2pi = 0.3989422804014   # (2 pi)^(-1/2)
#     mpshift  = -0.22278298       # Landau maximum location
#
#     # Control constants
#     np = 100.0      # number of convolution steps
#     sc =   5.0      # convolution extends to +-sc Gaussian sigmas
#
#     # Variables
#     sum = 0.0
#
#
#     # MP shift correction
#     mpc = par[1] - mpshift * par[0]
#
#     # Range of convolution integral
#     xlow = x - sc * par[3]
#     xupp = x + sc * par[3]
# #     xlow = x[0] - sc * par[3]
# #     xupp = x[0] + sc * par[3]
#
#     step = (xupp-xlow) / np
#
#     # Convolution integral of Landau and Gaussian by sum
#     #for(i=1.0; i<=np/2; i++) {
#     for i in range(1,np/2 + 1) :
#         xx = xlow + (i-.5) * step
#         fland = TMath.Landau(xx,mpc,par[0]) / par[0]
#         sum = sum + fland * TMath.Gaus(x[0],xx,par[3])
#
#         xx = xupp - (i-.5) * step;
#         fland = TMath.Landau(xx,mpc,par[0]) / par[0]
#         sum = sum + fland * TMath.Gaus(x[0],xx,par[3])
#
#
#     return (par[2] * step * sum * invsq2pi / par[3])
#
#
#
#
# def langaufit(his, fitrange,startvalues, parlimitslo, parlimitshi, fitparams, fiterrors, ChiSqr, NDF) :
#
#     # Once again, here are the Landau * Gaussian parameters:
#     #   par[0]=Width (scale) parameter of Landau density
#     #   par[1]=Most Probable (MP, location) parameter of Landau density
#     #   par[2]=Total area (integral -inf to inf, normalization constant)
#     #   par[3]=Width (sigma) of convoluted Gaussian function
#     #
#     # Variables for langaufit call:
#     #   his             histogram to fit
#     #   fitrange[2]     lo and hi boundaries of fit range
#     #   startvalues[4]  reasonable start values for the fit
#     #   parlimitslo[4]  lower parameter limits
#     #   parlimitshi[4]  upper parameter limits
#     #   fitparams[4]    returns the final fit parameters
#     #   fiterrors[4]    returns the final fit errors
#     #   ChiSqr          returns the chi square
#     #   NDF             returns ndf
#
#
# #     FunName = "Fitfcn_%s"%his.GetName()
#     # sprintf(FunName,"Fitfcn_%s",his->GetName());
#
# #     ffitold = gROOT.GetListOfFunctions().FindObject(FunName)
# #     if (ffitold) :
# #         ffitold.Delete
#
#     ffit = TF1("ffit",langaufun,fitrange[0],fitrange[1],4)
#     ffit.SetParameters(startvalues[0],startvalues[1],startvalues[2],startvalues[3])
#     ffit.SetParNames("Width","MP","Area","GSigma")
#
#     for i in range (0,4) :
#         ffit.SetParLimits(i, parlimitslo[i], parlimitshi[i])
#
#
#     his.Fit("ffit","RB0")   # fit within specified range, use ParLimits, do not plot
#
#     fitparams=ffit.GetParameters()    # obtain fit parameters
#     for i in range (0,4) :
#         fiterrors.append(ffit.GetParError(i))     # obtain fit parameter errors
# #         fiterrors[i] = ffit.GetParError(i)     # obtain fit parameter errors
#
#     ChiSqr = ffit.GetChisquare()   # obtain chi^2
#     NDF = ffit.GetNDF()           # obtain ndf
# #     ChiSqr[0] = ffit.GetChisquare()   # obtain chi^2
# #     NDF[0] = ffit.GetNDF()           # obtain ndf
#
#     return (ffit)              # return fit function
#
#
#
#
# def langaupro(params, maxx, FWHM) :
#
#     # Seaches for the location (x value) at the maximum of the
#     # Landau-Gaussian convolute and its full width at half-maximum.
#     #
#     # The search is probably not very efficient, but it's a first try.
#
#     i = 0
#     MAXCALLS = 10000
#
# #for test, to comment after !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#     params.append(1000)
#     params.append(1000)
#
#     # Search for maximum
#
#     p = params[1] - 0.1 * params[0]
#     step = 0.05 * params[0]
#     lold = -2.0
#     l    = -1.0
#
#
#     while ( (l != lold) and (i < MAXCALLS) ) :
#         i = i + 1
#
#         lold = l
#         x = p + step
#         l = langaufun(x,params)
#
#         if (l < lold) :
#             step = -step/10
#
#         p = p + step
#
#
#     if (i == MAXCALLS) :
#         return (-1)
#
#     maxx = x
#
#     fy = l/2
#
#
#     # Search for right x location of fy
#
#     p = maxx + params[0]
#     step = params[0]
#     lold = -2.0
#     l    = -1e300
#     i    = 0
#
#
#     while ( (l != lold) and (i < MAXCALLS) ) :
#         i = i + 1
#
#         lold = l
#         x = p + step
#         l = Abs(langaufun(x,params) - fy)
#
#         if (l > lold) :
#             step = -step/10
#
#         p = p + step
#
#
#     if (i == MAXCALLS) :
#         return (-2)
#
#     fxr = x
#
#
#     # Search for left x location of fy
#
#     p = maxx - 0.5 * params[0]
#     step = -params[0]
#     lold = -2.0
#     l    = -1e300
#     i    = 0
#
#     while ( (l != lold) and (i < MAXCALLS) ) :
#         i = i + 1
#
#         lold = l
#         x = p + step
#         l = Abs(langaufun(x,params) - fy)
#
#         if (l > lold) :
#             step = -step/10
#
#         p = p + step
#
#
#     if (i == MAXCALLS) :
#         return (-3)
#
#
#     fxl = x
#
#     FWHM = fxr - fxl
#     return (0)
