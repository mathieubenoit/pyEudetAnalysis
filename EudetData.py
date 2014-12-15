from ROOT import *
import ROOT
from math import fsum,fabs
from array import array
from scipy.cluster.hierarchy import fclusterdata

import pyximport; pyximport.install(pyimport=True)
from ROOT import TFile,TH1D,TH2D,TF1
from Cluster import *
from Track import  *
#gROOT.LoadMacro("Track.C+")
#from ROOT import Track
from Constant import *
from ToolBox import *
from PersistentList import *
from itertools import product
###############################################################################################################################
#
#        A container for TBTrack Data: contains all the informations about tracks and cluster from data
#
###############################################################################################################################

class EudetData:
    """A container for TBTrack Data """

    RunNumber = 0

    tbtrack_file = 0
    pixelTree = ROOT.TTree()
    TrackTree = ROOT.TTree()

    Chi2 = TH1D()
    Chi2ndof = TH1D()
    Chi2_Cut=10000000
    EnergyCut = 0.
    scale =1.

    p_nEntries = 0
    t_nEntries = 0
    entry = 0

    sigmaX=0.005
    sigmaY=0.005

    #Track Data Holders
    t_nTrackParams = 0
    t_euEv= 0
    t_posX = 0
    t_posY= 0
    t_dxdz= 0
    t_dydz= 0
    t_iden= 0
    t_trackNum= 0
    t_chi2= 0
    t_ndof = 0
#     t_posX = []
#     t_posY = []
#     t_dxdz = []
#     t_dydz = []
#     t_iden = []
#     t_trackNum = []
#     t_chi2 = []
#     t_ndof = []

    #Pixel Data holders
    p_nHits= 0
    p_col= 0
    p_row= 0
    p_tot= 0
    p_lv1= 0
    p_chip= 0
    p_iden= 0
    p_euEv = 0


    edge = 0
#     p_col = []
#     p_row = []
#     p_tot = []
#     p_lv1 = []
#     p_chip = []
#     p_iden = []
#     p_euEv = []

    #hotpixel firing matrix
    #hit_map = [[0]*npix_Y]*npix_X
    #frequency_map = [[0.0]*npix_Y]*npix_X

    hit_map =       [[0 for x in xrange(npix_Y)] for y in xrange(npix_X)]
    frequency_map = [[0 for x in xrange(npix_Y)] for y in xrange(npix_X)]
        
    hotpixels = []

    mode = ""

    def __init__(self,filename,ECut,edge=0,scale=1.0, Run = 0,mode="tbtrack"):

        self.RunNumber = Run
        self.edge=edge
        #self.AllClusters = PersistentList("cluster_%i"%self.RunNumber,250)
        #self.AllTracks = PersistentList("Track_%i"%self.RunNumber,250)
        self.AllClusters = []
        self.AllTracks = []

        self.mode=mode

        self.scale=scale


        self.tbtrack_file = TFile(filename)
        print "Opening %s"%filename
        if(self.mode=="tbtrack"):
            print "Reading in tbtrack mode"
            self.TrackTree = self.tbtrack_file.Get("eutracks")
            self.pixelTree = self.tbtrack_file.Get("zspix")
            self.p_nEntries = self.pixelTree.GetEntries()
            self.t_nEntries = self.TrackTree.GetEntries()
        elif(self.mode=="pyEudetNTuple") :
            print "Reading in pyEudetNTuple mode"
            self.TrackTree = self.tbtrack_file.Get("tracks")
            self.pixelTree = self.tbtrack_file.Get("clusters")
            self.TrackTree.Print()
            self.pixelTree.Print()
            self.p_nEntries = self.pixelTree.GetEntries()
            self.t_nEntries = self.TrackTree.GetEntries()
        else :
            print "Wrong mode. Exiting ..."
            exit()

        self.EnergyCut=ECut
        self.frequency_map = [[0.0]*npix_Y]*npix_X
        for i in range(len(self.hit_map)):
            for j in range(len(self.hit_map[0])):
                self.hit_map[i][j]=0

    def GetChi2Cut(self,reduction_factor=0.95,applyChiCut=True) :

        self.TrackTree.Draw("chi2 >> chi2plot(1000,0,1000)","","goff")
        self.Chi2 = gROOT.FindObject("chi2plot")

        self.TrackTree.Draw("chi2/ndof >> chi2ndofplot(1000,0,100)","","goff")
        self.Chi2ndof = gROOT.FindObject("chi2ndofplot")


        totIntegral = reduction_factor*self.Chi2.Integral()

        aBin = 0
        while(self.Chi2.Integral(0,aBin)<totIntegral) :
            aBin+=1


        self.Chi2_Cut = 1000000000
        if(applyChiCut==True):
            print "Cutting at Chi2 = %f"%(aBin*self.Chi2.GetBinWidth(0))
            self.Chi2_Cut = aBin*self.Chi2.GetBinWidth(0)
        return self.Chi2,self.Chi2ndof





    def FindHotPixel(self,threshold,Nevents=-1):

        # will calculate the frequency with which each pixel fires
        # threshold (0 -> 1) defines hot pixel cut

        n_max = 0
        prev_pixel_xhits = []
        unique_events = 0

        histo_nhits = TH1D("nhit","N Pixel Fires",40,0,39)
        histo_hitpixel = TH2D("hit","Hit Pixel Map",npix_Y,0,npix_Y-1,npix_X,0,npix_X-1)
        histo_frequency = TH1D("freq","Pixel Firing Frequency",10000,0,1)
        histo_hotpixel = TH2D("hot","Hot Pixel Map",npix_Y,0,npix_Y-1,npix_X,0,npix_X-1)

        if Nevents>self.p_nEntries or Nevents==-1:
            n_max = self.p_nEntries
        elif Nevents < 10000 and self.p_nEntries >= 10000:
            print "FindHotPixel over-riding requested nevents"
            print "FindHotPixel must be run on atleast 10000 events (for a threshold of 0.01) to be accurate"
            print "FindHotPixel will use 10000 events"
            n_max = 10000
        elif Nevents < 10000 and self.p_nEntries < 10000:
            print "FindHotPixel over-riding requested nevents"
            print "FindHotPixel must be run on atleast 10000 events (for a threshold of 0.01) to be accurate"
            print "FindHotPixel will use as many events as exist in this run"
            n_max = self.p_nEntries
        else :
            n_max = Nevents

        # loop through events to find unique events
        # for each fired pixel in each event, increment hit map
        for i in range(n_max) :
            self.getEvent(i)
            if i%10000 == 0 :
                print " [Hot Pixel Finder] Parsing event %i" %i

            # is this a new frame, or the next event in the same frame?
            npixels_hit = len(self.p_col)
            pixel_x_hits = []
            for k in xrange(npixels_hit):
                pixel_x_hits.append(self.p_col[k])

            if (pixel_x_hits == prev_pixel_xhits):
                # another track in the same event
                continue
            else:
                # this is a new event
                unique_events = unique_events + 1

            prev_pixel_xhits = pixel_x_hits

            for j in range(len(self.p_row)) :
                self.hit_map[self.p_row[j]][self.p_col[j]] += 1
                histo_hitpixel.Fill(self.p_col[j],self.p_row[j])

        # loop through hitmap
        # fill freq map with hits / nevents
        print "Ran over", n_max, "events, found", unique_events, "unique pixel maps"
        for i in range(npix_X):
            for j in range(npix_Y) :
		self.frequency_map[i][j]=self.hit_map[i][j]*(1.0/float(unique_events))
                histo_nhits.Fill(self.hit_map[i][j])
                histo_frequency.Fill(self.frequency_map[i][j])

                # if freq > threshold, make a hotpixel
                if(self.frequency_map[i][j] > threshold):
                    histo_hotpixel.Fill(i,j,self.frequency_map[i][j]) 
                    self.hotpixels.append([i,j])

        print "##### Hot Pixel Report #####"
        print " %i Hot pixel found at  : "%(len(self.hotpixels))
        print self.hotpixels
        print "############################"
        return histo_nhits,histo_hitpixel,histo_hotpixel,histo_frequency


    def getEvent(self,i):

        self.entry = self.TrackTree.GetEntry(i)

        self.t_nTrackParams = self.TrackTree.nTrackParams
        self.t_euEv= self.TrackTree.euEvt
        self.t_posX = self.TrackTree.xPos
        self.t_posY= self.TrackTree.yPos 
        self.t_dxdz= self.TrackTree.dxdz
        self.t_dydz= self.TrackTree.dydz
        self.t_iden= self.TrackTree.iden
        self.t_trackNum= self.TrackTree.trackNum
        self.t_chi2= self.TrackTree.chi2
        self.t_ndof = self.TrackTree.ndof

        for index,Xval in enumerate(self.t_posX) :
            self.t_posX[index]+=pitchX/2.
            self.t_posY[index]+=pitchY/2.

        self.entry = self.pixelTree.GetEntry(i)

        self.p_nHits= self.pixelTree.nPixHits
        self.p_col= self.pixelTree.col
        self.p_row= self.pixelTree.row
        self.p_tot= self.pixelTree.tot
        self.p_lv1= self.pixelTree.lv1
        #self.p_chip= self.pixelTree.chip
        self.p_iden= self.pixelTree.iden
        self.p_euEv = self.pixelTree.euEvt

 #       for index,totvalue in enumerate(self.p_tot) :
 #           self.p_tot[index]=float(totvalue)/self.scale
 
 
    def PlotFrame(self,i,c,n_pix_min=0) : 
        
	plot = TH2D("frame %i" %i, "frame %i" %i, npix_X, 0, npix_X, npix_Y, 0, npix_Y)
	self.getEvent(i)
	
	if(len(self.p_col) > n_pix_min):
	    for j in xrange(len(self.p_col)) : 
                plot.Fill(self.p_col[j],self.p_row[j],self.p_tot[j])
	
            plot.Draw("colz")
            c.Update()
	    print "press enter for next frame, ctrl-D to exit"
	    a=raw_input()
	else : 
	    print "Skipping event %i, does not have more than minimum number of hits (%i)" %(i,n_pix_min)

    def WriteReconstructedData(self,filename,dut=20) :
        outfile = TFile(filename,'recreate')
        self.DumpTrackTree(outfile,dut)
        self.DumpClusterTree(outfile,dut)
        outfile.Close()

    def ReadReconstructedData(self,NEvents=-1) :
        print "Reading data ..."
        self.ReadTrackTree()
        self.ReadClusterTree()

    def ReadTrackTree(self,NEvents=-1):
        event = 0
        nentries = self.TrackTree.GetEntriesFast()
        events = []
        for i in xrange(nentries) :
            self.entry = self.TrackTree.GetEntry(i)
            events.append(self.TrackTree.event)
        for i in range(max(events)+1):
            self.AllTracks.append([])
        print "%i events in track file"%max(events)

        for i in xrange(nentries) :
            self.entry = self.TrackTree.GetEntry(i)
            track = Track()
            for i in range(self.TrackTree.size) :
                track.trackX.append(self.TrackTree.trackX[i])
                track.trackY.append(self.TrackTree.trackY[i])
                track.chi2.append(self.TrackTree.chi2[i])
                track.ndof.append(self.TrackTree.ndof[i])
                track.trackNum.append(self.TrackTree.trackNum[i])
                track.dxdz.append(self.TrackTree.dxdz[i])
                track.dydz.append(self.TrackTree.dydz[i])
                track.iden.append(self.TrackTree.iden[i])
            event = self.TrackTree.event
            track.cluster = self.TrackTree.cluster
            self.AllTracks[event].append(track)

    def ReadClusterTree(self):
        event = 0
        events = []
        nentries = self.pixelTree.GetEntriesFast()
        for i in xrange(nentries) :
            self.entry = self.pixelTree.GetEntry(i)
            events.append(self.pixelTree.event)
        for i in range(max(events)+1):
            self.AllClusters.append([])
        print "%i events in cluster file"%max(events)
        for i in xrange(nentries) :
            self.entry = self.pixelTree.GetEntry(i)
            cluster = Cluster()
            for i in range(self.pixelTree.size) :
                cluster.col.append(self.pixelTree.col[i])
                cluster.row.append(self.pixelTree.row[i])
                cluster.tot.append(self.pixelTree.tot[i])
            event = self.pixelTree.event
            cluster.sizeX = self.pixelTree.sizeX
            cluster.sizeY =self.pixelTree.sizeY
            cluster.size  = self.pixelTree.size
            cluster.totalTOT =self.pixelTree.totalTOT
            cluster.aspectRatio = self.pixelTree.aspectRatio
            cluster.relX = self.pixelTree.relX
            cluster.relY = self.pixelTree.relY
            cluster.absX = self.pixelTree.absX
            cluster.absY =      self.pixelTree.absY
            cluster.resX = self.pixelTree.resX
            cluster.resY = self.pixelTree.resY
            cluster.id = self.pixelTree.id
            cluster.trackNum = self.pixelTree.trackNum
            self.AllClusters[event].append(cluster)



    def DumpTrackTree(self,outfile,dut=20):

        outfile.cd()
        trackTree = TTree('tracks','TestBeam track tree')
        nplanes = 7
        trackX=array( 'f', nplanes*[ 0. ] )
        trackY=array( 'f', nplanes*[ 0. ] )
        iden=array( 'i', nplanes*[ 0 ] )
        chi2 =array( 'f', nplanes*[ 0. ] )
        event=array( 'i', [ 0 ] )
        ndof =array( 'f', nplanes*[ 0. ] )
        trackNum =array( 'i', nplanes*[ 0 ] )
        dxdz =array( 'f', nplanes*[ 0. ] )
        dydz =array( 'f', nplanes*[ 0. ] )
        cluster= array( 'i', [ 0 ] )
        clusterX=array( 'f', [ 0. ] )
        clusterY=array( 'f', [ 0. ] )
        size= array( 'i', [ 0 ] )

        trackTree.Branch( 'size', size, 'size/I' )
        trackTree.Branch( 'trackX', trackX, 'trackX[size]/F' )
        trackTree.Branch( 'trackY', trackY, 'trackY[size]/F' )
        trackTree.Branch( 'iden', iden, 'iden[size]/I' )
        trackTree.Branch( 'chi2', chi2, 'chi2[size]/F' )
        trackTree.Branch( 'event', event, 'event/I' )
        trackTree.Branch( 'ndof', ndof, 'ndof[size]/F' )
        trackTree.Branch( 'trackNum', trackNum, 'trackNum[size]/I' )
        trackTree.Branch( 'dxdz', dxdz, 'dxdz[size]/F' )
        trackTree.Branch( 'dydz', dydz, 'dydz[size]/F' )
        trackTree.Branch( 'cluster', cluster, 'cluster/I' )
        trackTree.Branch( 'clusterX', clusterX, 'clusterX/F' )
        trackTree.Branch( 'clusterY', clusterY, 'clusterY/F' )


        for j,tracks in enumerate(self.AllTracks) :
            for track in tracks :
                size[0]=nplanes
                for index,t in enumerate(track.trackX):
                    trackX[index]=t
                for index,t in enumerate(track.trackY):
                    trackY[index]=t
                for index,t in enumerate(track.chi2):
                    chi2[index]=t
                for index,t in enumerate(track.ndof):
                    ndof[index]=t
                for index,t in enumerate(track.trackNum):
                    trackNum[index]=t
                for index,t in enumerate(track.dxdz):
                    dxdz[index]=t
                for index,t in enumerate(track.dydz):
                    dydz[index]=t
                for index,t in enumerate(track.iden):
                    iden[index]=t
                event[0]=j
                cluster[0]=track.cluster
                if track.cluster!=-11 and len(self.AllClusters[j])!=0 :
                    clusterX[0]=self.AllClusters[j][track.cluster].absX
                    clusterY[0]=self.AllClusters[j][track.cluster].absY
                else:
                    clusterX[0]=-1000
                    clusterY[0]=-1000
                trackTree.Fill()
        outfile.Write()


    def DumpClusterTree(self,outfile,dut=20):

        outfile.cd()

        clusterTree = TTree('clusters','Timepix cluster tree')

        maxn = 500
        col = array( 'i', maxn*[ 0 ] )
        row = array( 'i', maxn*[ 0 ] )
        tot = array( 'f', maxn*[ 0. ] )
        event = array( 'i', [ 0 ] )
        sizeX = array( 'i', [ 0 ] )
        sizeY = array( 'i', [ 0 ] )
        size  = array( 'i', [ 0 ] )
        totalTOT =array( 'f', [ 0. ] )
        aspectRatio = array( 'f', [ 0. ] )
        relX = array( 'f', [ 0. ] )
        relY = array( 'f', [ 0. ] )
        absX = array( 'f', [ 0. ] )
        absY = array( 'f', [ 0. ] )
        resX = array( 'f', [ 0. ] )
        resY = array( 'f', [ 0. ] )
        trackX = array( 'f', [ 0. ] )
        trackY = array( 'f', [ 0. ] )
        id = array( 'f', [ 0. ] )
        trackNum = array( 'i', [ 0 ] )

        clusterTree.Branch( 'event', event, 'event/I' )
        clusterTree.Branch( 'size', size, 'size/I' )
        clusterTree.Branch( 'sizeX', sizeX, 'sizeX/I' )
        clusterTree.Branch( 'sizeY', sizeY, 'sizeY/I' )
        clusterTree.Branch( 'totalTOT', totalTOT, 'totalTOT/F' )
        clusterTree.Branch( 'aspectRatio', aspectRatio, 'aspectRatio/F' )
        clusterTree.Branch( 'relX', relX, 'relX/F' )
        clusterTree.Branch( 'relY', relY, 'relY/F' )
        clusterTree.Branch( 'absX', absX, 'relX/F' )
        clusterTree.Branch( 'absY', absY, 'relY/F' )
        clusterTree.Branch( 'resX', resX, 'resX/F' )
        clusterTree.Branch( 'resY', resY, 'resY/F' )
        clusterTree.Branch( 'trackX', trackX, 'trackX/F' )
        clusterTree.Branch( 'trackY', trackY, 'trackY/F' )
        clusterTree.Branch( 'id', id, 'id/I' )
        clusterTree.Branch( 'trackNum', trackNum, 'trackNum/I' )
        clusterTree.Branch( 'col', col, 'col[size]/I' )
        clusterTree.Branch( 'row', row, 'row[size]/I' )
        clusterTree.Branch( 'tot', tot, 'tot[size]/F' )


        for j,clusters in enumerate(self.AllClusters) :
            for cluster in clusters :

                if(len(cluster.col)<maxn) :
		
	   	    for i in range(len(cluster.col)) :
                        col[i]=cluster.col[i]
                    for i in range(len(cluster.row)) :
                        row[i]=cluster.row[i]
                    for i in range(len(cluster.tot)) :
                        tot[i]=cluster.tot[i]
			
		else : 
	   	    for i in range(maxn) :
                        col[i]=cluster.col[i]
                    for i in range(maxn) :
                        row[i]=cluster.row[i]
                    for i in range(maxn) :
                        tot[i]=cluster.tot[i]
					

                sizeX[0]=cluster.sizeX
                sizeY[0]=cluster.sizeY
                size[0]=cluster.size
                totalTOT[0]=cluster.totalTOT
                aspectRatio[0]=cluster.aspectRatio
                relX[0]=cluster.relX
                relY[0]=cluster.relY
                resX[0]=cluster.resX
                resY[0]=cluster.resY
                absX[0]=cluster.absX
                absY[0]=cluster.absY
                id[0]=cluster.id
                event[0]=j
                trackNum[0]=cluster.tracknum
                if (cluster.tracknum!=-1):
                    try :
                        trackX[0]=self.AllTracks[j][cluster.tracknum].trackX[self.AllTracks[j][cluster.tracknum].iden.index(dut)]
                        trackY[0]=self.AllTracks[j][cluster.tracknum].trackY[self.AllTracks[j][cluster.tracknum].iden.index(dut)]
                    except :
                        trackX[0]=0
                        trackY[0]=0
                clusterTree.Fill()
        outfile.Write()

    def IsInEdges(self,track,dut=20):
        is_in = False
        
	if(fabs(track.trackX[track.iden.index(dut)])<=(halfChip_X+self.edge) and fabs(track.trackY[track.iden.index(dut)])<=(halfChip_Y+self.edge)):
            is_in = True
            if(fabs(track.trackX[track.iden.index(dut)])<=(halfChip_X) and fabs(track.trackY[track.iden.index(dut)])<=(halfChip_Y)):
                is_in=False
        return is_in


    def IsInGoodRegion(self,track,dut=20) : 
    
        if (track.trackX[track.iden.index(dut)]>=((pitchX*npix_X)/2-4*pitchX) and track.trackX[track.iden.index(dut)]<=((pitchX*npix_X)/2)) and (track.trackY[track.iden.index(dut)]>=(-(pitchY*npix_Y)/2.) and track.trackY[track.iden.index(dut)]<=((pitchY*npix_Y)/2.)) :
 	    return true
	else : 
	    return false

    def ComputeResiduals(self,i,dut=20) :

        nmatch_in_main = 0.
        nmatch_in_edge = 0.

        for track in self.AllTracks[i] :
            if(i<len(self.AllClusters)) :
                for cluster in self.AllClusters[i] :

                    if cluster.id == track.cluster  :
                        cluster.GetResiduals(track.trackX[track.iden.index(dut)],track.trackY[track.iden.index(dut)])

                        if(self.IsInEdges(track)):
                            nmatch_in_edge += 1.
                        else :
                            nmatch_in_main += 1.

        return nmatch_in_main, nmatch_in_edge



    def PrintResiduals(self,i) :
        print "###################### Event : %d ######################"%i
        self.getEvent(i)
        for cluster in self.AllClusters[i] :
            cluster.GetResiduals(self.t_posX[3],self.t_posY[3])
            print "resX = %f resY = %f"%(cluster.resX,cluster.resY)
        print "#######################################################"

    def PrintEvent(self,i):
        self.getEvent(i)
        print "###################### Event : %d ######################"%i
        outstr =""
        print "posX = %f posY = %f"%(self.t_posX[3],self.t_posY[3])


#    for j in self.t_posX :
#            outstr+="%.3f "%j
#        print "posX = %s"%outstr
#        outstr =""
#        for j in self.t_posY :
#            outstr+="%.3f "%j
#        print "posY = %s"%outstr
        #=======================================================================
        # outstr =""
        # for j in self.t_dxdz :
        #    outstr+="%.6e "%j
        # print "dxdz = %s"%outstr
        # outstr =""
        # for j in self.t_dydz :
        #    outstr+="%.6e "%j
        # print "dydz = %s"%outstr
        #=======================================================================
        print "#######################################################"


    def GetTrack(self,i) :

        self.getEvent(i)
        posX_tmp = []
        posY_tmp = []
        dxdz_tmp = []
        dydz_tmp = []
        iden_tmp = []
        chi2_tmp = []
        ndof_tmp = []
        trackNum_tmp = []
        nTrackParams_tmp = 0


        tracks = []

        #--------------------- self.t_nTrackParams = self.TrackTree.nTrackParams
        #------------------------------------- self.t_euEv= self.TrackTree.euEvt
        #------------------------------------- self.t_posX = self.TrackTree.xPos
        #-------------------------------------- self.t_posY= self.TrackTree.yPos
        #-------------------------------------- self.t_dxdz= self.TrackTree.dxdz
        #-------------------------------------- self.t_dydz= self.TrackTree.dydz
        #-------------------------------------- self.t_iden= self.TrackTree.iden
        #------------------------------ self.t_trackNum= self.TrackTree.trackNum
        #-------------------------------------- self.t_chi2= self.TrackTree.chi2
        #------------------------------------- self.t_ndof = self.TrackTree.ndof

        for data in self.t_posX :
            posX_tmp.append(data)
        for data in self.t_posY :
            posY_tmp.append(data)
        for data in self.t_iden :
            iden_tmp.append(data)
        for data in self.t_dxdz :
            dxdz_tmp.append(data)
        for data in self.t_dydz :
            dydz_tmp.append(data)
        for data in self.t_chi2 :
            chi2_tmp.append(data)
        for data in self.t_ndof :
            ndof_tmp.append(data)
        for data in self.t_trackNum :
            trackNum_tmp.append(data)

        nTrackParams_tmp=self.t_nTrackParams

        trackNum_tmp=[0]

        if len(trackNum_tmp)>0 :
            for j in range(max(trackNum_tmp)+1) :
                aTrack = Track()
                ndata = nTrackParams_tmp/(max(trackNum_tmp)+1)
                #print "nTrackParam : %i len(trackNum %i)"%(nTrackParams_tmp,max(trackNum_tmp)+1)
                aTrack.trackX = posX_tmp[j*ndata:j*ndata+ndata]
                aTrack.trackY = posY_tmp[j*ndata:j*ndata+ndata]
                
        for index,element in enumerate(aTrack.trackX) :
            aTrack.trackX[index] = aTrack.trackX[index]-npix_X*pitchX/2.-pitchX/2.
            aTrack.trackY[index] = aTrack.trackY[index]-npix_Y*pitchY/2.-pitchY/2.
                
            aTrack.iden = iden_tmp[j*ndata:j*ndata+ndata]
            aTrack.chi2 = chi2_tmp[j*ndata:j*ndata+ndata]
            aTrack.trackNum = trackNum_tmp[j*ndata:j*ndata+ndata]
            aTrack.ndof = ndof_tmp[j*ndata:j*ndata+ndata]
            aTrack.dxdz = dxdz_tmp[j*ndata:j*ndata+ndata]
            aTrack.dydz = dydz_tmp[j*ndata:j*ndata+ndata]
            
            #print aTrack.chi2
            #print self.Chi2_Cut
            if(aTrack.chi2[0]<self.Chi2_Cut):
                tracks.append(aTrack)
        self.AllTracks.append(tracks)




    #---------------------------------------------------------------- trackX=[]
    #----------------------------------------------------------------- trackY=[]
    #------------------------------------------------------------------- chi2=0.
    #------------------------------------------------------------------- event=0
    #-------------------------------------------------------------------- ndof=0
    #------------------------------------------------------------------- iden=[]
    #--------------------------------------------------------------- trackNum=[]
    #----------------------------------------------------------------- cluster=0

    def FindMatchedCluster(self,i,r_max,dut=20,distances_histo=None,filter_cluster=False,TrackingRes=0.03) :
        # find the clusters closest to the tracks in this event
        # clusters are matched to tracks using GetPixelResiduals
        # r_max: maximum radial distance allowed between track and any pixel of the cluster

        try : 
            clusters_tmp = self.AllClusters[i]
        except: 
            clusters_tmp = []
        matched_clusters = []
        good_count = 0

        for track in self.AllTracks[i] :
            if len(clusters_tmp)!=0 :
                dut_iden = track.iden.index(dut)
                distances = []

                for cluster in clusters_tmp :
                    
                    #cluster.Print()
                    #track.Print()
                    
                    mdr, mdx, mdy = cluster.GetPixelResiduals(track.trackX[dut_iden],track.trackY[dut_iden])
                    distances.append(mdr)

                    if distances_histo:
                        distances_histo.Fill(mdr)

                cluster = clusters_tmp[distances.index(min(distances))]

                if ((fabs(track.trackX[dut_iden])<=(halfChip_X+self.edge+TrackingRes))and(fabs(track.trackY[dut_iden])<=(halfChip_Y+self.edge+TrackingRes))):
                   #track.Print()
                    #cluster.Print()
                    
                    if (min(distances) < r_max) :
                        # matched cluster
                        cluster.id = good_count
                        track.cluster = cluster.id
                        cluster.tracknum = track.trackNum[dut_iden]
                        matched_clusters.append(cluster)
                        good_count += 1

                    else :
                        # unmatched cluster
                        track.cluster = -11

                else :
                    # track outside DUT
                    track.cluster = -11

            else :
                # no clusters
                track.cluster = -11


        if(filter_cluster) :
            if(i<len(self.AllClusters)):
	    	self.AllClusters[i] = matched_clusters
        else :
            self.AllClusters[i] = clusters_tmp



    def DoPatternRecognition(self,i,tolerance,scale=1) :



        trackDistX = []
        clusterDistX = []

        trackDistY = []
        clusterDistY = []

        tmp_track_X = []
        tmp_track_Y = []
        for ind in range(i,i+scaler):
            for track in self.AllTracks[ind] :
                tmp_track_X.append(track.trackX[3])
                tmp_track_Y.append(track.trackY[3])

        allTrack_tmp = []
        allCluster_tmp = []

        for tracks in self.AllTracks[i:i+scale] :
            for track in tracks:
                allTrack_tmp.append(track)

        allCluster_tmp = self.AllClusters[i]


        pattern = TH2D("","",14000,-npix_Y*pitchY/2,npix_Y*pitchY/2,len(allCluster_tmp)+len(allTrack_tmp),0,len(allCluster_tmp)+len(allTrack_tmp))

        count = 0

        for index,cluster in enumerate(allCluster_tmp) :
            pass



#        # Generate cluster distance lists
#        for index1,cluster1 in enumerate(allCluster_tmp) :
#            tmpx = []
#            tmpy = []
#            for index2,cluster2 in enumerate(allCluster_tmp) :
#                tmpx.append(cluster1.absX-cluster2.absX)
#                tmpy.append(cluster1.absY-cluster2.absY)
#            clusterDistX.append(tmpx)
#            clusterDistY.append(tmpy)
#            print "[PR] cluster %i"%index1
#
#            for value in tmpx :
#                pattern.Fill(value,count,1)
#            count=count+1
#
#            tmpx.sort()
#            tmpy.sort()
#
#            print tmpx
#            print tmpy
#
#        for index1,track1 in enumerate(allTrack_tmp) :
#            tmpx = []
#            tmpy = []
#            if(fabs(track1.trackX[3])<npix_X*pitchX/2 and fabs(track1.trackY[3])<npix_Y*pitchY/2 ):
#                for index2,track2 in enumerate(allTrack_tmp) :
#                    if(fabs(track2.trackX[3])<npix_X*pitchX/2 and fabs(track2.trackY[3])<npix_Y*pitchY/2 ):
#                        tmpx.append(track1.trackX[3]-track2.trackX[3])
#                        tmpy.append(track1.trackY[3]-track2.trackY[3])
#                trackDistX.append(tmpx)
#                trackDistY.append(tmpy)
#                print "[PR] Track %i"%index1
#                for value in tmpx :
#                    pattern.Fill(value,count,2)
#                count=count+1
#                tmpx.sort()
#                tmpy.sort()
#                print tmpx
#                print tmpy
#
#        can=TCanvas()
#        pattern.Draw("colz")
#        c=raw_input()




    def ComputePosition(self,i,method="QWeighted",sigma=0.003):

        
	if(i<len(self.AllClusters)):
    	    for cluster in self.AllClusters[i] :
                cluster.Statistics()
                if (method=="QWeighted"):
                    cluster.GetQWeightedCentroid()
                elif (method=="DigitalCentroid"):
                    cluster.GetDigitalCentroid()
                elif (method=="maxTOT"):
                    cluster.GetMaxTOTCentroid()
                elif (method=="EtaCorrection"):
                    cluster.GetEtaCorrectedQWeightedCentroid(sigma,sigma)



    def ClusterEvent(self,i,method="QWeighted",sigma=0.003):

        self.getEvent(i)

        row_tmp = [s for s in self.p_row]
        col_tmp = [s for s in self.p_col]
        tot_tmp = [s for s in self.p_tot]
	


        # remove hot pixels
        hpindex = 0
        if len(self.hotpixels)>0:
            while(hpindex<len(row_tmp)) :
                if([col_tmp[hpindex],row_tmp[hpindex]] in self.hotpixels):
                    col_tmp.pop(hpindex)
                    row_tmp.pop(hpindex)
                    tot_tmp.pop(hpindex)
                else :
                    hpindex+=1

        # set a maximum number of hit pixels to be clustered (skips large events)
        if len(col_tmp) < 5000:
            try : 
                clusters = self.SciPyClustering(col_tmp,row_tmp,tot_tmp)
            except : 
                clusters=[]
        else:
            print "Event", i, "not beng clustered,", len(col_tmp), "hit pixels"
            clusters=[]

	
        for cluster in clusters :
            cluster.Statistics()
	
        clusters = [cluster for cluster in clusters if cluster.totalTOT>0]

        clusterid = 0

        for cluster in clusters :
            if (method=="QWeighted"):
                cluster.GetQWeightedCentroid()
            elif (method=="DigitalCentroid"):
                cluster.GetDigitalCentroid()
            elif (method=="maxTOT"):
                cluster.GetMaxTOTCentroid()
            elif (method=="EtaCorrection"):
                cluster.GetEtaCorrectedQWeightedCentroid(sigma,sigma)

            cluster.id=clusterid
            clusterid+=1
            cluster=0

	
        self.AllClusters.append(clusters)
        del clusters



    def SciPyClustering(self,col,row,tot):

        pixels = [[col[i],row[i]] for i,x in enumerate(col)]
        if(len(pixels)>1):
            result=fclusterdata(pixels,sqrt(2.),criterion="distance")
            clusters=[Cluster() for i in range(max(result))]
            [clusters[x-1].addPixel(col[j],row[j],tot[j]) for j,x in enumerate(result)]
        else:
            if(len(pixels)==1):
                c=Cluster()
                c.addPixel(col[0],row[0],tot[0])
                clusters=[c]        
	
	return clusters

    def RecursiveClustering(self,row,col,tot) :
        clusters = []
        while(len(row)!=0) :

            cluster = Cluster()
            cluster.addPixel(col[0], row[0], tot[0])
            #print "[DEBUG] adding pixel col=%d row=%d as seed"%(col_tmp[0],row_tmp[0])
            row.pop(0)
            col.pop(0)
            tot.pop(0)
            while(self.addNeighbor(cluster, col,row, tot)>0):
                pass
            clusters.append(cluster)
        return clusters


    def FixedFrameClustering(self,X,Y,TOT):

        frame = [[0 for i in xrange(256)] for j in xrange(256)]
        totframe = [[0 for i in xrange(256)] for j in xrange(256)]
        for i,x in enumerate(X) :
            frame[X[i]][Y[i]]=-1
            totframe[X[i]][Y[i]]=TOT[i]

        cluster_number = 1

        for i,j in [[i,j] for i,j in product(xrange(256),xrange(256)) if frame[i][j]==-1] :
            for u,v in [[u,v] for u,v in product([-1,0,1],[-1,0,1]) if (((i+u>=0 and i+u<=255) and (j+v>=0 and j+v<=255)) and (u!=0 or v!=0)) ] :
                if(frame[i+u][j+v]==-1) :
                    frame[i][j]=cluster_number
                    frame[i+u][j+v]=cluster_number
                    cluster_number+=1
                elif (frame[i+u][j+v]>0) :
                    frame[i][j]= frame[i+u][j+v]

        clusters = {}
        for i,j in [[i,j] for i,j in product(xrange(256),xrange(256)) if frame[i][j]>0] :
            try :
                clusters[frame[i][j]].addPixel(i,j,totframe[i][j])
            except KeyError :
                clusters[frame[i][j]]=Cluster()
                clusters[frame[i][j]].addPixel(i,j,totframe[i][j])

        del frame
        del totframe
        #print clusters.ite
        return clusters.values()


    def addNeighbor(self,cluster,col,row,tot):
        counter =0
        i=0
        j=0
        len_col=len(col)
        len_clu_col=len(cluster.col)
        while(i<len_col):
            j=0
            while(j<len_clu_col):

                if((col[i]-cluster.col[j])**2>1) :
                    j+=1
                    continue

                if((row[i]-cluster.row[j])**2>1) :
                    j+=1
                    continue

                cluster.addPixel(col[i],row[i],tot[i])

            #print "[DEBUG] after adding pixel col=%d row=%d to existing cluster as neighbor to x=%d y=%d "%(col[i],row[i],cluster.col[j],cluster.row[j])

                col.pop(i)
                row.pop(i)
                tot.pop(i)
                counter+=1
                i+=-1
                len_col=len(col)
                len_clu_col=len(cluster.col)
                break
            i+=1
        return counter

    def PrintClusters(self,i):
        print "########## Event %d ##########"%i
        for j,c in enumerate(self.AllClusters[i]):
            if(c.totalTOT<self.EnergyCut):
                print "##### Cluster %d #####"%j
                c.Print()

    def PrintTBranchElement(self):
        for event_i in self.TrackTree :
            print "new event........................"
            print self.TrackTree.xPos
            print "number of entries : "
            print self.TrackTree.xPos.size()
            for entry_j in range (0,self.TrackTree.xPos.size()) :
                print self.TrackTree.xPos[entry_j]
            print "trackNum : "
            print self.TrackTree.trackNum.size()
            for entry_j in range (0,self.TrackTree.trackNum.size()) :
                print self.TrackTree.trackNum[entry_j]
            print "nTrackParams : "
            print self.TrackTree.nTrackParams
