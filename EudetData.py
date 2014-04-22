from ROOT import *
import ROOT
from math import fsum,fabs
from array import array
import pyximport; pyximport.install(pyimport=True)
from ROOT import TFile,TH1D,TH2D,TF1
from Cluster import *
from Track import  *
#gROOT.LoadMacro("Track.C+")
#from ROOT import Track
from Constant import *
from ToolBox import *
from PersistentList import *
from scipy.cluster.hierarchy import fclusterdata
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

    hit_map =[[0 for x in xrange(npix_X)] for x in xrange(npix_Y)]
    frequency_map = [[0 for x in xrange(npix_X)] for x in xrange(npix_Y)]
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





    def FilterHotPixel(self,threshold,Nevents=-1,scaler=1.):

        # Threshold between [0,1], cut firing frequency of pixels
        # Nevent is how many event for building the frequency Matrix

        nevent_tmp =0
        n_max = 0

#        last =0
#        now =0

#        counter_tmp=0

        histo_frequency = TH1D("freq","Firing Frequency",10000,0,1)
        histo_hotpixel = TH2D("hot","hot pixel map",256,0,255,256,0,255)

        if Nevents>self.p_nEntries or Nevents==-1:
            n_max = self.p_nEntries
        else :
            n_max = Nevents


        for ii in range(n_max) :
            self.getEvent(ii)
            if ii%1000==0 :
                print " [Hot Pixel Finder] Parsing event %i"%ii
            for jj in range(len(self.p_row)) :

#                last = self.hit_map[129][20]

                self.hit_map[self.p_col[jj]][self.p_row[jj]]+=1./scaler

#                now = self.hit_map[129][20]

#                if last != now :
#                    counter_tmp+=1
#                    print "adding %i %i , counter = %i %i"%(self.p_col[jj],self.p_row[jj],counter_tmp,self.hit_map[self.p_col[jj]][self.p_row[jj]])

            nevent_tmp+=1



        for ii in range(npix_X):
            for jj in range(npix_Y) :

#                if ii==129 and jj==20 :
#                    print "[DEDUG !!!!!!!!!] %f nevents : %i"%(self.hit_map[ii][jj],nevent_tmp)
                self.frequency_map[ii][jj]=self.hit_map[ii][jj]*(1.0/float(nevent_tmp))
#                if ii==129 and jj==20 :
#                    print "[DEDUG !!!!!!!!!] %f nevents : %i"%(self.frequency_map[ii][jj],nevent_tmp)
                histo_hotpixel.Fill(ii,jj,self.frequency_map[ii][jj])

                #if(self.frequency_map[ii][jj]>1):
                    #print ii,jj,self.frequency_map[ii][jj]

                histo_frequency.Fill(self.frequency_map[ii][jj])
                if(self.frequency_map[ii][jj]>threshold):
                    #print "hotpixel %i %i"%(ii,jj)
                    self.hotpixels.append([ii,jj])

        print "##### Hot Pixel Report #####"
        print " %i Hot pixel found at  : "%(len(self.hotpixels))
        print self.hotpixels
        print "############################"
        return histo_hotpixel,histo_frequency


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



        self.entry = self.pixelTree.GetEntry(i)

        self.p_nHits= self.pixelTree.nPixHits
        self.p_col= self.pixelTree.col
        self.p_row= self.pixelTree.row
        self.p_tot= self.pixelTree.tot
        self.p_lv1= self.pixelTree.lv1
        self.p_chip= self.pixelTree.chip
        self.p_iden= self.pixelTree.iden
        self.p_euEv = self.pixelTree.euEvt

        for index,totvalue in enumerate(self.p_tot) :
            self.p_tot[index]=float(totvalue)/self.scale


    def WriteReconstructedData(self,filename,dut=6) : 
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
            cluster.absY = 	self.pixelTree.absY
            cluster.resX = self.pixelTree.resX
            cluster.resY = self.pixelTree.resY
            cluster.id = self.pixelTree.id
            cluster.trackNum = self.pixelTree.trackNum
            self.AllClusters[event].append(cluster)

    
    
    def DumpTrackTree(self,outfile,dut=6):

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
        

    def DumpClusterTree(self,outfile,dut=6):

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
        
                for i in range(len(cluster.col)) :
                    col[i]=cluster.col[i]
                for i in range(len(cluster.row)) :
                    row[i]=cluster.row[i]
                for i in range(len(cluster.tot)) :
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
                    trackX[0]=self.AllTracks[j][cluster.tracknum].trackX[self.AllTracks[j][cluster.tracknum].iden.index(dut)]
                    trackY[0]=self.AllTracks[j][cluster.tracknum].trackY[self.AllTracks[j][cluster.tracknum].iden.index(dut)]
                clusterTree.Fill()
        outfile.Write()

    def IsInEdges(self,track,dut=6):
        is_in = False
        if(fabs(track.trackX[track.iden.index(dut)])<=(halfChip_X+self.edge) and fabs(track.trackY[track.iden.index(dut)])<=(halfChip_Y+self.edge)):
            is_in = True
            if(fabs(track.trackX[track.iden.index(dut)])<=(halfChip_X) and fabs(track.trackY[track.iden.index(dut)])<=(halfChip_Y)):
                is_in=False
        return is_in
        
        


    def ComputeResiduals(self,i,dut=6) :

        nmatch = 0.
        nmatch_edge = 0
        for track in self.AllTracks[i] :
            for index,cluster in enumerate(self.AllClusters[i]) :
                #print "looking at track with id %i and cluster with id %i"%(track.cluster,cluster.id)
                if cluster.id == track.cluster  :
                    cluster.GetResiduals(track.trackX[track.iden.index(dut)],track.trackY[track.iden.index(dut)])
                    #print "after match resX : %f resY : %f"%(cluster.resX,cluster.resY)
                    nmatch+=1
                    if(self.IsInEdges(track)):
                        nmatch_edge+=1
        return nmatch,nmatch_edge
#     def ComputeResiduals(self,i) :
#         for cluster in self.AllClusters[i] :
#             cluster.GetResiduals(self.t_posX[3],self.t_posY[3])

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
        for j in range(max(trackNum_tmp)+1) :
            aTrack = Track()
            ndata = nTrackParams_tmp/(max(trackNum_tmp)+1)
            #print "nTrackParam : %i len(trackNum %i)"%(nTrackParams_tmp,max(trackNum_tmp)+1)
            aTrack.trackX = posX_tmp[j*ndata:j*ndata+ndata]
            aTrack.trackY = posY_tmp[j*ndata:j*ndata+ndata]
            for index,element in enumerate(aTrack.trackX) :
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
                aTrack.trackX[index] = aTrack.trackX[index]-npix_X*pitchX/2.-pitchX/2.
                aTrack.trackY[index] = aTrack.trackY[index]-npix_Y*pitchY/2.-pitchY/2.
            aTrack.iden = iden_tmp[j*ndata:j*ndata+ndata]
            aTrack.chi2 = chi2_tmp[j*ndata:j*ndata+ndata]
            aTrack.trackNum = trackNum_tmp[j*ndata:j*ndata+ndata]
            aTrack.ndof = ndof_tmp[j*ndata:j*ndata+ndata]
            aTrack.dxdz = dxdz_tmp[j*ndata:j*ndata+ndata]
            aTrack.dydz = dydz_tmp[j*ndata:j*ndata+ndata]

            if(aTrack.chi2[0]<self.Chi2_Cut):
                tracks.append(aTrack)

            #print aTrack.chi2
        self.AllTracks.append(tracks)




    #---------------------------------------------------------------- trackX=[]
    #----------------------------------------------------------------- trackY=[]
    #------------------------------------------------------------------- chi2=0.
    #------------------------------------------------------------------- event=0
    #-------------------------------------------------------------------- ndof=0
    #------------------------------------------------------------------- iden=[]
    #--------------------------------------------------------------- trackNum=[]
    #----------------------------------------------------------------- cluster=0

    def FindMatchedCluster(self,i,r_max,dut=6,filter_cluster=True,TrackingRes = 0.003) :

        # i : event number
        # r_max_X,Y maximum distance in X,Y between track and cluster
        # dut = iden of the Device Under Test

        clusters_tmp = self.AllClusters[i]
        good_clusters = []
        good_cnt =0
        for track in self.AllTracks[i] :
            if len(clusters_tmp)!=0 :
                dut_iden = track.iden.index(dut)
                distances = []
                for cluster in clusters_tmp :
                    cluster.GetResiduals(track.trackX[dut_iden],track.trackY[dut_iden])
                    distances.append(sqrt(cluster.resX**2 + cluster.resY**2))   
#                if(i%250==0):
#                    print "example of distance vector"
#                    print [x for x in distances if x<1.]    
                cluster = clusters_tmp[distances.index(min(distances))]                     
                if((fabs(track.trackX[dut_iden])<=(halfChip_X+self.edge+TrackingRes))and(fabs(track.trackY[dut_iden])<=(halfChip_Y+self.edge+TrackingRes))):
                    if((cluster.resX**2 + cluster.resY**2)<r_max**2) :
                        cluster.id=good_cnt
                        track.cluster=cluster.id
                        cluster.tracknum=track.trackNum[dut_iden]
                        good_clusters.append(cluster)
                        good_cnt+=1
#                        print "resX : %f resY : %f"%(cluster.resX,cluster.resY)
#                        cluster.Print()
#                        track.Print()
                        #clusters_tmp.pop(index)
                        #break
                    else :
                        track.cluster=-11
#                        print "Found an unmatched "
#                        print "resX : %f resY : %f"%(cluster.resX,cluster.resY)
#                        cluster.Print()
#                        track.Print()
                else :
                    track.cluster=-11
#       for u,cl1 in enumerate(good_clusters) :
#               for v,cl2 in enumerate(good_clusters[u+1:]) :
#                       if(cl1.tracknum==cl2.tracknum) :
#                               if((cl1.resX**2 + cl2.resY**2)>=(cl2.resX**2 + cl2.resY**2)) :
#                                       good_clusters.pop(v)
#                                       break
#                               else :
#                                       good_clusters.pop(u)
#                                       break


        if(filter_cluster) :
            self.AllClusters[i]=good_clusters
        else :
            self.AllClusters[i]=clusters_tmp



#             if(track.cluster!=-11):
#
#                 print "#### a Match #####"
#                 self.AllClusters[i][track.cluster].Print()
#                 track.Print()
#                 print "#################"





#     def ClusterEvent(self,i,method="QWeighted"):


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




    def ClusterEvent(self,i,method="QWeighted",sigma=0.003, scaler=1):

        self.getEvent(i)

        clusterid=0

        row_tmp = [s for s in self.p_row]
        col_tmp = [s for s in self.p_col]
        tot_tmp = [s for s in self.p_tot]


#       for index in self.p_row:
#            row_tmp.append(index)
#        for index in self.p_col:
#            col_tmp.append(index)
#        for index in self.p_tot:
#            tot_tmp.append(index)

        hpindex=0

        if len(self.hotpixels)>0:
            while(hpindex<len(row_tmp)) :
                if([col_tmp[hpindex],row_tmp[hpindex]] in self.hotpixels):
                        #print "Removing hot pixel x: %i y:%i"%(col_tmp[hpindex],row_tmp[hpindex])
                    col_tmp.pop(hpindex)
                    row_tmp.pop(hpindex)
                    tot_tmp.pop(hpindex)
                else :
                    hpindex+=1



        try : 
            clusters = self.SciPyClustering(col_tmp,row_tmp,tot_tmp)
        except : 
            clusters=[]
        
        
        for cluster in clusters : 
            cluster.Statistics()
            if (method=="QWeighted"):
                cluster.GetQWeightedCentroid()
            elif (method=="DigitalCentroid"):
                cluster.GetDigitalCentroid()
            elif (method=="maxTOT"):
                cluster.GetMaxTOTCentroid()
            elif (method=="EtaCorrection"):
#                 cluster.GetEtaCorrectedQWeightedCentroid()
                cluster.GetEtaCorrectedQWeightedCentroid(sigma)
            # to be implemented in the future:
            # digital, maxTOT/maxQ, eta corrected

            cluster.id=clusterid
            clusterid+=1
            cluster=0

        for ind in range(i,i+scaler):
            self.AllClusters.append(clusters)
        del clusters	
	
	
	#while(len(row_tmp)!=0) :

            #cluster = Cluster()
            #cluster.addPixel(col_tmp[0], row_tmp[0], tot_tmp[0])
            ##print "[DEBUG] adding pixel col=%d row=%d as seed"%(col_tmp[0],row_tmp[0])
            #row_tmp.pop(0)
            #col_tmp.pop(0)
            #tot_tmp.pop(0)
            #while(self.addNeighbor(cluster, col_tmp, row_tmp, tot_tmp)>0):
                #pass


            #cluster.Statistics()
            #if (method=="QWeighted"):
                #cluster.GetQWeightedCentroid()
            #elif (method=="DigitalCentroid"):
                #cluster.GetDigitalCentroid()
            #elif (method=="maxTOT"):
                #cluster.GetMaxTOTCentroid()
            #elif (method=="EtaCorrection"):
##                 cluster.GetEtaCorrectedQWeightedCentroid()
                #cluster.GetEtaCorrectedQWeightedCentroid(sigma)
            ## to be implemented in the future:
            ## digital, maxTOT/maxQ, eta corrected

            #cluster.id=clusterid
            #clusters.append(cluster)
            #clusterid+=1
            #cluster=0

        #for ind in range(i,i+scaler):
            #self.AllClusters.append(clusters)
        #del clusters


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
