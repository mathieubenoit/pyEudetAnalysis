'''
Created on Oct 9, 2013

@author: mbenoit
'''
import time,os
from math  import sqrt
from numpy import matrix
from array import array
from ROOT import *
import pyximport; pyximport.install(pyimport=True)

from itertools import product

class CalibTreeMaker:
    '''
    classdocs
    '''
    last_time = time.time()
    AllClusters = []

    def __init__(self,filename,outfile):
        '''
        Constructor
        '''
        self.ReadFile(filename,outfile)
        
    def ReadFile(self,filename,outfile):
        
        
        data_file = open(filename,"r")
        lines = data_file.readlines()
        
        X = []
        Y = []
        TOT = []
        nFrames = 0
       
        size1cnt = 0
        self.last_time=time.time()
        
        outfile = TFile(outfile,'recreate')
        pixelTree = TTree("pixels","pixelstree")
        
        xt=array( 'i', [ 0 ] )
        yt=array( 'i', [ 0 ] )
        tott=array( 'i', [ 0 ])
        
        pixelTree.Branch( 'col', xt, 'col/I' )
        pixelTree.Branch( 'row', yt, 'row/I' )      
        pixelTree.Branch( 'tot', tott, 'tot/I' )       
        
        
        self.last_time = time.time()
        for line in lines : 
            
            if "#" in line : 
                nFrames+=1
                clusters = self.SimpleClustering(X,Y,TOT)
                for cluster in clusters : 
                    size1cnt+=1           
                    xt[0]=cluster[0]
                    yt[0]=cluster[1]
                    tott[0]=cluster[2]
                    pixelTree.Fill()              
                
                if(nFrames%1000==0):
                    print "Processed Frame %i (%.5fs/frame)"%(nFrames,(time.time()-self.last_time )/1000.)
                    self.last_time = time.time()
                X = []
                Y = []
                TOT = []
                
            
            else : 
                data = line.split()
                X.append(int(data[0]))
                Y.append(int(data[1]))
                TOT.append(int(data[2]))
                del data
        
        del lines 
        data_file.close()
            
        
        pixelTree.Write() 
        outfile.Close()
        
        print "found %i single pixel clusters"%size1cnt
 
    
    
   
    def SimpleClustering(self,X,Y,TOT):
        
        frame = [[0 for i in xrange(256)] for j in xrange(256)]     
        for i,x in enumerate(X) :
            frame[X[i]][Y[i]]=1
                
        for i,j in [[i,j] for i,j in product(xrange(256),xrange(256)) if frame[i][j]==1] :                     
            for u,v in [[u,v] for u,v in product([-1,0,1],[-1,0,1]) if (((i+u>=0 and i+u<=255) and (j+v>=0 and j+v<=255)) and (u!=0 or v!=0)) ] :
                if(frame[i+u][j+v]==1) : 
                    frame[i][j]=0
                    frame[i+u][j+v]=0
                                 
        pixels = [[X[i],Y[i],TOT[i]] for i in range(len(X)) if frame[X[i]][Y[i]]==1]
     
        del frame
        return pixels

    
    
    
    
#    def SciPyClustering(self,col,row,tot):
#        
#        pixels = [[col[i],row[i]] for i,x in enumerate(col)]
#        if(len(pixels)>1):
#            result=fclusterdata(pixels,sqrt(2.),criterion="distance",method="single")        
#            clusters=[Cluster() for i in range(max(result))]
#            [clusters[x-1].addPixel(col[j],row[j],tot[j]) for j,x in enumerate(result)]
#        else:
#            if(len(pixels)==1):
#                c=Cluster()
#                c.addPixel(col[0],row[0],tot[0]) 
#                clusters=[c]
#        
#        return clusters
#    
#    
#    def Clustering(self,row,col,tot) :
#        clusters = []
#        while(len(row)!=0) :
#
#            cluster = Cluster()
#            cluster.addPixel(col[0], row[0], tot[0])
#            #print "[DEBUG] adding pixel col=%d row=%d as seed"%(col_tmp[0],row_tmp[0])
#            row.pop(0)
#            col.pop(0)
#            tot.pop(0)
#            while(self.addNeighbor(cluster, col,row, tot)>0):
#                pass
#            clusters.append(cluster)
#        return clusters
# 
#    def addNeighbor(self,cluster,col,row,tot):
#        counter =0
#        i=0
#        j=0
#                        
#        len_col=len(col)
#        len_clu_col=len(cluster.col)
#        while(i<len_col):
#            j=0
#            while(j<len_clu_col):
#
#                if((col[i]-cluster.col[j])**2>1) :
#                    j+=1
#                    continue
#
#                if((row[i]-cluster.row[j])**2>1) :
#                    j+=1
#                    continue
#
#                cluster.addPixel(col[i],row[i],tot[i])
#
#            #print "[DEBUG] after adding pixel col=%d row=%d to existing cluster as neighbor to x=%d y=%d "%(col[i],row[i],cluster.col[j],cluster.row[j])
#
#                col.pop(i)
#                row.pop(i)
#                tot.pop(i)
#                counter+=1
#                i+=-1
#                len_col=len(col)
#                len_clu_col=len(cluster.col)
#                break
#            i+=1
#        return counter 
#aCalibDataSet=  CalibTreeMaker("/home/mbenoit/Calibration_Data_Buffer/A06-W0110/Fe55/Fe55_08-10-2013_A06-W0110_r0000_.dat","/home/mbenoit/Calibration_Data_Buffer/A06-W0110/Fe55/Fe55_A06-W0110_0.root")

#aCalibDataSet=  CalibTreeMaker("/home/mbenoit/Calibration_Data_Buffer/C06-W0110/Cd109/Cd109_10-10-2013_C06-W0110_dat_r0000_","/home/mbenoit/Calibration_Data_Buffer/C06-W0110/Cd109/Cd109_C06-W0110_0.root")
#aCalibDataSet=  CalibTreeMaker("/home/mbenoit/Calibration_Data_Buffer/C06-W0110/Cd109/Cd109_10-10-2013_C06-W0110_dat_r0001_","/home/mbenoit/Calibration_Data_Buffer/C06-W0110/Cd109/Cd109_C06-W0110_1.root")
#aCalibDataSet=  CalibTreeMaker("/home/mbenoit/Calibration_Data_Buffer/C06-W0110/Cd109/Cd109_10-10-2013_C06-W0110_dat_r0002_","/home/mbenoit/Calibration_Data_Buffer/C06-W0110/Cd109/Cd109_C06-W0110_2.root")

for i in range(14,43) : 
    aCalibDataSet=  CalibTreeMaker("/home/mbenoit/Calibration_Data_Buffer/C06-W0110/Cd109/Cd109_10-10-2013_C06-W0110_dat_r00%02i_"%i,"/home/mbenoit/Calibration_Data_Buffer/C06-W0110/Cd109/Cd109_C06-W0110_%i.root"%i)
   

#aCalibDataSet=  CalibTreeMaker("/home/mbenoit/Calibration_Data_Buffer/C06-W0110/Fe55/Fe55_09-10-2013_C06-W0110__r0000_.dat","/home/mbenoit/Calibration_Data_Buffer/C06-W0110/Fe55/Fe55_C06-W0110_0.root")
#aCalibDataSet=  CalibTreeMaker("/home/mbenoit/Calibration_Data_Buffer/C06-W0110/Fe55/Fe55_09-10-2013_C06-W0110__r0001_.dat","/home/mbenoit/Calibration_Data_Buffer/C06-W0110/Fe55/Fe55_C06-W0110_1.root")  
#aCalibDataSet=  CalibTreeMaker("/home/mbenoit/Calibration_Data_Buffer/C06-W0110/Fe55/Fe55_09-10-2013_C06-W0110__r0002_.dat","/home/mbenoit/Calibration_Data_Buffer/C06-W0110/Fe55/Fe55_C06-W0110_2.root")  
#aCalibDataSet=  CalibTreeMaker("/home/mbenoit/Calibration_Data_Buffer/C06-W0110/Fe55/Fe55_09-10-2013_C06-W0110__r0003_.dat","/home/mbenoit/Calibration_Data_Buffer/C06-W0110/Fe55/Fe55_C06-W0110_3.root")    