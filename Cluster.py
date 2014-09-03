'''
author: Anne-Laure pequegnot
'''

from math import fsum
from Constant import *
from ROOT import TMath
from ToolBox import *



#
#Compute the additional shift for the hit position du to the charge sharing (eta correction) when neighbor pixels are on the same raw or the same column (neighbor pixels aligned)
#first parameter:sigma of the eta correction (charge sharing)
#second parameter:relative charge i.e. Qrel = (charge of the pixel with the highest energy)/(total charge of the cluster)
#shiftLat
def shiftLat(sigma_tmp,Qrel_tmp):
    return sigma_tmp*TMath.ErfInverse(2.*Qrel_tmp-1.)

#
#Compute the additional shift for the hit position du to the charge sharing (eta correction) when neighbor pixels are on a diadonal
#first parameter:sigma of the eta correction (charge sharing)
#second parameter:relative charge i.e. Qrel = (charge of the pixel with the highest energy)/(total charge of the cluster)
#
def shiftDiag(sigmaX_tmp,sigmaY_tmp,Qrel_tmp):
    return sigmaX_tmp*TMath.ErfInverse(2.*Qrel_tmp-1.)*1./sqrt(2.),sigmaY_tmp*TMath.ErfInverse(2.*Qrel_tmp-1.)*1./sqrt(2.)



###############################################################################################################################
#
#        Class for the clusters and their properties
#
###############################################################################################################################

class Cluster:

    col = []
    row = []
    tot = []
    sizeX = 0
    sizeY = 0
    size  = 0
    totalTOT =0
    aspectRatio = 0

# local coordinates
    relX = 0.
    relY = 0.
    relZ = 0.

# telescope coordinates
    absX =-10000.
    absY =-10000.
    absZ =-10000.

    resX = -10000.
    resY = -10000.
    id = 0

    # track number
    tracknum = -1

    def __init__(self):
        self.sizeX = 0
        self.sizeY = 0
        self.size  = 0
        self.col = []
        self.row = []
        self.tot = []

    def addPixel(self,col,row,tot):
        self.col.append(col)
        self.row.append(row)
        self.tot.append(tot)

    def Print(self):
        for i in range(len(self.col)):
            print "x:%d y%d tot:%.3f"%(self.col[i],self.row[i],self.tot[i])
        print "Cluster Total size = %d , in X = %d Y = %d , Aspect Ratio = %.3f , Total Energy (keV) = %.1f ID = %i"%(self.size,self.sizeX,self.sizeY,self.aspectRatio,self.totalTOT,self.id)
        print "Position in sensor X = %.3f Y = %.3f"%(self.absX,self.absY)

    def Statistics(self) :
        self.totalTOT=fsum(self.tot)
        self.size=len(self.col)
        self.sizeX=max(self.col)-min(self.col)+1
        self.sizeY=max(self.row)-min(self.row)+1
        self.aspectRatio=float(self.sizeY)/self.sizeX

#
#compute the hit position as the mean of the fired pixels positions weighted by their deposited energy
#
    def GetQWeightedCentroid(self) :
        self.relX=0.
        self.relY=0.

        if(self.totalTOT==0):
            self.GetDigitalCentroid()
        else :
            for index,tot_tmp in enumerate(self.tot) :
                self.relX+=(self.col[index]*pitchX)*tot_tmp
                self.relY+=(self.row[index]*pitchY)*tot_tmp

            self.relX = self.relX/self.totalTOT + pitchX/2.
            self.relY = self.relY/self.totalTOT + pitchY/2.

            self.absX=self.relX - npix_X*pitchX/2.
            self.absY=self.relY - npix_Y*pitchY/2.
            self.absZ=0

#
#compute the hit position as the mean of the fired pixels positions (digital method)
#
    def GetDigitalCentroid(self) :
        self.relX=0.
        self.relY=0.
        for index,col_tmp in enumerate(self.col) :
            self.relX+=(self.col[index]*pitchX)
            self.relY+=(self.row[index]*pitchY)

        self.relX = self.relX/len(self.col) + pitchX/2.
        self.relY = self.relY/len(self.row) + pitchY/2.

        self.absX=self.relX - npix_X*pitchX/2.
        self.absY=self.relY - npix_Y*pitchY/2.
        self.absZ=0

#
#compute the hit position as the center of the fired pixel with the highest energy
#
    def GetMaxTOTCentroid(self) :
        maxTOTindex_tmp=0
        maxTOT_tmp=self.tot[0]
        for index,tot_tmp in enumerate(self.tot) :
            if self.tot[index]>maxTOT_tmp:
                maxTOT_tmp=self.tot[index]
                maxTOTindex_tmp=index

        self.relX=self.col[maxTOTindex_tmp]*pitchX + pitchX/2.
        self.relY=self.row[maxTOTindex_tmp]*pitchY + pitchY/2.

        self.absX=self.relX - npix_X*pitchX/2.
        self.absY=self.relY - npix_Y*pitchY/2.
        self.absZ=0

#
#compute the hit position as the Qweighted method but adding an eta correction du to the charge sharing between the fired pixels
#
    def GetEtaCorrectedQWeightedCentroid(self,sigmaX=0.003,sigmaY=0.003) :
#     def GetEtaCorrectedQWeightedCentroid(self) :
#         global n_sizeX2sizeY2
        maxTOTindex_tmp = 0
        minTOTindex_tmp = 0
        maxTOT_tmp=self.tot[0]

        maxTOTindex1_tmp = 0
        minTOTindex1_tmp = 0
        maxTOT1_tmp=self.tot[0]
        maxTOTindex2_tmp = 0
        minTOTindex2_tmp = 0
        index2_tmp = []
        maxTOT2_tmp=self.tot[0]

        self.relX=-1000
        self.relY=-1000
        self.absX=-1000
        self.absY=-1000
        #self.GetQWeightedCentroid()

        if(self.size==2) :
            for index,tot_tmp in enumerate(self.tot) :
#looking for the pixel with high energy deposition
                if self.tot[index]>maxTOT_tmp:
                    maxTOT_tmp=self.tot[index]
                    maxTOTindex_tmp=index
                if(maxTOTindex_tmp==1)  :
                    minTOTindex_tmp =0
                else :
                    minTOTindex_tmp = 1


#computing the relative charge i.e. Qrel = (charge of the pixel with the highest energy)/(total charge of the cluster)
            Qrel = (self.tot[maxTOTindex_tmp])/(self.totalTOT)

#check-cross
#                 self.Print()
#                 print "TMath.ErfInverse(2.*Qrel-1.) : %f"%(float(TMath.ErfInverse(2.*Qrel-1.)))
#                 print "sigma : %f"%(float(sigma))
#                 print sigma
#                 print "sigma*TMath.ErfInverse(2.*Qrel-1.) : %f"%(float(sigma)*TMath.ErfInverse(2.*Qrel-1.))
#                 print 'self.sizeX : %i'%self.sizeX
#                 print 'self.sizeY : %i'%self.sizeY
#                 print 'maxTOTindex_tmp : %i'%maxTOTindex_tmp
#                 print 'self.row[maxTOTindex_tmp]*pitchY : %f'%(float(self.row[maxTOTindex_tmp])*pitchY)
#                 print 'self.col[maxTOTindex_tmp]*pitchY : %f'%(float(self.col[maxTOTindex_tmp])*pitchY)

#distinguishing the 3 cases 'cluster size 2x1','cluster size 1x2' and 'cluster size 2x2'
            if(self.sizeX==2 and self.sizeY==1) :
            #cluster size 2x1
                if(self.col[maxTOTindex_tmp]>self.col[minTOTindex_tmp]) :
                #neighbor on the left side
                    self.relX = self.col[maxTOTindex_tmp]*pitchX + shiftLat(sigmaX,Qrel)
                    self.relY = self.row[maxTOTindex_tmp]*pitchY +pitchY/2.
                elif(self.col[maxTOTindex_tmp]<self.col[minTOTindex_tmp]) :
                #neighbor on the right side
                    self.relX = (self.col[maxTOTindex_tmp]+1.)*pitchX - shiftLat(sigmaX,Qrel)
                    self.relY = self.row[maxTOTindex_tmp]*pitchY +pitchY/2.
#                     print "relX : %f"%float(self.relX)
#                     print "relY : %f"%float(self.relY)
            elif(self.sizeX==1 and self.sizeY==2) :
            #cluster size 1x2
                if(self.row[maxTOTindex_tmp]>self.row[minTOTindex_tmp]) :
                #neighbor on the bottom side
                    self.relX = self.col[maxTOTindex_tmp]*pitchX + pitchX/2.
                    self.relY = self.row[maxTOTindex_tmp]*pitchY + shiftLat(sigmaY,Qrel)
                elif(self.row[maxTOTindex_tmp]<self.row[minTOTindex_tmp]) :
                #neighbor on the top side
                    self.relX = self.col[maxTOTindex_tmp]*pitchX + pitchX/2.
                    self.relY = (self.row[maxTOTindex_tmp]+1.)*pitchY - shiftLat(sigmaY,Qrel)
            elif(self.sizeX==2 and self.sizeY==2) :
            #cluster size 2 with sizeX = 2 and sizeY = 2 i.e. 2 pixels on a diagonal
                #print"cluster sizeX : 2 ; sizeY : 2 ; size : 2"
#                     n_sizeX2sizeY2 = n_sizeX2sizeY2 + 1
                shiftX,shiftY = shiftDiag(sigmaX,sigmaY,Qrel)
                if(self.col[maxTOTindex_tmp]>self.col[minTOTindex_tmp] and self.row[maxTOTindex_tmp]>self.row[minTOTindex_tmp]) :
                    self.relX = self.col[maxTOTindex_tmp]*pitchX + shiftX
                    self.relY = self.row[maxTOTindex_tmp]*pitchY + shiftY
                    #self.relX = 1111
                    #self.relY = 1111

                elif(self.col[maxTOTindex_tmp]<self.col[minTOTindex_tmp] and self.row[maxTOTindex_tmp]<self.row[minTOTindex_tmp]) :
                    self.relX = (self.col[maxTOTindex_tmp]+1.)*pitchX - shiftX
                    self.relY = (self.row[maxTOTindex_tmp]+1.)*pitchY - shiftY
                    #self.relX = 1111
                    #self.relY = 1111

                elif(self.col[maxTOTindex_tmp]>self.col[minTOTindex_tmp] and self.row[maxTOTindex_tmp]<self.row[minTOTindex_tmp]) :
                    self.relX = self.col[maxTOTindex_tmp]*pitchX + shiftX
                    self.relY = (self.row[maxTOTindex_tmp]+1.)*pitchY - shiftY
                    #self.relX = 1111
                    #self.relY = 1111


                elif(self.col[maxTOTindex_tmp]<self.col[minTOTindex_tmp] and self.row[maxTOTindex_tmp]>self.row[minTOTindex_tmp]) :
                    self.relX = (self.col[maxTOTindex_tmp]+1.)*pitchX - shiftX
                    self.relY = self.row[maxTOTindex_tmp]*pitchY + shiftY
                    #self.relX = 1111
                    #self.relY = 1111



        elif(self.size==4) :
            if(self.sizeX==2 and self.sizeY==2) :

#looking for the pixel with the highest energy
                for index,tot_tmp in enumerate(self.tot) :
                    if self.tot[index]>maxTOT_tmp:
                        maxTOT1_tmp=self.tot[index]
                        maxTOTindex1_tmp=index

#finding the second pixel to build a pixels pair on a diagonal
                for index,tot_tmp in enumerate(self.tot) :
                    if (abs(self.col[maxTOTindex1_tmp] - self.col[index]) == 1 and abs(self.row[maxTOTindex1_tmp] - self.row[index]) == 1) :
                        minTOTindex1_tmp = index

#computing the first hit position shift du to the eta correction for this first pixels couple
                Qrel1 = (self.tot[maxTOTindex1_tmp])/(self.tot[maxTOTindex1_tmp] + self.tot[minTOTindex1_tmp])
                shift1X,shift1Y = shiftDiag(sigmaX,sigmaY,Qrel1)

#finding the second pixels pair on a diagonal and looking for the one between the 2 which has the highest energy
                for index,tot_tmp in enumerate(self.tot) :
                    if (index != maxTOTindex1_tmp and index != minTOTindex1_tmp) :
                        index2_tmp.append(index)

                if(self.tot[index2_tmp[0]] > self.tot[index2_tmp[1]]) :
                    maxTOTindex2_tmp = index2_tmp[0]
                    minTOTindex2_tmp = index2_tmp[1]

                else :
                    maxTOTindex2_tmp = index2_tmp[1]
                    minTOTindex2_tmp = index2_tmp[0]

#computing the second hit position shift du to the eta correction for this second pixels couple
                Qrel2 = (self.tot[maxTOTindex2_tmp])/(self.tot[maxTOTindex2_tmp] + self.tot[minTOTindex2_tmp])
                shift2X,shift2Y = shiftDiag(sigmaX,sigmaY,Qrel2)

#computing the eta corrected position of the hit
                if(self.col[maxTOTindex1_tmp] < self.col[minTOTindex1_tmp] and self.row[maxTOTindex1_tmp] > self.row[minTOTindex1_tmp]) :
                #pixel from pair 1 with highest energy on the left top corner
                    if(self.col[maxTOTindex2_tmp] < self.col[minTOTindex2_tmp]) :
                    #pixel from pair 2 with highest energy on the left bottom corner
                        self.relX = (self.col[maxTOTindex1_tmp]+1.)*pitchX - shift1X - shift2X
                        self.relY = (self.row[maxTOTindex1_tmp])*pitchY + shift1Y - shift2Y
                    else :
                    #pixel from pair 2 with highest energy on the right top corner
                        self.relX = (self.col[maxTOTindex1_tmp]+1.)*pitchX - shift1X + shift2X
                        self.relY = (self.row[maxTOTindex1_tmp])*pitchY + shift1Y + shift2Y
                elif(self.col[maxTOTindex1_tmp] > self.col[minTOTindex1_tmp] and self.row[maxTOTindex1_tmp] < self.row[minTOTindex1_tmp]) :
                #pixel from pair 1 with highest energy on the right bottom corner
                    if(self.col[maxTOTindex2_tmp] < self.col[minTOTindex2_tmp]) :
                    #pixel from pair 2 with highest energy on the left bottom corner
                        self.relX = (self.col[maxTOTindex1_tmp])*pitchX + shift1X - shift2X
                        self.relY = (self.row[maxTOTindex1_tmp]+1.)*pitchY - shift1Y - shift2Y
                    else :
                    #pixel from pair 2 with highest energy on the right top corner
                        self.relX = (self.col[maxTOTindex1_tmp])*pitchX + shift1X + shift2X
                        self.relY = (self.row[maxTOTindex1_tmp]+1.)*pitchY - shift1Y + shift2Y
                elif(self.col[maxTOTindex1_tmp] < self.col[minTOTindex1_tmp] and self.row[maxTOTindex1_tmp] < self.row[minTOTindex1_tmp]) :
                #pixel from pair 1 with highest energy on the left bottom corner
                    if(self.col[maxTOTindex2_tmp] < self.col[minTOTindex2_tmp]) :
                        self.relX = (self.col[maxTOTindex1_tmp]+1.)*pitchX - shift1X - shift2X
                        self.relY = (self.row[maxTOTindex1_tmp]+1.)*pitchY - shift1Y + shift2Y
                    else :
                        self.relX = (self.col[maxTOTindex1_tmp]+1.)*pitchX - shift1X + shift2X
                        self.relY = (self.row[maxTOTindex1_tmp]+1.)*pitchY - shift1Y - shift2Y
                elif(self.col[maxTOTindex1_tmp] > self.col[minTOTindex1_tmp] and self.row[maxTOTindex1_tmp] > self.row[minTOTindex1_tmp]) :
                #pixel from pair 1 with highest energy on the right top corner
                    if(self.col[maxTOTindex2_tmp] < self.col[minTOTindex2_tmp]) :
                        self.relX = (self.col[maxTOTindex1_tmp])*pitchX + shift1X - shift2X
                        self.relY = (self.row[maxTOTindex1_tmp])*pitchY + shift1Y + shift2Y
                    else :
                        self.relX = (self.col[maxTOTindex1_tmp])*pitchX + shift1X + shift2X
                        self.relY = (self.row[maxTOTindex1_tmp])*pitchY + shift1Y - shift2Y

            else : # not 2x2 -> using the simple Qweighted centroid
                self.relX=0.
                self.relY=0.
                for index,tot_tmp in enumerate(self.tot) :
                    self.relX+=(self.col[index]*pitchX + pitchX/2.)*tot_tmp
                    self.relY+=(self.row[index]*pitchY + pitchY/2.)*tot_tmp
                self.relX/=self.totalTOT
                self.relY/=self.totalTOT


        elif(self.size==3) :
            if(self.sizeX==2 and self.sizeY==2) :
                # copy original row, col, tot to keep safe
                orig_row = list(self.row)
                orig_col = list(self.col)
                orig_tot = list(self.tot)

                # calculate and add missing pixel
                for i in xrange(len(self.row)):
                    if self.row.count(self.row[i]) == 1:
                        self.row.append(self.row[i])
                    if self.col.count(self.col[i]) == 1:
                        self.col.append(self.col[i])
                        self.tot.append(1.0)

                # proceed as 4 cluster case
                for i in xrange(len(self.row)):
                    if self.row[i] == min(self.row) and self.col[i] == min(self.col):
                        bottomlefti = i
                    if self.row[i] == max(self.row) and self.col[i] == max(self.col):
                        toprighti = i
                    if self.row[i] == min(self.row) and self.col[i] == max(self.col):
                        toplefti = i
                    if self.row[i] == max(self.row) and self.col[i] == min(self.col):
                        bottomrighti = i
                Qrel1 = self.tot[bottomlefti] / (self.tot[bottomlefti] + self.tot[toprighti])
                Qrel2 = self.tot[bottomrighti] / (self.tot[bottomrighti] + self.tot[toplefti])

                shift1X,shift1Y = shiftDiag(sigmaX,sigmaY,Qrel1)
                shift2X,shift2Y = shiftDiag(sigmaX,sigmaY,Qrel2)

                self.relX = max(self.col)*pitchX - shift1X - shift2X
                self.relY = max(self.row)*pitchY - shift1Y + shift2Y

                # put back original row, col, tmp
                self.row = orig_row
                self.col = orig_col
                self.tot = orig_tot

            else : # not 2x2 -> using the simple Qweighted centroid
                self.relX=0.
                self.relY=0.
                for index,tot_tmp in enumerate(self.tot) :
                    self.relX+=(self.col[index]*pitchX + pitchX/2.)*tot_tmp
                    self.relY+=(self.row[index]*pitchY + pitchY/2.)*tot_tmp
                self.relX/=self.totalTOT
                self.relY/=self.totalTOT


        else : #other cluster sizes->using the simple Qweighted centroid
            self.relX=0.
            self.relY=0.
            for index,tot_tmp in enumerate(self.tot) :
                self.relX+=(self.col[index]*pitchX + pitchX/2.)*tot_tmp
                self.relY+=(self.row[index]*pitchY + pitchY/2.)*tot_tmp

            if(self.totalTOT>0):
                self.relX/=self.totalTOT
                self.relY/=self.totalTOT
            else :
                self.GetDigitalCentroid()


        #print "relX : %f"%float(self.relX)
        #print "relY : %f"%float(self.relY)
        self.absX=self.relX  -npix_X*pitchX/2.
        self.absY=self.relY  -npix_Y*pitchY/2.
        self.absZ=0





    def GetResiduals(self,x,y) :
        self.resX = self.absX-(x)
        self.resY = self.absY-(y)


    def GetPixelResiduals(self,trackx,tracky) :
        # compute the x, y distances between a track 
        # and the centre of each pixel in the cluster
        # return the smallest combined

        dr = []
        dx = []
        dy = []
        for i in xrange(self.size):
            resX = self.col[i]*pitchX + pitchX/2. - halfChip_X - trackx
            resY = self.row[i]*pitchY + pitchY/2. - halfChip_Y - tracky

            dr.append(sqrt(resX**2 + resY**2))
            dx.append(resX)
            dy.append(resY)

        return min(dr), dx[dr.index(min(dr))], dy[dr.index(min(dr))]
