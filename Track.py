###############################################################################################################################
#
#        Class for the tracks and their properties
#
###############################################################################################################################

class Track :

    trackX=[]
    trackY=[]
    chi2 = []
    event=0
    ndof = []
    iden = []
    trackNum = []
    dxdz = []
    dydz = []
    cluster=-11

    def __init__(self):
        self.trackX=[]
        self.trackY=[]
        self.chi2= []
        self.event=0
        self.ndof= []
        self.iden=[]
        self.trackNum=[]

    def Fill(self,x,y,chi2,event,ndof,iden,trackNum):
        self.trackX=x
        self.trackY=y
        self.chi2=chi2
        self.event=event
        self.ndof=ndof
        self.iden=iden
        self.trackNum=trackNum

    def FindCluster(self,clusters) :
        resX_tmp = []
        resY_tmp = []

    def Print(self):

        print "##### Track #####"
        for i,x in enumerate(self.trackX) :
            print "Track X : %f Track Y : %f iden : %i cluster = %i"%(self.trackX[i],self.trackY[i],self.iden[i],self.cluster)
        print "#################"
