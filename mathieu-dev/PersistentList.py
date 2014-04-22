import shelve
import os

class PersistentList :

    # A Class acting as a list, but persistified on Disk to run on large data set whithout saturating RAM

    database  = 0
    index = 0
    sync_count = 0
    n_items = 0

    theBuffer = []
    buffer_length = 0
    def __init__(self):
        pass
    def __init__(self,name,sync_count_max=10000,buffer_length=10000):
        os.system("rm -fr %s"%name)
        self.database = shelve.open(name, writeback=True)
        self.index = 0
        self.maxsync=sync_count_max
        self.buffer_length = buffer_length
        self.theBuffer = {}

    def __iter__(self):
        return self

    def __getitem__(self,i) :
        #self.database.sync()
        #print "Length of DB %i"%len(self.database)
        self.sync_count+=1
        if self.sync_count == self.maxsync :
            self.database.sync()
            self.sync_count=0
        return self.database['%i'%i]

    def __setitem__(self,i,value) :
        if "%i"%i not in self.database :
            self.n_items+=1
        self.database["%i"%i]= value
        #self.sync_count+=1
#               if self.sync_count == self.maxsync :
#                       self.database.sync()
#                       self.sync_count=0

    def __len__(self):
        return len(self.database.keys())


    def next(self):

        if(self.index%self.buffer_length==0):
            for key in self.theBuffer :
                self.database[key]=self.theBuffer[key]

            self.theBuffer = {}
            for i in range(self.index,self.index+self.buffer_length):
                try :
                    self.theBuffer["%i"%i]=self.database["%i"%i]
                except :
                    pass

        if self.index >= len(self.database):
            self.index=0
            raise StopIteration
        self.index = self.index + 1
        return self.theBuffer["%i"%(self.index-1)]

    def append(self,value) :
        self.database['%i'%(self.index)]=value
        self.n_items+=1
        self.index+=1
        self.sync_count+=1
        if self.sync_count == self.maxsync :
            self.database.sync()
            self.sync_count=0

    def pop(self,i) :
        del self.database["%i"%i]
        self.n_items-=1
        self.sync_count+=1
        if self.sync_count == self.maxsync :
            self.database.sync()
            self.sync_count=0

#       def insert(self,i,value) :
#               if "%i"%i not in self.database :
#                       self.n_items+=1
#               self.database["%i"%i]=value
##              self.sync_count+=1
#               if self.sync_count == self.maxsync :
#                       self.database.sync()
#                       self.sync_count=0

    def sclose(self) :
        self.database.close()
