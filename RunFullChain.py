import time,os
from optparse import OptionParser

parser = OptionParser()

parser.add_option("-b", "--begin",
                  help="First Run", dest="BEGIN", type="int")

parser.add_option("-f", "--end",
                  help="Last Run", dest="END", type="int")
		  
parser.add_option("-n", "--nevent",
                  help="Number of events to process", dest="NEVENT")

parser.add_option("-m", "--method",
                  help="Position Reconstruction Method, QWeighted, DigitalCentroid, maxTOT, EtaCorrection", dest="METHOD", default="QWeighted")

parser.add_option("-d", "--data",
                  help="pyEudetNtuple Input Folder", dest="INPUT")

parser.add_option("-o", "--output",
                  help="Histograms and results output folder", dest="OUTPUT", default=".")

parser.add_option("-a", "--alignment",
                  help="alignement file", dest="ALIGNMENT", default="alignement.dat")

parser.add_option("-e", "--edge",
                  help="edge width", dest="EDGE", default=0.0, type="float")


parser.add_option("-t", "--target",
                  help="target : alignmnent, reconstruction, reconstruction+analysis, analysis, all", dest="TARGET")

(options, args) = parser.parse_args()

     
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
    elif(options.METHOD=="all"):
        method_name=options.METHOD
    else :
        print "Please provide a valid cluster position reconstruction method ( -m [method]  QWeighted,maxTOT,DigitalControid,EtaCorrection)"
        parser.print_help()
        exit()

else:
    print "Please provide a valid cluster position reconstruction method ( -m [method]  QWeighted,maxTOT,DigitalControid,EtaCorrection)"
    parser.print_help()
    exit()

if(options.INPUT):
    input_folder=options.INPUT
else :
    print "Please provide an input folder with pyEudetNtuple files (-d [PathToData] , put no / at the end )"
    parser.print_help()
    exit()

if(options.OUTPUT):
    PlotPath=options.OUTPUT
else :
    print "Please provide an output folder (-o [PathToOutput] , put no / at the end )"
    parser.print_help()
    exit()

if(options.ALIGNMENT):
    AlignementPath = "%s"%(options.ALIGNMENT)
else :
    print "Please provide an Alignment File (-a [PathToFile]  0 0 0 0 0 if no alignement needed )"
    parser.print_help()
    exit()


if (options.BEGIN) : 
	begin=int(options.BEGIN)
else : 
	print "provide a first run!"
	exit()
		
			
if (options.END) : 
	end=int(options.END)
else : 
	print "provide a last run!"
	exit()
		
	
if (options.TARGET) : 
	target=options.TARGET
else : 
	print "provide a target!"
	exit()
			

for run in range(begin,end+1) : 
	
	AlignementFile="%s/Alignment_run%i.txt"%(AlignementPath,run) 
	
	if (target=="alignment" or target == "all") :
		os.system("python ComputeAlignment.py -r %i -d %s -o %s -a %s -e %f -m %s -n 5000"%(run,input_folder,PlotPath,AlignementFile,edge_width,"QWeighted"))
	
	
	if ("reconstruction" in target or target == "all" ) :
	
		if method_name=="all" :
			os.system("python pyEudetReconstructionOnly.py -r %i -d %s -o %s -a %s -e %f -m %s"%(run,input_folder,PlotPath,AlignementFile,edge_width,"QWeighted"))
			#os.system("python pyEudetAnalysisOnly.py -r %i -d %s -o %s -a %s -e %f -m %s"%(run,input_folder,PlotPath,AlignementFile,edge_width,"QWeighted"))
			os.system("python pyEudetReconstructionOnly.py -r %i -d %s -o %s -a %s -e %f -m %s"%(run,input_folder,PlotPath,AlignementFile,edge_width,"EtaCorrection"))
			#os.system("python pyEudetAnalysisOnly.py -r %i -d %s -o %s -a %s -e %f -m %s "%(run,input_folder,PlotPath,AlignementFile,edge_width,"EtaCorrection"))
			os.system("python pyEudetReconstructionOnly.py -r %i -d %s -o %s -a %s -e %f -m %s"%(run,input_folder,PlotPath,AlignementFile,edge_width,"DigitalCentroid"))
			#os.system("python pyEudetAnalysisOnly.py -r %i -d %s -o %s -a %s -e %f -m %s "%(run,input_folder,PlotPath,AlignementFile,edge_width,"DigitalCentroid"))
		else : 
			os.system("python pyEudetReconstructionOnly.py -r %i -d %s -o %s -a %s -e %f -m %s "%(run,input_folder,PlotPath,AlignementFile,edge_width,method_name))
			#os.system("python pyEudetAnalysisOnly.py -r %i -d %s -o %s -a %s -e %f -m %s "%(run,input_folder,PlotPath,AlignementFile,edge_width,method_name))
	
	if ("analysis" in target or target == "all" ) :

	
		if method_name=="all" :
			#os.system("python pyEudetReconstructionOnly.py -r %i -d %s -o %s -a %s -e %f -m %s"%(run,input_folder,PlotPath,AlignementFile,edge_width,"QWeighted"))
			os.system("python pyEudetAnalysisOnly.py -r %i -d %s -o %s -a %s -e %f -m %s"%(run,input_folder,PlotPath,AlignementFile,edge_width,"QWeighted"))
			#os.system("python pyEudetReconstructionOnly.py -r %i -d %s -o %s -a %s -e %f -m %s"%(run,input_folder,PlotPath,AlignementFile,edge_width,"EtaCorrection"))
			os.system("python pyEudetAnalysisOnly.py -r %i -d %s -o %s -a %s -e %f -m %s "%(run,input_folder,PlotPath,AlignementFile,edge_width,"EtaCorrection"))
			#os.system("python pyEudetReconstructionOnly.py -r %i -d %s -o %s -a %s -e %f -m %s"%(run,input_folder,PlotPath,AlignementFile,edge_width,"DigitalCentroid"))
			os.system("python pyEudetAnalysisOnly.py -r %i -d %s -o %s -a %s -e %f -m %s "%(run,input_folder,PlotPath,AlignementFile,edge_width,"DigitalCentroid"))
		else : 
			#os.system("python pyEudetReconstructionOnly.py -r %i -d %s -o %s -a %s -e %f -m %s "%(run,input_folder,PlotPath,AlignementFile,edge_width,method_name))
			os.system("python pyEudetAnalysisOnly.py -r %i -d %s -o %s -a %s -e %f -m %s "%(run,input_folder,PlotPath,AlignementFile,edge_width,method_name))




