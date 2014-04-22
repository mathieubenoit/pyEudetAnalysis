import os
from os import listdir
from os import walk
import time,os.path
import re
from optparse import OptionParser

parser = OptionParser()
parser.add_option("-s", "--step",
                  help="step (Reco of Analysis)", dest="STEP")
(options, args) = parser.parse_args()

if(options.STEP) :
    if(options.STEP=="Reco"):
        print "Launching Batch for Reconstruction"
    elif(options.STEP=="Analysis"):
        print "Launching Batch for Analysis"   
else :
    print "Please provide a step to perform (Reco or Analysis)"
    parser.print_help()
    exit()


def FetchRuns(folder):
    runs = []
    for (dirpath, dirnames, filenames) in walk(folder):
        
        for file in filenames :
            if (file.find("tbtrack")>=0):
                print "found file %s"%file
                tmp = int(''.join(x for x in file if x.isdigit()))
                runs.append(tmp)
        break
    return runs
 
 
 
 # Important folders for the analysis (input, output, plots)

data_folder = "/afs/cern.ch/eng/clic/data/DESY_TB_DATA_October2013_results/tbtrack"
result_folder = "/afs/cern.ch/eng/clic/data/DESY_TB_DATA_October2013_results/pyEudetAnalysis_plots"
LogFolder = "/afs/cern.ch/eng/clic/data/DESY_TB_DATA_October2013_results/pyEudetAnalysis_plots/launch/Logs"
pyEudetFolder = "/afs/cern.ch/eng/clic/data/DESY_TB_DATA_October2013_results/pyEudetNtuples"
Alignement_file =  ""
 
 
 
 
# Establish list of processable runs in data folder folder   
runs=FetchRuns(data_folder)


# You Can filter the runs to process like this 
runs = [x for x in runs if ( (x in range(1005,1049))]

## 1) Generate range for your data set 
#run_B04_W0110_90deg_bias_scan = range(48,67)
#run_A06_W0110_90deg_bias_scan = range(96,154)
#run_L05_W0125_90deg_bias_scan = range(155,184)
#
##2) Filter only to keep run in the range you want
#runs = [x for x in runs if ((x in run_B04_W0110_90deg_bias_scan ) or (x in run_A06_W0110_90deg_bias_scan) or (x in run_L05_W0125_90deg_bias_scan) )]

print "launching batch for Runs : "
print runs


#Queue Selection 
queue = "8nh"
if(options.STEP=="Reco"):
    queue = "8nh"
elif(options.STEP=="Analysis"):
    queue = "1nh"
    
print "using queue %s"%queue







# Create space for your batch files and logs 
os.system("mkdir %s/launch"%result_folder)
launch_folder = "%s/launch"%result_folder
print "Writing Batch Files to %s"%launch_folder





# Build launch scripts 

logs = []
batch = []

if(options.STEP=="Reco"):
    method_names = ["QWeighted"]
    for method in method_names :
        for run in runs :
            filename="%s/Batch_%s_Run%i.sh"%(launch_folder,method,run)
            subLogFolder = "%s/Run%i_%s"%(LogFolder,run,method)
    
            f=open(filename,'w')
            f.write("cd /afs/cern.ch/user/m/mbenoit/workspace/pyEudetAna \n")
            f.write("source setup_CERN.sh \n")
            
            
            if run in range(48,67) :
                Alignement_file =  "/afs/cern.ch/eng/clic/data/DESY_TB_DATA_19-30-08-2013_results/pyEudetAnalysis_plots/Alignment/Alignment_Run48.txt" 
                f.write("python pyEudetReconstructionOnly.py -r %i -m %s -d %s -o %s -a %s -e %f\n"%(run,method,data_folder,pyEudetFolder,Alignement_file,0.05))
                batch.append(filename)
                logs.append(subLogFolder)
            elif run in range(96,154) : 
                Alignement_file =  "/afs/cern.ch/eng/clic/data/DESY_TB_DATA_19-30-08-2013_results/pyEudetAnalysis_plots/Alignment/Alignment_Run130.txt" 
                f.write("python pyEudetReconstructionOnly.py -r %i -m %s -d %s -o %s -a %s -e %f\n"%(run,method,data_folder,pyEudetFolder,Alignement_file,0.02))       
                batch.append(filename)
                logs.append(subLogFolder)
            elif run in range(155,184) : 
                Alignement_file =  "/afs/cern.ch/eng/clic/data/DESY_TB_DATA_19-30-08-2013_results/pyEudetAnalysis_plots/Alignment/Alignment_Run170.txt" 
                f.write("python pyEudetReconstructionOnly.py -r %i -m %s -d %s -o %s -a %s -e %f\n"%(run,method,data_folder,pyEudetFolder,Alignement_file,0.05))       
                batch.append(filename) 
                logs.append(subLogFolder)
            else : 
                print "No alignment for run %i"%run    
            
            
            os.system("chmod u+rwx %s"%filename)

if(options.STEP=="Analysis"):
    method_names = ["QWeighted","maxTOT","DigitalCentroid","EtaCorrection"]
    for method in method_names :
        for run in runs :
            filename="%s/Batch_%s_Run%i.sh"%(launch_folder,method,run)
            subLogFolder = "%s/Run%i_%s"%(LogFolder,run,method)
    
            f=open(filename,'w')
            f.write("cd /afs/cern.ch/user/m/mbenoit/workspace/pyEudetAna \n")
            f.write("source setup_CERN.sh \n")
            
            
            if run in range(48,67) :
                Alignement_file =  "/afs/cern.ch/eng/clic/data/DESY_TB_DATA_19-30-08-2013_results/pyEudetAnalysis_plots/Alignment/Alignment_Run48.txt" 
                f.write("python pyEudetAnalysisOnly.py -r %i -m %s -d %s -o %s -a %s -e %f \n"%(run,method,pyEudetFolder,result_folder,Alignement_file,0.1))
                batch.append(filename)
                logs.append(subLogFolder)
            elif run in range(96,154) : 
                Alignement_file =  "/afs/cern.ch/eng/clic/data/DESY_TB_DATA_19-30-08-2013_results/pyEudetAnalysis_plots/Alignment/Alignment_Run130.txt" 
                f.write("python pyEudetAnalysisOnly.py -r %i -m %s -d %s -o %s -a %s -e %f\n"%(run,method,pyEudetFolder,result_folder,Alignement_file,0.1))       
                batch.append(filename)
                logs.append(subLogFolder)
            elif run in range(155,184) : 
                Alignement_file =  "/afs/cern.ch/eng/clic/data/DESY_TB_DATA_19-30-08-2013_results/pyEudetAnalysis_plots/Alignment/Alignment_Run170.txt" 
                f.write("python pyEudetAnalysisOnly.py -r %i -m %s -d %s -o %s -a %s -e %f\n"%(run,method,pyEudetFolder,result_folder,Alignement_file,0.1))       
                batch.append(filename) 
                logs.append(subLogFolder)
            else : 
                print "No alignment for run %i"%run    
            
            
            os.system("chmod u+rwx %s"%filename)
            
        


# Launch the batch to lxbatch 

for i,job in enumerate(batch) :
    print "Launching %s, STDOUT at %s"%(job,logs[i])
    os.system("mkdir %s"%logs[i])
    os.system("bsub -o %s/STDOUT -q %s %s"%(logs[i],queue,job))
