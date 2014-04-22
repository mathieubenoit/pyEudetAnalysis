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
    if(options.STEP=="Alignment"):
        print "Launching Batch for Alignment"
    elif(options.STEP=="Reco"):
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

data_folder = "/afs/cern.ch/eng/clic/TBData/DESY_TB_DATA_August2013_results/tbtrack"
result_folder = "/afs/cern.ch/eng/clic/TBData/DESY_TB_DATA_August2013_results/pyEudetAnalysis_plots"
LogFolder = "/afs/cern.ch/eng/clic/TBData/DESY_TB_DATA_August2013_results/pyEudetAnalysis_plots/launch/Logs"
pyEudetFolder = "/afs/cern.ch/eng/clic/TBData/DESY_TB_DATA_August2013_results/pyEudetNtuples"
Alignement_file =  ""




# Establish list of processable runs in data folder folder
#runs=FetchRuns(data_folder)
#runs=range(1,2)


#runs = range(2032,2043)
runs = range(48,52) + range(97,99) + range(150,154) + range(314,333) +range(353,388) +range(390,399)+[465,466,473,474]

# You Can filter the runs to process like this
#runs = [x for x in runs if ( (x in range(1005,1049))]

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
if(options.STEP=="Alignment"):
    queue = "1nd"
elif(options.STEP=="Reco"):
    queue = "8nh"
elif(options.STEP=="Analysis"):
    queue = "1nd"

print "using queue %s"%queue







# Create space for your batch files and logs
os.system("mkdir %s/launch"%result_folder)
launch_folder = "%s/launch"%result_folder
print "Writing Batch Files to %s"%launch_folder





# Build launch scripts

logs = []
batch = []

if(options.STEP=="Alignment"):
    method_names = ["QWeighted"]
    for method in method_names :
        for run in runs :
            filename="%s/Batch_%s_Run%i.sh"%(launch_folder,method,run)
            subLogFolder = "%s/Run%i_%s"%(LogFolder,run,method)

            f=open(filename,'w')
            f.write("cd /afs/cern.ch/user/m/mbenoit/workspace/pyEudetAna \n")
            f.write("source setup_CERN.sh \n")

            Alignement_file =  "/afs/cern.ch/eng/clic/TBData/DESY_TB_DATA_August2013_results/Alignment/Alignment_Run%i.txt"%run
            f.write("python ComputeAlignment.py -r %i -m %s -d %s -o %s -a %s -e %f -n %i\n"%(run,method,data_folder,pyEudetFolder,Alignement_file,0.05,15000))
            batch.append(filename)
            logs.append(subLogFolder)
            os.system("chmod u+rwx %s"%filename)

elif(options.STEP=="Reco"):
    method_names = ["QWeighted"]
    for method in method_names :
        for run in runs :
            filename="%s/Batch_%s_Run%i.sh"%(launch_folder,method,run)
            subLogFolder = "%s/Run%i_%s"%(LogFolder,run,method)

            f=open(filename,'w')
            f.write("cd /afs/cern.ch/user/m/mbenoit/workspace/pyEudetAna \n")
            f.write("source setup_CERN.sh \n")


            Alignement_file =  "/afs/cern.ch/eng/clic/TBData/DESY_TB_DATA_August2013_results/Alignment/Alignment_Run%i.txt"%run
            f.write("python pyEudetReconstructionOnly.py -r %i -m %s -d %s -o %s -a %s -e %f\n"%(run,method,data_folder,pyEudetFolder,Alignement_file,0.05))
            batch.append(filename)
            logs.append(subLogFolder)


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


            Alignement_file =  "/afs/cern.ch/eng/clic/TBData/DESY_TB_DATA_August2013_results/Alignment/Alignment_Run%i.txt"%run	    
            Alignement_file =  "/afs/cern.ch/eng/clic/TBData/DESY_TB_DATA_August2013_results/Alignment/Alignment_Run3208.txt"
            f.write("python pyEudetAnalysisOnly.py -r %i -m %s -d %s -o %s -a %s -e %f \n"%(run,method,pyEudetFolder,result_folder,Alignement_file,0.1))
            batch.append(filename)
            logs.append(subLogFolder)


            os.system("chmod u+rwx %s"%filename)




# Launch the batch to lxbatch

for i,job in enumerate(batch) :
    print "Launching %s, STDOUT at %s"%(job,logs[i])
    os.system("mkdir %s"%logs[i])
    os.system("bsub -o %s/STDOUT -q %s %s"%(logs[i],queue,job))
