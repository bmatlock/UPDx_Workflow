 #!/bin/bash


#Use getopts code to establish arguments that the user will need to input

while getopts ":l:f:d:p:n:h" flag;
do
    case "${flag}" in 
    p) Pflag="${OPTARG}";;
    f) Fflag="${OPTARG}";;
    d) Dflag="${OPTARG}";;

    h) 
        echo "Certain Inputs are needed by the user in order to run this workflow"
    echo "          -p(Essential)           specify the path to where the workflow CDC-Complete-UPDx-Workflow script is located"  
    echo "          -f(Essential)           specify the path to where the fastq files are located"
    echo "          -d(Essential)           specify the analysis date in 'mmddyyyyy' format (just numbers no spaces; underscores; hyphens; or dashes)"         
    exit 1
    ;;

    \?)
        echo "Not a valid argument: -$OPTARG. WORKFLOW HALTED"
        exit 1
    ;;

    esac
done

#Send back error if the user did not input any arguments

while test $# -lt 1; do

echo "ERROR You have not supplied any argmuments"
echo "use -h argument to get help on what the arguments are and which are necessary"

exit
done

#Establish the arguments 

if [ -z "$Pflag" ]; then
    echo "ERROR: Please specify the full path to the folder where the UPDx Workflow Script is located using '-p' to run this workflow. WORKFLOW HALTED"
    exit 1
    elif [[ $Fflag = *"="* ]] || [[ $Fflag = *":"* ]] || [[ $Fflag = *";"* ]]|| [[ $Fflag = *","* ]] || [[ $Fflag = *"--"* ]]; then
    echo "ERROR: Please do not provide an '=' sign, colon ':' or any other separator in front of your p argument. WORKFLOW HALTED"
    exit 1
    fi

if [ -z "$Fflag" ]; then
    echo "ERROR: Please specify the full path to the folder where the fastq files are located using '-f' to run this workflow. WORKFLOW HALTED"
    exit 1
    elif [[ $Fflag = *"="* ]] || [[ $Fflag = *":"* ]] || [[ $Fflag = *";"* ]]|| [[ $Fflag = *","* ]] ||  [[ $Fflag = *"--"* ]]; then
    echo "ERROR: Please do not provide an '=' sign, colon ':' or any other separator in front of your f arguments. WORKFLOW HALTED"
    exit 1
    fi

if [ ! -d "$Fflag" ]; then
    echo "ERROR: the directory $Fflag does not exist. WORKFLOW HALTED"
    exit 1
fi

if [ -z "$Dflag" ]; then
    echo "ERROR: the date for this experiment has not been input using -d and is necessary for workflow to continue. WORKFLOW HALTED"
    exit 1
    elif [[ $Dflag = *"="* ]] || [[ $Dflag = *":"* ]] || [[ $Dflag = *";"* ]]|| [[ $Dflag = *","* ]] ||  [[ $Dflag = *"--"* ]]; then
    echo "ERROR: Please do not provide an '=' sign, colon ':' or any other separator in front of your d argument. WORKFLOW HALTED"
    exit 1
    elif [[ $Dflag != [[:digit:]]* ]]; then
    echo "ERROR: Please provide a numeric value for argument '-d' in mmddyyyy format. WORKFLOW HALTED"
    exit 1
fi

my_date_and_timestamp=$(date +"%m-%d-%y-%H_%M_%S")


my_timestamp=$(date +"%H_%M_%S")

input_reads=$Fflag

length_of_specimen_names=$Lflag

length_of_negative_names=$Nflag

date=$Dflag-$my_timestamp

PTW=$Pflag

mkdir Experiment_$my_date_and_timestamp


PTE=$PTW/Experiment_$my_date_and_timestamp

PTR=$PTW/18S_Results

py_scripts=$PTW/Python_Scripts 

#Establish Variables for Database and Primers

ref_database=$PTW/UPDx_Variables/UPDx_18s_Reference_Database/Ad_UPDx_18s_Blood_and_Tissue_Reference_Database_-_BM20230912.fasta
ref_database_for_tree=$PTW/UPDx_Variables/UPDx_18s_Reference_Database/Ad_UPDx_18s_Blood_and_Tissue_Reference_Database_For_Tree_-_BM20230912.fasta
exc_database=$PTW/UPDx_Variables/UPDx_18s_Exclusion_Database/UPDx_EXCLUSION_version_11Unpaired.fasta
f_primer=$PTW/UPDx_Variables/UPDx_F_primer/UPDx_18S_F_Primer.fasta
r_primer=$PTW/UPDx_Variables/UPDx_R_primer/UPDx_18S_R_Primer.fasta

#Create folders to store fastq reads for the experiment that need to be analyzed

cd $PTW || exit

if [ -d FQC_$date ]; then
    echo "The Previous FQC Folder is still present and is being deleted"
    rm -rf FQC_$date
fi

if [ -d Negative_Reads_$date ]; then
  echo "Warning previous Negative Reads folder is still present and is being deleted"
  rm -rf Negative_Reads_$date
fi


mkdir FQC_$date

cd FQC_$date || exit

mkdir R1_Files_$date
R1_file_location=$PTW/FQC_$date/R1_Files_$date

mkdir R2_Files_$date
R2_file_location=$PTW/FQC_$date/R2_Files_$date

cd $PTW/FQC_$date || exit

mkdir $PTW/FQC_$date/Negative_Reads_$date

mkdir $PTW/FQC_$date/Negative_Reads_$date/R1_Files
Neg_R1_file_location=$PTW/FQC_$date/Negative_Reads_$date/R1_Files

mkdir $PTW/FQC_$date/Negative_Reads_$date/R2_Files
Neg_R2_file_location=$PTW/FQC_$date/Negative_Reads_$date/R2_Files

#Copy fastq reads from path input by user to our newly created folders

cd $input_reads || exit

cp *_R1_* $R1_file_location
cp *_R2_* $R2_file_location

#Copy negative control fastq reads to separate folders to be processed separately
#Need to change from any letters to whatever the decide on naming the negative reads

cd $R1_file_location || exit
mv Neg*.fastq* $Neg_R1_file_location

cd $R2_file_location || exit
mv Neg*.fastq* $Neg_R2_file_location


#Rename fastq files to replace hyphen with underscore and use length input by user to leave just the sample name and the read number

cd $R1_file_location || exit

for file in *
    do 
    mv "$file" `echo "$file" | sed -e 's/-/_/g' | sed -e "s/\(.*\)_\(R[0-9]\).*.\(.\{3\}\)/\1.\2.fastq\3/"` 
done
   
cd $R2_file_location || exit

for file in *
    do 
    mv "$file" `echo "$file" | sed -e 's/-/_/g' | sed -e "s/\(.*\)_\(R[0-9]\).*.\(.\{3\}\)/\1.\2.fastq\3/"` 
done

PTF=$PTW/FQC_$date
cd $PTF

#FORNEGATIVES

cd $PTF/Negative_Reads_$date/R1_Files || exit

for file in *
  do
   mv "$file" `echo "$file" | sed -e 's/-/_/g' | sed -e "s/\(.*\)_\(R[0-9]\).*.\(.\{3\}\)/\1.\2.fastq\3/"` 
done

cd $PTF/Negative_Reads_$date/R2_Files || exit

for file in *
  do
   mv "$file" `echo "$file" | sed -e 's/-/_/g' | sed -e "s/\(.*\)_\(R[0-9]\).*.\(.\{3\}\)/\1.\2.fastq\3/"`  
done

###Start Sequence Processing the reads#####


#Activate the conda environment with the packages needed to perform the processing steps
#Going to be different based on the path to the users conda environment


# echo "Activating conda environment"
# eval "$($(which conda) 'shell.bash' 'hook')"
module load miniconda3
#ENVS=$(conda env list | awk '{print$1}')

# if [[ $ENVS = *updx_environment* ]]; then
# 	echo "UPDx environment exists and is being activated"
#   source activate $PTW/updx_environment
# else
# 	echo "UPDx environment does not exist. It will be created now"
# 	conda env create --prefix=$PTW/updx_environment --file=$PTW/linux_updx_environment.yml
#   echo "UPDx environment has been created and is being activated"
#   source activate $PTW/updx_environment
#   echo "Necessary pip packages now being installed"
#   pip install argparse
#   pip install cdhit-reader
#   pip install cutadapt
#   pip install pandas
#   pip install regex
#   echo "Pip packages are installed and environment is ready"
# fi;

source activate $PTW/updx_environment

echo "Now entered into Conda environment"

cd $PTF || exit

if [ -d failed_sequence_runs ]; then
    echo "Warning previous failed runs folder is still present and being deleted!"
    rm -rf failed_sequence_runs
fi

if [ -d tempfolder ]; then
    echo "Warning previous tempfolder is still present and being deleted!"
    rm -rf tempfolder
fi

if [ -d fqc_output ]; then
    echo "Warning previous fqc_output is still present and being deleted!"
    rm -rf fqc_output
fi

if [ ! -z shortnames.txt ]; then
    echo "Warning previous shortnames file is still present and being deleted!"
    rm -rf shortnames.txt
fi

if [ -d post_trimming_calculations ]; then
  echo "Warning previous Post-Trimming calculations is still present and is being deleted"
  rm -rf post_trimming_calculations
fi


#make directories needed for holding the reads that passed the quality check and for displaying the specimens that failed
mkdir $PTF/tempFolder
temp=$PTF/tempFolder
mkdir $PTF/failed_sequence_runs

#FORNEGATIVES
mkdir $PTF/Negative_Reads_$date/tempFolder
ntemp=$PTF/Negative_Reads_$date/tempFolder
mkdir $PTF/Negative_Reads_$date/failed_negative_runs

#Create a text file that holds the names of all the specimen

ls $PTF/R1_Files_$date | awk -F "/" '{print$NF}' | awk -F "R" '{print$1}'  >> $PTF/shortnames.txt

#FORNEGATIVES

ls $PTF/Negative_Reads_$date/R1_Files | awk -F "/" '{print$NF}' | awk -F "R" '{print$1}' >> $PTF/Negative_Reads_$date/negnames.txt 

#Create all the fqc_output folders

mkdir -p $PTF/fqc_output/raw_seqs_fqc_results
mkdir -p $PTF/fqc_output/trimmed_seqs_fqc_results
mkdir -p $PTF/fqc_output/BBMerge_Output
mkdir -p $PTF/fqc_output/BBMerge_Output/trimmed_merged_reads
mkdir -p $PTF/fqc_output/BBMerge_Output/trimmed_unmerged_reads
mkdir -p $PTF/fqc_output/phiX_output
mkdir -p $PTF/fqc_output/bbduk_output
mkdir -p $PTF/fqc_output/BBMerge_Output

#FORNEGATIVES
mkdir -p $PTF/Negative_Reads_$date/fqc_output/raw_seqs_fqc_results
mkdir -p $PTF/Negative_Reads_$date/fqc_output/trimmed_seqs_fqc_results
mkdir -p $PTF/fqc_output/BBMerge_Output
mkdir -p $PTF/Negative_Reads_$date/fqc_output/BBMerge_Output/trimmed_merged_reads
mkdir -p $PTF/Negative_Reads_$date/fqc_output/BBMerge_Output/trimmed_unmerged_reads
mkdir -p $PTF/Negative_Reads_$date/fqc_output/phiX_output
mkdir -p $PTF/Negative_Reads_$date/fqc_output/bbduk_output

#Path to Bowtie2 tools; needs to be saved in the same folder as the UPDx Script (This might not be necessary because of the conda environment)
PTB2=$PTW/PhiX_Illumina_RTA/PhiX/Illumina/RTA/Sequence/Bowtie2Index/genome
PTBBD=$PTW/BBTools

#Copy your files into a temporary folder for pre-trimming quality check

cat $PTF/shortnames.txt | while read specimen_name
do
    cp -a $PTF/R1_Files_$date/. $temp
    cp -a $PTF/R2_Files_$date/. $temp
    fastqc $PTF/tempFolder/$specimen_name*.gz -o $PTF/fqc_output/raw_seqs_fqc_results
    bowtie2 -p 5 -x $PTB2 -1 $PTF/tempFolder/$specimen_name\R1.fastq.gz -2 $PTF/tempFolder/$specimen_name\R2.fastq.gz -S $PTF/fqc_output/phiX_output/$specimen_name.phix.sam &> $PTF/fqc_output/phiX_output/$specimen_name.log
done

#FORNNEGATIVES
cat $PTF/Negative_Reads_$date/negnames.txt | while read negative_name
do
    cp -a $PTF/Negative_Reads_$date/R1_Files/. $ntemp
    cp -a $PTF/Negative_Reads_$date/R2_Files/. $ntemp
    fastqc $PTF/Negative_Reads_$date/tempFolder/$negative_name*.gz -o $PTF/Negative_Reads_$date/fqc_output/raw_seqs_fqc_results
    bowtie2 -p 5 -x $PTB2 -1 $PTF/Negative_Reads_$date/tempFolder/$negative_name\R1.fastq.gz -2 $PTF/Negative_Reads_$date/tempFolder/$negative_name\R2.fastq.gz -S $PTF/Negative_Reads_$date/fqc_output/phiX_output/$negative_name.phix.sam &> $PTF/Negative_Reads_$date/fqc_output/phiX_output/$negative_name.log
done

#Check the fastqc results and create a text file that contains the main pre-trimming quality check
cd $PTF/fqc_output/raw_seqs_fqc_results
for x in *.zip; do unzip "$x"; done 

cat $PTF/shortnames.txt | while read specimen_name;
do
    echo "Average Quality" `sed -n '/^#Base.*Mean/,/^>>END_MODULE/p' $specimen_name*R1*/fastqc_data.txt | grep -v "#Base" | grep -v ">>END_MODULE" | awk '{total += $2; count++} END {print total/count}'` > $specimen_name.rawR1_qualityCheck.txt
    echo "Average Quality" `sed -n '/^#Base.*Mean/,/^>>END_MODULE/p' $specimen_name*R2*/fastqc_data.txt | grep -v "#Base" | grep -v ">>END_MODULE" | awk '{total += $2; count++} END {print total/count}'` > $specimen_name.rawR2_qualityCheck.txt
    grep "Total Sequences" $specimen_name*R1*/fastqc_data.txt >> $specimen_name.rawR1_qualityCheck.txt
    grep "Total Sequences" $specimen_name*R2*/fastqc_data.txt >> $specimen_name.rawR2_qualityCheck.txt
    grep "Sequence length" $specimen_name*R1*/fastqc_data.txt >> $specimen_name.rawR1_qualityCheck.txt
    grep "Sequence length" $specimen_name*R2*/fastqc_data.txt >> $specimen_name.rawR2_qualityCheck.txt
    grep ">>Basic Statistics" $specimen_name*R1*/fastqc_data.txt >> $specimen_name.rawR1_qualityCheck.txt
    grep ">>Basic Statistics" $specimen_name*R2*/fastqc_data.txt >> $specimen_name.rawR2_qualityCheck.txt
    grep ">>Per base sequence quality" $specimen_name*R1*/fastqc_data.txt >> $specimen_name.rawR1_qualityCheck.txt
    grep ">>Per base sequence quality" $specimen_name*R2*/fastqc_data.txt >> $specimen_name.rawR2_qualityCheck.txt
    grep ">>Per sequence quality scores" $specimen_name*R1*/fastqc_data.txt >> $specimen_name.rawR1_qualityCheck.txt
    grep ">>Per sequence quality scores" $specimen_name*R2*/fastqc_data.txt >> $specimen_name.rawR2_qualityCheck.txt
    grep ">>Per base sequence content" $specimen_name*R1*/fastqc_data.txt >> $specimen_name.rawR1_qualityCheck.txt
    grep ">>Per base sequence content" $specimen_name*R2*/fastqc_data.txt >> $specimen_name.rawR2_qualityCheck.txt
    grep ">>Per sequence GC content" $specimen_name*R1*/fastqc_data.txt >> $specimen_name.rawR1_qualityCheck.txt
    grep ">>Per sequence GC content" $specimen_name*R2*/fastqc_data.txt >> $specimen_name.rawR2_qualityCheck.txt
    grep ">>Sequence Duplication Levels" $specimen_name*R1*/fastqc_data.txt >> $specimen_name.rawR1_qualityCheck.txt
    grep ">>Sequence Duplication Levels" $specimen_name*R2*/fastqc_data.txt >> $specimen_name.rawR2_qualityCheck.txt
    grep ">>Overrepresented sequences" $specimen_name*R1*/fastqc_data.txt >> $specimen_name.rawR1_qualityCheck.txt
    grep ">>Overrepresented sequences" $specimen_name*R2*/fastqc_data.txt >> $specimen_name.rawR2_qualityCheck.txt
done


#FORNEGATIVES
cd $PTF/Negative_Reads_$date/fqc_output/raw_seqs_fqc_results
for x in *.zip; do unzip "$x"; done

cat $PTF/Negative_Reads_$date/negnames.txt | while read specimen_name;
do
    echo "Average Quality" `sed -n '/^#Base.*Mean/,/^>>END_MODULE/p' $specimen_name*R1*/fastqc_data.txt | grep -v "#Base" | grep -v ">>END_MODULE" | awk '{total += $2; count++} END {print total/count}'` > $specimen_name.rawR1_qualityCheck.txt
    echo "Average Quality" `sed -n '/^#Base.*Mean/,/^>>END_MODULE/p' $specimen_name*R2*/fastqc_data.txt | grep -v "#Base" | grep -v ">>END_MODULE" | awk '{total += $2; count++} END {print total/count}'` > $specimen_name.rawR2_qualityCheck.txt
    grep "Total Sequences" $specimen_name*R1*/fastqc_data.txt >> $specimen_name.rawR1_qualityCheck.txt
    grep "Total Sequences" $specimen_name*R2*/fastqc_data.txt >> $specimen_name.rawR2_qualityCheck.txt
    grep "Sequence length" $specimen_name*R1*/fastqc_data.txt >> $specimen_name.rawR1_qualityCheck.txt
    grep "Sequence length" $specimen_name*R2*/fastqc_data.txt >> $specimen_name.rawR2_qualityCheck.txt
    grep ">>Basic Statistics" $specimen_name*R1*/fastqc_data.txt >> $specimen_name.rawR1_qualityCheck.txt
    grep ">>Basic Statistics" $specimen_name*R2*/fastqc_data.txt >> $specimen_name.rawR2_qualityCheck.txt
    grep ">>Per base sequence quality" $specimen_name*R1*/fastqc_data.txt >> $specimen_name.rawR1_qualityCheck.txt
    grep ">>Per base sequence quality" $specimen_name*R2*/fastqc_data.txt >> $specimen_name.rawR2_qualityCheck.txt
    grep ">>Per sequence quality scores" $specimen_name*R1*/fastqc_data.txt >> $specimen_name.rawR1_qualityCheck.txt
    grep ">>Per sequence quality scores" $specimen_name*R2*/fastqc_data.txt >> $specimen_name.rawR2_qualityCheck.txt
    grep ">>Per base sequence content" $specimen_name*R1*/fastqc_data.txt >> $specimen_name.rawR1_qualityCheck.txt
    grep ">>Per base sequence content" $specimen_name*R2*/fastqc_data.txt >> $specimen_name.rawR2_qualityCheck.txt
    grep ">>Per sequence GC content" $specimen_name*R1*/fastqc_data.txt >> $specimen_name.rawR1_qualityCheck.txt
    grep ">>Per sequence GC content" $specimen_name*R2*/fastqc_data.txt >> $specimen_name.rawR2_qualityCheck.txt
    grep ">>Sequence Duplication Levels" $specimen_name*R1*/fastqc_data.txt >> $specimen_name.rawR1_qualityCheck.txt
    grep ">>Sequence Duplication Levels" $specimen_name*R2*/fastqc_data.txt >> $specimen_name.rawR2_qualityCheck.txt
    grep ">>Overrepresented sequences" $specimen_name*R1*/fastqc_data.txt >> $specimen_name.rawR1_qualityCheck.txt
    grep ">>Overrepresented sequences" $specimen_name*R2*/fastqc_data.txt >> $specimen_name.rawR2_qualityCheck.txt
  done

#Run Quality trimming on sequences using BBDuk and run fastqc on newly trimmed sequences

cd $PTW

if [ -d sequence_analysis_$date ]; then
    echo "Warning previous sequence_analysis directory is still present and being deleted!"
    rm -rf sequence_analysis_$date
fi

mkdir $PTW/sequence_analysis_$date
mkdir $PTW/sequence_analysis_$date/Negative_Reads
PTA=$PTW/sequence_analysis_$date
PTAN=$PTW/sequence_analysis_$date/Negative_Reads

cat $PTF/shortnames.txt | while read specimen_name
do
    bbduk.sh in1=$PTF/tempFolder/$specimen_name\R1.fastq.gz in2=$PTF/tempFolder/$specimen_name\R2.fastq.gz out1=$PTF/fqc_output/bbduk_output/$specimen_name.clean1.fq out2=$PTF/fqc_output/bbduk_output/$specimen_name.clean2.fq minlen=50 \
  ref=$PTBBD/adapters.fa ktrim=r k=27 mink=11 hdist=1 -Xmx2g qtrim=r trimq=10
    fastqc $PTF/fqc_output/bbduk_output/$specimen_name.clean* -o $PTF/fqc_output/trimmed_seqs_fqc_results
done

#FORNEGATIVES

cat $PTF/Negative_Reads_$date/negnames.txt | while read negative_name
do
    bbduk.sh in1=$PTF/Negative_Reads_$date/tempFolder/$negative_name\R1.fastq.gz in2=$PTF/Negative_Reads_$date/tempFolder/$negative_name\R2.fastq.gz out1=$PTF/Negative_Reads_$date/fqc_output/bbduk_output/$negative_name.clean1.fq out2=$PTF/Negative_Reads_$date/fqc_output/bbduk_output/$negative_name.clean2.fq minlen=50 \
    ref=$PTBBD/adapters.fa ktrim=r k=27 mink=11 hdist=1 -Xmx2g qtrim=r trimq=10
    fastqc $PTF/Negative_Reads_$date/fqc_output/bbduk_output/$negative_name.clean* -o $PTF/Negative_Reads_$date/fqc_output/trimmed_seqs_fqc_results
done

#Create a new text file that contains the main post-trimming quality metrics

cd $PTF/fqc_output/trimmed_seqs_fqc_results
for x in *.zip; do unzip "$x"; done

cat $PTF/shortnames.txt | while read specimen_name;
do
    echo "Average Quality" `sed -n '/^#Base.*Mean/,/^>>END_MODULE/p' $specimen_name*clean1*/fastqc_data.txt | grep -v "#Base" | grep -v ">>END_MODULE" | awk '{total += $2; count++} END {print total/count}'` > $specimen_name.cleanR1_qualityCheck.txt
    echo "Average Quality" `sed -n '/^#Base.*Mean/,/^>>END_MODULE/p' $specimen_name*clean2*/fastqc_data.txt | grep -v "#Base" | grep -v ">>END_MODULE" | awk '{total += $2; count++} END {print total/count}'` > $specimen_name.cleanR2_qualityCheck.txt
    grep "Total Sequences" $specimen_name*clean1*/fastqc_data.txt >> $specimen_name.cleanR1_qualityCheck.txt
    grep "Total Sequences" $specimen_name*clean2*/fastqc_data.txt >> $specimen_name.cleanR2_qualityCheck.txt
    grep "Sequence length" $specimen_name*clean1*/fastqc_data.txt >> $specimen_name.cleanR1_qualityCheck.txt
    grep "Sequence length" $specimen_name*clean2*/fastqc_data.txt >> $specimen_name.cleanR2_qualityCheck.txt
    grep ">>Basic Statistics" $specimen_name*clean1*/fastqc_data.txt >> $specimen_name.cleanR1_qualityCheck.txt
    grep ">>Basic Statistics" $specimen_name*clean2*/fastqc_data.txt >> $specimen_name.cleanR2_qualityCheck.txt
    grep ">>Per base sequence quality" $specimen_name*clean1*/fastqc_data.txt >> $specimen_name.cleanR1_qualityCheck.txt
    grep ">>Per base sequence quality" $specimen_name*clean2*/fastqc_data.txt >> $specimen_name.cleanR2_qualityCheck.txt
    grep ">>Per sequence quality scores" $specimen_name*clean1*/fastqc_data.txt >> $specimen_name.cleanR1_qualityCheck.txt
    grep ">>Per sequence quality scores" $specimen_name*clean2*/fastqc_data.txt >> $specimen_name.cleanR2_qualityCheck.txt
    grep ">>Per base sequence content" $specimen_name*clean1*/fastqc_data.txt >> $specimen_name.cleanR1_qualityCheck.txt
    grep ">>Per base sequence content" $specimen_name*clean2*/fastqc_data.txt >> $specimen_name.cleanR2_qualityCheck.txt
    grep ">>Per sequence GC content" $specimen_name*clean1*/fastqc_data.txt >> $specimen_name.cleanR1_qualityCheck.txt
    grep ">>Per sequence GC content" $specimen_name*clean2*/fastqc_data.txt >> $specimen_name.cleanR2_qualityCheck.txt
    grep ">>Sequence Duplication Levels" $specimen_name*clean1*/fastqc_data.txt >> $specimen_name.cleanR1_qualityCheck.txt
    grep ">>Sequence Duplication Levels" $specimen_name*clean2*/fastqc_data.txt >> $specimen_name.cleanR2_qualityCheck.txt
    grep ">>Overrepresented sequences" $specimen_name*clean1*/fastqc_data.txt >> $specimen_name.cleanR1_qualityCheck.txt
    grep ">>Overrepresented sequences" $specimen_name*clean2*/fastqc_data.txt >> $specimen_name.cleanR2_qualityCheck.txt
done

#FORNEGATIVES

cd $PTF/Negative_Reads_$date/fqc_output/trimmed_seqs_fqc_results
for x in *.zip; do unzip "$x"; done

cat $PTF/Negative_Reads_$date/negnames.txt | while read specimen_name;
do
    echo "Average Quality" `sed -n '/^#Base.*Mean/,/^>>END_MODULE/p' $specimen_name*clean1*/fastqc_data.txt | grep -v "#Base" | grep -v ">>END_MODULE" | awk '{total += $2; count++} END {print total/count}'` > $specimen_name.cleanR1_qualityCheck.txt
    echo "Average Quality" `sed -n '/^#Base.*Mean/,/^>>END_MODULE/p' $specimen_name*clean2*/fastqc_data.txt | grep -v "#Base" | grep -v ">>END_MODULE" | awk '{total += $2; count++} END {print total/count}'` > $specimen_name.cleanR2_qualityCheck.txt
    grep "Total Sequences" $specimen_name*clean1*/fastqc_data.txt >> $specimen_name.cleanR1_qualityCheck.txt
    grep "Total Sequences" $specimen_name*clean2*/fastqc_data.txt >> $specimen_name.cleanR2_qualityCheck.txt
    grep "Sequence length" $specimen_name*clean1*/fastqc_data.txt >> $specimen_name.cleanR1_qualityCheck.txt
    grep "Sequence length" $specimen_name*clean2*/fastqc_data.txt >> $specimen_name.cleanR2_qualityCheck.txt
    grep ">>Basic Statistics" $specimen_name*clean1*/fastqc_data.txt >> $specimen_name.cleanR1_qualityCheck.txt
    grep ">>Basic Statistics" $specimen_name*clean2*/fastqc_data.txt >> $specimen_name.cleanR2_qualityCheck.txt
    grep ">>Per base sequence quality" $specimen_name*clean1*/fastqc_data.txt >> $specimen_name.cleanR1_qualityCheck.txt
    grep ">>Per base sequence quality" $specimen_name*clean2*/fastqc_data.txt >> $specimen_name.cleanR2_qualityCheck.txt
    grep ">>Per sequence quality scores" $specimen_name*clean1*/fastqc_data.txt >> $specimen_name.cleanR1_qualityCheck.txt
    grep ">>Per sequence quality scores" $specimen_name*clean2*/fastqc_data.txt >> $specimen_name.cleanR2_qualityCheck.txt
    grep ">>Per base sequence content" $specimen_name*clean1*/fastqc_data.txt >> $specimen_name.cleanR1_qualityCheck.txt
    grep ">>Per base sequence content" $specimen_name*clean2*/fastqc_data.txt >> $specimen_name.cleanR2_qualityCheck.txt
    grep ">>Per sequence GC content" $specimen_name*clean1*/fastqc_data.txt >> $specimen_name.cleanR1_qualityCheck.txt
    grep ">>Per sequence GC content" $specimen_name*clean2*/fastqc_data.txt >> $specimen_name.cleanR2_qualityCheck.txt
    grep ">>Sequence Duplication Levels" $specimen_name*clean1*/fastqc_data.txt >> $specimen_name.cleanR1_qualityCheck.txt
    grep ">>Sequence Duplication Levels" $specimen_name*clean2*/fastqc_data.txt >> $specimen_name.cleanR2_qualityCheck.txt
    grep ">>Overrepresented sequences" $specimen_name*clean1*/fastqc_data.txt >> $specimen_name.cleanR1_qualityCheck.txt
    grep ">>Overrepresented sequences" $specimen_name*clean2*/fastqc_data.txt >> $specimen_name.cleanR2_qualityCheck.txt
done

#Check to see if specimens pass the key quality metrics and remove the fastq files from the tempfolder if they do not
#Make directories so the User can see why a particular specimen failed
cd $PTF

mkdir $PTF/failed_sequence_runs/total_sequences_fail
mkdir $PTF/failed_sequence_runs/per_base_sequence_quality_fail
mkdir $PTF/failed_sequence_runs/basic_statistics_fail
mkdir $PTF/failed_sequence_runs/per_sequence_quality_scores_fail
mkdir $PTF/fqc_output/post_trimming_calculations
mkdir $PTF/failed_sequence_runs/specimen_total_sequences_fail

mkdir $PTF/Negative_Reads_$date/failed_negative_runs/total_sequences_fail
mkdir $PTF/Negative_Reads_$date/failed_negative_runs/per_base_sequence_quality_fail
mkdir $PTF/Negative_Reads_$date/failed_negative_runs/basic_statistics_fail
mkdir $PTF/Negative_Reads_$date/failed_negative_runs/per_sequence_quality_scores_fail
mkdir $PTF/Negative_Reads_$date/fqc_output/post_trimming_calculations
mkdir $PTF/Negative_Reads_$date/failed_negative_runs/negative_total_sequences_fail

#Check that the 3 main key quality matrics all passed

#PER BASE SEQUENCE QUALITY

cat $PTF/shortnames.txt | while read specimen_name;
do
  cd $PTF/fqc_output/trimmed_seqs_fqc_results
  PBSQ=`grep ">>Per base sequence quality" $specimen_name.cleanR1_qualityCheck.txt | awk '{print$5}'`
    if [[ $PBSQ == fail ]];
    then
        cd $PTF/failed_sequence_runs/per_base_sequence_quality_fail
        echo "$specimen_name failed per base sequence quality and can't continue in the workflow" > $specimen_name.txt
        cd $PTF/fqc_output/bbduk_output
        rm -rf $specimen_name.clean1.fq
        rm -rf $specimen_name.clean2.fq
    fi
done

#FORNEGATIVES

cat $PTF/Negative_Reads_$date/negnames.txt | while read negative_name;
do
  cd $PTF/Negative_Reads_$date/fqc_output/trimmed_seqs_fqc_results
  NPBSQ=`grep ">>Per base sequennce quality" $negative_name.cleanR1_qualityCheck.txt | awk '{print$5}'`
    if [[ $NPBSQ == fail ]];
    then
      cd $PTF/Negative_Reads_$date/failed_negative_runs/basic_statistics_fail
      echo "$negative_name failed per base sequence quality annd can't continue in the workflow" > $negative_name.txt
      cd $PTF/Negative_Reads_$date/fqc_output/bbduk_output
      rm -rf $specimen_name.clean1.fq
      rm -rf $specimen_name.clean2.fq
    fi
done

#BASIC STATISTICS

cat $PTF/shortnames.txt | while read specimen_name;
do
  cd $PTF/fqc_output/trimmed_seqs_fqc_results
  BS=`grep ">>Basic Statistics" $specimen_name.cleanR1_qualityCheck.txt | awk '{print$3}'`
    if [[ $BS == fail ]];
    then
        cd $PTF/failed_sequence_runs/basic_statistics_fail
        echo "$specimen_name failed per basic statistics and can't continue in the workflow" > $specimen_name.txt
        cd $PTF/fqc_output/bbduk_output
        rm -rf $specimen_name.clean1.fq
        rm -rf $specimen_name.clean2.fq
    fi
done

#FORNEGATIVES

cat $PTF/Negative_Reads_$date/negnames.txt | while read negative_name;
do
  cd $PTF/Negative_Reads_$date/fqc_output/trimmed_seqs_fqc_results
  NBS=`grep ">>Basic Statistics" $negative_name.cleanR1_qualityCheck.txt | awk '{print$3}'`
    if [[ $NBS == fail ]];
    then
      cd $PTF/Negative_Reads_$date/failed_negative_runs/basic_statistics_fail
      echo "$negative_name failed per basic statistics and can't continue in the workflow" > $negative_name.txt
      cd $PTF/Negative_Reads_$date/fqc_output/bbduk_output
      rm -rf $negative_name.clean1.fq
      rm -rf $negative_name.clean2.fq
    fi
done

#PER SEQUENCE QUALITY SCORES

cat $PTF/shortnames.txt | while read specimen_name;
do
  cd $PTF/fqc_output/trimmed_seqs_fqc_results
  PSQS=`grep ">>Per sequence quality scores" $specimen_name.cleanR1_qualityCheck.txt | awk '{print$5}'`
    if [[ $PSQS == fail ]];
    then
        cd $PTF/failed_sequence_runs/per_sequence_quality_scores_fail/
        echo "$specimen_name failed per base sequence quality and can't continue in the workflow" > $specimen_name.txt
        cd $PTF/fqc_output/bbduk_output
        rm -rf $specimen_name.clean1.fq
        rm -rf $specimen_name.clean2.fq
    fi
done

#FORNEGATIVES

cat $PTF/Negative_Reads_$date/negnames.txt | while read negative_name;
do
  cd $PTF/Negative_Reads_$date/fqc_output/trimmed_seqs_fqc_results
  NPSQS=`grep ">>Per sequence quality scores" $negative_name.cleanR1_qualityCheck.txt | awk '{print$5}'`
  if [[ $NPSQS == fail ]];
  then
    cd $PTF/Negative_Reads_$date/failed_negative_runs/per_sequence_quality_scores_fail
    echo "$negative_name failed per base sequence quality and can't continue in the workflow" > $negative_name.txt
    cd $PTF/Negative_Reads_$date/bbduk_output
    rm -rf $negative_name.clean1.fq
    rm -rf $negative_name.clean2.fq
  fi
done


#Check that the average of the total reads for all specimen is above at least 17,000

cat $PTF/shortnames.txt | while read specimen_name;
do
  cd $PTF/fqc_output/trimmed_seqs_fqc_results
  TS1=`grep "Total Sequences" $specimen_name.cleanR1_qualityCheck.txt | awk '{print$3}'`
  TS2=`grep "Total Sequences" $specimen_name.cleanR2_qualityCheck.txt | awk '{print$3}'`
  TS=`expr $TS1 + $TS2`
  echo $TS >> $PTF/fqc_output/post_trimming_calculations/total_sample_sequences.txt
done

#FORTHENEGATIVES

cat $PTF/Negative_Reads_$date/negnames.txt | while read negative_name;
do
  cd $PTF/Negative_Reads_$date/fqc_output/trimmed_seqs_fqc_results
  NTS1=`grep "Total Sequences" $negative_name.cleanR1_qualityCheck.txt | awk '{print$3}'`
  NTS2=`grep "Total Sequences" $negative_name.cleanR2_qualityCheck.txt | awk '{print$3}'`
  NTS=`expr $NTS1 + $NTS2`
  echo $NTS >> $PTF/Negative_Reads_$date/fqc_output/post_trimming_calculations/total_negative_sequences.txt
done

##Combine total sequences for negatives with total sequences for sample before checking the average and standard deviation## 
cat $PTF/Negative_Reads_$date/fqc_output/post_trimming_calculations/total_negative_sequences.txt $PTF/fqc_output/post_trimming_calculations/total_sample_sequences.txt > $PTF/fqc_output/post_trimming_calculations/total_sequences.txt

cd $PTF/fqc_output/post_trimming_calculations

avg=`cat total_sequences.txt | datamash mean 1`
avg=`printf "%.0f\n" "$avg"`
echo $avg > $PTF/fqc_output/post_trimming_calculations/ts_avg.txt


if [ $avg -lt 17000 ];
then
    echo "The average number of sequences for this run is too low and therfore can not be further analyzed by this workflow. WORKFLOW HALTED" 
    mkdir $PTF/failed_sequence_runs/$date.experiment_failed_average_sequences
    exit 0
fi

#Check that the total reads for each specimen is above the standard deviation of the total reads for all specimens

stdev=`cat total_sequences.txt | datamash sstdev 1`
stdev=`printf "%.0f\n" "$stdev"`
echo $stdev > $PTF/fqc_output/post_trimming_calculations/ts_stdev.txt

cd $PTF/fqc_output/post_trimming_calculations

cat $PTF/shortnames.txt | while read specimen_name;
do
  cd $PTF/fqc_output/trimmed_seqs_fqc_results
  TS1=`grep "Total Sequences" $specimen_name.cleanR1_qualityCheck.txt | awk '{print$3}'`
  TS2=`grep "Total Sequences" $specimen_name.cleanR2_qualityCheck.txt | awk '{print$3}'`
  TS=`expr $TS1 + $TS2`
  cd $PTF/fqc_output/post_trimming_calculations
    if [[ $TS -lt $stdev ]];
    then
      cd $PTF/failed_sequence_runs/specimen_total_sequences_fail
      echo "The total sequences for $specimen_name is less than the overall STDEV for all sequences and can't continue in this workflow" > $specimen_name.txt
      cd $PTF/fqc_output/bbduk_output
      rm -rf $specimen_name.clean1.fq
      rm -rf $specimen_name.clean2.fq
    fi
done

#FORTHENEGATIVES

cat $PTF/Negative_Reads_$date/negnames.txt | while read negative_name;
do
  cd $PTF/Negative_Reads_$date/fqc_output/trimmed_seqs_fqc_results
  NTS1=`grep "Total Sequences" $negative_name.cleanR1_qualityCheck.txt | awk '{print$3}'`
  NTS2=`grep "Total Sequences" $negative_name.cleanR2_qualityCheck.txt | awk '{print$3}'`
  NTS=`expr $NTS1 + $NTS2`
  cd $PTF/fqc_output/post_trimming_calculations
    if [[ $NTS -lt $stdev ]];
    then
      cd $PTF/Negative_Reads_$date/failed_negative_runs/negative_total_sequences_fail
      echo "The total sequences for $negative_name is less than the overall STDEV for all sequences and can't continue in this workflow" > $negative_name.txt
      cd $PTF/Negative_Reads_$date/fqc_output/bbduk_output
      rm -rf $negative_name.clean1.fq
      rm -rf $negative_name.clean2.fq
    fi
done


#Ensure that there are still at least 3 negatives that are continuing in the workflow for our cut-off threshold calculation and halt the workflow if there is not

cd $PTF/Negative_Reads_$date/fqc_output

if [ -d $PTF/Negative_Reads_$date/fqc_output/tempFolder ];
then
    echo "Previous Negative Read temporary folder is still present and is being deleted"
    rm -rf $PTF/Negative_Reads_$date/fqc_output/tempFolder
fi

if [ ! -z $PTF/Negative_Reads_$date/fqc_output/negatives_remaining.txt ];
then
    echo "Previous negative read text file is still present and is being deleted"
    rm -rf $PTF/Negative_Reads_$date/fqc_output/negatives_remaining.txt
fi

mkdir $PTF/Negative_Reads_$date/fqc_output/tempFolder

cp $PTF/Negative_Reads_$date/fqc_output/bbduk_output/*.clean1.fq $PTF/Negative_Reads_$date/fqc_output/tempFolder

ls $PTF/Negative_Reads_$date/fqc_output/tempFolder | awk -F "/" '{print$NF}' | awk -F "R" '{print$1}'  >> $PTF/Negative_Reads_$date/fqc_output/negatives_remaining.txt

cd $PTF/Negative_Reads_$date/fqc_output

negatives_remaining=`wc -l < negatives_remaining.txt`
echo $negatives_remaining

if [ $negatives_remaining -lt 3 ];
then
    echo "The amount of negatives remaining in the workflow is too few to calculate a reliable cut-off threshold. WORKFLOW HALTED!"
    exit 0
fi

rm -rf $PTF/Negative_Reads_$date/tempFolder
rm -rf $PTF/Negative_Reads_$date/fqc_output/negatives_remaining.txt
rm -rf $PTF/tempFolder

#Trim F and R UPDx Primers

cd $PTF/fqc_output/bbduk_output
mkdir $PTF/fqc_output/bbduk_output/primers_trimmed

cat $PTF/shortnames.txt | while read specimen_name;
do
  cutadapt -g file:$r_primer \
  -g file:$f_primer \
  -G file:$r_primer \
  -G file:$f_primer \
  -j 0 \
  -o $PTF/fqc_output/bbduk_output/primers_trimmed/$specimen_name.Ptrimmed.r1.fq \
  -p $PTF/fqc_output/bbduk_output/primers_trimmed/$specimen_name.Ptrimmed.r2.fq \
  $PTF/fqc_output/bbduk_output/$specimen_name.*clean1*.fq $PTF/fqc_output/bbduk_output/$specimen_name.*clean2*.fq \
  >> $PTF/fqc_output/bbduk_output/primers_trimmed/cutadapt_primer_trimming_stats.txt
done

#FORNEGATIVES

cd $PTF/Negative_Reads_$date/fqc_output/bbduk_output || exit
mkdir $PTF/Negative_Reads_$date/fqc_output/bbduk_output/primers_trimmed

cd $PTF/Negative_Reads_$date/fqc_output/bbduk_output || exit


cat $PTF/Negative_Reads_$date/negnames.txt | while read negative_name;
do
  cutadapt -g file:$r_primer \
  -g file:$f_primer \
  -G file:$r_primer \
  -G file:$f_primer \
  -j 0 \
  -o $PTF/Negative_Reads_$date/fqc_output/bbduk_output/primers_trimmed/$negative_name.Ptrimmed.r1.fq \
  -p $PTF/Negative_Reads_$date/fqc_output/bbduk_output/primers_trimmed/$negative_name.Ptrimmed.r2.fq \
  $PTF/Negative_Reads_$date/fqc_output/bbduk_output/$negative_name.*clean1*.fq $PTF/Negative_Reads_$date/fqc_output/bbduk_output/$negative_name.*clean2*.fq
  >> $PTF/Negative_Reads_$date/fqc_output/bbduk_output/primers_trimmed/cutadapt_primer_trimming_stats.txt
done

#Take the clean reads that passed the quality check and merge the paired reads using bbmerge
mkdir $PTF/fqc_output/BBMerge_Output
mkdir $PTA/fqc_processed_reads 
cd $PTF/fqc_output/bbduk_output

cat $PTF/shortnames.txt | while read specimen_name;
do
    bbmerge-auto.sh \
    in1=$PTF/fqc_output/bbduk_output/primers_trimmed/$specimen_name.Ptrimmed.r1.fq \
    in2=$PTF/fqc_output/bbduk_output/primers_trimmed/$specimen_name.Ptrimmed.r2.fq \
    out=$PTF/fqc_output/BBMerge_Output/trimmed_merged_reads/$specimen_name.merged.fq \
    outu=$PTF/fqc_output/BBMerge_Output/trimmed_unmerged_reads/$specimen_name.unmerged.fq \
    -Xmx2g \
    rem k=62 \
    extend2=50 ecct 
    cp $PTF/fqc_output/BBMerge_Output/trimmed_merged_reads/$specimen_name.merged.fq $PTA/fqc_processed_reads
done

#FORNEGATIVES

mkdir $PTAN/fqc_processed_reads

cd $PTF/Negative_Reads_$date/fqc_output/bbduk_output

cat $PTF/Negative_Reads_$date/negnames.txt | while read negative_name;
do
  bbmerge-auto.sh \
  in1=$PTF/Negative_Reads_$date/fqc_output/bbduk_output/primers_trimmed/$negative_name.Ptrimmed.r1.fq \
  in2=$PTF/Negative_Reads_$date/fqc_output/bbduk_output/primers_trimmed/$negative_name.Ptrimmed.r2.fq \
  out=$PTF/Negative_Reads_$date/fqc_output/BBMerge_Output/trimmed_merged_reads/$negative_name.merged.fq \
  outu=$PTF/Negative_Reads_$date/fqc_output/BBMerge_Output/trimmed_unmerged_reads/$negative_name.unmerged.fq \
  -Xmx2g \
  rem k=62 \
  extend2=50 ecct
  cp $PTF/Negative_Reads_$date/fqc_output/BBMerge_Output/trimmed_merged_reads/$negative_name.merged.fq $PTAN/fqc_processed_reads
done
#Assemble clusters using CDHIT and blast to the UPDx exclusion database filtering out the clusters that have hits

cd $PTA

if [ -d $PTA/cd_hit_output ]; then
  echo "previous clusters folder was still present and is being deleted"
  rm -rf cd_hit_output
fi

if [ -d $PTA/exc_blast ]; then
  echo "previous exclusion blast folder was still present and is being deleted"
  rm -rf exc_blast
fi



mkdir $PTA/cd_hit_output
mkdir $PTA/exc_blast



#FORNEGATIVES 

cd $PTAN || exit

if [ -d $PTAN/cd_hit_output ]; then
  echo "previous clusters folder was still present and is being deleted"
  rm -rf cd_hit_output
fi

if [ -d $PTAN/exc_blast ]; then
  echo "previous exclusion blast folder was still present and is being deleted"
  rm -rf exc_blast
fi


mkdir $PTAN/cd_hit_output
mkdir $PTAN/exc_blast
mkdir $PTA/tempFolder
mkdir $PTAN/tempFolder



cd $PTA/fqc_processed_reads || exit

cat $PTF/shortnames.txt | while read specimen_name;
do
  seqtk seq -a $PTA/fqc_processed_reads/$specimen_name.merged.fq > $PTA/fqc_processed_reads/$specimen_name.merged.fa
  cd-hit-est -i $PTA/fqc_processed_reads/$specimen_name.merged.fa \
            -o $PTA/cd_hit_output/$specimen_name.clusters.fa \
            -d 0 \
            -c 0.99 \
            -s 1.00 \
            -g 1.00
####(added in order to make the count script work)#### 
  cd $PTA/cd_hit_output
  new_clstr=`sed 's/>Cluster />Cluster/g' $specimen_name.clusters.fa.clstr |  sed 's/nt, /nt,/g' | sed 's/,>M0/,M0/g' > $PTA/cd_hit_output/$specimen_name.clusters.updated.clstr`
  ####                                              ####
#Call Python script to create a csv file that contains each cluster and the total reads count
  cd $py_scripts
  python3 cluster_counter.py \
  $PTA/cd_hit_output/$specimen_name.clusters.updated.clstr \
  $PTA/cd_hit_output/count_summary_$specimen_name.csv 
#Call Python script to find Clusters with more than 20 sequences
  python3 cluster_filter.py \
  --count_df $PTA/cd_hit_output/count_summary_$specimen_name.csv  \
  --clusters_filtered $PTA/cd_hit_output/$specimen_name.clusters_greater_than_20_temp.txt
#Copy documents into a temporary folder to get the sequence ID of clusters with more than 20 sequences
  cp $PTA/cd_hit_output/$specimen_name.clusters_greater_than_20_temp.txt $PTA/tempFolder
  rm -rf $PTA/cd_hit_output/$specimen_name.clusters_greater_than_20_temp.txt
  cp $PTA/cd_hit_output/$specimen_name.clusters.updated.clstr $PTA/tempFolder
  cp $PTA/cd_hit_output/$specimen_name.clusters.fa $PTA/tempFolder
  cd $PTA/tempFolder
#Parse the .clstr file to match the seqID with the cluster number
  awk -v OFS="\n" '/^>/ {getline seq; print $0, seq}' $PTA/tempFolder/$specimen_name.clusters.updated.clstr > $PTA/tempFolder/$specimen_name.clstr_file_simplified_temp.txt
  sed 's/nt,/nt, /g' < $PTA/tempFolder/$specimen_name.clstr_file_simplified_temp.txt | sed 's/\.\.\.//g' > $PTA/tempFolder/$specimen_name.clstr_file_simplified.txt
  rm -rf $PTA/tempFolder/$specimen_name.clstr_file_simplified_temp.txt
  awk '$1~">Cluster[0-9]*"{print$1,$2}' < $PTA/tempFolder/$specimen_name.clstr_file_simplified.txt > $PTA/tempFolder/$specimen_name.cluster_numbers.txt
  awk '$1!~">Cluster[0-9]*"{print$3}' < $PTA/tempFolder/$specimen_name.clstr_file_simplified.txt > $PTA/tempFolder/$specimen_name.seq_IDs.txt
  paste $PTA/tempFolder/$specimen_name.cluster_numbers.txt $PTA/tempFolder/$specimen_name.seq_IDs.txt > $PTA/tempFolder/$specimen_name.cluster_numbers_matched_to_seq_IDs.txt
  sed 's/Cluster/>Cluster/g' $PTA/tempFolder/$specimen_name.clusters_greater_than_20_temp.txt > $PTA/tempFolder/$specimen_name.clusters_greater_than_20_temp2.txt
  sed 's/ >Cluster/>Cluster/g' $PTA/tempFolder/$specimen_name.clusters_greater_than_20_temp2.txt > $PTA/tempFolder/$specimen_name.clusters_greater_than_20.txt
  rm -rf $PTA/tempFolder/$specimen_name.clusters_greater_than_20_temp.txt
  rm -rf $PTA/tempFolder/$specimen_name.clusters_greater_than_20_temp2.txt
#Filter the cluster_numbers_matched_to_seq_IDs for only cluster number that are in the clusters_greater_than_20 file
  grep -f $PTA/tempFolder/$specimen_name.clusters_greater_than_20.txt $PTA/tempFolder/$specimen_name.cluster_numbers_matched_to_seq_IDs.txt > $PTA/tempFolder/$specimen_name.cluster_numbers_matched_to_seq_IDs_filtered.txt
#Get just the seqIDs from newly created filtered text file to match with fasta file
  awk '{print$2}' $PTA/tempFolder/$specimen_name.cluster_numbers_matched_to_seq_IDs_filtered.txt > $PTA/tempFolder/$specimen_name.seq_IDs_filtered_temp.txt
  sed -r 's/^/>/g' $PTA/tempFolder/$specimen_name.seq_IDs_filtered_temp.txt > $PTA/tempFolder/$specimen_name.seq_Ids_filtered.txt
  rm -rf $PTA/tempFolder/$specimen_name.seq_IDs_filtered_temp.txt
#Get the sequence on the same line as the ID in the fasta file
  paste -s -d' \n' $PTA/tempFolder/$specimen_name.clusters.fa > $PTA/tempFolder/$specimen_name.clusters_filtered_temp1.fa
#Filter the fasta file for only the clusters with greater than 20 sequences
  grep -f $PTA/tempFolder/$specimen_name.seq_Ids_filtered.txt $PTA/tempFolder/$specimen_name.clusters_filtered_temp1.fa > $PTA/tempFolder/$specimen_name.clusters_filtered_temp2.fa
#Get the filtered fasta file back into the correct format
  sed -E 's/ ([A-Z]*)$/\n\1/g' $PTA/tempFolder/$specimen_name.clusters_filtered_temp2.fa > $PTA/tempFolder/$specimen_name.clusters_filtered1.fa
  cp $PTA/tempFolder/$specimen_name.clusters_filtered1.fa $PTA/cd_hit_output
  cd $PTW/UPDx_Variables/UPDx_18s_Exclusion_Database
  makeblastdb -in $exc_database -dbtype nucl -parse_seqids
  blastn -db $exc_database -query $PTA/cd_hit_output/$specimen_name.clusters_filtered1.fa -outfmt 6 -out $PTA/exc_blast/$specimen_name.exclusion_blast.fa -evalue 0.00001 -perc_identity 96 -max_target_seqs 1
  awk '{ if ($0 ~ /_/) { printf ">"; } print $0; }' $PTA/exc_blast/$specimen_name.exclusion_blast.fa > $PTA/exc_blast/$specimen_name.exclusion_blast_modified.fa
#Using the output of the exclusion BLAST filter out representative sequences that had hits towards the exclusion database
  filterbyname.sh in=$PTA/cd_hit_output/$specimen_name.clusters_filtered1.fa out=$PTA/cd_hit_output/$specimen_name.clusters.filtered.fa names=$PTA/exc_blast/$specimen_name.exclusion_blast_modified.fa include=f substring=header
done
#rm -rf $PTA/tempFolder

#FORNEGATIVES

cd $PTAN/fqc_processed_reads || exit

cat $PTF/Negative_Reads_$date/negnames.txt | while read negative_name;
do
  seqtk seq -a $PTAN/fqc_processed_reads/$negative_name.merged.fq > $PTAN/fqc_processed_reads/$negative_name.merged.fa
  cd-hit-est -i $PTAN/fqc_processed_reads/$negative_name.merged.fa \
            -o $PTAN/cd_hit_output/$negative_name.clusters.fa \
            -d 0 \
            -c 0.99 \
            -s 1.00 \
            -g 1.00
####(added in order to make the count script work)#### 
  cd $PTAN/cd_hit_output
  new_clstr=`sed 's/>Cluster />Cluster/g' $negative_name.clusters.fa.clstr |  sed 's/nt, /nt,/g' | sed 's/,>M0/,M0/g' > $PTAN/cd_hit_output/$negative_name.clusters.updated.clstr`
####                                              ####
#Call Python script to create a csv file that contains each cluster and the total reads count
  cd $py_scripts
  python3 cluster_counter.py \
  $PTAN/cd_hit_output/$negative_name.clusters.updated.clstr \
  $PTAN/cd_hit_output/count_summary_$negative_name.csv
  python3 cluster_filter.py \
  --count_df $PTAN/cd_hit_output/count_summary_$negative_name.csv  \
  --clusters_filtered $PTAN/cd_hit_output/$negative_name.clusters_greater_than_20_temp.txt
  python3 cluster_filter.py \
  --count_df $PTAN/cd_hit_output/count_summary_$negative_name.csv  \
  --clusters_filtered $PTAN/cd_hit_output/$negative_name.clusters_greater_than_20_temp.txt
  cp $PTAN/cd_hit_output/$negative_name.clusters_greater_than_20_temp.txt $PTAN/tempFolder
  rm -rf $PTAN/cd_hit_output/$negative_name.clusters_greater_than_20_temp.txt
  cp $PTAN/cd_hit_output/$negative_name.clusters.updated.clstr $PTAN/tempFolder
  cp $PTAN/cd_hit_output/$negative_name.clusters.fa $PTAN/tempFolder
  cd $PTAN/tempFolder
  awk -v OFS="\n" '/^>/ {getline seq; print $0, seq}' $PTAN/tempFolder/$negative_name.clusters.updated.clstr > $PTAN/tempFolder/$negative_name.clstr_file_simplified_temp.txt
  sed 's/nt,/nt, /g' < $PTAN/tempFolder/$negative_name.clstr_file_simplified_temp.txt | sed 's/\.\.\.//g' > $PTAN/tempFolder/$negative_name.clstr_file_simplified.txt
  rm -rf $PTAN/tempFolder/$negative_name.clstr_file_simplified_temp.txt
  awk '$1~">Cluster[0-9]*"{print$1,$2}' < $PTAN/tempFolder/$negative_name.clstr_file_simplified.txt > $PTAN/tempFolder/$negative_name.cluster_numbers.txt
  awk '$1!~">Cluster[0-9]*"{print$3}' < $PTAN/tempFolder/$negative_name.clstr_file_simplified.txt > $PTAN/tempFolder/$negative_name.seq_IDs.txt
  paste $PTAN/tempFolder/$negative_name.cluster_numbers.txt $PTAN/tempFolder/$negative_name.seq_IDs.txt > $PTAN/tempFolder/$negative_name.cluster_numbers_matched_to_seq_IDs.txt
  sed 's/Cluster/>Cluster/g' $PTAN/tempFolder/$negative_name.clusters_greater_than_20_temp.txt > $PTAN/tempFolder/$negative_name.clusters_greater_than_20_temp2.txt
  sed 's/ >Cluster/>Cluster/g' $PTAN/tempFolder/$negative_name.clusters_greater_than_20_temp2.txt > $PTAN/tempFolder/$negative_name.clusters_greater_than_20.txt
  rm -rf $PTAN/tempFolder/$negative_name.clusters_greater_than_20_temp.txt
  rm -rf $PTAN/tempFolder/$negative_name.clusters_greater_than_20_temp2.txt
  grep -f $PTAN/tempFolder/$negative_name.clusters_greater_than_20.txt $PTAN/tempFolder/$negative_name.cluster_numbers_matched_to_seq_IDs.txt > $PTAN/tempFolder/$negative_name.cluster_numbers_matched_to_seq_IDs_filtered.txt
  awk '{print$2}' $PTAN/tempFolder/$negative_name.cluster_numbers_matched_to_seq_IDs_filtered.txt > $PTAN/tempFolder/$negative_name.seq_IDs_filtered_temp.txt
  sed -r 's/^/>/g' $PTAN/tempFolder/$negative_name.seq_IDs_filtered_temp.txt > $PTAN/tempFolder/$negative_name.seq_Ids_filtered.txt
  rm -rf $PTAN/tempFolder/$negative_name.seq_IDs_filtered_temp.txt
  paste -s -d' \n' $PTAN/tempFolder/$negative_name.clusters.fa > $PTAN/tempFolder/$negative_name.clusters_filtered_temp1.fa
  grep -f $PTAN/tempFolder/$negative_name.seq_Ids_filtered.txt $PTAN/tempFolder/$negative_name.clusters_filtered_temp1.fa > $PTAN/tempFolder/$negative_name.clusters_filtered_temp2.fa
  sed -E 's/ ([A-Z]*)$/\n\1/g' $PTAN/tempFolder/$negative_name.clusters_filtered_temp2.fa > $PTAN/tempFolder/$negative_name.clusters_filtered1.fa
  cp $PTAN/tempFolder/$negative_name.clusters_filtered1.fa $PTAN/cd_hit_output
  cd $PTW/UPDx_Variables/UPDx_18s_Exclusion_Database
  makeblastdb -in $exc_database -dbtype nucl -parse_seqids
  blastn -db $exc_database -query $PTAN/cd_hit_output/$negative_name.clusters_filtered1.fa -outfmt 6 -out $PTAN/exc_blast/$negative_name.exclusion_blast.fa -evalue 0.00001 -perc_identity 96 -max_target_seqs 1
  awk '{ if ($0 ~ /_/) { printf ">"; } print $0; }' $PTAN/exc_blast/$negative_name.exclusion_blast.fa > $PTAN/exc_blast/$negative_name.exclusion_blast_modified.fa
  filterbyname.sh in=$PTAN/cd_hit_output/$negative_name.clusters_filtered1.fa out=$PTAN/cd_hit_output/$negative_name.clusters.filtered.fa names=$PTAN/exc_blast/$negative_name.exclusion_blast_modified.fa include=f substring=header
done
#rm -rf $PTAN/tempFolder


#Take the representative sequences that were not filtered out and BLAST them to the reference database

cd $PTA

if [ -d $PTA/ref_blast ]; then
    echo "Previous reference BLAST folder is still present and is being deleted!"
    rm -rf $PTA/ref_blast
fi
mkdir $PTA/ref_blast
nucl
PTP=$PTA/ref_blast
cd $PTW/UPDx_Variables/UPDx_18s_Reference_Database


cat $PTF/shortnames.txt | while read specimen_name;
do
    cd $PTP
    makeblastdb -in $ref_database -dbtype nucl
    blastn -query $PTA/cd_hit_output/$specimen_name.clusters.filtered.fa \
            -db $ref_database \
            -outfmt "6 qseqid pident sseqid qseq evalue"  \
            -evalue 0.001 \
            -max_target_seqs 1 \
            -out $PTP/$specimen_name.parasite_BLAST_results
    awk '{ if ($0 ~ /_/) { printf ">"; } print $0; }' $PTP/$specimen_name.parasite_BLAST_results > $PTP/$specimen_name.parasite_BLAST_results.modified
done



#FORNEGATIVES

cd $PTAN

if [ -d $PTAN/ref_blast ]; then
    echo "Previous negative reference BLAST folder is still present and is being deleted!"
    rm -rf $PTAN/ref_blast
fi

mkdir $PTAN/ref_blast

PTPN=$PTAN/ref_blast


cat $PTF/Negative_Reads_$date/negnames.txt | while read negative_name;
do
    cd $PTPN 
    makeblastdb -in $ref_database -dbtype nucl
    blastn -query  $PTAN/cd_hit_output/$negative_name.clusters.filtered.fa \
            -db $ref_database \
            -outfmt "6 qseqid pident sseqid qseq evalue"  \
            -evalue 0.001 \
            -max_target_seqs 1 \
            -out $PTPN/$negative_name.parasite_BLAST_results
    awk '{ if ($0 ~ /_/) { printf ">"; } print $0; }' $PTPN/$negative_name.parasite_BLAST_results > $PTPN/$negative_name.parasite_BLAST_results.modified
done

Final_Blast_Time=$(date +"%H_%M_%S")


if [ -d $PTW/Final_run_summary_$date ]; then
    echo "Previous Summary folder is still present and is being deleted!"
    rm -rf $PTW/Final_run_summary_$date
fi

mkdir $PTW/FINAL_run_summary_$date
PTS=$PTW/FINAL_run_summary_$date
mkdir $PTS/Samples
mkdir $PTS/Negatives


cat $PTF/shortnames.txt | while read specimen_name;
do
  if [ -d $PTS/Samples/$specimen_name.summary ]; then
    echo "Previous Summary folder is still present and is being deleted!"
    rm -rf $PTS/Samples/$specimen_name.summary
  fi
  cd $PTS/Samples
  mkdir $PTS/Samples/$specimen_name.summary
  cd $PTA/cd_hit_output
  cp $specimen_name.clusters.updated.clstr $PTS/Samples/$specimen_name.summary
  rm -rf $PTA/cd_hit_output/$specimen_name.clusters.updated.clstr

  #Call Python script to create a library that connects each cluster to its representative sequence
  cd $py_scripts
  python3 cluster_library.py \
  --clstr_file $PTS/Samples/$specimen_name.summary/$specimen_name.clusters.updated.clstr \
  --fasta_file $PTA/cd_hit_output/$specimen_name.clusters.fa \
  --summary_table $PTS/Samples/$specimen_name.summary/summary.$specimen_name.csv 

  #Call Python script to create a csv file that contains each cluster and the total reads count
  cd $PTA/cd_hit_output
  cp count_summary_$specimen_name.csv $PTS/Samples/$specimen_name.summary
  #Call Python script to merge the two dataframes together and create one summary file
  cd $py_scripts
  python3 column_merge.py \
  --count_summary $PTS/Samples/$specimen_name.summary/count_summary_$specimen_name.csv \
  --rep_seq_summary $PTS/Samples/$specimen_name.summary/summary.$specimen_name.csv \
  --full_summary $PTS/Samples/$specimen_name.summary/temp1_Summary_Table.csv

  #Call Python script to pull out the sequences from the fasta file and add it to the summary file
  cd $py_scripts
  python3 sequence_identifier.py \
  --fasta_file $PTA/cd_hit_output/$specimen_name.clusters.fa \
  --summary_table $PTS/Samples/$specimen_name.summary/temp1_Summary_Table.csv \
  --new_summary_table $PTS/Samples/$specimen_name.summary/Cluster_Summary_Table.csv

  #Call Python script to turn the exc BLAST output into an overall BLAST Summary Table
  cd $py_scripts
  python3 exc_BLAST_Summary_Table.py \
  --blast_output $PTA/exc_blast/$specimen_name.exclusion_blast_modified.fa \
  --BLAST_summary_table $PTS/Samples/$specimen_name.summary/exc_BLAST_Summary_Table.csv

  #Call Python script to turn the ref BLAST output into an overall BLAST Summary Table
  cd $py_scripts
  python3 BLAST_Summary_Table.py \
  --blast_output $PTP/$specimen_name.parasite_BLAST_results.modified \
  --cluster_summary $PTS/Samples/$specimen_name.summary/Cluster_Summary_Table.csv \
  --Complete_Summary_Table $PTS/Samples/$specimen_name.summary/temp2_Summary_Table.csv \
  --BLAST_summary_table $PTS/Samples/$specimen_name.summary/ref_BLAST_Summary_Table.csv
done


#FORNEGATIVES

PTSN=$PTW/FINAL_run_summary_$date/Negatives

cat $PTF/Negative_Reads_$date/negnames.txt | while read negative_name;
do
  cd $PTSN
  mkdir $PTSN/$negative_name.summary
  cd $PTAN/cd_hit_output
  ####(added in order to make the count script work)### 
  #new_clstr=`sed 's/>Cluster />Cluster/g' $negative_name.clusters.fa.clstr |  sed 's/nt, /nt,/g' | sed 's/,>M0/,M0/g' > $PTAN/cd_hit_output/$negative_name.clusters.updated.clstr`
  cp $negative_name.clusters.updated.clstr /$PTSN/$negative_name.summary
  rm -rf $PTAN/cd_hit_output/$negative_name.clusters.updated.clstr

  #Call Python script to create a library that connects each cluster to its representative sequence
  cd $py_scripts
  python3 cluster_library.py \
  --clstr_file $PTSN/$negative_name.summary/$negative_name.clusters.updated.clstr \
  --fasta_file $PTAN/cd_hit_output/$negative_name.clusters.fa \
  --summary_table $PTSN/$negative_name.summary/summary.$negative_name.csv 

  #Call Python script to create a csv file that contains each cluster and the total reads count
  # cd $py_scripts
  # python3 cluster_counter.py \
  # $PTSN/$negative_name.summary/$negative_name.clusters.updated.clstr \
  # $PTSN/$negative_name.summary/count_summary_$negative_name.csv
  cp $PTAN/cd_hit_output/count_summary_$negative_name.csv $PTSN/$negative_name.summary
  #Call Python script to merge the two dataframes together and create one summary file
  cd $py_scripts
  python3 column_merge.py \
  --count_summary $PTSN/$negative_name.summary/count_summary_$negative_name.csv \
  --rep_seq_summary $PTSN/$negative_name.summary/summary.$negative_name.csv \
  --full_summary $PTSN/$negative_name.summary/temp1_Summary_Table.csv

  #Call Python script to pull out the sequences from the fasta file and add it to the summary file
  cd $py_scripts
  python3 sequence_identifier.py \
  --fasta_file $PTAN/cd_hit_output/$negative_name.clusters.fa \
  --summary_table $PTSN/$negative_name.summary/temp1_Summary_Table.csv \
  --new_summary_table $PTSN/$negative_name.summary/Cluster_Summary_Table.csv

  #Call Python script to turn the exc BLAST output into an overall BLAST Summary Table
  cd $py_scripts
  python3 exc_BLAST_Summary_Table.py \
  --blast_output $PTA/Negative_Reads/exc_blast/$negative_name.exclusion_blast_modified.fa \
  --BLAST_summary_table $PTSN/$negative_name.summary/exc_BLAST_Summary_Table.csv

  #Call Python script to turn the BLAST output into an overall BLAST Summary Table
  cd $py_scripts 
  python3 neg_BLAST_Summary_Table.py \
  --blast_output $PTPN/$negative_name.parasite_BLAST_results.modified \
  --cluster_summary $PTSN/$negative_name.summary/Cluster_Summary_Table.csv \
  --Complete_Summary_Table $PTSN/$negative_name.summary/COMPLETE_Summary_Table.csv \
  --BLAST_summary_table $PTSN/$negative_name.summary/ref_BLAST_Summary_Table.csv 

  #Call Python script to calculate the proportion of contaminated reads for each negative
  cd $py_scripts
  python3 Proportion_Calculator.py \
  --Complete_Summary_Table $PTSN/$negative_name.summary/COMPLETE_Summary_Table.csv \
  --proportion_txt $PTSN/$negative_name.proportion.txt \
  --list_of_parasites $PTSN/$negative_name.parasites_used_for_cutoff.txt \
  --parasite_table $PTSN/$negative_name.parasites_used_for_cutoff.csv

  #Modify the text file to get it in the correct formula for further calculations
  cd $PTSN
  paste -d, *.proportion.txt > $PTS/Samples/proportion_of_contam_reads.txt
  #sed -r 's/(0.[0-9]+)(0.[0-9]+)(0.[0-9]+)(.*)/\1,\2,\3,\4/g' $PTS/Samples/proportion_of_contam_reads.txt > $PTS/Samples/proportion_of_contam_reads.modified.txt 

  #Call Python script to calculate the average, standard deviation, and cut-off multiplier
  cd $py_scripts
  python3 Cut_off_calculator.py \
  --proportions_txt $PTS/Samples/proportion_of_contam_reads.txt \
  --cutoff_multiplier $PTS/Samples/cutoff_multiplier.txt
done

#Call Python Script to calculate the cutoff threshold for each sample and create the COMPLETE summary table and the parasite hits table for each sample
#Remove all the temporary tables created so you are left with just the Complete Summary Table, Cluster Summary Table, and BLAST summary Table
mkdir $PTS/Samples/Parasite_Hits
mkdir $PTS/Cutoff_Stats
cat $PTF/shortnames.txt | while read specimen_name;
do
  cd $py_scripts
  python3 Threshold_Calculator.py \
  --cutoff_multiplier $PTS/Samples/cutoff_multiplier.txt \
  --Complete_Summary_Table $PTS/Samples/$specimen_name.summary/temp2_Summary_Table.csv \
  --Cutoff_Summary_Table $PTS/Samples/$specimen_name.summary/COMPLETE_Summary_Table.csv \
  --Sample_Parasite_Table $PTS/Samples/Parasite_Hits/$specimen_name.parasite_hits.csv \
  --markdown_summary_table $PTS/Samples/$specimen_name.summary/markdown_summary_table.csv \
  --threshold_summary_text $PTS/Samples/$specimen_name.summary/threshold_int.txt \
  --results_table $PTS/Samples/$specimen_name.summary/results_table.csv \
  --specimen_name $specimen_name \
  --unclear_seqs_fasta_table $PTS/Samples/$specimen_name.summary/unclear_seqs_fasta_table.csv
  awk -F , '{print ">"$1"\n"$2}' $PTS/Samples/$specimen_name.summary/unclear_seqs_fasta_table.csv > $PTS/Samples/$specimen_name.summary/unclear_hit_sequences.fasta
  cd $PTS/Samples/$specimen_name.summary
  rm -rf summary.$specimen_name.csv
  rm -rf count_summary_$specimen_name.csv
  rm -rf temp1_Summary_Table.csv
  rm -rf temp2_Summary_Table.csv
  rm -rf $specimen_name.clusters.updated.clstr
done

cd $PTS/Samples/Parasite_Hits
cat * > $PTS/Samples/Parasite_Hits/Parasite_Summary_Table.csv


cd $PTS
cp $PTS/Samples/proportion_of_contam_reads.txt $PTS/Cutoff_Stats
cp $PTS/Samples/cutoff_multiplier.txt $PTS/Cutoff_Stats 
rm -rf $PTS/Samples/proportion_of_contam_reads.txt
rm -rf $PTS/Samples/cutoff_multiplier.txt

mkdir $PTS/Cutoff_Stats/temp_folder

cat $PTF/Negative_Reads_$date/negnames.txt | while read negative_name;
do
  cd $PTSN/$negative_name.summary
  rm -rf summary.$negative_name.csv
  rm -rf count_summary_$negative_name.csv
  rm -rf temp1_Summary_Table.csv
  rm -rf temp2_Summary_Table.csv
  rm -rf $negative_name.clusters.updated.clstr
  cd $PTSN
  cp $PTSN/$negative_name.proportion.txt $PTS/Cutoff_Stats
  cp $PTSN/$negative_name.parasites_used_for_cutoff.txt $PTS/Cutoff_Stats
  cp $PTSN/$negative_name.parasites_used_for_cutoff.csv $PTS/Cutoff_Stats/temp_folder
  rm -rf $PTSN/$negative_name.proportion.txt
  rm -rf $PTSN/$negative_name.parasites_used_for_cutoff.txt
  rm -rf $PTSN/$negative_name.parasites_used_for_cutoff.csv
done

cd $PTS/Cutoff_Stats/temp_folder
cat * > parasites_used_for_cutoff.csv

cd $PTW

cp -r FQC_$date $PTE
rm -rf FQC_$date

cp -r $PTW/sequence_analysis_$date $PTE
rm -rf $PTW/sequence_analysis_$date

#Generate a tree with the positive sequences compared to the reference database


# cat $PTE/FQC_$date/shortnames.txt | while read specimen_name;
# do
#   mkdir $PTS/Samples/$specimen_name.summary/Unclear_Hits
#   mkdir $PTS/Samples/$specimen_name.summary/Unclear_Hits/BLAST
#   mkdir $PTS/Samples/$specimen_name.summary/Unclear_Hits/BLAST/tempFolder
#   mkdir $PTS/Samples/$specimen_name.summary/Unclear_Hits/BLAST/Results
#   mkdir $PTS/Samples/$specimen_name.summary/Unclear_Hits/Clustering
#   mkdir $PTS/Samples/$specimen_name.summary/Unclear_Hits/Orientation_Check
#   cd $PTS/Samples/$specimen_name.summary
#   cp $PTS/Samples/$specimen_name.summary/unclear_hit_sequences.fasta $PTS/Samples/$specimen_name.summary/Unclear_Hits/Orientation_Check
#   cp $ref_database_for_tree $PTS/Samples/$specimen_name.summary/Unclear_Hits/BLAST
#   cp $ref_database $PTS/Samples/$specimen_name.summary/Unclear_Hits/Orientation_Check
# done

 

# #Check orientation of samples with unclear hits before collecting top 40 BLAST hits


# cat $PTE/FQC_$date/shortnames.txt | while read specimen_name;
# do
# #perform blast step make sure to add sstrand argument to see the orientation of the query sequence
#   cd $PTS/Samples/$specimen_name.summary/Unclear_Hits/Orientation_Check
#   makeblastdb -in $ref_database -dbtype nucl
#   blastn -db $ref_database -query unclear_hit_sequences.fasta -max_target_seqs 1 \
#   -out $PTS/Samples/$specimen_name.summary/Unclear_Hits/Orientation_Check/Unclear_sequences_orientation_blast_results -outfmt "6 qseqid pident sseqid qseq evalue sstrand"
#   awk '{ if ($0 ~ /_/) { printf ">"; } print $0; }' $PTS/Samples/$specimen_name.summary/Unclear_Hits/Orientation_Check/Unclear_sequences_orientation_blast_results > $PTS/Samples/$specimen_name.summary/Unclear_Hits/Orientation_Check/Unclear_sequences_orientation_blast_results.modified
#   mkdir $PTS/Samples/$specimen_name.summary/Unclear_Hits/Orientation_Check/fasta_split
#   cp $PTS/Samples/$specimen_name.summary/Unclear_Hits/Orientation_Check/Unclear_sequences_orientation_blast_results.modified $PTS/Samples/$specimen_name.summary/Unclear_Hits/Orientation_Check/fasta_split
# #Use awk script to split fasta file into multiple files each with one entry
#   cd $PTS/Samples/$specimen_name.summary/Unclear_Hits/Orientation_Check/fasta_split
#   awk 'BEGIN {n=0;} /^>/ {if(n%1==0){file=sprintf("chunk%d.fa",n);} print >> file; n++; next;} { print >> file; }' < Unclear_sequences_orientation_blast_results.modified
#   rm -rf Unclear_sequences_orientation_blast_results.modified
# #Establish a text file with the names of all the individual fasta files to use in a while read loop
#   cd $PTS/Samples/$specimen_name.summary/Unclear_Hits/Orientation_Check/fasta_split
#   ls | awk -F ".fa" '{print$1}' >> temp_fasta_names.txt
#   sed -e 's/temp//g' $PTS/Samples/$specimen_name.summary/Unclear_Hits/Orientation_Check/fasta_split/temp_fasta_names.txt > $PTS/Samples/$specimen_name.summary/Unclear_Hits/Orientation_Check/fasta_split/fasta_names.txt
#   rm -rf $PTS/Samples/$specimen_name.summary/Unclear_Hits/Orientation_Check/fasta_split/temp_fasta_names.txt
#   sed -i '.bak' '/^[[:space:]]*$/d' $PTS/Samples/$specimen_name.summary/Unclear_Hits/Orientation_Check/fasta_split/fasta_names.txt
#   mkdir $PTS/Samples/$specimen_name.summary/Unclear_Hits/Orientation_Check/fasta_split/correct_orientation
#   mkdir $PTS/Samples/$specimen_name.summary/Unclear_Hits/Orientation_Check/fasta_split/reversed_orientation
# #Find the orientation of each sequence and copy the reversed sequences to a new folder
#   cd $PTS/Samples/$specimen_name.summary/Unclear_Hits/Orientation_Check/fasta_split
#   cat fasta_names.txt | while read fasta_names;
#   do
#     orientation=`grep '>' $fasta_names.fa | awk '{print$6}'`
#     if [[ $orientation == "minus" ]]
#     then
#         cp $fasta_names.fa reversed_orientation
#         rm -rf $fasta_names.fa
#     fi
#   done
# #move remaining sequences into the correct orientation folder
#   if [ ! -z "$(ls -A $PTS/Samples/$specimen_name.summary/Unclear_Hits/Orientation_Check/fasta_split)" ]; then
#     cp $PTS/Samples/$specimen_name.summary/Unclear_Hits/Orientation_Check/fasta_split/chunk*.fa $PTS/Samples/$specimen_name.summary/Unclear_Hits/Orientation_Check/fasta_split/correct_orientation
#     rm -rf $PTS/Samples/$specimen_name.summary/Unclear_Hits/Orientation_Check/fasta_split/chunk*.fa
#   fi
# #Turn both the reverse orientation and correct orientation into fasta format, use seqtk to reverse the sequence of the necessary entries
#   mkdir $PTS/Samples/$specimen_name.summary/Unclear_Hits/Orientation_Check/fasta_split/fixed_orientation
#   cd $PTS/Samples/$specimen_name.summary/Unclear_Hits/Orientation_Check/fasta_split/reversed_orientation
#   cat $PTS/Samples/$specimen_name.summary/Unclear_Hits/Orientation_Check/fasta_split/fasta_names.txt | while read fasta_names;
#   do
#     awk 'BEGIN { OFS = "\n" } { print ">"$1, $4 }' < $fasta_names.fa > $fasta_names.fasta
#     seqtk seq -r $fasta_names.fasta > $fasta_names.reversed.fasta
#   done
# #Remove unncessary files created due to fasta_names text file
#   cp $PTS/Samples/$specimen_name.summary/Unclear_Hits/Orientation_Check/fasta_split/reversed_orientation/*.reversed.fasta $PTS/Samples/$specimen_name.summary/Unclear_Hits/Orientation_Check/fasta_split/fixed_orientation
#   cd $PTS/Samples/$specimen_name.summary/Unclear_Hits/Orientation_Check/fasta_split/fixed_orientation
#   cat *.fasta > $specimen_name.unclear_parasite_hits_fixed_orientation.fasta
#   cp $specimen_name.unclear_parasite_hits_fixed_orientation.fasta $PTS/Samples/$specimen_name.summary/Unclear_Hits/BLAST
#   cd $PTS/Samples/$specimen_name.summary
#   rm -rf unclear_hit_sequences.fasta
#   rm -rf unclear_seqs_fasta_table.csv
# #Split unclear fasta into individual fasta files
#   cd $PTS/Samples/$specimen_name.summary/Unclear_Hits/BLAST
#   awk 'BEGIN {n=0;} /^>/ {if(n%1==0){file=sprintf("chunk%d.fa",n);} print >> file; n++; next;} { print >> file; }' < $specimen_name.unclear_parasite_hits_fixed_orientation.fasta
#   rm -rf $PTS/Samples/$specimen_name.summary/Unclear_Hits/BLAST/$specimen_name.unclear_parasite_hits_fixed_orientation.fasta
#   for chunk in $PTS/Samples/$specimen_name.summary/Unclear_Hits/BLAST/*.fa
#     do
#     sed -i '.bak' 's/>>/> /g' $chunk
#     sample_name=`grep ">" $chunk | awk '{print$2}'`
#     sed 's/> />/g' $chunk > $PTS/Samples/$specimen_name.summary/Unclear_Hits/BLAST/$sample_name.fa
#     rm -rf $chunk.bak
#     rm -rf $chunk
#   done
#   cp $PTS/Samples/$specimen_name.summary/Unclear_Hits/BLAST/*.fa $PTS/Samples/$specimen_name.summary/Unclear_Hits/BLAST/tempFolder
#   cd $PTS/Samples/$specimen_name.summary/Unclear_Hits/BLAST/tempFolder
#   ls | awk -F ".fa" '{print$1}' >> $PTS/Samples/$specimen_name.summary/Unclear_Hits/BLAST/unclear_fasta_names.txt
# ##BLAST unclear fasta to find the best 30 hits by bit score to the reference database##
#   cd $PTS/Samples/$specimen_name.summary/Unclear_Hits/BLAST
#   path_to_clustering_blast=$PTS/Samples/$specimen_name.summary/Unclear_Hits/BLAST
#   path_to_clustering_tree=$PTS/Samples/$specimen_name.summary/Unclear_Hits/Clustering
#   cat $PTS/Samples/$specimen_name.summary/Unclear_Hits/BLAST/unclear_fasta_names.txt | while read unclear_sequence;
#   do
#     makeblastdb -in $ref_database_for_tree -dbtype nucl
#     blastn -query  $path_to_clustering_blast/$unclear_sequence.fa \
#     -db $ref_database_for_tree \
#     -outfmt "6 qseqid bitscore pident sseqid qseq evalue"  \
#     -evalue 0.001 \
#     -max_target_seqs 50 \
#     -out $path_to_clustering_blast/Results/$unclear_sequence.results
#     awk '{ if ($0 ~ /_/) { printf ">"; } print $0; }' $path_to_clustering_blast/Results/$unclear_sequence.results > $path_to_clustering_blast/Results/$unclear_sequence.results.modified
# #Call python script to build summary table and fasta table from BLAST results
#     cd $py_scripts
#     python3 Clustering.py \
#     --blast_output $path_to_clustering_blast/Results/$unclear_sequence.results.modified \
#     --unclear_sequence_BLAST_summary_table $path_to_clustering_blast/Results/$unclear_sequence.BLAST_Summary_Table.csv \
#     --unclear_sequence_fasta_table $path_to_clustering_blast/Results/$unclear_sequence.BLAST_Fasta_Table.csv
#     awk -F , '{print ">"$1"\n"$2}' $path_to_clustering_blast/Results/$unclear_sequence.BLAST_Fasta_Table.csv > $path_to_clustering_blast/Results/$unclear_sequence.BLAST.fasta
#     mkdir $path_to_clustering_tree/$unclear_sequence.clusters
#     cp $path_to_clustering_blast/Results/$unclear_sequence.BLAST.fasta $path_to_clustering_tree/$unclear_sequence.clusters
#     cp $path_to_clustering_blast/$unclear_sequence.fa $path_to_clustering_tree/$unclear_sequence.clusters
#     cd $path_to_clustering_tree/$unclear_sequence.clusters
#     cat * > $unclear_sequence.with_top_BLAST_results.fasta
#     cd $path_to_clustering_tree/$unclear_sequence.clusters
#     rm -rf $path_to_clustering_tree/$unclear_sequence.clusters/$unclear_sequence.BLAST.fasta
#     rm -rf $path_to_clustering_tree/$unclear_sequence.clusters/$unclear_sequence.fa
#   done
# done


echo "Activating conda R environment"
eval "$($(which conda) 'shell.bash' 'hook')"
#module load miniconda3
ENVS=$(conda env list | awk '{print$1}')

if [[ $ENVS = *updx_r_environment* ]]; then
	echo "UPDx R environment exists and is being activated"
  source activate $PTW/updx_r_environment
else
	echo "UPDx R environment does not exist. It will be created now"
	conda env create --prefix=$PTW/updx_r_environment --file=$PTW/linux_updx_r_environment.yml
  echo "UPDx R environment has been created and is being activated"
  conda activate $PTW/updx_r_environment
fi;

echo "Now entered into Conda R environment"




#Generate R-markdown summary report for each sample

mkdir $PTE/Important_Run_Information
cd $PTS
mkdir summary_reports

cat $PTE/FQC_$date/shortnames.txt | while read specimen_name;
do
  threshold_int=`cat $PTS/Samples/$specimen_name.summary/threshold_int.txt`
  echo $@ > $PTE/Important_Run_Information/Command_line_input.txt
  Rscript -e "rmarkdown::render('$PTW/summary_reports.Rmd',params=list(args = myarg), output_file = paste('$PTS/Samples/$specimen_name.summary/$specimen_name.summary_report.html'))" \
  $specimen_name \
  $PTS/Samples/$specimen_name.summary/results_table.csv \
  $PTS/Samples/$specimen_name.summary/markdown_summary_table.csv \
  $threshold_int \
  $PTE/Important_Run_Information/Command_line_input.txt \
  $PTS/Cutoff_Stats/temp_folder/parasites_used_for_cutoff.csv \
  $PTE/FQC_$date/fqc_output/raw_seqs_fqc_results/$specimen_name.rawR1_qualityCheck.txt \
  $PTE/FQC_$date/fqc_output/raw_seqs_fqc_results/$specimen_name.rawR2_qualityCheck.txt \
  $PTE/FQC_$date/fqc_output/trimmed_seqs_fqc_results/$specimen_name.cleanR1_qualityCheck.txt \
  $PTE/FQC_$date/fqc_output/trimmed_seqs_fqc_results/$specimen_name.cleanR2_qualityCheck.txt \
  $PTE/FQC_$date/fqc_output/trimmed_seqs_fqc_results/$specimen_name.clean1_fastqc.html \
  $PTE/FQC_$date/fqc_output/trimmed_seqs_fqc_results/$specimen_name.clean2_fastqc.html \
  $PTE/FQC_$date/fqc_output/post_trimming_calculations/ts_avg.txt \
  $PTE/FQC_$date/fqc_output/post_trimming_calculations/ts_stdev.txt \
  $PTS/Samples/$specimen_name.summary/Cluster_Summary_Table.csv \
  $PTS/Samples/$specimen_name.summary/exc_BLAST_Summary_Table.csv \
  $PTS/Samples/$specimen_name.summary/ref_BLAST_Summary_Table.csv \
  $PTS/Samples/$specimen_name.summary/Tree_Building/parasite_hits/Final_Tree/$specimen_name.sequences_with_reference_database_tree.jpg \
  $PTS/Samples/$specimen_name.summary/Tree_Building/parasite_hits/Final_Tree/$specimen_name.sequences_with_reference_database_tree.pdf
  cd $PTS/Samples/$specimen_name.summary
  cp $specimen_name.summary_report.html $PTS/summary_reports
  rm -rf $PTS/Samples/$specimen_name.summary/$specimen_name.summary_report.html
done

#Generate a r-markdown report for the positive parasite hits
cat $PTE/FQC_$date/shortnames.txt | while read specimen_name;
do
  Rscript -e "rmarkdown::render('$PTW/parasite_reports.Rmd',params=list(args = myarg), output_file = paste('$PTS/Samples/Parasite_Hits/parasite.summary_report.html'))" \
  $PTS/Samples/Parasite_Hits/Parasite_Summary_Table.csv
  cd $PTS/Samples/Parasite_Hits
  cp parasite.summary_report.html $PTS/summary_reports
  rm -rf $PTS/Samples/Parasite_Hits/parasite.summary_report.html
done




cd $PTS/Cutoff_Stats
rm -rf temp_folder

cd $PTW 

cp -r FINAL_run_summary_$date $PTE
rm -rf FINAL_run_summary_$date



Final_Time=$(date +"%H_%M_%S")

echo $Final_Time > $PTE/Important_Run_Information/End_of_run_timestamp.txt
echo $my_timestamp > $PTE/Important_Run_Information/Start_of_run_timestamp.txt
echo $Final_Blast_Time > $PTE/Important_Run_Information/End_of_Final_BLAST_timestamp.txt
echo "$USER" > $PTE/Important_Run_Information/User_for_this_run.txt

cp -r $PTE $PTR
rm -rf $PTE


  