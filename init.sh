#! bin/bash

#loads functions into session
. helper.sh
. functions.sh

ORI=$(pwd)

#Get into project directory
echo "Please input project directory."
getDir P_DIR

# ask if they have a starGenome already; if so ask for stargenome directory
while true; do
    read -p "Does STAR Genome already exist? " yn
    case $yn in
        [Yy1]* ) echo "Please Input STARgenome directory."; getDir starGenome; break;;
        [Nn0]* ) starGenome=0;break;;
        * ) echo "Please answer yes or no. ";;
    esac
done

# get reference genome path and gtf file path
if [ ! -d $starGenome ]; then 
    echo "Please input path to reference genome .fasta file."; 
    getFile REF_GENOME; REF_GENOME=$(readlink -e $REF_GENOME);
fi;
echo "Please input path to GTF file."; 
getFile GTF; GTF=$(readlink -e $GTF);
# get gene associations file path
echo "Please input path to gene associations file."; 
getFile GA; GA=$(readlink -e $GA);

# get fastq folder path
echo "Please input path to fastq folder. This folder should contain all the fastq files for both groups. "; getDir FASTQ_DIR;
# get group 1 and 2 fastq filenames
cd $FASTQ_DIR
echo "Please input the first filename within group A. ";
getFileArray GROUP_A; GROUP_A=($GROUP_A);
echo "Please input the first filename within group B. ";
getFileArray GROUP_B; GROUP_B=($GROUP_B);

GROUP_ALL=(${GROUP_A[@]} ${GROUP_B[@]})

cd $ORI