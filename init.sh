#! bin/bash

# Find program folder
DIR="$( readlink -e "$( dirname "${BASH_SOURCE[0]}" )" )"

#loads functions into session
. "$DIR"/helper.sh
. "$DIR"/functions.sh

ORI=$(pwd)

#Get into project directory
echo
echo "Please input project directory."
echo "The tools necessary for this workflow should be installed within the scope of this directory"
echo "A results folder and a log folder will be created in this directory"
getDir P_DIR

if [[ ! -d $P_DIR/bin ]]; then
  echo "Please make sure the tools for this workflow are installed in the scope of this the project directory";
  exit;
fi

# checks if logs and results folders already exist and if so asks whether or not to replace
dirExists "$P_DIR" logs
dirExists "$P_DIR" results

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
