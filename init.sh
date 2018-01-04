#! bin/bash

# TODO: use getopts to allow this whole workflow to be run with a one line command


#Get into project directory
echo "Please input project directory."
echo "The tools necessary for this workflow should be installed within the scope of this directory"
echo "A results folder and a log folder will be created in this directory"
getDir P_DIR

if [[ ! -d $P_DIR/bin ]]; then
  echo "bin folder not found. Please make sure the tools for this workflow are installed in the scope of this the project directory";
  exit;
fi

# ask whether we wanna keep the intermediate files after the workflow is done
while true; do
    read -p "Should intermediate files be deleted after the workflow is done? " yn
    case $yn in
        [Yy1]* ) RESET=true; break;;
        [Nn0]* ) RESET=false; break;;
        * ) echo "Please answer yes or no. ";;
    esac
done

# ask if they have a starGenome already; if so ask for stargenome directory
while true; do
    read -p "Does STAR Genome already exist? " yn
    case $yn in
        [Yy1]* ) echo "Please Input STAR Genome directory."; getDir starGenome; break;;
        [Nn0]* ) starGenome="$P_DIR"/STARgenome; break;;
        * ) echo "Please answer yes or no. ";;
    esac
done

# get reference genome path and gtf file path
if [ ! -d $starGenome ]; then
    echo "STAR Genome will be created in project directory."
    echo "Please input path to reference genome .fasta file.";
    getFile REF_GENOME; REF_GENOME=$(readlink -e $REF_GENOME);
fi;
echo "Please input path to GTF file.";
getFile GTF; GTF=$(readlink -e $GTF);
# get gene associations file path
echo "Please input path to gene associations file.";
getFile GA; GA=$(readlink -e $GA);

# allow user to save configuration
while true; do
    read -p "Would you like to save the current configuration? " yn
    case $yn in
        [Yy1]* )
          printf "P_DIR=\"$P_DIR\"\nstarGenome=\"$starGenome\"\nREF_GENOME=\"$REF_GENOME\"\nGTF=\"$GTF\"\nGA=\"$GA\"\nRESET=\"$RESET\"\n" > "$DIR"/.config;
          echo "Saved ";
          break;;
        [Nn0]* ) break;;
        * ) echo "Please answer yes or no. ";;
    esac
done


# ## get fastq folder path
# echo "Please input path to fastq folder. This folder should contain all the fastq files for both groups. "; getDir FASTQ_DIR;
# # get group 1 and 2 fastq filenames
# cd $FASTQ_DIR
# read -p "Please input the name of group A. " GROUP_A_NAME;
# echo "Please input the first filename within group $GROUP_A_NAME. ";
# getFileArray GROUP_A; GROUP_A=($GROUP_A);
# read -p "Please input the name of group B. " GROUP_B_NAME;
# echo "Please input the first filename within group $GROUP_B_NAME. ";
# getFileArray GROUP_B; GROUP_B=($GROUP_B);
#
# GROUP_ALL=(${GROUP_A[@]} ${GROUP_B[@]})
