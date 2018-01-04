#! /bin/bash
clear

# welcome message
cat << EOF
+----------------------------------------------------------------+
|                                                                |
|                                                                |
|          Differential Expression Analysis Workflow             |
|                  By Sam Chen & Nikhil Gowda                    |
|                                                                |
|                                                                |
+----------------------------------------------------------------+

EOF

# checks if the file was sourced
if [[ ${BASH_SOURCE[0]} != $0 ]]; then

cat << EOF
This script is being sourced, so it may not run properly.
Please quit and run this workflow within a subshell.
Use the command: 'bash path/to/run.sh'
EOF

  read -n 1 -s -r -p "Press CTRL+C to quit or press any key to continue sourcing."
  echo
fi

# Find directory program was started from
ORI=$(pwd)
# Find program folder
DIR="$( readlink -e "$( dirname "${BASH_SOURCE[0]}" )" )"
# Sourcing helpers and Functions
. "$DIR"/helper.sh
. "$DIR"/functions.sh

# Ask whether we wanna use default configs or set anew
# P_DIR, starGenome, REF_GENOME, GTF, GA, RESET
if [[ -f $DIR/.config ]]; then
  while true; do
      read -p "Use saved configuration? " yn
      case $yn in
          [Yy1]* ) . "$DIR"/.config; break;;
          [Nn0]* ) . "$DIR"/init.sh; break;;
          * ) echo "Please answer yes or no. ";;
      esac
  done
else
  . "$DIR"/init.sh
fi

# get group A and group B files by folder
echo
read -p "Please input the name of group A. " GROUP_A_NAME;
read -p "Please input the name of group B. " GROUP_B_NAME;
echo "Please make sure $GROUP_A_NAME and $GROUP_B_NAME fastq files are separated into different folders."
echo "Please input group $GROUP_A_NAME directory path"
getDir GROUP_A_DIR;
IFS=$'\n'
GROUP_A=($(readlink -e "$GROUP_A_DIR"/* | grep .fastq$));
unset IFS
echo "Found ${#GROUP_A[@]} fastq files."
echo
echo "Please input group $GROUP_B_NAME directory path"
getDir GROUP_B_DIR;
IFS=$'\n'
GROUP_B=($(readlink -e "$GROUP_B_DIR"/* | grep .fastq$));
unset IFS
echo "Found ${#GROUP_B[@]} fastq files."
GROUP_ALL=(${GROUP_A[@]} ${GROUP_B[@]})

# make results and log folder names particular to this run
# c=$(echo $b | tr -cd [a-z0-9\-\_])
LOGS="$( echo "$GROUP_A_NAME-$GROUP_B_NAME"_logs | tr -cd [A-Za-z0-9\-\_])"
RESULTS="$( echo "$GROUP_A_NAME-$GROUP_B_NAME"_results | tr -cd [A-Za-z0-9\-\_])"

# checks if logs and results folders already exist and if so asks whether or not to replace
dirExists "$P_DIR" "$LOGS"
dirExists "$P_DIR" "$RESULTS"


mkdir "$P_DIR"/"$LOGS"
cat << EOF
Workflow starting in background with PID: ${$}
Logs will be stored in $LOGS/log.txt
If this script was started using the bash function
then it will continue running after log out
EOF

. "$P_DIR"/workflow.sh

cd $ORI
