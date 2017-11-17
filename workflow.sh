#! /bin/bash

# clears the console
clear

# checks if the file was sourced
if [[ ${BASH_SOURCE[0]} != $0 ]]; then
  echo "This script is being sourced, so it may not run properly.";
  echo
  echo "Please quit and run this workflow within a subshell."
  echo
  echo "Use the command: 'bash path/to/workflow.sh'"
  echo
  read -n 1 -s -r -p "Press CTRL+C to quit or press any key to continue sourcing."
  echo
fi

# Find program folder
DIR="$( readlink -e "$( dirname "${BASH_SOURCE[0]}" )" )"
# Initializing Params and Helper Functions
. "$DIR"/init.sh

mkdir "$P_DIR"/logs
echo
echo "Workflow starting in background with PID: ${$}";
echo 'Logs will be stored in logs/log.txt';
echo 'If this script was started using the bash function';
echo 'then it will continue running after log out';
######################### START WORKFLOW ############################
{
# Make sure we're in the project directory before we actually start working...
cd "$P_DIR"

echo $(date): 'STARTING PRETRIMQC'
mkdir results
mkdir results/preTrimQC
loopThru GROUP_ALL pretrim > logs/preTrimQC.log.out
echo $(date): 'FINISHED PRETRIMQC'

echo $(date): "STARTING TRIMMOMATIC WITH COMMAND HEADCROP 13"
loopThru GROUP_ALL trim > logs/trim.log.out
echo $(date): "FINISHED TRIMMING"

echo $(date): 'STARTING POSTTRIMQC'
mkdir results/postTrimQC
loopThru GROUP_ALL posttrim > logs/postTrimQC.log.out
echo $(date): 'FINISHED POSTTRIMQC'

if [[ ! -d $starGenome ]]; then
    echo $(date): 'MAKING STARGENOME'
    mkdir STARgenome
    starGenome="$P_DIR"/STARgenome
    genStarGenome > logs/genStarGenome.log.out
else
    echo $(date): "STARGENOME REFERENCE PROVIDED, SKIPPING STARGENOME CREATION"
fi

echo $(date): 'RUNNING STAR PASS 1'
mkdir results/STARp1
loopThru GROUP_ALL star1 > logs/STARp1.log.out
echo $(date): 'FINISHED STAR P1'

#makes array of all sjdb files
formatStringArray sjdbFileString GROUP_ALL "results/STARp1/" ".trim.SJ.out.tab"
sjdbFiles=($sjdbFileString)

echo $(date): 'RUNNING STAR P2'
mkdir results/STARp2
loopThru GROUP_ALL star2 > logs/STARp2.log.out
echo $(date): 'FINISHED STAR P2'

echo $(date): 'RUNNING PICARD'
mkdir results/bam_drem
loopThru GROUP_ALL picard > logs/picard.log.out
echo $(date): 'FINISHED PICARD'

echo $(date): 'SUBREAD FEATURECOUNTS'
mkdir results/counts
loopThru GROUP_ALL subread > logs/subread.log.out
echo $(date): 'FINISHED FEATURECOUNTS'

echo $(date): 'LIMMA VOOM ANALYSIS'
mkdir results/voom
module add R
export R_LIBS="$P_DIR/bin/R_libs/" #adds R_LIBS path as an environment variable

formatStringArray A_COUNTS GROUP_A "$P_DIR/results/counts/" ".count.txt"; A_COUNTS=($A_COUNTS)
formatStringArray B_COUNTS GROUP_B "$P_DIR/results/counts/" ".count.txt"; B_COUNTS=($B_COUNTS)
A_COUNTS_PATHS=$(joinBy , "${A_COUNTS[@]}")
B_COUNTS_PATHS=$(joinBy , "${B_COUNTS[@]}")

Rscript --vanilla "$DIR"/limmavoom.R "$P_DIR/results/voom" "$A_COUNTS_PATHS" "$B_COUNTS_PATHS" "$GA" > logs/voom.log.out
echo $(date): 'ALL DONE!'
} &> logs/log.out &
