#! /bin/bash

# Initializing Params and Helper Functions
. init.sh

mkdir $P_DIR/logs
echo '';
echo 'Workflow starting in background.';
echo 'Logs will be stored in logs/log.txt';
echo 'If this script was started using the bash function';
echo 'then it will continue running after log out';
######################### START WORKFLOW ############################
{
# Make sure we're in the project directory before we actually start working...
cd $P_DIR

echo 'STARTING PRETRIMQC'
mkdir results
mkdir results/preTrimQC
loopThru GROUP_ALL pretrim 
echo 'FINISHED PRETRIMQC'

echo "STARTING TRIMMOMATIC WITH COMMAND HEADCROP 13"
loopThru GROUP_ALL trim
echo "FINISHED TRIMMING"

echo 'STARTING POSTTRIMQC'
mkdir results/postTrimQC
loopThru GROUP_ALL posttrim 
echo 'FINISHED POSTTRIMQC'

if [ ! -d $starGenome ]; then 
    echo 'MAKING STARGENOME'
    mkdir STARgenome
    starGenome=$P_DIR/STARgenome
    genStarGenome 
else
    echo "STARGENOME REFERENCE PROVIDED, SKIPPING STARGENOME CREATION"
fi

echo 'RUNNING STAR PASS 1'
mkdir results/STARp1
loopThru GROUP_ALL star1 
echo 'FINISHED STAR P1'

#makes array of all sjdb files
formatStringArray sjdbFileString GROUP_ALL "results/STARp1/" ".trim.SJ.out.tab"
sjdbFiles=($sjdbFileString)

echo 'RUNNING STAR P2'
mkdir results/STARp2
loopThru GROUP_ALL star2 
echo 'FINISHED STAR P2'

echo 'RUNNING PICARD' 
mkdir results/bam_drem
loopThru GROUP_ALL picard 
echo 'FINISHED PICARD'

echo 'SUBREAD FEATURECOUNTS'
mkdir results/counts
loopThru GROUP_ALL subread 
echo 'FINISHED FEATURECOUNTS'

echo 'LIMMA VOOM ANALYSIS'
mkdir results/voom
module add R
export R_LIBS=$P_DIR"/bin/R_libs/" #adds R_LIBS path as an environment variable

formatStringArray A_COUNTS GROUP_A $P_DIR"/results/counts/" ".count.txt"; A_COUNTS=($A_COUNTS)
formatStringArray B_COUNTS GROUP_B $P_DIR"/results/counts/" ".count.txt"; B_COUNTS=($B_COUNTS)
A_COUNTS_PATHS=$(joinBy , "${A_COUNTS[@]}")
B_COUNTS_PATHS=$(joinBy , "${B_COUNTS[@]}")

Rscript --vanilla limmavoom.R "$P_DIR/results/voom" "$A_COUNTS_PATHS" "$B_COUNTS_PATHS" "$GA"
} &> logs/log.txt &