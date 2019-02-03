#! /usr/bin/env bash

log 'STARTING PRETRIMQC'
mkdir "$RESULTS"
mkdir "$RESULTS"/preTrimQC
loopThru GROUP_ALL pretrim >> "$LOGS"/preTrimQC.log.out 2>&1
log 'FINISHED PRETRIMQC'

log "STARTING TRIMMOMATIC WITH COMMAND HEADCROP 13"
mkdir "$RESULTS"/trim
loopThru GROUP_ALL trim >> "$LOGS"/trim.log.out 2>&1
log "FINISHED TRIMMING"

log 'STARTING POSTTRIMQC'
mkdir "$RESULTS"/postTrimQC
loopThru GROUP_ALL posttrim >> "$LOGS"/postTrimQC.log.out 2>&1
log 'FINISHED POSTTRIMQC'

if [[ ! -d $starGenome ]]; then
    log 'MAKING STARGENOME'
    mkdir "$P_DIR"/STARgenome
    genStarGenome >> "$LOGS"/genStarGenome.log.out 2>&1
    log 'FINISHED MAKING STARGENOME'
else
    log "STARGENOME REFERENCE PROVIDED, SKIPPING STARGENOME CREATION"
fi

log 'RUNNING STAR PASS 1'
mkdir "$RESULTS"/STARp1
loopThru GROUP_ALL star1 >> "$LOGS"/STARp1.log.out 2>&1
log 'FINISHED STAR PASS 1'

#makes array of all sjdb files
#formatStringArray sjdbFileString GROUP_ALL "$RESULTS/STARp1/" ".trim.SJ.out.tab"
sjdbFileString=($(for file in "${GROUP_ALL[@]}"; do basename "$file"| sed -e "s,.*, $RESULTS/STARp1/&.trim.SJ.out.tab,"; done;))
sjdbFiles=("${sjdbFileString[@]}")

log 'RUNNING STAR PASS 2'
mkdir "$RESULTS"/STARp2
loopThru GROUP_ALL star2 >> "$LOGS"/STARp2.log.out 2>&1
log 'FINISHED STAR PASS 2'

log 'RUNNING PICARD'
mkdir "$RESULTS"/bam_drem
loopThru GROUP_ALL picard >> "$LOGS"/picard.log.out 2>&1
log 'FINISHED PICARD'

log 'SUBREAD FEATURECOUNTS'
mkdir "$RESULTS"/counts
loopThru GROUP_ALL subread >> "$LOGS"/subread.log.out 2>&1
log 'FINISHED FEATURECOUNTS'

log 'LIMMA VOOM ANALYSIS'
mkdir "$RESULTS"/voom
module add R
export R_LIBS="$P_DIR/bin/R_libs/" #adds R_LIBS path as an environment variable

A_COUNTS=($(for file in "${GROUP_A[@]}"; do basename "$file"| sed -e "s,.*, $P_DIR/$RESULTS/counts/&.count.txt,"; done;))
B_COUNTS=($(for file in "${GROUP_B[@]}"; do basename "$file"| sed -e "s,.*, $P_DIR/$RESULTS/counts/&.count.txt,"; done;))
# formatStringArray A_COUNTS GROUP_A "$P_DIR/$RESULTS/counts/" ".count.txt"; A_COUNTS=($A_COUNTS)
# formatStringArray B_COUNTS GROUP_B "$P_DIR/$RESULTS/counts/" ".count.txt"; B_COUNTS=($B_COUNTS)
A_COUNTS_PATHS=$(joinBy , "${A_COUNTS[@]}")
B_COUNTS_PATHS=$(joinBy , "${B_COUNTS[@]}")

Rscript --vanilla "$DIR"/limmavoom.R "$P_DIR/$RESULTS/voom" "$A_COUNTS_PATHS" "$B_COUNTS_PATHS" "$GA" "$GROUP_A_NAME" "$GROUP_B_NAME" >> "$LOGS"/voom.log.out 2>&1
log "Final results in $P_DIR/$RESULTS/voom"

# removes intermediate results files if we want to reset
if [ "$RESET" = true ]; then
  log "Deleting intermediate files and folders in $P_DIR/$RESULTS"
  # IDEA: make a new make directory function within helper and keep track of the creation of each new results folder in an array so its more elegant
  rm -rf "$RESULTS"/trim "$RESULTS"/STARp1 "$RESULTS"/STARp2 "$RESULTS"/bam_drem "$RESULTS"/counts
fi

log "ALL DONE!"
