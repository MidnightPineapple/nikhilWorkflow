#! /usr/bin/env bash

uses "pretrim" "trim" "posttrim" "generateStarGenome" "star1" "star2" "removeDuplicates" "countGenes" "limma"

depend "groups"

########### PRETRIM ############
loopThru pretrim groups 

########### TRIM #############
loopThru trim groups 

########### POSTTRIM ###########
loopThru posttrim groups 

########### STAR GENOME GENERATE ############
depend "starGenome" "outputDirectory"
if [[ ! -d $starGenome ]]; then
    generateStarGenome 
else
    log "STAR genome provided. Skipping STAR genomeGenerate"
fi

############ STAR PASS 1 ###############
loopThru star1 groups 

########### STAR PASS 2 ##############
loopThru star2 groups 

########## REMOVE DUPLICATES ############
loopThru removeDuplicates groups 

######### COUNT FEATURES ###############
loopThru countGenes groups

######### LIMMA VOOM ##############
limma

log "All Done! Cheers" $'\U1F37B'
