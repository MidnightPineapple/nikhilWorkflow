#! /usr/bin/env bash

uses "pretrim" "trim" "posttrim" "generateStarGenome" "star1" "star2" "removeDuplicates" "countGenes" "limma"

########### PRETRIM ############
loopThru pretrim __groups 

########### TRIM #############
loopThru trim __groups 

########### POSTTRIM ###########
loopThru posttrim __groups 

########### STAR GENOME GENERATE ############
depend "starGenome" "outputDirectory"
if [[ ! -d $starGenome ]]; then
    generateStarGenome 
else
    log "STAR genome provided. Skipping STAR genomeGenerate"
fi

############ STAR PASS 1 ###############
loopThru star1 __groups 

########### STAR PASS 2 ##############
loopThru star2 __groups 

########## REMOVE DUPLICATES ############
loopThru removeDuplicates __groups 

######### COUNT FEATURES ###############
loopThru countGenes __groups

######### LIMMA VOOM ##############
limma