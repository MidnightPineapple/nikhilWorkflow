#! /usr/bin/env bash

usePretrim() {
    pretrim() {
        load fastqc

        makeDirectoryIfNotExists "$__results_directory"/preTrimQC

        local currentFilePath="$1"
        local currentFile="$(basename $currentFilePath)"

        log "Starting PretrimQC for $currentFile"
        iFastqc                                       \
        "$currentFilePath"                            \
        --outdir="$__results_directory"/preTrimQC/    ;
        log "Finished PretrimQC for $currentFile"
    }
}

useTrim() {    
    trim() {
        load trimmomatic

        makeDirectoryIfNotExists "$__results_directory"/trim

        local currentFilePath="$1";
        local currentFile="$(basename $currentFilePath)"
        log "Starting Trim SE HEADCROP:13 for $currentFile"
        iTrimmomatic                                            \
        SE                                                      \
        "$currentFilePath"                                      \
        "$__results_directory/trim/$currentFile.trim"           \
        HEADCROP:13                                             ;
        log "Finished Trim for $currentFile"
    }
}

usePosttrim() {
    posttrim() {
        load fastqc

        makeDirectoryIfNotExists "$__results_directory"/postTrimQC

        local currentFilePath="$1";
        local currentFile="$(basename $currentFilePath)"
        log "Starting PosttrimQC for $currentFile"
        iFastqc                                         \
        "$__results_directory/trim/$currentFile.trim"   \
        --outdir="$__results_directory"/postTrimQC/     ;
        log "Finished PosttrimQC for $currentFile"
    }
}

useGenerateStarGenome() {
    depend "referenceGenome" "featureAnnotationsFile" "outputDirectory"
    generateStarGenome() {
        load star

        makeDirectoryIfNotExists "$outputDirectory/StarGenome"

        log "Starting STAR genomeGenerate"
        iStar                                       \
        --runMode genomeGenerate                    \
        --genomeFastaFiles "$referenceGenome"       \
        --genomeDir "$outputDirectory/StarGenome"   \
        --sjdbGTFfile "$featureAnnotationsFile"     ;
        log "Finished STAR genomeGenerate"

        starGenome="$outputDirectory/StarGenome"
    }
    #    --sjdbOverhang 33                       ;
}

useStar1() {
    depend "starGenome"
    star1() {
        load star

        makeDirectoryIfNotExists "$__results_directory"/STARp1

        local currentFilePath="$1"
        local currentFile="$(basename $currentFilePath)"
        log "Starting STAR first pass for $currentFile"
        iStar                                                                   \
        --genomeDir "$starGenome"                                               \
        --alignIntronMax 5000                                                   \
        --readFilesIn "$__results_directory/trim/$currentFile.trim"             \
        --outFileNamePrefix "$__results_directory/STARp1/$currentFile.trim."    \
        --outSAMtype BAM Unsorted                                               ;
        log "Finished STAR first pass for $currentFile"
    }
}

useStar2() {
    depend "starGenome" "featureAnnotationsFile"
    star2() {
        need "groups"

        load star

        makeDirectoryIfNotExists "$__results_directory"/STARp2

        local currentFilePath="$1"
        local currentFile="$(basename $currentFilePath)"

        local IFS=$'\n'
        local sjdbFiles=( $(for file in "${groups[@]}"; do basename "$file"| sed -e "s,.*, $__results_directory/STARp1/&.trim.SJ.out.tab,"; done;) )
        unset IFS

        log "Starting STAR second pass for $currentFile"

        # ! Why tf does this fail if sjdbFiles has quotes around it...

        iStar                                                                   \
        --genomeDir "$starGenome"                                               \
        --alignIntronMax 5000                                                   \
        --readFilesIn "$__results_directory/trim/$currentFile.trim"             \
        --outFileNamePrefix "$__results_directory/STARp2/$currentFile.trim."    \
        --outSAMunmapped Within                                                 \
        --outSAMtype BAM SortedByCoordinate                                     \
        --sjdbFileChrStartEnd ${sjdbFiles[@]}                                   \
        --sjdbGTFfile "$featureAnnotationsFile"                                 \
        --outFilterType BySJout                                                 ;
        log "Finished STAR second pass for $currentFile"
    }
}

useRemoveDuplicates() {
    removeDuplicates() {
        load picard

        makeDirectoryIfNotExists "$__results_directory"/bam_drem

        local currentFilePath="$1"
        local currentFile="$(basename $currentFilePath)"
        log "Removing duplicates for $currentFile"
        iPicard                                                                 \
        MarkDuplicates                                                          \
        I="$__results_directory"/STARp2/$currentFile.trim.Aligned.sortedByCoord.out.bam     \
        O="$__results_directory"/bam_drem/$currentFile.bam                                  \
        M="$__results_directory"/bam_drem/$currentFile.metrics.txt                          \
        REMOVE_DUPLICATES=true                                                  \
        CREATE_INDEX=true                                                       ;
        log "Finished removing duplicates for $currentFile"
    }
}

useCountGenes() {
    depend "featureAnnotationsFile"
    countGenes() {
        load featureCounts

        makeDirectoryIfNotExists "$__results_directory"/counts

        local currentFilePath="$1"
        local currentFile="$(basename $currentFilePath)"
        log "Counting features for $currentFile"
        iFeatureCounts                              \
        -R BAM                                      \
        -g gene_id                                  \
        -s 0                                        \
        -t exon                                     \
        -a "$featureAnnotationsFile"                \
        -o "$__results_directory"/counts/$currentFile.count.txt    \
        "$__results_directory"/bam_drem/$currentFile.bam           ;
        log "Finished counting features for $currentFile"
    }
}

useLimma() {
    depend "goAnnotationsFile" "group1Name" "group2Name"
    limma() {
        need "group1" "group2"

        load r

        makeDirectoryIfNotExists "$__results_directory"/voom

        log "Starting Limma Voom Analysis"
        
        local IFS=$'\n'
        local group1CountsArray=($(for file in "${group1[@]}"; do basename "$file"| sed -e "s,.*,$__results_directory/counts/&.count.txt\n,"; done;))
        local group2CountsArray=($(for file in "${group2[@]}"; do basename "$file"| sed -e "s,.*,$__results_directory/counts/&.count.txt\n,"; done;))
        unset IFS

        local group1Counts=$(joinBy , "${group1CountsArray[@]}")
        local group2Counts=$(joinBy , "${group2CountsArray[@]}")

        iRscript "$__dirname"/scripts/limmavoom.R "$__results_directory/voom" "$group1Counts" "$group2Counts" "$goAnnotationsFile" "$group1Name" "$group2Name"
        
        log "Final results in $__results_directory/voom"

    }
}