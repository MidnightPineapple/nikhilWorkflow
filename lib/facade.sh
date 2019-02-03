#! usr/bin/env bash

pretrim() {
    load fastqc

    local currentFilePath="$1"
    iFastqc                                       \
    "$currentFilePath"                            \
    --outdir="$__results_directory"/preTrimQC/    ;
}

trim() {
    load trimmomatic

    local currentFilePath="$1";
    local currentFile="$(basename $currentFilePath)"
    iTrimmomatic                                            \
    SE                                                      \
    "$currentFilePath"                                      \
    "$__results_directory/trim/$currentFile.trim"           \
    HEADCROP:13                                             ;
}

posttrim() {
    load fastqc

    local currentFilePath="$1";
    local currentFile="$(basename $currentFilePath)"
    iFastqc                                         \
    "$__results_directory/trim/$currentFile.trim"   \
    --outdir="$__results_directory"/postTrimQC/     ;
}

depend "referenceGenome" "featureAnnotationsFile" "outputDirectory"
generateStarGenome() {
    load star
    iStar                                       \
    --runMode genomeGenerate                    \
    --genomeFastaFiles "$referenceGenome"       \
    --genomeDir "$outputDirectory/StarGenome"   \
    --sjdbGTFfile "$featureAnnotationsFile"     ;
}
#    --sjdbOverhang 33                       ;

depend "starGenome"
star1() {
    load star

    local currentFilePath="$1"
    local currentFile="$(basename $currentFilePath)"
    iStar                                                                   \
    --genomeDir "$starGenome"                                               \
    --alignIntronMax 5000                                                   \
    --readFilesIn "$__results_directory/trim/$currentFile.trim"             \
    --outFileNamePrefix "$__results_directory/STARp1/$currentFile.trim."    \
    --outSAMtype BAM Unsorted                                               ;
}

depend "starGenome" "featureAnnotationsFile"
star2() {
    load star
    need "sjdbFiles" # ! currently we dont support testing if array
    # TODO: I should compute the filepaths for sjdbFiles inside here
    # means I'll need a standard global var for holding all sample files we're interested in 

    local currentFilePath="$1"
    local currentFile="$(basename $currentFilePath)"
    iStar                                                                   \
    --genomeDir "$starGenome"                                               \
    --alignIntronMax 5000                                                   \
    --readFilesIn "$__results_directory/trim/$currentFile.trim"             \
    --outFileNamePrefix "$__results_directory/STARp2/$currentFile.trim."    \
    --outSAMunmapped Within                                                 \
    --outSAMtype BAM SortedByCoordinate                                     \
    --sjdbFileChrStartEnd "${sjdbFiles[@]}"                                 \
    --sjdbGTFfile "$featureAnnotationsFile"                                 \
    --outFilterType BySJout                                                 ;
}

removeDuplicates() {
    load picard

    local currentFilePath="$1"
    local currentFile="$(basename $currentFilePath)"
    iPicard                                                                 \
    MarkDuplicates                                                          \
    I="$__results_directory"/STARp2/$currentFile.trim.Aligned.sortedByCoord.out.bam     \
    O="$__results_directory"/bam_drem/$currentFile.bam                                  \
    M="$__results_directory"/bam_drem/$currentFile.metrics.txt                          \
    REMOVE_DUPLICATES=true                                                  \
    CREATE_INDEX=true                                                       ;
}

depend "featureAnnotationsFile"
countGenes() {
    load subread

    local currentFilePath="$1"
    local currentFile="$(basename $currentFilePath)"
    iFeatureCounts                              \
    -R                                          \
    -g gene_id                                  \
    -s 0                                        \
    -t exon                                     \
    -a "$featureAnnotationsFile"                \
    -o "$__results_directory"/counts/$currentFile.count.txt    \
    "$__results_directory"/bam_drem/$currentFile.bam           ;
}
