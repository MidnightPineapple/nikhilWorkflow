#! bin/bash

#############DATA PROCESSING FUNCTIONS DOWN HERE####################

pretrim() {
    #req $RESULTS
    local currentFilePath=$1;
    local currentFile=$(basename $currentFilePath);
    bin/FastQC/fastqc              \
    $currentFilePath               \
    --outdir="$RESULTS"/preTrimQC/ ;
}

trim() {
    #req $RESULTS
    local currentFilePath=$1;
    local currentFile=$(basename $currentFilePath);
    java -jar bin/Trimmomatic-0.36/trimmomatic-0.36.jar     \
    SE                                                      \
    $currentFilePath                                        \
    "$RESULTS"/trim/$currentFile.trim                       \
    HEADCROP:13                                             ;
}

posttrim() {
    #req $RESULTS
    local currentFilePath=$1;
    local currentFile=$(basename $currentFilePath);
    bin/FastQC/fastqc                    \
    "$RESULTS"/trim/$currentFile.trim    \
    --outdir="$RESULTS"/postTrimQC/    ;
}

genStarGenome() {
    #req $REF_GENOME $GTF
    bin/STAR-2.5.2b/bin/Linux_x86_64/STAR   \
    --runMode genomeGenerate                \
    --genomeFastaFiles $REF_GENOME          \
    --genomeDir STARgenome                  \
    --sjdbGTFfile $GTF                      ;
}
#    --sjdbOverhang 33                       ;

star1() {
    #requires $starGenome $RESULTS defined
    local currentFilePath=$1;
    local currentFile=$(basename $currentFilePath);
    bin/STAR-2.5.2b/bin/Linux_x86_64/STAR                      \
    --genomeDir $starGenome                                    \
    --alignIntronMax 10000                                     \
    --readFilesIn "$RESULTS"/trim/$currentFile.trim            \
    --outFileNamePrefix "$RESULTS"/STARp1/$currentFile.trim.   \
    --outSAMtype BAM Unsorted                                  ;
}

star2() {
    #requires $starGenome $RESULTS $sjdbfiles $GTF defined
    local currentFilePath=$1;
    local currentFile=$(basename $currentFilePath);
    bin/STAR-2.5.2b/bin/Linux_x86_64/STAR                      \
    --genomeDir $starGenome                                    \
    --alignIntronMax 10000                                     \
    --readFilesIn "$RESULTS"/trim/$currentFile.trim            \
    --outFileNamePrefix "$RESULTS"/STARp2/$currentFile.trim.   \
    --outSAMunmapped Within                                    \
    --outSAMtype BAM SortedByCoordinate                        \
    --sjdbFileChrStartEnd ${sjdbFiles[@]}                      \
    --sjdbGTFfile $GTF                                         \
    --outFilterType BySJout                                    ;
}

picard() {
    #requires $RESULTS
    local currentFilePath=$1;
    local currentFile=$(basename $currentFilePath);
    java -jar bin/picard.jar                                                \
    MarkDuplicates                                                          \
    I="$RESULTS"/STARp2/$currentFile.trim.Aligned.sortedByCoord.out.bam     \
    O="$RESULTS"/bam_drem/$currentFile.bam                                  \
    M="$RESULTS"/bam_drem/$currentFile.metrics.txt                          \
    REMOVE_DUPLICATES=true                                                  \
    CREATE_INDEX=true                                                       ;
}

subread() {
    #req $GTF $RESULTS
    local currentFilePath=$1;
    local currentFile=$(basename $currentFilePath);
    bin/subread-1.5.1-source/bin/featureCounts  \
    -R                                          \
    -g gene_id                                  \
    -s 0                                        \
    -t exon                                     \
    -a $GTF                                     \
    -o "$RESULTS"/counts/$currentFile.count.txt    \
    "$RESULTS"/bam_drem/$currentFile.bam           ;
}
