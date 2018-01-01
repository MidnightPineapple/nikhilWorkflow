#! bin/bash

#############DATA PROCESSING FUNCTIONS DOWN HERE####################

pretrim() {
    #req $FASTQ_DIR
    local currentFile=$1;
    bin/FastQC/fastqc           \
    $FASTQ_DIR/$currentFile     \
    --outdir="$RESULTS"/preTrimQC/ ;
}

trim() {
    #req $FASTQ_DIR
    local currentFile=$1;
    java -jar bin/Trimmomatic-0.36/trimmomatic-0.36.jar     \
    SE                                                      \
    $FASTQ_DIR/$currentFile                                 \
    "$RESULTS"/trim/$currentFile.trim                            \
    HEADCROP:13                                             ;
}

posttrim() {
    #req $FASTQ_DIR
    local currentFile=$1;
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
#    --sjdbOverhang 33

star1() {
    #requires $starGenome $FASTQ_DIR defined
    local currentFile=$1;
    bin/STAR-2.5.2b/bin/Linux_x86_64/STAR                   \
    --genomeDir $starGenome                                 \
    --alignIntronMax 10000                                  \
    --readFilesIn "$RESULTS"/trim/$currentFile.trim         \
    --outFileNamePrefix "$RESULTS"/STARp1/$currentFile.trim.  \
    --outSAMtype BAM Unsorted                               ;
}

star2() {
    #requires $starGenome $FASTQ_DIR $sjdbfiles $GTF defined
    local currentFile=$1;
    bin/STAR-2.5.2b/bin/Linux_x86_64/STAR                   \
    --genomeDir $starGenome                                 \
    --alignIntronMax 10000                                  \
    --readFilesIn "$RESULTS"/trim/$currentFile.trim         \
    --outFileNamePrefix "$RESULTS"/STARp2/$currentFile.trim.   \
    --outSAMunmapped Within                                 \
    --outSAMtype BAM SortedByCoordinate                     \
    --sjdbFileChrStartEnd ${sjdbFiles[@]}                   \
    --sjdbGTFfile $GTF                                      \
    --outFilterType BySJout                                 ;
}

picard() {
    local currentFile=$1;
    java -jar bin/picard.jar                                            \
    MarkDuplicates                                                      \
    I="$RESULTS"/STARp2/$currentFile.trim.Aligned.sortedByCoord.out.bam    \
    O="$RESULTS"/bam_drem/$currentFile.bam                                 \
    M="$RESULTS"/bam_drem/$currentFile.metrics.txt                         \
    REMOVE_DUPLICATES=true                                              \
    CREATE_INDEX=true                                                   ;
}

subread() {
    #req $GTF
    local currentFile=$1;
    bin/subread-1.5.1-source/bin/featureCounts  \
    -R                                          \
    -g gene_id                                  \
    -s 0                                        \
    -t exon                                     \
    -a $GTF                                     \
    -o "$RESULTS"/counts/$currentFile.count.txt    \
    "$RESULTS"/bam_drem/$currentFile.bam           ;
}
