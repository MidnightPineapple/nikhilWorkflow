#! /usr/bin/env bash

# A facade for the loading and aliasing of modules

usageLoad() { echo "Usage: $0 {fastqc|trimmomatic|star|picard|featureCounts|r|salmon}"; }

if [[ $# -ne 1 ]]; then 
    echo "Invalid arguments: $@"
    usage
fi

__module="$1"

case "$__module" in 
fastqc)
    module load FastQC/0.11.7-Java-1.8.0_162
    iFastqc() { fastqc "$@" ; }
    ;;
[Tt]rimmomatic)
    module load Trimmomatic/0.38-Java-1.8.0_162
    iTrimmomatic() { java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.38.jar "$@" ; }
    ;;
[Ss][Tt][Aa][Rr])
    module load GNU/6.4.0-2.28  OpenMPI/2.1.2 STAR/2.6.0c
    iStar() { STAR "$@" ; }
    ;;
picard) 
    module load picard/2.18.1-Java-1.8.0_152
    iPicard() { java -jar $EBROOTPICARD/picard.jar "$@" ; }
    ;;
feature[Cc]ounts)
    module load Subread/1.6.2
    iFeatureCounts() { featureCounts "$@" ; }
    ;;
[Rr])
    module load GNU/7.3.0-2.30  OpenMPI/3.1.1 R/3.5.1-X11-20180604
    iR() { R "$@" ; }
    iRscript() { Rscript --vanilla "$@" ; }
    ;;
salmon)
    module load Salmon/0.11.3
    iSalmon() { salmon "$@" ; }
    ;;
*)
    usageLoad
    echo "$__module not found."
    return 1
    ;;
esac

unset __module

return 0

