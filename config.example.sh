outputDirectory=
workflowFile=
referenceGenome=
featureAnnotationsFile=
goAnnotationsFile=
starGenome=

group1Directory=
group2Directory=

group1DirectoryAbsolute="$(readlink -e $group1Directory)"
group2DirectoryAbsolute="$(readlink -e $group2Directory)"
group1Name="$(basename $group1DirectoryAbsolute)"
group2Name="$(basename $group2DirectoryAbsolute)"
IFS=$'\n'
group1=("$(readlink -e "$__group_1_path/*" | grep .fastq)")
group2=("$(readlink -e "$__group_2_path/*" | grep .fastq)")
unset IFS
groups=( "${__group_1[@]}" "${__group_2[@]}" )