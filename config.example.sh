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
group1=( $(readlink -e "$group1DirectoryAbsolute"/* | grep .fastq) )
group2=( $(readlink -e "$group2DirectoryAbsolute"/* | grep .fastq) )
unset IFS
groups=( "${group1[@]}" "${group2[@]}" )