#! /usr/bin/env bash

need "__config_filepath" "__workflow_filepath"

if [[ -f "$__config_filepath" ]]; then
    source "$__config_filepath"
else 
    error "Config file not found at $__config_filepath" 
fi

if [[ -f "$__workflow_filepath" ]]; then
    __workflow_filepath="$(readlink -e "$__workflow_filepath")"
else 
    error "Workflow file not found at $__workflow_filepath"
fi

need "outputDirectory" 

makeDirectoryIfNotExists "$outputDirectory"
outputDirectory="$(readlink -e "$outputDirectory")"

__log_file="$outputDirectory/log.txt"
__results_directory="$outputDirectory/results"
makeDirectoryIfNotExists "$__results_directory"

IFS=$'\n'
for line in $(awk '$0 ~ "^(depend|uses)" { print }' "$__workflow_filepath"); do
    unset IFS
    eval $line
    IFS=$'\n'
done
unset IFS

need "group1Directory" "group2Directory"

IFS=$'\n'
group1=( $(readlink -e "$group1Directory"/* | grep .fastq) )
group2=( $(readlink -e "$group2Directory"/* | grep .fastq) )
unset IFS
groups=( "${group1[@]}" "${group2[@]}" )

finalize 
__configure=true