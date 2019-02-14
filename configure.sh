need "__config_filepath"

if [[ -f "$__config_filepath" ]]; then
    source "$__config_filepath"
else 
    error "Config file not found at $__config_filepath" 
fi

need "outputDirectory" "workflowFile"

# TODO: find abs path here

__log_directory="$outputDirectory/logs"
__results_directory="$outputDirectory/results"
makeDirectoryIfNotExists "$__results_directory"
makeDirectoryIfNotExists "$__log_directory"

IFS=$'\n'
for line in $(awk '$0 ~ "^(depend|uses)" { print }' $workflowFile); do
    unset IFS
    eval $line
    IFS=$'\n'
done
unset IFS

validate

__configure=true

