need "__config_filepath"
depend "outputDirectory" "workflowFile"

if [[ -f "$__config_filepath" ]]; then
    source "$__config_filepath"
else 
    error "Config file not found at $__config_filepath" 
fi

__log_directory="$outputDirectory/logs"
__results_directory="$outputDirectory/results"

makeDirectoryIfNotExists "$__results_directory"
makeDirectoryIfNotExists "$__log_directory"

validate

__configure=true

