#! /usr/bin/env bash

loopThru() {
    local operation="$1"; shift
    local arrayName="$1[@]"; shift
    local array=("${!arrayName}")
    local counter=0
    while [ $counter -lt ${#array[@]} ]; do
        log "Performing $operation iteration $(expr $counter + 1)"
        "$operation" "${array[$counter]}" "${@}"
        let counter=counter+1
    done
}

forEach() {
    local arrayName="$1[@]"; shift
    local operation="$1"; shift
    local array=("${!arrayName}")
    local counter=0
    while [ $counter -lt ${#array[@]} ]; do
        log "Performing $operation iteration $(expr $counter + 1)"
        "$operation" "${array[$counter]}" "${@}"
        let counter=counter+1
    done
}

joinBy() {
    local IFS="$1";
    shift;
    echo "$*";
}

makeDirectoryIfNotExists() {
    local dirPath="$1"
    if [[ ! -d "$dirPath" ]]; then
        mkdir -p "$dirPath"
    fi
}

error() {
    local message="$*"
    echo -e "\e[2;34m$(date) \e[22;1;91mERROR\e[0m: $message" 1>&2
    exit 1
}

log() {
    echo -e "\e[2;34m$(date) \e[22;1;94mLOG\e[0m: $@"
}

warn() {
    echo -e "\e[2;34m$(date) \e[22;1;93mWARN\e[0m: $@"
}

clink() {
    echo $'\U1F37B'
}

load() {
    if [[ "$#" -ne 1 ]]; then 
        error "load takes one argument"
    fi

    source $__dirname/lib/load.bash "$1"

    if [[ "$?" -ne 0 ]]; then
        error "Unable to load $1"
    fi
}

cleanup() {
    cd "$__wd"
    unset __should_run
}