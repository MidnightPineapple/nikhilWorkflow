#! /usr/bin/env bash

loopThru() {
    local operation="$1"; shift
    local arrayName="$1[@]"; shift
    local array=("${!arrayName}")
    local counter=0
    while [ $counter -lt ${#array[@]} ]; do
        echo "Performing operation iteration number $(expr $counter + 1)"
        "$operation" "${array[$counter]}" "${@}"
        let counter=counter+1
    done
}

makeDirectoryIfNotExists() {
    local dirPath="$1"
    if [[ ! -d "$dirPath" ]]; then
        mkdir -p "$dirPath"
    fi
}

error() {
    local message="$1"
    echo -e "\e[2;34m$(date) \e[22;1;91mERROR\e[0m: $message" 1>&2
    exit 1
}

log() {
    echo -e "\e[2;34m$(date) \e[22;1;94mLOG\e[0m: $@"
}

warn() {
    echo -e "\e[2;34m$(date) \e[22;1;93mWARN\e[0m: $@"
}