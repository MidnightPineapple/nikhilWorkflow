#! /usr/bin/env bash

if [[ -z "$__necessary_globals" ]]; then
    __necessary_globals=( )
fi

depend() {
    __necessary_globals=( "${__necessary_globals[@]}" "$@")
}

need() {
    for __varname in "$@"; do
        if [[ -z "${!__varname+x}" ]]; then
            error "Variable $__varname not defined" 
        fi
    done
}

validate() {
    for __name in "${__necessary_globals[@]}"; do
        need "$__name"
    done
}

absolute() {
    for __name in "${__necessary_globals[@]}"; do
        if [[ -d "${!__name}" ]] || [[ -f "${!__name}" ]]; then 
            eval "$__name"="\"$(readlink -e "${!__name}")\""
            log "$__name is set to ${!__name}"
        else 
            warn "$__name is not a file or directory"
        fi
    done
}

finalize() {
    validate
    absolute
}