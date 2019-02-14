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
        eval "$__name"="\"$(readlink -e "${!__name}")\""
    done
}

finalize() {
    validate
    absolute
}