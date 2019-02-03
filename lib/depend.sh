#! usr/bin/env bash

if [[ -z "$__necessary_globals" ]]; then
    __necessary_globals=( )
fi

depend() {
    __necessary_globals=( "${__necessary_globals[@]}" "$@")
}

need() {
    for __varname in "$@"; do
        if [[ -z "${!__varname+x}" ]]; then
            echo "Variable $__varname not defined" 1>&2
            exit 1;
        fi
    done
}

validate() {
    for __name in "${__necessary_globals[@]}"; do
        need "$__name"
    done
}


