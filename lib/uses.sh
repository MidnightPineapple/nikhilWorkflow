uses() {
    for __facade in "$@"; do
        local __useFacade="use${__facade^}"
        if [[ -n "$(type -t $__useFacade)" ]] && [[ "$(type -t $__useFacade)" = "function" ]]; then
            $__useFacade
        else 
            error "$__facade is not a facade"
        fi
    done
}