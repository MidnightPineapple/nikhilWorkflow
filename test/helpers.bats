#! /usr/bin/env bats

load ../lib/helpers

@test "log function prints a log" {

    run log "Hello World"
    [[ $status -eq 0 ]]
    [[ $output =~ Log.+Hello World ]]

}
