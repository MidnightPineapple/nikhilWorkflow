#! /usr/bin/env bats

load ../lib/helpers

@test "log function prints a log" {

    run log Hello World
    [[ $status -eq 0 ]]
    [[ $output =~ "LOG" ]]
    [[ $output =~ "Hello World" ]]

}

@test "warn function prints a warning" {

    run warn Hello World
    [[ $status -eq 0 ]]
    [[ $output =~ "WARN" ]]
    [[ $output =~ "Hello World" ]]

}

@test "error function prints an error and exits 1" {

    run error Hello World
    [[ $status -eq 1 ]]
    [[ $output =~ "ERROR" ]]
    [[ $output =~ "Hello World" ]]

}

