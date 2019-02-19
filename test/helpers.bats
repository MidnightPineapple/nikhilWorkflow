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

@test "joinby can join a list of parameters with some char" {

    run joinBy "," "wa" "ha" "ha"
    [[ $output = "wa,ha,ha" ]]

    run joinBy ":" "wa" "ha" "ha"
    [[ $output = "wa:ha:ha" ]]

}

@test "joinby can join an array expansion with some char" {

    local arr=( "wa" "ha" "ha" )
    run joinBy "," "${arr[@]}"    
    [[ $output = "wa,ha,ha" ]]

    run joinBy ":" "${arr[@]}"    
    [[ $output = "wa:ha:ha" ]]

}

@test "loopThru goes through each element of the array and performs some operation" {

    local arr=( 1 2 3 )
    run loopThru echo arr "foo" "bar" "baz"
    [[ "${lines[0]}" =~ "iteration 1" ]]
    [[ "${lines[1]}" = "1 foo bar baz" ]]
    [[ "${lines[2]}" =~ "iteration 2" ]]
    [[ "${lines[3]}" = "2 foo bar baz" ]]
    [[ "${lines[4]}" =~ "iteration 3" ]]
    [[ "${lines[5]}" = "3 foo bar baz" ]]

}

@test "forEach goes through each element of the array and performs some operation" {

    local arr=( 1 2 3 )
    run forEach arr echo "foo" "bar" "baz"
    [[ "${lines[0]}" =~ "iteration 1" ]]
    [[ "${lines[1]}" = "1 foo bar baz" ]]
    [[ "${lines[2]}" =~ "iteration 2" ]]
    [[ "${lines[3]}" = "2 foo bar baz" ]]
    [[ "${lines[4]}" =~ "iteration 3" ]]
    [[ "${lines[5]}" = "3 foo bar baz" ]]

}

@test "clink prints a beer mug" {
    
    local mug=$'\U1F37B'
    run clink
    [[ $output = $mug ]]

}