
load "../lib/uses"

@test "uses finds camelcased useFacade function and runs it" {

    function useTest() {
        test=1
        echo "ran test"
    }

    run uses test

    [[ "$output" = "ran test" ]]

}