#! bin/bash

####################HELPER FUNCTIONS####################3

getDir() {
    local __outVar=$1;
    local DIR=0;
    while [ ! -d $DIR ]; do
        read -p "" DIR;
        if [ ! -d $DIR ]; then echo "Directory does not exist. Please try again. "; fi;
    done;
    DIR=$(readlink -e $DIR)
    eval $__outVar="'$DIR'"
}

getFile() {
    local __outVar=$1;
    local FILE=0;
    while [ ! -f $FILE ]; do
        read -p "" FILE;
        if [ ! -f $FILE ]; then echo "File does not exist. Please try again. "; fi;
    done;
    eval $__outVar="'$FILE'"
}
#FIXME: If user submits more than 1 file at once, it won't care if the file actually exists!
#FIXME: If user puts in something other than y or n, it will treat it as a yes and no more chances to answer yn prob b/c it encounters getFile first.
getFileArray() {
    local __outVar=$1;
    local ARRAY=();
    while true; do
        getFile newFile
        ARRAY=("${ARRAY[@]}" $newFile)
        read -p "Do you have more files? " yn
        case $yn in
            [Yy1]* ) echo "Please input next filename. ";;
            [Nn0]* ) break;;
            * ) echo "Please answer yes or no. ";;
        esac
    done
    eval $__outVar="'${ARRAY[@]}'"
}

formatStringArray() {
    local __outVar=$1
    local ARRAY_NAME=$2[@]
    local INPUT_ARRAY=("${!ARRAY_NAME}")
    local FRONT=$3
    local BACK=$4
    local ARRAY=();
    for i in "${INPUT_ARRAY[@]}"; do
        ARRAY=("${ARRAY[@]}" $FRONT$i$BACK);
    done;
    eval $__outVar="'${ARRAY[@]}'"
}

joinBy() {
    local IFS="$1";
    shift;
    echo "$*";
}

loopThru() {
    local ARRAY_NAME=$1[@];
    local ARRAY=("${!ARRAY_NAME}")
    local COUNTER=0;
    while [ $COUNTER -lt ${#ARRAY[@]} ]; do
        echo "Performing operation iteration number `expr $COUNTER + 1`";
        "$2" "${ARRAY[$COUNTER]}";
        let COUNTER=COUNTER+1;
    done;
}

dirExists() {
  local DIR=$1
  local CHK=$2
  while true; do
    if [[ -d $DIR/$CHK ]]; then
      read -p "The $CHK folder already exists. Would you like to overwrite it? " yn
      case $yn in
        [Yy1]* ) rm -rf "$DIR"/"$CHK"; break;;
        [Nn0]* ) echo "Files will be added to the $CHK folder"; break;;
        * ) echo "Please answer yes or no. ";;
      esac
    else break;
    fi
  done
}
