#!/bin/bash

# remove_lib_path.sh

# set -x 

# Function to remove a library path from R environment
remove_home_R_library_path() {
    # Current library paths
    # current_paths=$(R -q -e '.libPaths()')

    # # # Library path to remove
    # path_to_remove="$HOME/R-3.5.1/library"

    # # Remove the path from the library paths using awk
    # new_paths=$(echo "$current_paths" | awk -v path="$path_to_remove" '{ if ($0 != path) { print $0 } }')


    # echo "new Paths"
    # echo $new_paths
    #R -q -e ".libPaths(c('$new_paths'))"

    #echo $current_paths
    # R -q -e ".libPaths(.libPaths()[2:3])"

    export R_LIBS="/usr/lib64/R/library"
    export R_LIBS_USER=""
    export R_LIBS_SITE=""
}

remove_home_R_library_path