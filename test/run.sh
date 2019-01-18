#!/bin/bash

# This script executes all the binaries under the bin directory.
# It expects the following bin directory structure.
#
# bin
#  |
#  +--- subdir_1
#         |
#         +--- subsubdir_1
#         |
#         +--- subsubdir_2
#  |
#  +--- ...
#  |
#  +--- subdir_n
#         |
#         +--- subsubdir_1
#
#
# A subdir can have any name. They refer to the general method the binary uses.
# A subsubdir's name should start with "sequential" or "parallel" followed by
# "_" if the name continues.


# Miscellaneous
## Color codes
RED='\033[0;31m'
GREEN='\033[0;32m'
NO_COLOR='\033[0m'

# Full path to base directory, where the script is found
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null && pwd )"

# Full path to bin directory inside the base directory
BIN_DIR="${DIR}/bin"

ExecuteAndCompare() {
    MODE=$1
    N=$2
    DATA_FILE=$3
    RESULTS_DIR=$4
    EXPECTED_RESULTS_DIR=$5

    if [[ "${MODE}" == "sequential" ]]; then
        echo -e "\tExecuting: ./prog -n ${N} -i "${DATA_FILE}" -d "${RESULTS_DIR}""
        ./prog -n ${N} -i "${DATA_FILE}" -d "${RESULTS_DIR}"
    elif [[ "${MODE}" == "parallel" ]]; then
        echo -e "\tExecuting: mpirun -np 2 ./prog -n ${N} -i "${DATA_FILE}" -d "${RESULTS_DIR}""
        mpirun -np 2 ./prog -n ${N} -i "${DATA_FILE}" -d "${RESULTS_DIR}"
    else
        echo This should not execute ever!
    fi

    FILENAME=$(basename -- "${DATA_FILE}")
    DIFF="$( diff \
             <(awk '{print $1, $2, $3}' "${RESULTS_DIR}/${FILENAME}.out") \
             <(awk '{print $1, $2, $3}' "${EXPECTED_RESULTS_DIR}/${FILENAME}.out") \
           )"

    if [[ "${DIFF}" == "" ]]; then
        echo -e "\t${GREEN}OK${NO_COLOR}"
    else
        echo -e "\t${RED}NOT OK${NO_COLOR}"
    fi
}

# For each file inside the bin dir
for f1 in ${BIN_DIR}/*; do
    # If it is a dir
    if [[ -d ${f1} ]]; then
        # For each of its files
        for f2 in ${f1}/*; do
            # If it is a dir
            if [[ -d ${f2} ]]; then
                # Change to that directory
                cd ${f2}
                # Get directory base name
                CUR_DIR="$( basename "${f2}" )"

                # Prefix of directory's name
                # Should be 'sequential' or 'parallel'
                MODE=${CUR_DIR%%_*}

                echo Current dir: "${f2}"
                echo Executing in "${MODE}"

                DATA_DIR="resources/testsets"
                EXPECTED_RESULTS_DIR="resources/expected"
                RESULTS_DIR="results"

                for data_file in ${DATA_DIR}/*; do
                    # Number of particles 
                    N="$( wc -l < "${data_file}" )"

                    ExecuteAndCompare "${MODE}" "${N}" "${data_file}" "${RESULTS_DIR}" "${EXPECTED_RESULTS_DIR}"
                done

                # Print new line
                echo
            fi
        done
    fi
done
