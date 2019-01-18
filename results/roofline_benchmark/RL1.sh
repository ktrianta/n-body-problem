#!/bin/bash

## declare an array variable
declare -a arr1=("../src/naive/parallel_one_sided/main.cpp" "../src/barnes/parallel_one_sided/main.cpp" "../src/barnes/parallel_balanced/main.cpp")
declare -a arr2=("../src/naive/parallel_one_sided/main_RF.cpp" "../src/barnes/parallel_one_sided/main_RF.cpp" "../src/barnes/parallel_balanced/main_RF.cpp")
declare -a arr3=("../src/naive/parallel_one_sided/CMakeLists.txt" "../src/barnes/parallel_one_sided/CMakeLists.txt" "../src/barnes/parallel_balanced/CMakeLists.txt")
declare -a arr4=("../src/naive/parallel_one_sided/CMakeLists_papi.txt" "../src/barnes/parallel_one_sided/CMakeLists_papi.txt" "../src/barnes/parallel_balanced/CMakeLists_papi.txt")


## now loop through the above array
for j in $(seq 1 1 ${#arr2[@]})
do
    cp "${arr2[$(($j-1))]}" "${arr1[$(($j-1))]}" 
    echo copied "${arr2[$(($j-1))]}" to "${arr1[$(($j-1))]}" 
    cp "${arr4[$(($j-1))]}" "${arr3[$(($j-1))]}" 
    echo copied "${arr4[$(($j-1))]}" to "${arr3[$(($j-1))]}" 
done

