#!/bin/bash

## declare an array variable
declare -a arr1=("src/naive/parallel_one_sided/main.cpp" "src/barnes/parallel_one_sided/main.cpp" "src/barnes/parallel_balanced/main.cpp")
declare -a arr2=("src/naive/parallel_one_sided/main_RF.cpp" "src/barnes/parallel_one_sided/main_RF.cpp" "src/barnes/parallel_balanced/main_RF.cpp")


## now loop through the above array
for j in $(seq 1 1 ${#arr2[@]})
do
    echo "${arr2[$(($j-1))]}" "${arr1[$(($j-1))]}" 
    cp "${arr2[$(($j-1))]}" "${arr1[$(($j-1))]}" 
done

