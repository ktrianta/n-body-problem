#!/bin/bash

## declare an array variable
declare -a arr1=("./bin/naive/parallel/prog" "./bin/naive/parallel_one_sided/prog" "./bin/naive/sequential/prog" "./bin/barnes/parallel/prog" "./bin/barnes/parallel_one_sided/prog" "./bin/barnes/parallel_balanced/prog" "./bin/barnes/sequential/prog")
declare -a arr2=("test/resources/testsets/galaxy128_ball.txt" "test/resources/testsets/galaxy512_ball.txt" "test/resources/testsets/galaxy1024_ball.txt" "test/resources/testsets/galaxy65536_ball.txt" "test/resources/testsets/galaxy131072_ball.txt ")
declare -a arr3=("128" "512" "1024" "65536" "131072")
declare -a arr4=("naive-parallel.data" "naive-parallel-one-sided.data" "naive-sequential.data" "barnes-parallel.data" "barnes-parallel-one-sided.data" "barnes-parallel-balanced.data" "barnes-sequential.data")

## now loop through the above array
for j in $(seq 1 1 ${#arr2[@]})
do

    for i in $(seq 1 1 ${#arr1[@]})
    do
        echo "${arr2[$(($j-1))]}" "${arr3[$(($j-1))]}" "${arr1[$(($i-1))]}" 
        for k in $(seq 1 1 1)
        do
            echo "$k"
            mpirun -n 1 "${arr1[$(($i-1))]}" -i ${arr2[$(($j-1))]} -n "${arr3[$(($j-1))]}"  
            # >> "${arr4[$(($i-1))]}"       
        done
   done
done

