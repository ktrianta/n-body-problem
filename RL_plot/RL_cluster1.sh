#BSUB -q bigmem.24h
#BSUB -o 4proc.out
#BSUB -e 4proc.err
#BSUB -n 1
#BSUB -R "rusage[mem=2500]"
#BSUB -R "select[model==XeonE5_2680v3]"


cat /proc/cpuinfo
module load new
module load gcc/6.3.0
module load mvapich2/2.2



declare -a arr1=("./../bin/naive/parallel_one_sided/prog" "./../bin/barnes/parallel_one_sided/prog -h 0.25" "./../bin/barnes/parallel_balanced/prog -h 0.25")

declare -a arr2=("../test/resources/testsets/galaxy128_ball.txt" "../test/resources/testsets/galaxy256_ball.txt" "../test/resources/testsets/galaxy512_ball.txt" "../test/resources/testsets/galaxy1024_ball.txt" "../test/resources/testsets/galaxy2048_ball.txt" "../test/resources/testsets/galaxy4096_ball.txt" "../test/resources/testsets/galaxy8192_ball.txt" "../test/resources/testsets/galaxy16384_ball.txt" "../test/resources/testsets/galaxy32768_ball.txt" "../test/resources/testsets/galaxy65536_ball.txt")

declare -a arr3=("128" "256" "512" "1024" "2048" "4096" "8192" "16384" "32768" "65536")

declare -a arr4=("naive-parallel-one-sided.data" "barnes-parallel-one-sided.data" "barnes-parallel-balanced.data")

## now loop through the above array
for j in $(seq 1 1 ${#arr2[@]})
do

    for i in $(seq 1 1 ${#arr1[@]})
    do
        echo "${arr2[$(($j-1))]}" "${arr3[$(($j-1))]}" "${arr1[$(($i-1))]}" 
        for k in $(seq 1 1 10)
        do
            echo "$k"
            mpirun -n 1 ${arr1[$(($i-1))]} -i ${arr2[$(($j-1))]} -n "${arr3[$(($j-1))]}"  >> "${arr4[$(($i-1))]}"
        done
   done
done


