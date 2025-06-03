#!/bin/bash
#SBATCH --job-name=scpobaf       # Job name   
#SBATCH --output=scpobaf-%j.out  # Output file
#SBATCH --error=scpobaf-%j.err   # Error file 
#SBATCH --ntasks=32              # Number of tasks
#SBATCH --cpus-per-task=1        # Number of CPUs per task
#SBATCH --mem-per-cpu=1024       # Memory per CPU
#SBATCH --time=03:00:00          # Wall clock time limit

# Load modules
# module load gcc openmpi python

# Open venv
# conda activate serinv_cpu

# create file
path=results
file_name=scpobaf_64_F.txt
> $path/$file_name

echo "run,id,n,bandwidth,arrowhead_blocksize,effective_bandwidth,diagonal_blocksize,n_offdiags,n_t,time,numpy_time,error"  | tee -a $path/$file_name

for ((j=9; j<13; j++)) do

    i=16
    inside_n=$((2**i))
    bandwidth=$((2**j+1))
    arrowhead_blocksize=$((64))
    n=$((inside_n+arrowhead_blocksize))

    n_runs=4
    
    echo "Running scpobaf with matrix size $n, bandwidth $bandwidth, j $j, arrowhead_blocksize $arrowhead_blocksize"
    for ((r=0; r<n_runs; r++)) do
        echo -n "$r,$j," | tee -a $path/$file_name
        python ../../dev/scaling_scpobaf.py --n=$n --bandwidth=$bandwidth --arrowhead_blocksize=$arrowhead_blocksize |  tee -a $path/$file_name
    done
done
