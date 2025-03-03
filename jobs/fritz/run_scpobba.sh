#!/bin/bash
#SBATCH --job-name=scpobbaf       # Job name   
#SBATCH --output=scpobbaf-%j.out  # Output file
#SBATCH --error=scpobbaf-%j.err   # Error file 
#SBATCH --ntasks=1                # Number of tasks
#SBATCH --cpus-per-task=32        # Number of CPUs per task
#SBATCH --mem-per-cpu=1024        # Memory per CPU
#SBATCH --time=03:00:00           # Wall clock time limit

# Load modules
# module load gcc openmpi python

# Open venv
# conda activate serinv_cpu

# Parameters
arrowhead_blocksize=$((64))
output_file="scpobbaf_$arrowhead_blocksize.txt"
script="scaling_scpobbaf.py"

# Create output files
> results/$output_file

echo "id,n_runs,n,bandwidth,arrowhead_blocksize,effective_bandwidth,diagonal_blocksize,n_offdiags,n_t,time,numpy_time,error,FLOPS"  | tee -a results/$output_file

i=16
inside_n=$((2**i))
n=$((inside_n+arrowhead_blocksize)) # total matrix size

for ((j=i-7; j<i-2; j++)) do

    bandwidth=$((2**j+1)) # must be odd
    n_offdiags_blk=1

    numpy_compare=0
    overwrite=1
    n_runs=6
    
    echo "Running $script with matrix size $n, bandwidth $bandwidth, j $j, diagonal_blocksize $diagonal_blocksize"
    echo -n "$j," | tee -a results/$output_file
    python ../../scaling/$script --n=$n --bandwidth=$bandwidth --arrowhead_blocksize=$arrowhead_blocksize --n_offdiags_blk=$n_offdiags_blk --overwrite=$overwrite --n_runs=$n_runs| tee -a results/$output_file
done