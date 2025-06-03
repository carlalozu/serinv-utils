#!/bin/bash
#SBATCH --job-name=scpobbaf_blocks      # Job name   
#SBATCH --output=scpobbaf_blocks-%j.out # Output file
#SBATCH --error=scpobbaf_blocks-%j.err  # Error file 
#SBATCH --ntasks=72               # Number of tasks
#SBATCH --cpus-per-task=1         # Number of CPUs per task
#SBATCH --mem-per-cpu=1024        # Memory per CPU
#SBATCH --time=02:00:00           # Wall clock time limit
#SBATCH -N 1                     # One node
#SBATCH --exclusive               # Exclusive access

# Load modules
# module load gcc openmpi python

# Open venv
# conda activate serinv_cpu

# Parameters
script="scaling_scpobbaf.py"
output_file=scpobbaf_blocks_64.txt

# Create output files
> results/$output_file

echo "run,id,n,bandwidth,arrowhead_blocksize,effective_bandwidth,diagonal_blocksize,n_offdiags,n_t,time,numpy_time,error,FLOPS"  | tee -a results/$output_file

i=16
inside_n=$((2**i))

for ((j=i-7; j<i-3; j++)) do

    bandwidth=$((2**j+1)) # must be odd
    arrowhead_blocksize=64
    n=$((inside_n+arrowhead_blocksize)) # total matrix size

    numpy_compare=0
    overwrite=1
    n_runs=6
    
    for ((k=0; k<6; k++)) do

        # Fixed number of off diagonal blocks
        n_offdiags_blk=$((2**k))

        echo "Running $script with matrix size $n, bandwidth $bandwidth, j $j, diagonal_blocksize $diagonal_blocksize"
        for ((r=0; r<n_runs; r++)) do
            echo -n "$r,$j," | tee -a results/$output_file
            python ../../dev/$script --n=$n --bandwidth=$bandwidth --arrowhead_blocksize=$arrowhead_blocksize --n_offdiags_blk=$n_offdiags_blk --numpy_compare=$numpy_compare --overwrite=$overwrite| tee -a results/$output_file
        done

    done

done