#!/bin/bash -l
#SBATCH --job-name=scpobbasi_blocks      # Job name   
#SBATCH --output=output/scpobbasi_blocks-%j.out # Output file
#SBATCH --error=output/scpobbasi_blocks-%j.err  # Error file 
#SBATCH --ntasks=72               # Number of tasks
#SBATCH --cpus-per-task=1         # Number of CPUs per task
#SBATCH --mem-per-cpu=1024        # Memory per CPU
#SBATCH --time=02:00:00           # Wall clock time limit
#SBATCH -N 1                     # One node


# Load modules
module load gcc openmpi python

# Open venv
conda activate serinv_cpu

# Parameters
script_path=../../scaling
script=scaling_scpobbasi.py
output_file=scpobbasi_blocks_64_16.csv


# Create output files
> results/$output_file


echo "n_runs,n,bandwidth,arrowhead_blocksize,effective_bandwidth,diagonal_blocksize,n_offdiags,n_t,time_c_mean,time_c_std,time_si_mean,time_si_std,flops_c,flops_si" | tee -a results/$output_file

i=16
inside_n=$((2**i))

for ((j=i-7; j<i-3; j+=2)) do

    bandwidth=$((2**j+1)) # must be odd
    arrowhead_blocksize=64
    n=$((inside_n+arrowhead_blocksize)) # total matrix size

    overwrite=1
    n_runs=10
    
    for ((k=0; k<6; k++)) do

        # Fixed number of off diagonal blocks
        n_offdiags_blk=$((2**k))

        echo "Running $script with matrix size $n, bandwidth $bandwidth, n_offdiags_blk $n_offdiags_blk, n_runs $n_runs"

        python $script_path/$script --n=$n --bandwidth=$bandwidth --arrowhead_blocksize=$arrowhead_blocksize --n_offdiags_blk=$n_offdiags_blk --n_runs=$n_runs --overwrite=$overwrite| tee -a results/$output_file

    done

done