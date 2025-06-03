#!/bin/bash -l
#SBATCH --job-name=compare        # Job name   
#SBATCH --output=compare-%j.out   # Output file
#SBATCH --error=compare-%j.err    # Error file 
#SBATCH --ntasks=32               # Number of tasks
#SBATCH --cpus-per-task=1         # Number of CPUs per task
#SBATCH --mem-per-cpu=1024        # Memory per CPU
#SBATCH --time=03:00:00           # Wall clock time limit

# Load modules
# module load gcc openmpi python

# Open venv
# conda activate serinv_cpu

# Create output file
> results_combined.txt

echo "run,n,bandwidth,n_offdiags,diagonal_blocksize,arrowhead_blocksize,time,numpy_time,error" | tee -a results_combined.txt

# Run scaling script
for ((j=10; j<13; j++)) do

    n=$((2**16))          # Size of matrix
    bandwidth=$((2**j*3)) # Full bandwidth of matrix

    n_runs=4

    echo "Running with matrix size $n, bandwidth $bandwidth"

    echo "Running scpobbasi"
    python scaling_scpobbasi.py --n=$n --n_offdiags_blk=$(((bandwidth/1024 - 1)/2))  --diagonal_blocksize=$((1024)) --arrowhead_blocksize=$((bandwidth/8)) |  tee -a results_combined.txt

    echo "Running pobtasi"
    python scaling_pobtasi.py --n=$n --diagonal_blocksize=$((bandwidth/3)) --arrowhead_blocksize=$((bandwidth/8)) |  tee -a results_combined.txt

    echo "Running scpobasi"
    python scaling_scpobasi.py --n=$n --n_offdiags=$((bandwidth/2-1)) --arrowhead_size=$((bandwidth/8)) |  tee -a results_combined.txt

done