#!/bin/bash -l
#SBATCH --job-name=pobtasi      # Job name   
#SBATCH --output=output/pobtasi-%j.out # Output file
#SBATCH --error=output/pobtasi-%j.err  # Error file 
#SBATCH --ntasks=1               # Number of tasks
#SBATCH --cpus-per-task=16       # Number of CPUs per task
#SBATCH --gres=gpu:a40:1         # Number of GPUs
#SBATCH --time=01:00:00          # Wall clock time limit
#SBATCH -N 1                     # One node
#SBATCH --exclusive               # Exclusive access


# Load modules
module load gcc openmpi python

# Open venv
conda activate serinv_env

# Create output file
script_path=../../scaling
output_file=pobtasi_64_16.csv

> results/$output_file
echo "n_runs,n,bandwidth,arrowhead_blocksize,effective_bandwidth,diagonal_blocksize,n_offdiags,n_t,time_f_median,time_f_std,time_si_median,time_si_std,flops_f,flops_si" | tee -a results/$output_file

i=16
for ((j=i-7; j<i-2; j+=2)) do

    inside_n=$((2**i))
    bandwidth=$((2**j+1)) # must be odd
    arrowhead_blocksize=$((64))
    n=$((inside_n+arrowhead_blocksize)) # total matrix size

    n_runs=10
    
    echo "Running pobtaf with matrix size $n, bandwidth $bandwidth, arrowhead_blocksize $arrowhead_blocksize, n_runs $n_runs"
    python $script_path/scaling_pobtasi.py --n=$n --bandwidth=$bandwidth --arrowhead_blocksize=$arrowhead_blocksize --n_runs=$n_runs |  tee -a results/$output_file

done