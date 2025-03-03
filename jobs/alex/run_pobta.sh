#!/bin/bash
#SBATCH --job-name=pobtasi      # Job name   
#SBATCH --output=outputs/pobtasi-%j.out # Output file
#SBATCH --error=outputs/pobtasi-%j.err  # Error file 
#SBATCH --ntasks=1               # Number of tasks
#SBATCH --cpus-per-task=16       # Number of CPUs per task
#SBATCH --gres=gpu:a100:1         # Number of GPUs
#SBATCH --time=01:00:00          # Wall clock time limit
#SBATCH -N 1                     # One node


# # Load modules
# module load gcc python cuda

# # Open venv
# conda activate serinv_env


# Parameters
arrowhead_blocksize=$((64))
output_file="pobtasi_$arrowhead_blocksize.txt"
script="scaling_pobtasi.py"

# Create output file
> results/$output_file

echo "id,n_runs,n,bandwidth,arrowhead_blocksize,effective_bandwidth,diagonal_blocksize,n_offdiags,n_t,pobtaf_time,pobtasi_time,pobtaf_FLOPS,pobtasi_FLOPS"  | tee -a results/$output_file

i=16
inside_n=$((2**i))
n=$((inside_n+arrowhead_blocksize)) # total matrix size

for ((j=i-7; j<i-2; j++)) do

    bandwidth=$((2**j+1)) # must be odd

    n_runs=6
    
    echo "Running $script with matrix size $n, bandwidth $bandwidth, j $j, diagonal_blocksize $diagonal_blocksize"
    echo -n "$j," | tee -a results/$output_file
    python ../../scaling/$script --n=$n --bandwidth=$bandwidth --arrowhead_blocksize=$arrowhead_blocksize --n_runs=$n_runs| tee -a results/$output_file
done