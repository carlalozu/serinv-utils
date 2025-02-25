#!/bin/bash
#SBATCH --job-name=scpobasi      # Job name   
#SBATCH --output=outputs/scpobasi-%j.out # Output file
#SBATCH --error=outputs/scpobasi-%j.err  # Error file 
#SBATCH --ntasks=1               # Number of tasks
#SBATCH --cpus-per-task=16       # Number of CPUs per task
#SBATCH --gres=gpu:a40:1         # Number of GPUs
#SBATCH --time=03:00:00          # Wall clock time limit
#SBATCH -N 1                     # One node


# # Load modules
# module load gcc python cuda

# # Open venv
# conda activate serinv_env

# Parameters
arrowhead_blocksize=$((64))
output_file="scpobasi_$arrowhead_blocksize.txt"
script="scaling_scpobasi.py"

# Create output files
> results/$output_file

echo "run,id,n,bandwidth,arrowhead_blocksize,effective_bandwidth,diagonal_blocksize,n_offdiags,n_t,scpobaf_time,scpobasi_time,scpobaf_FLOPS,scpobasi_FLOPS"  | tee -a results/$output_file

i=16
inside_n=$((2**i))
n=$((inside_n+arrowhead_blocksize)) # total matrix size

for ((j=i-7; j<i-2; j++)) do

    bandwidth=$((2**j+1)) # must be odd

    overwrite=1
    n_runs=4
    
    echo "Running $script with matrix size $n, bandwidth $bandwidth, j $j, diagonal_blocksize $diagonal_blocksize"
    for ((r=0; r<n_runs; r++)) do
        echo -n "$r,$j," | tee -a results/$output_file
        python ../../scaling/$script --n=$n --bandwidth=$bandwidth --arrowhead_blocksize=$arrowhead_blocksize --overwrite=$overwrite| tee -a results/$output_file
    done
done