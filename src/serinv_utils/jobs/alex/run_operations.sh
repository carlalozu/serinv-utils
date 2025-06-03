#!/bin/bash
#SBATCH --job-name=operations      # Job name   
#SBATCH --output=output/operations-%j.out # Output file
#SBATCH --error=output/operations-%j.err  # Error file 
#SBATCH --ntasks=1               # Number of tasks
#SBATCH --cpus-per-task=16       # Number of CPUs per task
#SBATCH --gres=gpu:a40:1         # Number of GPUs
#SBATCH --time=03:00:00          # Wall clock time limit
#SBATCH -N 1                     # One node

# Load modules
module load gcc openmpi python

# Open venv
conda activate serinv_cpu

output_file=operations.csv

# # Create output files
> results/$output_file

echo "id,n_runs,diag_blocksize,arrowhead_blocksize,cholesky_ns3,triang_solve_ns3,triang_solve_ns2nb,dgemm_ns3,dgemm_ns2nb,dgemm_nb2ns,dgemm_nsnb2,matrix_vector_nsns,matrix_vector_nsnb,matrix_vector_nbns,dot_prod_ns,scale_ns" | tee -a results/$output_file

for ((j=4; j<13; j++)) do

    diag_blocksize=$((2**j))
    arrowhead_blocksize=$((64))
    n_runs=1000

    echo -n "$j," | tee -a results/$output_file
    python ../../scaling/scaling_operations.py --diag_blocksize=$diag_blocksize --arrowhead_blocksize=$arrowhead_blocksize --n_runs=$n_runs |  tee -a results/$output_file
done