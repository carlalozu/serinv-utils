#!/bin/bash
#SBATCH --job-name=operations      # Job name   
#SBATCH --output=operations-%j.out # Output file
#SBATCH --error=operations-%j.err  # Error file 
#SBATCH --ntasks=1               # Number of tasks
#SBATCH --cpus-per-task=16       # Number of CPUs per task
#SBATCH --gres=gpu:a40:1         # Number of GPUs
#SBATCH --time=03:00:00          # Wall clock time limit
#SBATCH -N 1                     # One node

# Load modules
# module load gcc openmpi python

# Open venv
# conda activate serinv_cpu

# Create output files
# > results/operations.txt

# echo "id,repetitions,diag_blocksize,arrowhead_blocksize,cholesky_ns3,triang_solve_ns3,triang_solve_ns2nb,dgemm_ns3,dgemm_ns2nb,dgemm_nb2ns,dgemm_nsnb2,matrix_vector_nsns,matrix_vector_nsnb,matrix_vector_nbns,dot_prod_ns,scale_ns" | tee -a results/operations.txt

for ((j=10; j<15; j++)) do

    diag_blocksize=$((2**j))
    arrowhead_blocksize=$((64))

    repetitions=1000

    echo -n "$j," | tee -a results/operations.txt
    python ../../dev/scaling_operations.py --diag_blocksize=$diag_blocksize --arrowhead_blocksize=$arrowhead_blocksize --repetitions=$repetitions |  tee -a results/operations.txt
done