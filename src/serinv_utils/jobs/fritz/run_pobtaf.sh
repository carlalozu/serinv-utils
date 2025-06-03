#!/bin/bash
#SBATCH --job-name=pobtaf      # Job name   
#SBATCH --output=pobtaf-%j.out # Output file
#SBATCH --error=pobtaf-%j.err  # Error file 
#SBATCH --ntasks=72               # Number of tasks
#SBATCH --cpus-per-task=1         # Number of CPUs per task
#SBATCH --mem-per-cpu=1024        # Memory per CPU
#SBATCH --time=03:00:00           # Wall clock time limit

# Load modules
# module load gcc openmpi python

# Open venv
# conda activate serinv_cpu

# Create output files
> results/pobtaf_128_new.txt

echo "run,id,n,bandwidth,arrowhead_blocksize,effective_bandwidth,diagonal_blocksize,n_offdiags,n_t,time,numpy_time,error"  | tee -a results/pobtaf_128_new.txt

for ((j=8; j<14; j++)) do

    i=16
    inside_n=$((2**i))
    bandwidth=$((2**j+1)) # must be odd
    arrowhead_blocksize=$((128))
    n=$((inside_n+arrowhead_blocksize)) # total matrix size

    n_runs=6
    
    echo "Running pobtaf with matrix size $n, bandwidth $bandwidth, j $j, arrowhead_blocksize $arrowhead_blocksize"
    for ((r=0; r<n_runs; r++)) do
        echo -n "$r,$j," | tee -a results/pobtaf_128_new.txt
        python ../../dev/scaling_pobtaf.py --n=$n --bandwidth=$bandwidth --arrowhead_blocksize=$arrowhead_blocksize |  tee -a results/pobtaf_128_new.txt
    done
done