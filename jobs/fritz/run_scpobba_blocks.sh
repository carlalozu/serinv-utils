#!/bin/bash
#SBATCH --job-name=scpobbasi_blocks      # Job name   
#SBATCH --output=outputs/scpobbasi_blocks-%j.out # Output file
#SBATCH --error=outputs/scpobbasi_blocks-%j.err  # Error file 
#SBATCH --ntasks=72               # Number of tasks
#SBATCH --cpus-per-task=1         # Number of CPUs per task
#SBATCH --mem-per-cpu=1024        # Memory per CPU
#SBATCH --time=02:00:00           # Wall clock time limit
#SBATCH -N 1                     # One node

# Load modules
# module load gcc openmpi python

# Open venv
# conda activate serinv_cpu

i=16
inside_n=$((2**i))
arrowhead_blocksize=128
n=$((inside_n+arrowhead_blocksize)) # total matrix size

# Parameters
script="scaling_scpobbasi.py"
output_file=scpobbasi_blocks_$arrowhead_blocksize.txt

# Create output files
> results/$output_file

echo "id,n_runs,n,bandwidth,arrowhead_blocksize,effective_bandwidth,diagonal_blocksize,n_offdiags,n_t,scpobbaf_time,scpobbasi_time,scpobbaf_FLOPS,scpobbasi_FLOPS"  | tee -a results/$output_file


for ((j=i-7; j<i-2; j++)) do

    bandwidth=$((2**j+1)) # must be odd

    overwrite=1
    n_runs=6
    
    for ((k=0; k<6; k++)) do

        # Fixed number of off diagonal blocks
        n_offdiags_blk=$((2**k))

        echo "Running $script with matrix size $n, bandwidth $bandwidth, j $j, diagonal_blocksize $diagonal_blocksize"
        echo -n "$j," | tee -a results/$output_file
        python ../../scaling/$script --n=$n --bandwidth=$bandwidth --arrowhead_blocksize=$arrowhead_blocksize --n_offdiags_blk=$n_offdiags_blk --overwrite=$overwrite --n_runs=$n_runs| tee -a results/$output_file
        done

    done

done