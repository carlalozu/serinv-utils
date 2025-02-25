#!/bin/bash
#SBATCH --job-name=scpobasi      # Job name   
#SBATCH --output=outputs/scpobasi-%j.out # Output file
#SBATCH --error=outputs/scpobasi-%j.err  # Error file 
#SBATCH --ntasks=32              # Number of tasks
#SBATCH --cpus-per-task=1        # Number of CPUs per task
#SBATCH --mem=64G
#SBATCH --time=03:00:00          # Wall clock time limit

# Load modules
# module load gcc openmpi python

# Open venv
# conda activate serinv_cpu

arrowhead_blocksize=$((64))
output_file="scpobasi_$arrowhead_blocksize.txt"

# create file
> results/$output_file


echo "run,id,n,bandwidth,arrowhead_blocksize,effective_bandwidth,diagonal_blocksize,n_offdiags,n_t,scpobaf_time,scpobasi_time,scpobaf_FLOPS,scpobasi_FLOPS"  | tee -a results/$output_file


for ((j=9; j<13; j++)) do

    i=16
    inside_n=$((2**i))
    bandwidth=$((2**j+1))
    n=$((inside_n+arrowhead_blocksize))

    n_runs=4
    
    echo "Running scpobasi with matrix size $n, bandwidth $bandwidth, j $j, arrowhead_blocksize $arrowhead_blocksize"
    for ((r=0; r<n_runs; r++)) do
        echo -n "$r,$j," | tee -a results/$output_file
        python ../../scaling/scaling_scpobasi.py --n=$n --bandwidth=$bandwidth --arrowhead_blocksize=$arrowhead_blocksize |  tee -a results/$output_file
    done

done


# python plot.py