#!/bin/bash -l
#SBATCH --job-name=scpobasi      # Job name   
#SBATCH --output=output/scpobasi-%j.out # Output file
#SBATCH --error=output/scpobasi-%j.err  # Error file 
#SBATCH --ntasks=32              # Number of tasks
#SBATCH --cpus-per-task=1        # Number of CPUs per task
#SBATCH --mem=64G
#SBATCH --time=03:00:00          # Wall clock time limit

# Load modules
module load gcc openmpi python

# Open venv
conda activate serinv_cpu

# create file
script_path=../../scaling

output_file=scpobasi_64.txt
> results/$output_file

echo "n_runs,n,bandwidth,arrowhead_blocksize,effective_bandwidth,diagonal_blocksize,n_offdiags,n_t,time_f_median,time_f_std,time_si_median,time_si_std,scpobaf_FLOPS,scpobasi_FLOPS" | tee -a results/scpobasi_64.txt

i=16
for ((j=i-7; j<i-2; j+=2)) do

    inside_n=$((2**i))
    bandwidth=$((2**j+1))
    arrowhead_blocksize=$((64))
    n=$((inside_n+arrowhead_blocksize))

    n_runs=4
    
    echo "Running scpobasi with matrix size $n, bandwidth $bandwidth, j $j, arrowhead_blocksize $arrowhead_blocksize"
    python $script_path/scaling_scpobasi.py --n=$n --bandwidth=$bandwidth --arrowhead_blocksize=$arrowhead_blocksize --n_runs=$n_runs |  tee -a results/$output_file

done


# python plot.py