#!/bin/bash
#SBATCH --job-name=game_of_life
#SBATCH --output=out_%j.out
#SBATCH --error=err_%j.err
#SBATCH --partition=thin_course
#SBATCH --nodes=1
#SBATCH --time=00:20  # 10 minutes time limit

module load 2023
module load MPICH/4.1.2-GCC-12.3.0

# Array of processes
processes=(1 2 4 8 16 32 64 128)

# Array of threads
threads=(1 2 4 8 16)

make

# Loop over processes
for p in "${processes[@]}"; do
    # Loop over threads
    for t in "${threads[@]}"; do
        echo "Running with $p processes and $t threads"
        mpirun -np $p --bind-to none ./game_of_life -omp $t
        echo "----------------------"
    done
done
