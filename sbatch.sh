#!/bin/bash
#SBATCH --job-name=game_of_life
#SBATCH --output=out_%j.out
#SBATCH --error=err_%j.err
#SBATCH --partition=thin_course
#SBATCH --nodes=1
#SBATCH --time=10:00  # 10 minutes time limit

# Array of processes
processes=(1 2 4 8 16 32 64 128)

# Array of threads
threads=(1 2 4 8 16)

# Loop over processes
for p in "${processes[@]}"; do
    # Loop over threads
    for t in "${threads[@]}"; do
        echo "Running with $p processes and $t threads"
        mpirun -np $p --bind-to none ./your_executable -omp $t
        echo "----------------------"
    done
done
