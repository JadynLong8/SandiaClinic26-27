#!/bin/bash
#SBATCH --job-name=2x2x2_BTO_Slab
#SBATCH -p RM-shared
#SBATCH -N 1
#SBATCH --ntasks-per-node=32
#SBATCH -t 12:00:00
#SBATCH --output=slurm-%j.out
#SBATCH --error=slurm-%j.err

set -e

#load the modules

module purge
# module load intel/20.4
# module load intelmpi/20.4-intel20.4
# module load mkl/2020.4.304
# module load allocations/1.0

# Load required modules
# module purge
# module load intel
module load openmpi/4.0.5-gcc10.2.0

# Add MKL library path so runtime linker can find .so files
# export LD_LIBRARY_PATH=$MKLROOT/lib/intel64:$LD_LIBRARY_PATH

cd "$PROJECT/BTO_slab_2x2x2_2023-24" || exit 1

echo "Working directory: $PWD"
echo "MPI launcher: $(which mpirun)"
mpirun --version

#mpirun "$PROJECT/software/seqquest_parallel/quest_265a.x"

ldd "$PROJECT/software/seqquest_parallel/quest_265a.x"

mpirun -np 32 "$PROJECT/software/seqquest_parallel/quest_265a.x"

