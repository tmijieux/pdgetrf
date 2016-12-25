#!/usr/bin/env bash
#SBATCH --job-name=mijieux-moussaid
#SBATCH --output=out.1
#SBATCH --error=err.1
#SBATCH -p mistral
#SBATCH --exclusive
#SBATCH --time=02:00:00
#SBATCH --nodes=9
#SBATCH --ntasks-per-node=1

# 9 noeud / 9 proc MPI / mkl parallele 20 threads

WORKDIR=${WORKDIR:-${HOME}/pdgetrf}

cd ${WORKDIR}
. ./.module.load

time mpirun -n 9 ./pdgetrf 20000 200
