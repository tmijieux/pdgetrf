#!/usr/bin/env bash
#SBATCH --job-name=thomas-mijieux
#SBATCH --output=job_test.out
#SBATCH --error=job_test.err
#SBATCH -p mistral
#SBATCH --time=00:60:00
#SBATCH --exclusive
#SBATCH --nodes=1 --ntasks-per-node=1

module load intel/mkl/64/11.2/2016.0.0
module load compiler/gcc/5.1.0
module load slurm/14.11.11
module load hardware/hwloc/1.11.0
module load mpi/openmpi/gcc/1.10.1-tm

cd /home/prcd2016-mijieux/matprodmpi

do_job() {
    mpiexec -n 1 --output-filename output/output ./matprod input/big_20K.txt input/big_20K.txt
}

do_job

