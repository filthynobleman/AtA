#!/bin/bash
#PBS -l nodes=3:ppn=5
# #PBS -q XXL
#PBS -l walltime=00:30:00
#PBS -o jobAtA.out
#PBS -e jobAtA.err
#PBS -N jobATA

echo "workdir=${PBS_O_WORKDIR}"
echo "nodefile=${PBS_NODEFILE}"

cd ${PBS_O_WORKDIR}
cd ../intel/oneapi
# source setvars.sh
cd ../../ATAMPI

source init.sh
mpirun --hostfile ${PBS_NODEFILE} ./main 20 20000 20000 results0016.txt
