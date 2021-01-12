#!/bin/bash
#PBS -l nodes=3:ppn=12
# #PBS -q XXL
#PBS -l walltime=02:00:00
#PBS -o jobSca0036.out
#PBS -e jobSca0036.err
#PBS -N jobSca0036

echo "workdir=${PBS_O_WORKDIR}"
echo "nodefile=${PBS_NODEFILE}"


cd ${PBS_O_WORKDIR}
cd ../intel/oneapi
# source setvars.sh
cd ../../ATAMPI

source init.sh
mpirun --hostfile ${PBS_NODEFILE} ./scaMain 20 60000 5000 scaResults/res/60000x5000/scaRes0036.txt 
