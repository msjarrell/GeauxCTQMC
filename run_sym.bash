#!/bin/bash
#PBS -q workq
#PBS -l walltime=12:00:00
#PBS -l nodes=1:ppn=20
#PBS -V
#PBS -N sym_test
#PBS -o /home/sfeng/mpi_wrapper/test.out
#PBS -e /home/sfeng/mpi_wrapper/test.err


module load impi/4.1.3.048/intel64  # load impi

export TASKS_PER_HOST=1 # number of MPI tasks per host
export THREADS_HOST=40    # number of OpenMP threads spawned by each task on the host
export TASKS_PER_MIC=1  # number of MPI tasks per MIC
export THREADS_MIC=240     # number of OpenMP threads spawned by each task on the MIC

#cd $PBS_O_WORKDIR        # go to where your PBS job is submitted if necessary
#micrun.sym -c /home/sfeng/mpi_wrapper/a.host -m /home/sfeng/mpi_wrapper/a.mic  # run micrun.sym

cd /home/sfeng/GeauxCTQMC/

micrun.sym -c ./ctqmc_cpu -m ./ctqmc_mic