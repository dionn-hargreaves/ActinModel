#!/bin/bash

# Example batch script for ShARC

#Request memory
#SBATCH --mem=2G

#Request CPU
#SBATCH -c 2

#Request max time
#SBATCH -t 10

#Email notifications - change EMAIL to your email address
#SBATCH --mail-user=d.hargreaves@sheffield.ac.uk

# Email notifications if the job fails
#SBATCH --mail-type=FAIL

# Change the name of the output log file.
#SBATCH --output=output.%j.test.out
# Rename the job's name
#SBATCH --job-name=actinModelling


#Load anaconda python
#module load apps/python/anaconda3-4.2.0

#Load compiler and boost
module load libs/boost/1.64.0/gcc-8.2-cmake-3.17.1

# Activate python environment
#source activate forActin

#Set OpenMP_NUM_THREADS
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK


#Change to directory - change USERNAME with your username on ShARC
#SBATCH --workdir=/mnt/parscratch/users/md1dhar/actinmodel/ActinModelling

# run make
make

#run test
bin/model config/Stanage_test_phago.cfg
