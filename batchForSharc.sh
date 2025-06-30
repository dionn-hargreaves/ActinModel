#!/bin/bash

# Example batch script for ShARC

#Request memory
#$ -l rmem=2G

#Request CPU
#$ -pe openmp 2

#Request max time
#$ -l h_rt=96:00:00

#Email notifications - change EMAIL to your email address
#$ -M EMAIL
#$ -m bea

# Output streams merge
#$ -j y


#Load anaconda python
module load apps/python/anaconda3-4.2.0

#Load compiler and boost
module load libs/boost/1.64.0/gcc-8.2-cmake-3.17.1

# Activate python environment
source activate forActin

#Set OpenMP_NUM_THREADS to 2
export OMP_NUM_THREADS=2

#Need to add python to library path
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$CONDA_PREFIX/lib

#Change to directory - change USERNAME with your username on ShARC
#$ -wd /data/USERNAME/ActinModelling/actinmodel

# run make
make

#run test
bin/model configExamples/SingleFila_ShARC.cfg
