#!/bin/bash
### run this file with command "sbatch slurm_taylor.sh"
#SBATCH --account=phangs
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=999G
#SBATCH --cpus-per-task=96

### CHANGE THESE!
#SBATCH --chdir=/gscratch/aweinbec/Pa_processing/Pa_scripts/
### Here is where the output of your program goes
#SBATCH --output=/gscratch/aweinbec/Pa_processing/Pa_scripts/logs/Pa_5-27_all.out
### Here is where the errors get sent
#SBATCH  --error=/gscratch/aweinbec/Pa_processing/Pa_scripts/logs/Pa_5-27_all.err

### you can continuously monitor these by typing:
###   watch -n 1 tail logs/Pa_5-27_all.out

### Don't change these!
echo "SLURM_JOB_ID:" $SLURM_JOB_ID
echo "SLURM_JOB_NAME:" $SLURM_JOB_NAME
echo "SLURM_JOB_NODELIST:" $SLURM_JOB_NODELIST

### Here's the name of the python script you want to run
python Pa_destripe_corrxxx.py

