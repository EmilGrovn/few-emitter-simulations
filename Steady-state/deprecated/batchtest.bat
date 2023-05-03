#!/bin/bash
#SBATCH --mail-type=NONE
#SBATCH --mail-user=s183774@student.dtu.dk  # The default value is the submitting user.
#SBATCH --partition=xeon56
#SBATCH -N 1      # Minimum of 1 nodes
#SBATCH -n 40    # 24 MPI processes per node, 48 tasks in total, appropriate for xeon24 nodes
#SBATCH --time=1-02:00:00
#SBATCH --output=mpi_job_slurm_output.log
#SBATCH --error=mpi_job_slurm_errors.log

module purge
module add intel

#~/bin/