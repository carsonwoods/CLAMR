#!/bin/bash
#SBATCH --job-name=CLAMR-unm-l7-mpi
#SBATCH --nodes=1
#SBATCH --tasks-per-node=36
#SBATCH --cpus-per-task=1

#SBATCH --time=1:00:00              # Time limit hrs:min:sec
#SBATCH --output=%x.%j.log # Standard output and error log

#SBATCH --sockets-per-node=2 
#SBATCH --cores-per-socket=18 
#SBATCH --threads-per-core=1
#SBATCH --partition pbatch

echo "Setting up ${SLURM_JOB_NAME} with ${SLURM_JOB_NUM_NODES} nodes"
DIR=/p/lustre1/bridges7/${SLURM_JOB_NAME}/${SLURM_JOB_NUM_NODES}
mkdir -p ${DIR}
rm -f ${DIR}/clamr_mpionly
cp ../build-quartz-nbr/clamr_mpionly ${DIR}/clamr_mpionly

JOBSIZE=$(echo "scale=1;2048*sqrt(${SLURM_JOB_NUM_NODES})" | bc | cut -f 1 -d.)
echo "Running ${SLURM_JOB_NAME} with N=${JOBSIZE} on ${SLURM_NTASKS} tasks on ${SLURM_JOB_NUM_NODES} nodes."
cd ${DIR}
srun ./clamr_mpionly -n ${JOBSIZE} -l 2 -t 500 -i 100 > ${SLURM_JOB_NAME}-${SLURM_JOB_NUM_NODES}.out
echo "Finished MPI Run"
