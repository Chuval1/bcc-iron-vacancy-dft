#!/bin/bash
#SBATCH --job-name=fe_conv_16a_k888
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --ntasks-per-node=32
#SBATCH --cpus-per-task=1
#SBATCH --time=01:30:00
#SBATCH --partition=regular
#SBATCH --output=fe_bulk_2x2x2_k888_%j.slurm.out
#SBATCH --error=fe_bulk_2x2x2_k888_%j.slurm.err
#SBATCH --mem=500G
# Using 32 MPI tasks for 16 atoms.

cd $SLURM_SUBMIT_DIR

module purge
module load QuantumESPRESSO/7.3.1-foss-2023a

export OMP_NUM_THREADS=1

echo "Job started at $(date)"
echo "Running on node: $(hostname)"
echo "SLURM_JOB_ID: $SLURM_JOB_ID"
echo "SLURM_NTASKS: $SLURM_NTASKS, OMP_NUM_THREADS: $OMP_NUM_THREADS"
echo "Input file: $1, Output file: $2"

if [ -z "$1" ] || [ -z "$2" ]; then
    echo "Error: SLURM script expects input_file and output_file as arguments."
    exit 1
fi

QE_INPUT_FILE=$1
QE_OUTPUT_FILE=$2

mkdir -p out

PSEUDO_PATH="./pseudo/Fe.pbe-spn-kjpaw_psl.1.0.0.UPF"
if [ ! -f "$PSEUDO_PATH" ]; then
    echo "ERROR: Pseudopotential not found at $PSEUDO_PATH (relative to $PWD)"
    exit 1
fi

echo "Running QE calculation with input=$QE_INPUT_FILE, output=$QE_OUTPUT_FILE"
echo "Command: mpirun -np $SLURM_NTASKS pw.x -nk 4 < $QE_INPUT_FILE > $QE_OUTPUT_FILE"
mpirun -np $SLURM_NTASKS pw.x -nk 4 < $QE_INPUT_FILE > $QE_OUTPUT_FILE

echo "Job finished at $(date)"
