#!/bin/bash
#SBATCH --job-name=fe_vacancy_nm
#SBATCH --nodes=1
#SBATCH --ntasks=128
#SBATCH --ntasks-per-node=128
#SBATCH --cpus-per-task=1 # For OpenMP, but we force OMP_NUM_THREADS=1
#SBATCH --time=05:30:00  # Increased time for relaxation
#SBATCH --partition=regular
#SBATCH --output=fe_vacancy_nm_3x3x3_%j.slurm.out
#SBATCH --error=fe_vacancy_nm_3x3x3_%j.slurm.err
#SBATCH --mem=500G
# Using 128 MPI tasks for 53 atoms (54 - 1 vacancy).

# Make sure we're in the correct directory (where sbatch was called from)
cd $SLURM_SUBMIT_DIR

# Load required modules
module purge # Good practice to start clean
module load QuantumESPRESSO/7.3.1-foss-2023a

# Set OMP_NUM_THREADS to 1 for pure MPI parallelization
export OMP_NUM_THREADS=1

echo "Job started at $(date)"
echo "Running on node: $(hostname)"
echo "SLURM_JOB_ID: $SLURM_JOB_ID"
echo "SLURM_NTASKS: $SLURM_NTASKS, OMP_NUM_THREADS: $OMP_NUM_THREADS"
echo "Input file: $1, Output file: $2"
echo "Running NON-MAGNETIC vacancy calculation"

# Check if proper input file is provided (passed as arg to sbatch script)
if [ -z "$1" ] || [ -z "$2" ]; then
echo "Error: SLURM script expects input_file and output_file as arguments."
echo "Usage: sbatch submit_vacancy_nm.sh vacancy_nm_input.in vacancy_nm_output.out"
exit 1
fi

QE_INPUT_FILE=$1
QE_OUTPUT_FILE=$2

# Create QE output directory if it doesn't exist (inside job dir)
mkdir -p out

# Verify pseudopotential is accessible (locally within the job directory)
PSEUDO_PATH="./pseudo/Fe.pbe-spn-kjpaw_psl.1.0.0.UPF"
if [ ! -f "$PSEUDO_PATH" ]; then
echo "ERROR: Pseudopotential not found at $PSEUDO_PATH (relative to $PWD)"
echo "Please ensure a 'pseudo' subdirectory with the UPF file exists here."
exit 1
fi

# Run the calculation with explicit file names and parallelization flags
echo "Running QE NM vacancy calculation with input=$QE_INPUT_FILE, output=$QE_OUTPUT_FILE"
echo "Command: mpirun -np $SLURM_NTASKS pw.x -nk 8 < $QE_INPUT_FILE > $QE_OUTPUT_FILE"
mpirun -np $SLURM_NTASKS pw.x -nk 8 < $QE_INPUT_FILE > $QE_OUTPUT_FILE

echo "Job finished at $(date)"