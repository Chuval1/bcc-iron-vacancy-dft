#!/bin/bash
#SBATCH --job-name=fe_conv
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --ntasks-per-node=4
#SBATCH --time=168:00:00
#SBATCH --partition=regular
#SBATCH --output=fe_conv_%j.out
#SBATCH --error=fe_conv_%j.err
#SBATCH --mem=16G
# Using 4 CPUs (2 per atom) for 2 atoms

# Make sure we're in the correct directory
cd $SLURM_SUBMIT_DIR

# Load required modules
module load QuantumESPRESSO/7.3.1-foss-2023a

# Echo start time and node info for debugging
echo "Job started at $(date)"
echo "Running on node: $(hostname)"
echo "Number of tasks: $SLURM_NTASKS"

# Check if proper input file is provided
if [ -z "$1" ] || [ -z "$2" ]; then
    echo "Error: Input and output files not specified"
    echo "Usage: sbatch submit.sh input_file output_file"
    exit 1
fi

INPUT_FILE=$1
OUTPUT_FILE=$2

# Create output directory
mkdir -p out

# Verify pseudopotential is accessible
if [ ! -f pseudo/Fe.pbe-spn-kjpaw_psl.1.0.0.UPF ]; then
    echo "ERROR: Pseudopotential not found in pseudo/ directory"
    exit 1
fi

# Run the calculation with explicit file names
echo "Running calculation with input=$INPUT_FILE, output=$OUTPUT_FILE"
mpirun -np $SLURM_NTASKS pw.x < $INPUT_FILE > $OUTPUT_FILE

# Echo completion time
echo "Job finished at $(date)"
