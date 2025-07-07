#!/bin/bash
#SBATCH --job-name=fe_fm_vacancy_paw_scf
#SBATCH --nodes=1
#SBATCH --ntasks=128
#SBATCH --ntasks-per-node=128
#SBATCH --cpus-per-task=1
#SBATCH --time=4:00:00
#SBATCH --partition=regular
#SBATCH --output=fe_fm_vacancy_paw_scf_%j.slurm.out
#SBATCH --error=fe_fm_vacancy_paw_scf_%j.slurm.err
#SBATCH --mem=500G

cd $SLURM_SUBMIT_DIR

module purge
module load QuantumESPRESSO/7.3.1-foss-2023a

export OMP_NUM_THREADS=1

echo "Job started at $(date)"
echo "Running SCF calculation on relaxed geometry"
echo "Input file: fe_fm_vac_3x3x3_k101010_paw_scf.in"
echo "Output file: fe_fm_vac_3x3x3_k101010_paw_scf.out"

# Check if pseudopotential exists
PSEUDO_PATH="./pseudo/Fe.pbe-spn-kjpaw_psl.1.0.0.UPF"
if [ ! -f $PSEUDO_PATH ]; then
    echo "ERROR: Pseudopotential not found at $PSEUDO_PATH"
    exit 1
fi

mkdir -p out

echo "Running QE SCF calculation..."
mpirun -np $SLURM_NTASKS pw.x -nk 8 < fe_fm_vac_3x3x3_k101010_paw_scf.in > fe_fm_vac_3x3x3_k101010_paw_scf.out

echo "Job finished at $(date)"
