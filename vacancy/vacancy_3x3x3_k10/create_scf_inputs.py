#!/usr/bin/env python3
"""
Extract final atomic positions from QE relaxation output and create SCF inputs
"""

import re
import os
import shutil

def extract_final_positions(output_file):
    """Extract the final atomic positions from QE relaxation output"""
    with open(output_file, 'r') as f:
        content = f.read()
    
    # Look for the final coordinates section
    pattern = r'Begin final coordinates\s*\n\s*ATOMIC_POSITIONS.*?\n(.*?)\nEnd final coordinates'
    match = re.search(pattern, content, re.DOTALL)
    
    if match:
        positions = match.group(1).strip()
        return positions
    else:
        print(f"Could not find final coordinates in {output_file}")
        return None

def create_scf_input(relax_input_file, scf_input_file, final_positions):
    """Create SCF input file from relaxation input with final positions"""
    
    with open(relax_input_file, 'r') as f:
        lines = f.readlines()
    
    scf_lines = []
    skip_ions = False
    in_atomic_positions = False
    
    for line in lines:
        # Change calculation type
        if 'calculation' in line and 'relax' in line:
            scf_lines.append("    calculation = 'scf'\n")
        # Skip the &ions section entirely
        elif line.strip().startswith('&ions'):
            skip_ions = True
            continue
        elif skip_ions and line.strip().startswith('/'):
            skip_ions = False
            continue
        elif skip_ions:
            continue
        # Replace atomic positions with final relaxed ones
        elif line.strip().startswith('ATOMIC_POSITIONS'):
            scf_lines.append(line)
            scf_lines.append(final_positions + '\n')
            in_atomic_positions = True
        elif in_atomic_positions and (line.strip().startswith('K_POINTS') or 
                                      line.strip().startswith('CELL_PARAMETERS') or
                                      line.strip() == ''):
            in_atomic_positions = False
            scf_lines.append(line)
        elif not in_atomic_positions:
            scf_lines.append(line)
    
    # Write SCF input file
    with open(scf_input_file, 'w') as f:
        f.writelines(scf_lines)

def process_calculation(calc_dir):
    """Process one calculation directory"""
    print(f"\nProcessing {calc_dir}...")
    
    # Find input and output files
    input_files = [f for f in os.listdir(calc_dir) if f.endswith('_relax.in')]
    output_files = [f for f in os.listdir(calc_dir) if f.endswith('_relax.out')]
    
    if not input_files or not output_files:
        print(f"Could not find input/output files in {calc_dir}")
        return
    
    relax_input = os.path.join(calc_dir, input_files[0])
    relax_output = os.path.join(calc_dir, output_files[0])
    
    # Create SCF input filename
    scf_input = relax_input.replace('_relax.in', '_scf.in')
    
    print(f"  Extracting positions from: {relax_output}")
    final_positions = extract_final_positions(relax_output)
    
    if final_positions:
        print(f"  Creating SCF input: {scf_input}")
        create_scf_input(relax_input, scf_input, final_positions)
        
        # Also create a submission script for SCF
        create_scf_submit_script(calc_dir, os.path.basename(scf_input))
        
        print(f"  ✓ Created {scf_input}")
    else:
        print(f"  ✗ Failed to extract positions from {relax_output}")

def create_scf_submit_script(calc_dir, scf_input_file):
    """Create submission script for SCF calculation"""
    
    scf_output_file = scf_input_file.replace('.in', '.out')
    
    # Determine job name and pseudopotential from directory
    if 'FM_PAW' in calc_dir:
        job_name = 'fe_fm_vacancy_paw_scf'
        pseudo_check = 'Fe.pbe-spn-kjpaw_psl.1.0.0.UPF'
    elif 'FM_USPP' in calc_dir:
        job_name = 'fe_fm_vacancy_uspp_scf'
        pseudo_check = 'Fe.pbe-spn-rrkjus_psl.0.2.1.UPF'
    elif 'NM_PAW' in calc_dir:
        job_name = 'fe_nm_vacancy_paw_scf'
        pseudo_check = 'Fe.pbe-spn-kjpaw_psl.1.0.0.UPF'
    elif 'NM_USPP' in calc_dir:
        job_name = 'fe_nm_vacancy_uspp_scf'
        pseudo_check = 'Fe.pbe-spn-rrkjus_psl.0.2.1.UPF'
    else:
        job_name = 'fe_vacancy_scf'
        pseudo_check = 'Fe.pbe-spn-*.UPF'
    
    submit_script_content = f'''#!/bin/bash
#SBATCH --job-name={job_name}
#SBATCH --nodes=1
#SBATCH --ntasks=128
#SBATCH --ntasks-per-node=128
#SBATCH --cpus-per-task=1
#SBATCH --time=4:00:00
#SBATCH --partition=regular
#SBATCH --output={job_name}_%j.slurm.out
#SBATCH --error={job_name}_%j.slurm.err
#SBATCH --mem=500G

cd $SLURM_SUBMIT_DIR

module purge
module load QuantumESPRESSO/7.3.1-foss-2023a

export OMP_NUM_THREADS=1

echo "Job started at $(date)"
echo "Running SCF calculation on relaxed geometry"
echo "Input file: {scf_input_file}"
echo "Output file: {scf_output_file}"

# Check if pseudopotential exists
PSEUDO_PATH="./pseudo/{pseudo_check}"
if [ ! -f $PSEUDO_PATH ]; then
    echo "ERROR: Pseudopotential not found at $PSEUDO_PATH"
    exit 1
fi

mkdir -p out

echo "Running QE SCF calculation..."
mpirun -np $SLURM_NTASKS pw.x -nk 8 < {scf_input_file} > {scf_output_file}

echo "Job finished at $(date)"
'''
    
    submit_script_path = os.path.join(calc_dir, 'submit_scf.sh')
    with open(submit_script_path, 'w') as f:
        f.write(submit_script_content)
    
    # Make it executable
    os.chmod(submit_script_path, 0o755)
    print(f"  Created submission script: {submit_script_path}")

# Main execution
if __name__ == "__main__":
    # Process all calculation directories
    calc_dirs = ['FM_PAW', 'FM_USPP', 'NM_PAW', 'NM_USPP']
    
    for calc_dir in calc_dirs:
        if os.path.exists(calc_dir):
            process_calculation(calc_dir)
        else:
            print(f"Directory {calc_dir} not found")
    
    print("\n" + "="*50)
    print("SCF input files and submission scripts created!")
    print("To submit the SCF jobs, run:")
    print("  cd FM_PAW && sbatch submit_scf.sh")
    print("  cd FM_USPP && sbatch submit_scf.sh") 
    print("  cd NM_PAW && sbatch submit_scf.sh")
    print("  cd NM_USPP && sbatch submit_scf.sh")