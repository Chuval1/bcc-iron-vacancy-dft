#!/usr/bin/env python3
"""
Script to create SCF input files using relaxed coordinates from relaxation calculations

For bulk relaxations with fixed cell (ion_dynamics only):
- Final coordinates are extracted from relaxation output 
- Cell parameters are taken from relaxation input (since they're fixed)
"""

import os
import re
from pathlib import Path

def extract_final_coordinates(relax_output_file):
    """
    Extract final coordinates from relaxation output
    Cell parameters are fixed in bulk relaxations, so we'll read them from input
    """
    with open(relax_output_file, 'r') as f:
        content = f.read()
    
    # Find the final coordinates section
    coord_match = re.search(r'Begin final coordinates\s*\n(.*?)\nEnd final coordinates', 
                           content, re.DOTALL)
    if not coord_match:
        print(f"ERROR: Could not find final coordinates in {relax_output_file}")
        return None
    
    coordinates = coord_match.group(1).strip()
    return coordinates

def extract_cell_from_input(relax_input_file):
    """
    Extract cell parameters from the original relaxation input file
    (since they're fixed during relaxation)
    """
    with open(relax_input_file, 'r') as f:
        content = f.read()
    
    # Find cell parameters in input file
    cell_match = re.search(r'CELL_PARAMETERS\s+angstrom\s*\n((?:.*?\n){3})', content, re.DOTALL)
    if not cell_match:
        print(f"ERROR: Could not find CELL_PARAMETERS in {relax_input_file}")
        return None
    
    cell_params = cell_match.group(1).strip()
    return cell_params

def create_scf_input(directory, coordinates, cell_params, is_magnetic, is_paw):
    """
    Create SCF input file with the correct parameters
    """
    # Determine file names and pseudopotential
    if is_magnetic:
        mag_type = "fm"
        nspin = 2
        mag_line = "    starting_magnetization(1) = 0.4"
    else:
        mag_type = "nm"
        nspin = 1
        mag_line = ""
    
    if is_paw:
        pp_type = "paw"
        pseudopotential = "Fe.pbe-spn-kjpaw_psl.1.0.0.UPF"
    else:
        pp_type = "uspp"
        pseudopotential = "Fe.pbe-spn-rrkjus_psl.0.2.1.UPF"
    
    # Create filename
    scf_filename = f"fe_{mag_type}_bulk_3x3x3_k101010_{pp_type}_scf.in"
    scf_filepath = os.path.join(directory, scf_filename)
    
    # SCF input template
    scf_input = f"""&control
    calculation = 'scf'
    verbosity = 'high'
    restart_mode = 'from_scratch'
    prefix = 'fe_{mag_type}_bulk_3x3x3_k101010_{pp_type}_scf'
    outdir = './out'
    pseudo_dir = './pseudo'
    tprnfor = .true.
    tstress = .true.
    disk_io = 'none'
    wf_collect = .false.
    max_seconds = 19800
/
&system
    ibrav = 0
    nat = 54
    ntyp = 1
    ecutwfc = 90
    ecutrho = 1080
    occupations = 'smearing'
    smearing = 'mv'
    degauss = 0.01
    nspin = {nspin}
{mag_line}
/
&electrons
    electron_maxstep = 300
    conv_thr = 1.0d-8
    mixing_beta = 0.2
/

ATOMIC_SPECIES
    Fe  55.845  {pseudopotential}

CELL_PARAMETERS angstrom
{cell_params}

{coordinates}

K_POINTS automatic
    10 10 10 0 0 0
"""
    
    # Write the SCF input file
    with open(scf_filepath, 'w') as f:
        f.write(scf_input)
    
    print(f"✅ Created: {scf_filepath}")

def main():
    """
    Main function to process all directories and create SCF inputs
    """
    # Define the calculation directories and their properties
    calculations = [
        {
            'dir': 'FM_PAW',
            'relax_output': 'fe_fm_bulk_3x3x3_k101010_paw_relax.out',
            'relax_input': 'fe_fm_bulk_3x3x3_k101010_paw_relax.in',
            'is_magnetic': True,
            'is_paw': True
        },
        {
            'dir': 'FM_USPP',
            'relax_output': 'fe_fm_bulk_3x3x3_k101010_uspp_relax.out',
            'relax_input': 'fe_fm_bulk_3x3x3_k101010_uspp_relax.in',
            'is_magnetic': True,
            'is_paw': False
        },
        {
            'dir': 'NM_PAW',
            'relax_output': 'fe_nm_bulk_3x3x3_k101010_paw_relax.out',
            'relax_input': 'fe_nm_bulk_3x3x3_k101010_paw_relax.in',
            'is_magnetic': False,
            'is_paw': True
        },
        {
            'dir': 'NM_USPP',
            'relax_output': 'fe_nm_bulk_3x3x3_k101010_uspp_relax.out',
            'relax_input': 'fe_nm_bulk_3x3x3_k101010_uspp_relax.in',
            'is_magnetic': False,
            'is_paw': False
        }
    ]
    
    print("🔧 Creating SCF input files with relaxed coordinates...")
    print("="*60)
    
    for calc in calculations:
        directory = calc['dir']
        relax_output = calc['relax_output']
        relax_input = calc['relax_input']
        is_magnetic = calc['is_magnetic']
        is_paw = calc['is_paw']
        
        # Check if directory exists
        if not os.path.exists(directory):
            print(f"❌ Directory {directory} does not exist")
            continue
            
        relax_output_path = os.path.join(directory, relax_output)
        relax_input_path = os.path.join(directory, relax_input)
        
        # Check if files exist
        if not os.path.exists(relax_output_path):
            print(f"❌ Relaxation output {relax_output_path} does not exist")
            continue
            
        if not os.path.exists(relax_input_path):
            print(f"❌ Relaxation input {relax_input_path} does not exist")
            continue
        
        print(f"\n📁 Processing {directory}...")
        print(f"   Reading coordinates from: {relax_output}")
        print(f"   Reading cell params from: {relax_input}")
        
        # Extract coordinates from output and cell from input
        coordinates = extract_final_coordinates(relax_output_path)
        cell_params = extract_cell_from_input(relax_input_path)
        
        if coordinates is None or cell_params is None:
            print(f"❌ Failed to extract data from {directory}")
            continue
        
        # Create SCF input file
        create_scf_input(directory, coordinates, cell_params, is_magnetic, is_paw)
    
    print("\n" + "="*60)
    print("✅ SCF input file generation complete!")
    print("\nWhat the script did:")
    print("- Extracted final coordinates from relaxation outputs")
    print("- Used original cell parameters from relaxation inputs (fixed during relax)")
    print("- Created SCF inputs with proper magnetic/non-magnetic settings")
    print("\nNext steps:")
    print("1. Review the generated input files")
    print("2. Submit the SCF calculations")
    print("3. Compare energies to calculate formation energies")

if __name__ == "__main__":
    main()