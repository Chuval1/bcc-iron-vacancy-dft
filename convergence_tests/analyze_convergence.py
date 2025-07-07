#!/usr/bin/env python3
"""
Analyzes Quantum ESPRESSO k-point convergence test output files.
This script processes output files with the naming convention:
 fe_bulk_NxNxN_kXXX.out
Examples:
 - fe_bulk_1x1x1_k161616.out (1x1x1 supercell, 16x16x16 k-mesh)
 - fe_bulk_2x2x2_k888.out (2x2x2 supercell, 8x8x8 k-mesh)
 - fe_bulk_3x3x3_k666.out (3x3x3 supercell, 6x6x6 k-mesh)
The script extracts energy, magnetization, and timing data, calculates differences,
and generates plots for k-point convergence analysis.
Configure the SUPERCELL_CONFIGS dictionary below to match your calculations.
"""
import os
import re
import glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter, EngFormatter
from pathlib import Path

# --- Configuration ---
# Directory containing the supercell size folders (e.g., '1x1x1', '2x2x2')
BASE_OUTPUT_DIR = Path("~/fe_defects/convergence/convergence_tests").expanduser()

# Supercell sizes and k-meshes generated
SUPERCELL_CONFIGS = {
    (1, 1, 1): [(16, 16, 16), (20, 20, 20), (24, 24, 24)],
    (2, 2, 2): [(8, 8, 8), (10, 10, 10), (12, 12, 12)],
    (3, 3, 3): [(6, 6, 6), (8, 8, 8), (10, 10, 10)],
    (4, 4, 4): [(4, 4, 4), (6, 6, 6), (8, 8, 8)]
}

# Output directory for plots and CSV files
ANALYSIS_RESULTS_DIR = BASE_OUTPUT_DIR / "analysis_results2"
ANALYSIS_RESULTS_DIR.mkdir(exist_ok=True)

# Conversion factor
RY_TO_EV = 13.605693122990

def count_atoms_from_forces(content):
    """
    Count atoms from the forces section as a backup method.
    """
    force_pattern = r'atom\s+(\d+)\s+type\s+\d+\s+force\s*='
    force_matches = re.findall(force_pattern, content)
    
    if force_matches:
        return max(int(match) for match in force_matches)
    return None

def extract_number_of_atoms(content):
    """
    Extract number of atoms using multiple methods with debugging.
    """
    # Method 1: Direct extraction from "number of atoms/cell"
    n_atoms_match = re.search(r'number of atoms/cell\s+=\s+(\d+)', content)
    if n_atoms_match:
        n_atoms = int(n_atoms_match.group(1))
        print(f"  Method 1 (atoms/cell): {n_atoms}")
        return n_atoms
    
    # Method 2: From nat parameter in &system section
    nat_match = re.search(r'nat\s*=\s*(\d+)', content, re.IGNORECASE)
    if nat_match:
        n_atoms = int(nat_match.group(1))
        print(f"  Method 2 (nat parameter): {n_atoms}")
        return n_atoms
    
    # Method 3: Count from forces section
    n_atoms_forces = count_atoms_from_forces(content)
    if n_atoms_forces:
        print(f"  Method 3 (forces count): {n_atoms_forces}")
        return n_atoms_forces
    
    # Method 4: Count from atomic positions section
    atomic_pos_section = re.search(r'ATOMIC_POSITIONS.*?\n(.*?)(?=\n\s*\n|\nK_POINTS|\nCELL_PARAMETERS|$)', 
                                  content, re.DOTALL | re.IGNORECASE)
    if atomic_pos_section:
        lines = atomic_pos_section.group(1).strip().split('\n')
        valid_lines = [line for line in lines if line.strip() and not line.strip().startswith('#')]
        if valid_lines:
            n_atoms = len(valid_lines)
            print(f"  Method 4 (atomic positions): {n_atoms}")
            return n_atoms
    
    print(f"  All methods failed to find number of atoms")
    return None

def parse_qe_output(filepath: Path):
    """
    Extracts data from a Quantum ESPRESSO output file.
    Expects files in format: fe_bulk_NxNxN_kXXX.out
    """
    try:
        with open(filepath, 'r') as f:
            content = f.read()
        
        filename = filepath.name
        
        # Parse filename in expected format: fe_bulk_2x2x2_k888.out
        filename_pattern = r'fe_bulk_(\d+x\d+x\d+)_k(\d+)\.out'
        match = re.search(filename_pattern, filename)
        if not match:
            print(f"Warning: Filename doesn't match expected pattern 'fe_bulk_NxNxN_kXXX.out': {filename}")
            return None
        
        supercell_str = match.group(1)  # e.g., "2x2x2"
        k_str = match.group(2)  # e.g., "888"
        
        # Parse supercell dimensions
        sc_dims = [int(x) for x in supercell_str.split('x')]
        
        # Parse k-points from the compact string
        # For most cases, assume symmetric k-grid (kx=ky=kz)
        if len(k_str) == 1:
            # Single digit: k=8 means 8x8x8
            k_x = k_y = k_z = int(k_str)
        elif len(k_str) == 2:
            # Two digits: k=12 means 12x12x12
            k_x = k_y = k_z = int(k_str)
        elif len(k_str) == 3:
            # Three digits: k=888 means 8x8x8
            k_x = k_y = k_z = int(k_str[0])
        elif len(k_str) == 4:
            # Four digits: k=1212 means 12x12x12
            k_x = k_y = k_z = int(k_str[:2])
        elif len(k_str) == 6:
            # Six digits: k=121212 means 12x12x12
            k_x = k_y = k_z = int(k_str[:2])
        else:
            print(f"Warning: Cannot parse k-points from: {k_str}")
            return None
        
        k_density_metric = k_x  # Using kx as representative
        k_label = f"{k_x}x{k_y}x{k_z}"
        
        # Extract number of atoms
        n_atoms = extract_number_of_atoms(content)
        if n_atoms is None:
            print(f"Warning: Could not determine number of atoms from {filename}")
            return None
        else:
            print(f"  Detected {n_atoms} atoms")
        
        # Extract energy
        energy_match = re.search(r'!\s+total energy\s+=\s+([-\d.]+)\s+Ry', content)
        energy_ry = float(energy_match.group(1)) if energy_match else None
        
        # Extract magnetization - be more specific about the pattern
        # Look for the final total magnetization value (sometimes there are multiple)
        total_mag_matches = re.findall(r'total magnetization\s+=\s+([-\d.]+)\s+Bohr mag/cell', content)
        if total_mag_matches:
            # Take the last occurrence (final converged value)
            total_mag_bohr = float(total_mag_matches[-1])
            print(f"  Found total magnetization: {total_mag_bohr} Bohr mag/cell")
        else:
            total_mag_bohr = None
            print(f"  Warning: No magnetization found in {filename}")
        
        # Extract CPU time - try multiple patterns
        cpu_time_s = None
        
        # Pattern 1: From JOB DONE section
        cpu_time_match_jobdone = re.search(r'JOB DONE\..*cpu time\s*:\s*([\d.]+)s', content, re.DOTALL | re.IGNORECASE)
        if cpu_time_match_jobdone:
            cpu_time_s = float(cpu_time_match_jobdone.group(1))
        else:
            # Pattern 2: From PWSCF line at the end
            cpu_time_match_pwscf = re.search(r'PWSCF\s+:\s+([\d.]+)s\s+CPU', content)
            if cpu_time_match_pwscf:
                cpu_time_s = float(cpu_time_match_pwscf.group(1))
            else:
                # Pattern 3: Alternative PWSCF format
                cpu_time_match_alt = re.search(r'PWSCF.*?:\s*([\d.]+)s CPU', content)
                if cpu_time_match_alt:
                    cpu_time_s = float(cpu_time_match_alt.group(1))
        
        if energy_ry is None:
            print(f"Warning: Could not extract energy from {filename}. Skipping.")
            return None
        
        energy_per_atom_ry = energy_ry / n_atoms
        energy_per_atom_ev = energy_per_atom_ry * RY_TO_EV
        mag_per_atom_bohr = total_mag_bohr / n_atoms if total_mag_bohr is not None else None
        
        # Validation check for magnetization
        if mag_per_atom_bohr is not None:
            print(f"  Calculated magnetization: {mag_per_atom_bohr:.3f} μB/atom")
            if mag_per_atom_bohr > 4.0:
                print(f"  WARNING: Magnetization seems too high for Fe (expected ~2.2 μB/atom)")
                print(f"  Check: total_mag={total_mag_bohr} / n_atoms={n_atoms}")
            elif mag_per_atom_bohr < 1.0:
                print(f"  WARNING: Magnetization seems too low for ferromagnetic Fe")
        
        return {
            'filepath': str(filepath),
            'filename': filename,
            'supercell': supercell_str,
            'n_atoms': n_atoms,
            'k_x': k_x,
            'k_y': k_y,
            'k_z': k_z,
            'k_label': k_label,
            'k_density_metric': k_density_metric,
            'total_energy_ry': energy_ry,
            'energy_per_atom_ry': energy_per_atom_ry,
            'energy_per_atom_ev': energy_per_atom_ev,
            'total_mag_bohr': total_mag_bohr,
            'mag_per_atom_bohr': mag_per_atom_bohr,
            'cpu_time_s': cpu_time_s,
        }
    
    except Exception as e:
        print(f"Error processing file {filepath}: {e}")
        return None

def find_output_files():
    """
    Finds all relevant QE output files based on SUPERCELL_CONFIGS.
    Constructs expected filenames and checks if they exist.
    """
    all_files = []
    expected_files = []
    
    for sc_dims, k_meshes in SUPERCELL_CONFIGS.items():
        sc_str = f"{sc_dims[0]}x{sc_dims[1]}x{sc_dims[2]}"
        for k_tuple in k_meshes:
            # Create compact k-string (e.g., (8,8,8) -> "888")
            k_compact_str = "".join(map(str, k_tuple))
            
            # Construct expected filename
            expected_filename = f"fe_bulk_{sc_str}_k{k_compact_str}.out"
            expected_files.append(expected_filename)
            
            # Try different possible directory structures
            possible_paths = [
                # Direct in supercell directory
                BASE_OUTPUT_DIR / sc_str / expected_filename,
                # In k-point subdirectory
                BASE_OUTPUT_DIR / sc_str / f"k{'x'.join(map(str, k_tuple))}" / expected_filename,
                # Direct in base directory
                BASE_OUTPUT_DIR / expected_filename,
                # In a subdirectory named after the k-points
                BASE_OUTPUT_DIR / f"k{k_compact_str}" / expected_filename,
            ]
            
            for filepath in possible_paths:
                if filepath.is_file():
                    all_files.append(filepath)
                    print(f"Found: {filepath}")
                    break
            else:
                print(f"Missing: {expected_filename}")
                # Try to find it with a broader search
                found_alternative = False
                for search_path in [BASE_OUTPUT_DIR]:
                    for alt_file in search_path.rglob(expected_filename):
                        all_files.append(alt_file)
                        print(f"Found (alternative location): {alt_file}")
                        found_alternative = True
                        break
                    if found_alternative:
                        break
    
    print(f"\nSummary:")
    print(f"Expected {len(expected_files)} files")
    print(f"Found {len(all_files)} files")
    
    if len(all_files) < len(expected_files):
        print(f"\nMissing files:")
        found_names = {f.name for f in all_files}
        for expected in expected_files:
            if expected not in found_names:
                print(f"  {expected}")
    
    return all_files

def plot_energy_convergence_combined(df, df_with_diff):
    """Plot combined energy convergence figure with two subplots side by side."""
    plt.style.use('default')
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 6))
    
    supercell_sizes = sorted(df['supercell'].unique())
    colors = plt.cm.viridis(np.linspace(0, 1, len(supercell_sizes)))
    
    # Subplot (a): Total Energy per Atom
    for i, sc_size in enumerate(supercell_sizes):
        sc_df = df[df['supercell'] == sc_size].sort_values('k_density_metric')
        if not sc_df.empty:
            ax1.plot(sc_df['k_density_metric'], sc_df['energy_per_atom_ev'], 'o-',
                    label=f"Supercell {sc_size}", color=colors[i], markersize=7, linewidth=2)
    
    ax1.set_xlabel('K-point mesh (k×k×k grid)', fontsize=12)
    ax1.set_ylabel('Energy per Atom (eV)', fontsize=12)
    ax1.legend(fontsize=10, title="System Size", title_fontsize=11)
    ax1.grid(True, linestyle=':', alpha=0.7)
    ax1.tick_params(axis='both', which='major', labelsize=11)
    ax1.ticklabel_format(useOffset=False, style='plain', axis='y')  # Remove offset
    ax1.text(0.02, 0.98, '(a)', transform=ax1.transAxes, fontsize=14, 
             verticalalignment='top', weight='bold')
    
    # Subplot (b): Absolute Energy Difference (LINEAR SCALE)
    if not df_with_diff.empty:
        for i, sc_size in enumerate(supercell_sizes):
            sc_df = df_with_diff[df_with_diff['supercell'] == sc_size].sort_values('k_density_metric')
            if not sc_df.empty:
                ax2.plot(sc_df['k_density_metric'], sc_df['energy_diff_abs_mev'], 'o-',
                        label=f"Supercell {sc_size}", color=colors[i], markersize=7, linewidth=2)
        
        # Add 0.5 meV threshold line
        ax2.axhline(y=0.5, color='red', linestyle='--', linewidth=1.5,
                   label='0.5 meV/atom threshold')
        
        # Set y-axis limits and format for linear scale
        ax2.set_ylim(0, 0.8)  # Set appropriate limits for your data
        ax2.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
    
    ax2.set_xlabel('K-point mesh (k×k×k grid)', fontsize=12)
    ax2.set_ylabel('|ΔE| per Atom (meV)', fontsize=12)
    ax2.legend(fontsize=10, title="System Size", title_fontsize=11)
    ax2.grid(True, linestyle=':', alpha=0.7)
    ax2.tick_params(axis='both', which='major', labelsize=11)
    ax2.text(0.02, 0.98, '(b)', transform=ax2.transAxes, fontsize=14, 
             verticalalignment='top', weight='bold')
    
    plt.tight_layout()
    plot_path = ANALYSIS_RESULTS_DIR / "convergence_energy_combined.png"
    plt.savefig(plot_path, dpi=300, bbox_inches='tight')
    print(f"Saved combined energy plot: {plot_path}")
    plt.close(fig)

def plot_convergence(df, property_key, ylabel, title_prefix, filename_suffix, log_scale=False):
    """Generic plotting function for convergence tests."""
    plt.style.use('default')
    fig, ax = plt.subplots(figsize=(10, 6.5))
    
    supercell_sizes = sorted(df['supercell'].unique())
    colors = plt.cm.viridis(np.linspace(0, 1, len(supercell_sizes)))
    
    for i, sc_size in enumerate(supercell_sizes):
        sc_df = df[df['supercell'] == sc_size].sort_values('k_density_metric')
        if not sc_df.empty:
            ax.plot(sc_df['k_density_metric'], sc_df[property_key], 'o-',
                   label=f"Supercell {sc_size}", color=colors[i], markersize=7, linewidth=2)
    
    ax.set_xlabel('K-point mesh size (k×k×k grid)', fontsize=12)
    ax.set_ylabel(ylabel, fontsize=12)
    
    ax.legend(fontsize=10, title="System Size", title_fontsize=11)
    ax.grid(True, linestyle=':', alpha=0.7)
    ax.tick_params(axis='both', which='major', labelsize=11)
    
    if log_scale:
        ax.set_yscale('log')
    else:
        # Only remove offset for non-log scale plots
        ax.ticklabel_format(useOffset=False, style='plain')
    
    if property_key == 'energy_diff_abs_mev':
        ax.yaxis.set_major_formatter(FormatStrFormatter('%.3f'))
    
    plt.tight_layout()
    plot_path = ANALYSIS_RESULTS_DIR / f"convergence_{filename_suffix}.png"
    plt.savefig(plot_path, dpi=300, bbox_inches='tight')
    print(f"Saved plot: {plot_path}")
    plt.close(fig)

def main():
    """Main function to perform convergence analysis."""
    print(f"Starting convergence analysis")
    print(f"Looking for files in: {BASE_OUTPUT_DIR}")
    print(f"Results will be saved to: {ANALYSIS_RESULTS_DIR}")
    
    # Show what files we expect to find
    print(f"\nExpected file format: fe_bulk_NxNxN_kXXX.out")
    print(f"Based on SUPERCELL_CONFIGS, looking for:")
    for sc_dims, k_meshes in SUPERCELL_CONFIGS.items():
        sc_str = f"{sc_dims[0]}x{sc_dims[1]}x{sc_dims[2]}"
        for k_tuple in k_meshes:
            k_compact_str = "".join(map(str, k_tuple))
            expected_filename = f"fe_bulk_{sc_str}_k{k_compact_str}.out"
            print(f"  {expected_filename}")
    
    if not BASE_OUTPUT_DIR.exists():
        print(f"Error: BASE_OUTPUT_DIR does not exist: {BASE_OUTPUT_DIR}")
        return
    
    output_filepaths = find_output_files()
    
    if not output_filepaths:
        print("\nNo output files found. Please check:")
        print("1. BASE_OUTPUT_DIR path is correct")
        print("2. Files follow the expected naming convention")
        print("3. SUPERCELL_CONFIGS matches your actual calculations")
        return
    
    results_data = []
    for fp in output_filepaths:
        print(f"\nProcessing: {fp}")
        data = parse_qe_output(fp)
        if data:
            results_data.append(data)
            print(f"  Successfully parsed: {data['supercell']}, k={data['k_label']}, "
                  f"E={data['energy_per_atom_ev']:.6f} eV/atom, N={data['n_atoms']} atoms")
        else:
            print(f"  Failed to parse")
    
    if not results_data:
        print("No data could be parsed from the found output files.")
        return
    
    df = pd.DataFrame(results_data)
    df = df.sort_values(['supercell', 'k_density_metric']).reset_index(drop=True)
    
    # Calculate energy differences
    df['energy_diff_mev'] = pd.NA
    df['energy_diff_abs_mev'] = pd.NA
    
    for sc_size in df['supercell'].unique():
        sc_mask = df['supercell'] == sc_size
        sc_df = df[sc_mask].sort_values('k_density_metric').copy()
        
        if len(sc_df) > 1:
            # Calculate differences between consecutive k-point densities
            energies = sc_df['energy_per_atom_ry'].values
            differences_ry = np.diff(energies)  # E[i+1] - E[i]
            differences_mev = differences_ry * RY_TO_EV * 1000  # Convert to meV/atom
            
            # Assign differences to the higher k-point density entry
            sc_indices = sc_df.index[1:]  # Skip first index (no difference for first point)
            df.loc[sc_indices, 'energy_diff_mev'] = differences_mev
            df.loc[sc_indices, 'energy_diff_abs_mev'] = np.abs(differences_mev)
    
    # Save detailed results to CSV
    csv_path = ANALYSIS_RESULTS_DIR / "convergence_analysis_all_data.csv"
    df.to_csv(csv_path, index=False, float_format='%.8f')
    print(f"\nDetailed analysis data saved to: {csv_path}")
    
    # Print Data Summary
    print(f"\n--- Data Summary ---")
    
    print(f"\n--- Magnetization Summary ---")
    print("Supercell K-mesh    Mag/atom (μB)")
    print("--------  --------  -------------")
    for idx, row in df.iterrows():
        if pd.notna(row['mag_per_atom_bohr']):
            print(f"{row['supercell']:8s}  {row['k_label']:8s}  {row['mag_per_atom_bohr']:8.3f}")
        else:
            print(f"{row['supercell']:8s}  {row['k_label']:8s}  {'N/A':8s}")
    
    for sc_size in sorted(df['supercell'].unique()):
        sc_df = df[df['supercell'] == sc_size].copy()
        sc_df = sc_df.sort_values('k_density_metric')
        
        print(f"\nSupercell {sc_size}:")
        print("  K-mesh    Energy/atom (eV)  |ΔE| (meV/atom)")
        print("  -------   ----------------  --------------")
        
        for idx, (_, row) in enumerate(sc_df.iterrows()):
            energy_str = f"{row['energy_per_atom_ev']:.6f}"
            
            if pd.notna(row['energy_diff_abs_mev']):
                diff_str = f"{row['energy_diff_abs_mev']:.3f}"
            else:
                diff_str = "N/A"
            
            print(f"  {row['k_label']:8s}  {energy_str:16s}  {diff_str:14s}")
    
    # Generate Plots
    print("\n--- Generating Plots ---")
    
    # Combined energy plot
    df_with_diff = df.dropna(subset=['energy_diff_abs_mev'])
    plot_energy_convergence_combined(df, df_with_diff)
    
    # Plot 3: Magnetization per atom vs. K-point density
    df_with_mag = df.dropna(subset=['mag_per_atom_bohr'])
    if not df_with_mag.empty:
        plot_convergence(df_with_mag, 'mag_per_atom_bohr', 'Magnetization per Atom (μ$_B$)',
                        'Magnetization per Atom', 'mag_vs_kpoints')
    else:
        print("Skipping magnetization plot as no data was found.")
    
    # Plot 4: CPU Time vs. K-point density
    df_with_time = df.dropna(subset=['cpu_time_s'])
    if not df_with_time.empty:
        plot_convergence(df_with_time, 'cpu_time_s', 'CPU Time (s)',
                        'CPU Time', 'cpu_time_vs_kpoints', log_scale=True)
    else:
        print("Skipping CPU time plot as no data was found.")
    
    print(f"\nAnalysis complete.")
    print(f"All plots and CSV data saved in: {ANALYSIS_RESULTS_DIR}")

if __name__ == "__main__":
    main()