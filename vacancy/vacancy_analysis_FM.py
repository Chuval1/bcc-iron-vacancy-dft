#!/usr/bin/env python3
"""
Complete vacancy analysis for BCC iron from DFT calculations.
Analyzes both PAW and USPP pseudopotential results.

This script combines:
1. Vacancy formation energy analysis
2. Structural relaxation analysis
3. Magnetic moment analysis
"""

import os
import re
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator
import pandas as pd
from pathlib import Path
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.colors as mcolors
from matplotlib.patches import Rectangle, Circle
from matplotlib.gridspec import GridSpec
from matplotlib.legend_handler import HandlerPatch

# Define base paths
SCRIPT_DIR = Path(__file__).parent.absolute()
RESULTS_DIR = SCRIPT_DIR / "complete_analysis_results"
RESULTS_DIR.mkdir(exist_ok=True)

# Define all file paths for both PAW and USPP
PATHS = {
    'PAW': {
        'bulk_scf': SCRIPT_DIR / "bulk_3x3x3_k10/FM_PAW/fe_fm_bulk_3x3x3_k101010_paw_scf.out",
        'bulk_relax': SCRIPT_DIR / "bulk_3x3x3_k10/FM_PAW/fe_fm_bulk_3x3x3_k101010_paw_relax.out",
        'vac_scf': SCRIPT_DIR / "vacancy_3x3x3_k10/FM_PAW/fe_fm_vac_3x3x3_k101010_paw_scf.out",
        'vac_relax': SCRIPT_DIR / "vacancy_3x3x3_k10/FM_PAW/fe_fm_vac_3x3x3_k101010_paw_relax.out"
    },
    'USPP': {
        'bulk_scf': SCRIPT_DIR / "bulk_3x3x3_k10/FM_USPP/fe_fm_bulk_3x3x3_k101010_uspp_scf.out",
        'bulk_relax': SCRIPT_DIR / "bulk_3x3x3_k10/FM_USPP/fe_fm_bulk_3x3x3_k101010_uspp_relax.out",
        'vac_scf': SCRIPT_DIR / "vacancy_3x3x3_k10/FM_USPP/fe_fm_vac_3x3x3_k101010_uspp_scf.out",
        'vac_relax': SCRIPT_DIR / "vacancy_3x3x3_k10/FM_USPP/fe_fm_vac_3x3x3_k101010_uspp_relax.out"
    }
}

# Constants
LATTICE_CONSTANT = 2.868  # Angstrom - BCC iron lattice constant
RY_TO_EV = 13.6057039763  # Conversion factor

# Literature values for comparison (eV)
LITERATURE_VALUES = {
    "DFT-GGA (Lavrentiev et al. 2018)": 2.215,
    "DFT-GGA (Nandi et al. 2010)": 2.18,
    "DFT-GGA (Domain & Becquart 2001)": 1.95,
    "Experiment (De Schepper et al. 1983)": 2.0,
    "Experiment (Schaefer et al. 1977)": 1.6,
}

ERROR_BARS = {
    "Experiment (De Schepper et al. 1983)": 0.2,
    "Experiment (Schaefer et al. 1977)": 0.15,
    "This study - PAW": 0.0504,  # Add error bar for PAW
    "This study - USPP": 0.0504,  # Add error bar for USPP
}

# ================= PARSING FUNCTIONS =================

def extract_energy(filename):
    """Extract the total energy from a Quantum ESPRESSO output file."""
    try:
        with open(filename, 'r') as file:
            content = file.read()
        
        # Look for the final total energy line
        match = re.search(r'!\s+total energy\s+=\s+([-\d.]+)\s+Ry', content)
        if match:
            energy_ry = float(match.group(1))
            return energy_ry
        else:
            raise ValueError(f"Could not find total energy in {filename}")
    except FileNotFoundError:
        print(f"Warning: File not found: {filename}")
        return None

def parse_atomic_positions_and_moments(filename):
    """Extract atomic positions and magnetic moments from a QE output file."""
    try:
        with open(filename, 'r') as file:
            content = file.read()
    except FileNotFoundError:
        print(f"Warning: File not found: {filename}")
        return None, None, None, None
    
    # Extract cell parameters
    cell_match = re.search(r'CELL_PARAMETERS.*?\n(.*?)\n(.*?)\n(.*?)\n', 
                          content, re.DOTALL)
    if cell_match:
        cell = []
        for i in range(1, 4):
            vector = cell_match.group(i).strip().split()
            cell.append([float(x) for x in vector])
        cell = np.array(cell)
    else:
        # Default cubic cell for 3x3x3 supercell
        cell = np.eye(3) * 8.604
    
    # Extract alat if present
    alat_angstrom = None
    alat_match = re.search(r'lattice parameter \(alat\)\s*=\s*([\d.]+)\s*a\.u\.', content)
    if alat_match:
        alat_bohr = float(alat_match.group(1))
        alat_angstrom = alat_bohr * 0.529177  # Bohr to Angstrom
    
    # Find the final atomic positions
    positions_pattern = r'ATOMIC_POSITIONS.*?\n(.*?)(?:\n\n|\nEnd)'
    all_positions = re.findall(positions_pattern, content, re.DOTALL)
    
    if not all_positions:
        print(f"Warning: Could not find atomic positions in {filename}")
        return None, None, None, None
    
    # Use the last set of positions (final relaxed structure)
    positions_text = all_positions[-1]
    position_lines = positions_text.strip().split('\n')
    
    atoms = []
    positions = []
    
    for line in position_lines:
        parts = line.split()
        if len(parts) >= 4 and parts[0] in ['Fe', 'Fe1', 'Fe2']:
            atoms.append('Fe')
            x, y, z = float(parts[1]), float(parts[2]), float(parts[3])
            positions.append([x, y, z])
    
    positions = np.array(positions)
    
    # Convert positions based on format
    if 'positions (alat units)' in content and alat_angstrom:
        # Positions are already in Angstroms when given in alat units
        # No conversion needed!
        print(f"  -> Positions in alat units, alat = {alat_angstrom:.3f} Å")
    elif 'positions (crystal)' in content or (positions.size > 0 and np.max(positions) <= 1.1):
        # Crystal coordinates - need to convert
        print("  -> Converting from crystal to Cartesian coordinates")
        positions = positions @ cell
    elif 'positions (angstrom)' in content:
        # Already in Angstroms
        print("  -> Positions already in Angstroms")
    else:
        # Check if values look like they're already in Angstroms
        if positions.size > 0 and np.max(positions) > 2.0:
            print("  -> Positions appear to be in Angstroms")
        else:
            print("  -> Warning: Unknown position format")
    
    # Find magnetic moments
    moment_pattern = r'Magnetic moment per site.*?\n(.*?)(?:End of self-consistent calculation|\n\n)'
    moment_matches = re.findall(moment_pattern, content, re.DOTALL)
    
    moments = []
    if moment_matches:
        moment_text = moment_matches[-1]
        moment_lines = moment_text.strip().split('\n')
        
        for line in moment_lines:
            if "atom" in line and "magn" in line:
                parts = line.split()
                if len(parts) >= 5:
                    moment = float(parts[-1])
                    moments.append(moment)
    
    if not moments:
        moments = [0.0] * len(positions)
    
    return np.array(atoms), positions, np.array(moments), cell

# ================= ANALYSIS FUNCTIONS =================

def calculate_formation_energy(e_vac, e_bulk, n_bulk=54):
    """Calculate vacancy formation energy."""
    if e_vac is None or e_bulk is None:
        return None
    e_form_ry = e_vac - ((n_bulk - 1) / n_bulk) * e_bulk
    e_form_ev = e_form_ry * RY_TO_EV
    return e_form_ev

def identify_vacancy_position(bulk_positions, vac_positions, lattice_constant):
    """Identify the vacancy position by finding the missing atom."""
    if bulk_positions is None or vac_positions is None:
        return None
    
    # Find which atom is missing by checking distances
    vacancy_candidates = []
    
    for bulk_pos in bulk_positions:
        min_dist = float('inf')
        for vac_pos in vac_positions:
            dist = np.linalg.norm(bulk_pos - vac_pos)
            if dist < min_dist:
                min_dist = dist
        
        # If no atom within 1 Angstrom, this is likely the vacancy position
        if min_dist > 1.0:
            vacancy_candidates.append(bulk_pos)
    
    if vacancy_candidates:
        # Return the first vacancy candidate
        return vacancy_candidates[0]
    else:
        # Fallback: center of the cell
        return np.mean(bulk_positions, axis=0)

def match_atoms_between_structures(bulk_positions, vac_positions):
    """Match atoms between bulk and vacancy structures."""
    matches = []
    
    for i, bulk_pos in enumerate(bulk_positions):
        best_match = -1
        min_dist = float('inf')
        
        for j, vac_pos in enumerate(vac_positions):
            dist = np.linalg.norm(bulk_pos - vac_pos)
            if dist < min_dist:
                min_dist = dist
                best_match = j
        
        if min_dist < 1.0:  # Only match if within reasonable distance
            matches.append({
                'bulk_index': i,
                'vac_index': best_match,
                'distance': min_dist
            })
    
    return pd.DataFrame(matches)

def calculate_displacement_vectors(bulk_positions, vac_positions):
    """Calculate displacement vectors for matched atoms."""
    matches = match_atoms_between_structures(bulk_positions, vac_positions)
    
    print(f"\nDEBUG: Displacement calculation")
    print(f"Bulk atoms: {len(bulk_positions)}, Vacancy atoms: {len(vac_positions)}")
    print(f"Matched atoms: {len(matches)}")
    
    displacements = []
    for _, row in matches.iterrows():
        bulk_idx = int(row['bulk_index'])
        vac_idx = int(row['vac_index'])
        
        bulk_pos = bulk_positions[bulk_idx]
        vac_pos = vac_positions[vac_idx]
        displacement = vac_pos - bulk_pos
        
        displacements.append({
            'bulk_index': bulk_idx,
            'vac_index': vac_idx,
            'bulk_position': bulk_pos,
            'vac_position': vac_pos,
            'displacement': displacement,
            'displacement_magnitude': np.linalg.norm(displacement)
        })
    
    df = pd.DataFrame(displacements)
    
    # Debug: Show displacement distribution
    if len(df) > 0:
        print(f"Displacement magnitudes - Min: {df['displacement_magnitude'].min():.4f}, Max: {df['displacement_magnitude'].max():.4f}")
        print(f"Number of atoms with displacement > 0.001 Å: {len(df[df['displacement_magnitude'] > 0.001])}")
    
    return df

# ================= PLOTTING FUNCTIONS =================

def plot_formation_energy_comparison(paw_energy, uspp_energy):
    """Create comparison plot of formation energies with improved error bars and legend."""
    # Update with our calculated values
    all_values = LITERATURE_VALUES.copy()
    if paw_energy is not None:
        all_values["This study - PAW\n54 atoms, 10×10×10 k-mesh"] = paw_energy
    if uspp_energy is not None:
        all_values["This study - USPP\n54 atoms, 10×10×10 k-mesh"] = uspp_energy
    
    # Create dataframe
    df = pd.DataFrame({
        'Source': list(all_values.keys()),
        'Formation Energy (eV)': list(all_values.values())
    })
    
    # Add error bars - updated mapping logic
    def get_error(source):
        # Check for exact match first
        if source in ERROR_BARS:
            return ERROR_BARS[source]
        # Then check for partial matches (for multiline sources)
        for key, value in ERROR_BARS.items():
            if key in source:
                return value
        return 0
    
    df['Error'] = df['Source'].map(get_error)
    
    # Sort by formation energy
    df = df.sort_values('Formation Energy (eV)')
    
    # Create figure
    fig, ax = plt.subplots(figsize=(12, 8), dpi=300)
    
    # Color scheme
    colors = []
    labels = []
    for source in df['Source']:
        if 'This study - PAW' in source:
            colors.append('#ff7f0e')
            labels.append('PAW')
        elif 'This study - USPP' in source:
            colors.append('#d62728')
            labels.append('USPP')
        elif 'Experiment' in source:
            colors.append('#2ca02c')
            labels.append('Experiment')
        else:
            colors.append('#1f77b4')
            labels.append('DFT Literature')
    
    # Plot bars
    bars = ax.barh(df['Source'], df['Formation Energy (eV)'], color=colors, alpha=0.8)
    
    # Add error bars with offset to avoid obscuring values
    for i, (_, row) in enumerate(df.iterrows()):
        if row['Error'] > 0:
            # Position error bars slightly above the bar to avoid obscuring text
            ax.errorbar(row['Formation Energy (eV)'], i + 0.1, xerr=row['Error'], 
                       fmt='none', ecolor='black', capsize=5, capthick=2)
    
    # Annotations - positioned to avoid error bars
    for i, bar in enumerate(bars):
        value = df['Formation Energy (eV)'].iloc[i]
        error = df['Error'].iloc[i]
        if pd.notna(value):
            text = f'{value:.3f} eV'
            # Position text slightly below the bar if there's an error bar
            y_offset = -0.1 if error > 0 else 0
            ax.text(value + 0.02, bar.get_y() + bar.get_height()/2 + y_offset, 
                   text, va='center', fontsize=10, fontweight='bold')
    
    # Create legend
    unique_labels = []
    unique_colors = []
    unique_names = []
    for color, label in zip(colors, labels):
        if label not in unique_names:
            unique_names.append(label)
            unique_colors.append(color)
            if label == 'PAW':
                unique_labels.append('This study - PAW')
            elif label == 'USPP':
                unique_labels.append('This study - USPP')
            elif label == 'DFT Literature':
                unique_labels.append('DFT Literature')
            else:
                unique_labels.append(label)
    
    # Create custom legend handles
    legend_elements = [Rectangle((0, 0), 1, 1, facecolor=color, alpha=0.8, edgecolor='none') 
                      for color in unique_colors]
    
    # Add legend outside plot area
    ax.legend(legend_elements, unique_labels, loc='center left', bbox_to_anchor=(1, 0.5),
             frameon=True, fancybox=True, shadow=True, fontsize=11)
    
    ax.set_xlabel('Vacancy Formation Energy (eV)', fontsize=12)
    ax.set_xlim(1.4, 2.5)
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.grid(True, linestyle='--', alpha=0.7, which='major')
    ax.grid(True, linestyle=':', alpha=0.4, which='minor')
    
    plt.tight_layout()
    plt.savefig(RESULTS_DIR / 'formation_energy_comparison.png', dpi=300, bbox_inches='tight')
    plt.savefig(RESULTS_DIR / 'formation_energy_comparison.pdf', bbox_inches='tight')
    plt.close()

def plot_radial_displacement_comparison(paw_data, uspp_data):
    """Plot radial displacement for both PAW and USPP."""
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6), dpi=300)
    
    for ax, (label, data) in zip([ax1, ax2], [('PAW', paw_data), ('USPP', uspp_data)]):
        if data is None:
            ax.text(0.5, 0.5, f'No {label} data available', 
                   transform=ax.transAxes, ha='center', va='center')
            continue
        
        displacements, vacancy_pos = data
        
        # Calculate distances from vacancy
        distances = []
        for _, row in displacements.iterrows():
            bulk_pos = row['bulk_position']
            dist = np.linalg.norm(bulk_pos - vacancy_pos)
            distances.append(dist)
        
        # Plot without colormap
        ax.scatter(distances, displacements['displacement_magnitude'],
                  color='#1f77b4', alpha=0.7, s=50, edgecolor='k', linewidth=0.5)
        
        # Fit curve
        from scipy.optimize import curve_fit
        
        def inverse_power(x, a, b, c):
            return a * np.power(x + 0.1, -b) + c  # Add small offset to avoid division by zero
        
        try:
            mask = displacements['displacement_magnitude'] > 0.001
            if np.sum(mask) > 3:  # Need at least 3 points to fit
                x_data = np.array(distances)[mask]
                y_data = displacements.loc[mask, 'displacement_magnitude'].values
                
                popt, _ = curve_fit(inverse_power, x_data, y_data, 
                                   p0=[0.1, 1.0, 0.0], maxfev=10000)
                
                x_fit = np.linspace(np.min(x_data), np.max(x_data), 100)
                y_fit = inverse_power(x_fit, *popt)
                ax.plot(x_fit, y_fit, 'r--', linewidth=2, 
                       label=f'Fit: {popt[0]:.3f}/r^{popt[1]:.2f} + {popt[2]:.3f}')
                ax.legend(fontsize=10)
        except:
            pass
        
        ax.set_xlabel('Distance from Vacancy (Å)', fontsize=12)
        ax.set_ylabel('Displacement Magnitude (Å)', fontsize=12)
        ax.set_title(f'{label} - Radial Atomic Displacement', fontsize=13)
        ax.grid(True, linestyle='--', alpha=0.7)
    
    plt.tight_layout()
    plt.savefig(RESULTS_DIR / 'radial_displacement_comparison.png', dpi=300, bbox_inches='tight')
    plt.savefig(RESULTS_DIR / 'radial_displacement_comparison.pdf', bbox_inches='tight')
    plt.close()

def plot_magnetic_moment_3d(paw_data, uspp_data, bulk_data=None, bulk_moment_ref=None):
    """Plot 3D magnetic moment visualization with neighbor shells."""
    # Create figure with 3 subplots
    fig = plt.figure(figsize=(18, 6), dpi=300)
    
    # Calculate neighbor shell radii based on BCC structure
    a = LATTICE_CONSTANT
    neighbor_radii = {
        0: a * np.sqrt(3) / 2,      # 1st NN: ~2.48 Å
        1: a,                        # 2nd NN: ~2.87 Å
        2: a * np.sqrt(2),          # 3rd NN: ~4.06 Å
    }
    
    # First plot: Bulk structure
    ax1 = fig.add_subplot(1, 3, 1, projection='3d')
    
    if bulk_data is not None:
        bulk_positions, bulk_moments = bulk_data
        
        # Check if bulk moments are essentially uniform
        moment_range = np.max(bulk_moments) - np.min(bulk_moments)
        moment_mean = np.mean(bulk_moments)
        
        if moment_range < 0.001:  # Essentially uniform
            # Use middle of colormap for all atoms
            norm = mcolors.Normalize(vmin=2.2, vmax=2.4)  # Fixed range around bulk value
            cmap = plt.cm.plasma
            color = cmap(norm(moment_mean))
            
            # Plot atoms with uniform color
            for i, (pos, moment) in enumerate(zip(bulk_positions, bulk_moments)):
                size = 20 + moment * 30  # Keep original size scaling
                ax1.scatter(pos[0], pos[1], pos[2], s=size, c=[color], alpha=0.8)
            
            # Still add colorbar but with fixed range
            sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
            sm.set_array([])
            cbar = plt.colorbar(sm, ax=ax1, pad=0.1)
            cbar.set_label('Magnetic Moment (μ$_B$)', fontsize=10)
            
            # Add text showing uniform moment value
            ax1.text2D(0.02, 0.98, f'Uniform: {moment_mean:.3f} μ$_B$', 
                      transform=ax1.transAxes, fontsize=10,
                      bbox=dict(boxstyle='round', facecolor='white', alpha=0.8),
                      verticalalignment='top')
        else:
            # Original code for non-uniform moments
            norm = mcolors.Normalize(vmin=np.min(bulk_moments), vmax=np.max(bulk_moments))
            cmap = plt.cm.plasma
            
            for i, (pos, moment) in enumerate(zip(bulk_positions, bulk_moments)):
                size = 20 + moment * 30
                color = cmap(norm(moment))
                ax1.scatter(pos[0], pos[1], pos[2], s=size, c=[color], alpha=0.8)
            
            sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
            sm.set_array([])
            cbar = plt.colorbar(sm, ax=ax1, pad=0.1)
            cbar.set_label('Magnetic Moment (μ$_B$)', fontsize=10)
        
        # Set axis labels
        ax1.set_xlabel('X (Å)', fontsize=12)
        ax1.set_ylabel('Y (Å)', fontsize=12)
        ax1.set_zlabel('Z (Å)', fontsize=12)
        
        # Set aspect ratio to be equal
        ax1.set_box_aspect([1, 1, 1])
        
        # Set viewing angle
        ax1.view_init(elev=20, azim=45)
    else:
        ax1.text2D(0.5, 0.5, 'No bulk data available', 
                   transform=ax1.transAxes, ha='center', va='center')
    
    # Second and third plots: Vacancy structures
    for idx, (label, data) in enumerate([('PAW', paw_data), ('USPP', uspp_data)]):
        ax = fig.add_subplot(1, 3, idx+2, projection='3d')
        
        if data is None:
            ax.text2D(0.5, 0.5, f'No {label} data available', 
                     transform=ax.transAxes, ha='center', va='center')
            continue
        
        vac_positions, vac_moments, vacancy_pos = data
        
        # For vacancy, we expect more variation, so use full colorbar
        # Set colorbar range to emphasize the variations
        vmin = max(np.min(vac_moments), 2.20)  # Don't go below bulk baseline
        vmax = min(np.max(vac_moments), 2.50)  # Reasonable upper limit
        
        norm = mcolors.Normalize(vmin=vmin, vmax=vmax)
        cmap = plt.cm.plasma
        
        # Plot atoms with sizes and colors based on magnetic moments
        for i, (pos, moment) in enumerate(zip(vac_positions, vac_moments)):
            size = 30 + (moment - 2.2) * 100  # Scale relative to bulk moment
            color = cmap(norm(moment))
            ax.scatter(pos[0], pos[1], pos[2], s=size, c=[color], 
                      alpha=0.8, edgecolor='k', linewidth=0.5)
        
        # Plot vacancy position
        ax.scatter(vacancy_pos[0], vacancy_pos[1], vacancy_pos[2],
                  color='red', s=150, marker='x', linewidth=3, label='Vacancy Site')
        
        # Plot shell boundaries as wire spheres (make them more subtle)
        for shell_num, radius in neighbor_radii.items():
            u, v = np.mgrid[0:2*np.pi:20j, 0:np.pi:10j]
            x = radius * np.cos(u) * np.sin(v) + vacancy_pos[0]
            y = radius * np.sin(u) * np.sin(v) + vacancy_pos[1]
            z = radius * np.cos(v) + vacancy_pos[2]
            
            # Use different colors for each shell with very low alpha
            shell_colors = ['cyan', 'yellow', 'magenta']
            ax.plot_surface(x, y, z, color=shell_colors[shell_num], alpha=0.05)
        
        # Add a colorbar
        sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
        sm.set_array([])
        cbar = plt.colorbar(sm, ax=ax, pad=0.1)
        cbar.set_label('Magnetic Moment (μ$_B$)', fontsize=10)
        
        # Calculate first neighbor enhancement
        # First neighbors in BCC are at distance a*sqrt(3)/2
        first_nn_dist = LATTICE_CONSTANT * np.sqrt(3) / 2
        tolerance = 0.3  # Angstrom tolerance for identifying neighbors
        
        # Find first neighbor atoms
        first_neighbor_moments = []
        for i, pos in enumerate(vac_positions):
            dist = np.linalg.norm(pos - vacancy_pos)
            if abs(dist - first_nn_dist) < tolerance:
                first_neighbor_moments.append(vac_moments[i])
        
        if first_neighbor_moments:
            avg_first_nn_moment = np.mean(first_neighbor_moments)
            # Use the reference bulk moment if provided, otherwise use default
            ref_moment = bulk_moment_ref if bulk_moment_ref is not None else 2.2
            enhancement = ((avg_first_nn_moment - ref_moment) / ref_moment) * 100
            text = f'1st NN enhancement: {enhancement:+.1f}%\n({len(first_neighbor_moments)} atoms)'
        else:
            text = '1st NN: No atoms found'
            
        ax.text2D(0.02, 0.98, text, 
                 transform=ax.transAxes, fontsize=10,
                 bbox=dict(boxstyle='round', facecolor='white', alpha=0.8),
                 verticalalignment='top')
        
        # Set axis labels
        ax.set_xlabel('X (Å)', fontsize=12)
        ax.set_ylabel('Y (Å)', fontsize=12)
        ax.set_zlabel('Z (Å)', fontsize=12)
        
        # Set aspect ratio to be equal
        ax.set_box_aspect([1, 1, 1])
        
        # Set viewing angle
        ax.view_init(elev=20, azim=45)
        
        # Add legend for vacancy marker
        ax.legend(loc='upper right', fontsize=10)
    
    plt.tight_layout()
    plt.savefig(RESULTS_DIR / 'magnetic_moment_3d.png', dpi=300, bbox_inches='tight')
    plt.savefig(RESULTS_DIR / 'magnetic_moment_3d.pdf', bbox_inches='tight')
    plt.close()

def plot_magnetic_moment_radial(paw_data, uspp_data, bulk_moment_ref=None):
    """Plot radial distribution of magnetic moments with baseline reference."""
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6), dpi=300)
    
    # Use the reference bulk moment if provided
    baseline_moment = bulk_moment_ref if bulk_moment_ref is not None else 2.2
    
    for ax, (label, data) in zip([ax1, ax2], [('PAW', paw_data), ('USPP', uspp_data)]):
        if data is None:
            ax.text(0.5, 0.5, f'No {label} data available', 
                   transform=ax.transAxes, ha='center', va='center')
            continue
        
        vac_positions, vac_moments, vacancy_pos = data
        
        # Calculate distances from vacancy
        distances = [np.linalg.norm(pos - vacancy_pos) for pos in vac_positions]
        
        # Create scatter plot without colormap
        ax.scatter(distances, vac_moments, color='#1f77b4', 
                  alpha=0.7, s=80, edgecolor='k', linewidth=0.5)
        
        # Add baseline magnetic moment line
        ax.axhline(y=baseline_moment, color='red', linestyle='--', linewidth=2,
                  label=f'Bulk Fe: {baseline_moment:.3f} μ$_B$')
        
        # Add neighbor shell boundaries with more visible styling
        a = LATTICE_CONSTANT
        neighbor_boundaries = {
            '1st-2nd': (a * np.sqrt(3) / 2 + a) / 2,           # ~2.68 Å
            '2nd-3rd': (a + a * np.sqrt(2)) / 2,               # ~3.47 Å
            '3rd-4th': (a * np.sqrt(2) + a * np.sqrt(11) / 2) / 2  # ~4.41 Å
        }
        
        for boundary_name, boundary_pos in neighbor_boundaries.items():
            ax.axvline(x=boundary_pos, color='gray', linestyle=':', alpha=0.7, linewidth=1.5)
        

        
        ax.set_xlabel('Distance from Vacancy (Å)', fontsize=12)
        ax.set_ylabel('Magnetic Moment (μ$_B$)', fontsize=12)
        ax.grid(True, linestyle='--', alpha=0.7)
        ax.legend(fontsize=11, loc='lower right')
        
        # Set axis limits for better visualization
        ax.set_xlim(2.0, 8.0)  # Start at 2.0 to show relaxed first neighbors
        ax.set_ylim(2.15, 2.55)
    
    plt.tight_layout()
    plt.savefig(RESULTS_DIR / 'magnetic_moment_radial.png', dpi=300, bbox_inches='tight')
    plt.savefig(RESULTS_DIR / 'magnetic_moment_radial.pdf', bbox_inches='tight')
    plt.close()

def create_summary_table(results):
    """Create a comprehensive summary table."""
    summary_data = []
    
    for pp_type in ['PAW', 'USPP']:
        if pp_type in results:
            data = results[pp_type]
            row = {
                'Pseudopotential': pp_type,
                'Formation Energy (eV)': f"{data['formation_energy']:.3f}" if data['formation_energy'] else 'N/A',
                'Bulk Energy (Ry)': f"{data['bulk_energy']:.6f}" if data['bulk_energy'] else 'N/A',
                'Vacancy Energy (Ry)': f"{data['vac_energy']:.6f}" if data['vac_energy'] else 'N/A',
                'Avg Bulk Moment (μ_B)': f"{data['avg_bulk_moment']:.3f}" if data['avg_bulk_moment'] else 'N/A',
                'Avg Vac Moment (μ_B)': f"{data['avg_vac_moment']:.3f}" if data['avg_vac_moment'] else 'N/A',
                'Max Displacement (Å)': f"{data['max_displacement']:.4f}" if data['max_displacement'] else 'N/A',
                'Mean Displacement (Å)': f"{data['mean_displacement']:.4f}" if data['mean_displacement'] else 'N/A'
            }
            summary_data.append(row)
    
    df = pd.DataFrame(summary_data)
    
    # Save to CSV
    df.to_csv(RESULTS_DIR / 'summary_table.csv', index=False)
    
    # Create a figure with the table
    fig, ax = plt.subplots(figsize=(12, 4), dpi=300)
    ax.axis('tight')
    ax.axis('off')
    
    table = ax.table(cellText=df.values, colLabels=df.columns, 
                    cellLoc='center', loc='center')
    table.auto_set_font_size(False)
    table.set_fontsize(10)
    table.scale(1.2, 1.5)
    
    # Style the header
    for i in range(len(df.columns)):
        table[(0, i)].set_facecolor('#4CAF50')
        table[(0, i)].set_text_props(weight='bold', color='white')
    
    # Alternate row colors
    for i in range(1, len(df) + 1):
        for j in range(len(df.columns)):
            if i % 2 == 0:
                table[(i, j)].set_facecolor('#f0f0f0')
    
    plt.savefig(RESULTS_DIR / 'summary_table.png', dpi=300, bbox_inches='tight')
    plt.savefig(RESULTS_DIR / 'summary_table.pdf', bbox_inches='tight')
    plt.close()
    
    return df

def main():
    """Main analysis function."""
    print("Starting comprehensive vacancy analysis...")
    print(f"Results will be saved to: {RESULTS_DIR}")
    
    results = {}
    bulk_moment_reference = None  # To store the bulk PAW moment for reference
    
    # Process both PAW and USPP
    for pp_type in ['PAW', 'USPP']:
        print(f"\n{'='*50}")
        print(f"Processing {pp_type} calculations...")
        print(f"{'='*50}")
        
        paths = PATHS[pp_type]
        
        # Extract energies from SCF calculations
        bulk_energy = extract_energy(paths['bulk_scf'])
        vac_energy = extract_energy(paths['vac_scf'])
        
        # Calculate formation energy
        formation_energy = calculate_formation_energy(vac_energy, bulk_energy)
        
        if formation_energy is not None:
            print(f"{pp_type} - Vacancy Formation Energy: {formation_energy:.3f} eV")
        
        # Parse atomic positions and moments from relaxed structures
        bulk_data = parse_atomic_positions_and_moments(paths['bulk_relax'])
        vac_data = parse_atomic_positions_and_moments(paths['vac_relax'])
        
        if bulk_data[0] is not None and vac_data[0] is not None:
            bulk_atoms, bulk_positions, bulk_moments, bulk_cell = bulk_data
            vac_atoms, vac_positions, vac_moments, vac_cell = vac_data
            
            # Calculate structural properties
            vacancy_pos = identify_vacancy_position(bulk_positions, vac_positions, LATTICE_CONSTANT)
            displacements = calculate_displacement_vectors(bulk_positions, vac_positions)
            
            # Store results
            results[pp_type] = {
                'formation_energy': formation_energy,
                'bulk_energy': bulk_energy,
                'vac_energy': vac_energy,
                'bulk_positions': bulk_positions,
                'vac_positions': vac_positions,
                'bulk_moments': bulk_moments,
                'vac_moments': vac_moments,
                'vacancy_pos': vacancy_pos,
                'displacements': displacements,
                'avg_bulk_moment': np.mean(bulk_moments) if len(bulk_moments) > 0 else None,
                'avg_vac_moment': np.mean(vac_moments) if len(vac_moments) > 0 else None,
                'max_displacement': displacements['displacement_magnitude'].max() if not displacements.empty else None,
                'mean_displacement': displacements['displacement_magnitude'].mean() if not displacements.empty else None
            }
            
            # Store PAW bulk moment as reference
            if pp_type == 'PAW' and results[pp_type]['avg_bulk_moment'] is not None:
                bulk_moment_reference = results[pp_type]['avg_bulk_moment']
            
            print(f"{pp_type} - Avg bulk moment: {results[pp_type]['avg_bulk_moment']:.3f} μ$_B$")
            print(f"{pp_type} - Avg vacancy moment: {results[pp_type]['avg_vac_moment']:.3f} μ$_B$")
            print(f"{pp_type} - Max displacement: {results[pp_type]['max_displacement']:.4f} Å")
        else:
            print(f"Warning: Could not parse atomic data for {pp_type}")
            results[pp_type] = {
                'formation_energy': formation_energy,
                'bulk_energy': bulk_energy,
                'vac_energy': vac_energy,
                'bulk_positions': None,
                'vac_positions': None,
                'bulk_moments': None,
                'vac_moments': None,
                'vacancy_pos': None,
                'displacements': None,
                'avg_bulk_moment': None,
                'avg_vac_moment': None,
                'max_displacement': None,
                'mean_displacement': None
            }
    
    # Create comparison plots
    print("\nGenerating comparison plots...")
    
    # DEBUG: Let's understand the data structure
    print("\n" + "="*50)
    print("DEBUG: Data structure analysis")
    print("="*50)
    
    if 'PAW' in results:
        print("\nPAW Data:")
        if results['PAW']['bulk_positions'] is not None:
            print(f"  Bulk positions shape: {results['PAW']['bulk_positions'].shape}")
        if results['PAW']['vac_positions'] is not None:
            print(f"  Vacancy positions shape: {results['PAW']['vac_positions'].shape}")
        if results['PAW']['displacements'] is not None:
            print(f"  Displacements dataframe shape: {results['PAW']['displacements'].shape}")
            print(f"  Displacements dataframe columns: {list(results['PAW']['displacements'].columns)}")
            
            # Show displacement distribution
            disp_df = results['PAW']['displacements']
            if 'displacement_magnitude' in disp_df.columns:
                print(f"\n  Displacement distribution:")
                print(f"    Min: {disp_df['displacement_magnitude'].min():.4f} Å")
                print(f"    Max: {disp_df['displacement_magnitude'].max():.4f} Å")
                print(f"    Mean: {disp_df['displacement_magnitude'].mean():.4f} Å")
                print(f"    Atoms with disp > 0.001 Å: {len(disp_df[disp_df['displacement_magnitude'] > 0.001])}")
                print(f"    Atoms with disp > 0.01 Å: {len(disp_df[disp_df['displacement_magnitude'] > 0.01])}")
                
                # Show distances from vacancy for displaced atoms
                if results['PAW']['vacancy_pos'] is not None:
                    print(f"\n  Distance from vacancy for atoms with disp > 0.01 Å:")
                    significant_disp = disp_df[disp_df['displacement_magnitude'] > 0.01]
                    for idx, row in significant_disp.iterrows():
                        dist = np.linalg.norm(row['bulk_position'] - results['PAW']['vacancy_pos'])
                        print(f"    Atom {idx}: distance = {dist:.3f} Å, displacement = {row['displacement_magnitude']:.4f} Å")
    
    # Formation energy comparison
    paw_fe = results.get('PAW', {}).get('formation_energy')
    uspp_fe = results.get('USPP', {}).get('formation_energy')
    plot_formation_energy_comparison(paw_fe, uspp_fe)
    
    # Structural analysis plots
    paw_struct_data = None
    uspp_struct_data = None
    
    if 'PAW' in results and results['PAW']['displacements'] is not None:
        paw_struct_data = (results['PAW']['displacements'], results['PAW']['vacancy_pos'])
    
    if 'USPP' in results and results['USPP']['displacements'] is not None:
        uspp_struct_data = (results['USPP']['displacements'], results['USPP']['vacancy_pos'])
    
    plot_radial_displacement_comparison(paw_struct_data, uspp_struct_data)
    
    # Magnetic analysis plots - now separated
    paw_mag_data = None
    uspp_mag_data = None
    bulk_mag_data = None
    
    # Debug: Check what data we have
    print("\n" + "="*50)
    print("DEBUG: Checking magnetic data before plotting")
    print("="*50)
    
    # Prepare bulk data for 3D plot (using PAW bulk as reference)
    if 'PAW' in results and results['PAW']['bulk_positions'] is not None:
        bulk_mag_data = (results['PAW']['bulk_positions'], results['PAW']['bulk_moments'])
        print(f"PAW bulk data: {len(results['PAW']['bulk_positions'])} atoms")
    
    if 'PAW' in results and results['PAW']['vac_positions'] is not None:
        paw_mag_data = (results['PAW']['vac_positions'], 
                       results['PAW']['vac_moments'],
                       results['PAW']['vacancy_pos'])
        print(f"PAW vacancy data: {len(results['PAW']['vac_positions'])} atoms")
        print(f"PAW vacancy position: {results['PAW']['vacancy_pos']}")
        
        # Show first few distances from vacancy
        distances = []
        for pos in results['PAW']['vac_positions'][:10]:
            dist = np.linalg.norm(pos - results['PAW']['vacancy_pos'])
            distances.append(dist)
        print(f"First 10 distances from vacancy: {[f'{d:.3f}' for d in sorted(distances)]}")
    
    if 'USPP' in results and results['USPP']['vac_positions'] is not None:
        uspp_mag_data = (results['USPP']['vac_positions'], 
                        results['USPP']['vac_moments'],
                        results['USPP']['vacancy_pos'])
        print(f"USPP vacancy data: {len(results['USPP']['vac_positions'])} atoms")
    
    # Create separate 3D and radial plots with bulk moment reference
    plot_magnetic_moment_3d(paw_mag_data, uspp_mag_data, bulk_data=bulk_mag_data, 
                           bulk_moment_ref=bulk_moment_reference)
    plot_magnetic_moment_radial(paw_mag_data, uspp_mag_data, 
                               bulk_moment_ref=bulk_moment_reference)
    
    # Create summary table
    summary_df = create_summary_table(results)
    
    print("\n" + "="*50)
    print("ANALYSIS COMPLETE")
    print("="*50)
    print("\nSummary of Results:")
    print(summary_df.to_string(index=False))
    
    # Save key values for manuscript
    key_values = {
        'PAW_formation_energy_eV': paw_fe,
        'USPP_formation_energy_eV': uspp_fe,
        'formation_energy_difference_eV': abs(paw_fe - uspp_fe) if paw_fe and uspp_fe else None,
        'PAW_avg_bulk_moment_muB': results.get('PAW', {}).get('avg_bulk_moment'),
        'USPP_avg_bulk_moment_muB': results.get('USPP', {}).get('avg_bulk_moment'),
        'PAW_avg_vac_moment_muB': results.get('PAW', {}).get('avg_vac_moment'),
        'USPP_avg_vac_moment_muB': results.get('USPP', {}).get('avg_vac_moment'),
        'PAW_max_displacement_A': results.get('PAW', {}).get('max_displacement'),
        'USPP_max_displacement_A': results.get('USPP', {}).get('max_displacement')
    }
    
    # Save to JSON for easy loading
    import json
    with open(RESULTS_DIR / 'key_manuscript_values.json', 'w') as f:
        json.dump(key_values, f, indent=2)
    
    print(f"\nAll results saved to: {RESULTS_DIR}")
    print("\nFiles generated:")
    for file in RESULTS_DIR.glob('*'):
        print(f"  - {file.name}")

if __name__ == "__main__":
    main()

# Quick test to verify positions
atoms, positions, moments, cell = parse_atomic_positions_and_moments(
    "bulk_3x3x3_k10/FM_PAW/fe_fm_bulk_3x3x3_k101010_paw_relax.out"
)

if positions is not None:
    print(f"Number of atoms: {len(positions)}")
    print(f"Min position: {np.min(positions):.3f} Å")
    print(f"Max position: {np.max(positions):.3f} Å")
    print(f"Expected max for 3×3×3: ~8.604 Å")
    
    # Check if it's BCC structure
    # First atom should be at origin, second at ~(1.434, 1.434, 1.434)
    print(f"\nFirst atom: {positions[0]}")
    print(f"Second atom: {positions[1]}")
    print(f"Expected for BCC: [0,0,0] and [1.434,1.434,1.434]")