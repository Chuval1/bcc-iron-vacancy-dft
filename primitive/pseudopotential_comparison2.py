#!/usr/bin/env python3
"""
Pseudopotential comparison for BCC iron: Analyzes and visualizes DFT results across
three different pseudopotentials and compares them with established literature values.

This script extracts data from Quantum ESPRESSO output files for different pseudopotentials
and creates horizontal bar charts comparing key properties with literature references.
"""

import re
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator
import pandas as pd
from pathlib import Path

# Define paths - you'll need to adjust these
ANALYSIS_DIR = Path(__file__).parent.absolute()
RESULTS_DIR = ANALYSIS_DIR / "results3"
RESULTS_DIR.mkdir(exist_ok=True)

# Output files to be processed - REPLACE THESE PATHS WITH YOUR ACTUAL FILES
PP1_FM_OUTPUT = ANALYSIS_DIR / "FM3/fe_bulk_fm_paw.out"  # kjpaw_psl.0.2.1
PP1_NM_OUTPUT = ANALYSIS_DIR / "NM3/fe_bulk_nm_paw.out"
PP2_FM_OUTPUT = ANALYSIS_DIR / "FM1/fe_bulk_fm.out"  # kjpaw_psl.1.0.0
PP2_NM_OUTPUT = ANALYSIS_DIR / "NM1/fe_bulk_nm.out"
PP3_FM_OUTPUT = ANALYSIS_DIR / "FM2/fe_bulk_fm_restart.out"  # rrkjus_psl.0.2.1
PP3_NM_OUTPUT = ANALYSIS_DIR / "NM2/fe_bulk_nm_restart.out"

# Pseudopotential labels for plots
PP1_LABEL = "Fe.pbe-spn-kjpaw_psl.0.2.1.UPF"
PP2_LABEL = "Fe.pbe-spn-kjpaw_psl.1.0.0.UPF"
PP3_LABEL = "Fe.pbe-spn-rrkjus_psl.0.2.1.UPF"

# Representative literature values for ferromagnetic BCC iron (α-Fe)
LATTICE_CONST_FM = {
    # Experimental
    "CRC Handbook (2023)\nCompiled experimental data": 2.866,
    
    # Computational
    "Materials Project (MP-13)\n(2023) DFT, GGA-PBE": 2.83,
}

MAGNETIC_MOMENT = {
    # Experimental
    "Crangle & Goodman\n(1950) Magnetization at 0 K": 2.22,
    
    # Computational
    "Herper et al.\n(1999) DFT, GGA": 2.22,
}

ENERGY_DIFF = {
    # Computational
    "Herper et al.\n(1999) DFT, GGA": 0.56,
    "Materials Project (MP-13)\n(2023) DFT, GGA-PBE": 0.53,
}

# Error bars for experimental values
ERROR_BARS = {
    "CRC Handbook\n(2023) Compiled experimental data": 0.001,
    "Crangle & Goodman\n(1950) Magnetization at 0 K": 0.02,
}

def extract_data_from_output(filename):
    """Extract key parameters from a Quantum ESPRESSO output file."""
    data = {}
    
    try:
        with open(filename, 'r') as file:
            content = file.read()
        
        # Extract lattice parameter (from celldm(1))
        match = re.search(r'celldm\(1\)\s*=\s*(\d+\.\d+)', content)
        if match:
            data['celldm1'] = float(match.group(1))
            data['lattice'] = data['celldm1'] * 0.529177  # convert from bohr to Å
        
        # Check if calculation is magnetic
        nspin_match = re.search(r'nspin\s*=\s*2', content)
        mag_results = re.search(r'total magnetization\s*=', content)
        data['is_magnetic'] = bool(nspin_match or mag_results)
        
        # Extract total energy - get the FINAL energy
        energy_entries = re.findall(r'!\s+total energy\s+=\s+([-\d.]+)\s+Ry', content)
        if energy_entries:
            data['total_energy'] = float(energy_entries[-1])  # Use the last entry (final SCF)
            data['total_energy_ev'] = data['total_energy'] * 13.6057  # convert Ry to eV
        
        # Extract pressure
        match = re.search(r'P=\s*([-+]?\d*\.\d+)', content)
        if match:
            data['pressure'] = float(match.group(1))  # in kbar
        
        # Extract magnetic moment - IMPROVED: get the final value after convergence
        mag_entries = re.findall(r'total magnetization\s*=\s*([-+]?\d*\.\d+)', content)
        if mag_entries:
            # Use the last (final) magnetization value
            data['total_magnetization'] = float(mag_entries[-1])  # in Bohr mag/cell
            data['is_magnetic'] = True
            
            # Print diagnostic info
            print(f"File: {Path(filename).name}, Found magnetization: {data['total_magnetization']} μB/cell")
        
        # Extract number of atoms
        match = re.search(r'number of atoms/cell\s*=\s*(\d+)', content)
        if match:
            data['num_atoms'] = int(match.group(1))
            
            # If magnetization is available, calculate magnetic moment per atom
            if 'total_magnetization' in data:
                data['magnetic_moment_per_atom'] = data['total_magnetization'] / data['num_atoms']
                print(f"  Per atom: {data['magnetic_moment_per_atom']} μB/atom (with {data['num_atoms']} atoms)")
        
        return data
    
    except Exception as e:
        print(f"Error extracting data from {filename}: {e}")
        return {}

def plot_pseudopotential_comparison(property_dict, calc_values, title, xlabel, filename, error_dict=None):
    """Create a horizontal bar chart comparing literature values with calculated values for multiple pseudopotentials."""
    # Combine literature and calculated values
    combined_dict = property_dict.copy()
    for label, value in calc_values.items():
        if pd.notna(value):
            combined_dict[label] = value
    
    # Create dataframe for plotting
    df = pd.DataFrame({
        'Source': list(combined_dict.keys()),
        'Value': list(combined_dict.values())
    })
    
    # Add error data if provided
    if error_dict:
        df['Error'] = df['Source'].map(lambda x: error_dict.get(x, 0))
    else:
        df['Error'] = 0
    
    # Sort by values, handling NaN
    df = df.sort_values('Value', na_position='first')
    
    # Set up the figure with publication quality - reduced height for smaller bars
    fig, ax = plt.subplots(figsize=(10, 5), dpi=300)
    
    # Plot horizontal bars
    colors = []
    for source in df['Source']:
        if PP1_LABEL in source or PP2_LABEL in source or PP3_LABEL in source:
            colors.append('#ff7f0e')  # Orange for our calculations
        elif 'CRC Handbook' in source or 'Crangle' in source:  # More specific check for experimental
            colors.append('#2ca02c')  # Green for experimental
        else:
            colors.append('#1f77b4')  # Blue for other computational
    
    # Filter out NaN values for plotting
    plot_df = df.dropna(subset=['Value'])
    plot_colors = [colors[i] for i in plot_df.index]
    
    # Smaller bar height
    bars = ax.barh(plot_df['Source'], plot_df['Value'], color=plot_colors, height=0.6)
    
    # Add error bars where applicable
    for i, (idx, row) in enumerate(plot_df.iterrows()):
        if row['Error'] > 0:
            ax.errorbar(row['Value'], i, xerr=row['Error'], 
                       fmt='none', ecolor='black', capsize=5)
    
    # Customize the plot
    ax.set_xlabel(xlabel, fontsize=14)
    ax.set_title(title, fontsize=16)
    
    # Set reasonable x-limits based on data (not starting from 0)
    values = plot_df['Value'].values
    if len(values) > 0:
        min_val, max_val = np.min(values), np.max(values)
        range_val = max_val - min_val
        ax.set_xlim(min_val - range_val * 0.1, max_val + range_val * 0.25)
    
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.tick_params(axis='both', which='major', labelsize=12)
    ax.tick_params(axis='y', labelsize=11)
    ax.grid(True, linestyle='--', alpha=0.7, which='major')
    ax.grid(True, linestyle=':', alpha=0.4, which='minor')
    
    # Add value annotations with larger font
    for i, bar in enumerate(bars):
        value = plot_df['Value'].iloc[i]
        error = plot_df['Error'].iloc[i]
        if pd.notna(value):
            if 'Lattice' in xlabel:
                text = f'{value:.4f} Å'
            elif 'Magnetic' in xlabel:
                text = f'{value:.2f} μ$_B$'
            else:  # Energy values
                text = f'{value:.3f} eV'
                
            # Calculate text position to avoid error bar overlap
            if error > 0:
                text += f' ±{error:.3f}'
                # Place text after the error bar with additional padding
                text_x = value + error + range_val * 0.02
            else:
                # Place text with standard padding
                text_x = value + range_val * 0.01
                
            ax.text(text_x, bar.get_y() + bar.get_height()/2, 
                   text, va='center', fontsize=12, fontweight='normal')
    
    # Add a legend
    handles = [
        plt.Rectangle((0,0),1,1, color='#2ca02c'),
        plt.Rectangle((0,0),1,1, color='#1f77b4'),
        plt.Rectangle((0,0),1,1, color='#ff7f0e')
    ]
    labels = ['Experimental', 'Computational (Literature)', 'This Study']
    ax.legend(handles, labels, loc='best', fontsize=12)
    
    fig.tight_layout()
    
    # Save the figure
    plt.savefig(RESULTS_DIR / f'{filename}.png', dpi=300, bbox_inches='tight')
    plt.savefig(RESULTS_DIR / f'{filename}.pdf', bbox_inches='tight')
    
    return fig

def plot_combined_comparison(all_data, pp_data):
    """Create a combined figure with all three comparisons as subplots."""
    # Optimized dimensions
    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(15, 5.5), dpi=300)
    
    # Helper function to plot on a specific axis
    def plot_on_axis(ax, property_dict, calc_values, title, xlabel, error_dict=None, subplot_label=''):
        # Combine literature and calculated values
        combined_dict = property_dict.copy()
        for label, value in calc_values.items():
            if pd.notna(value):
                combined_dict[label] = value
        
        # Create dataframe for plotting
        df = pd.DataFrame({
            'Source': list(combined_dict.keys()),
            'Value': list(combined_dict.values())
        })
        
        # Add error data if provided
        if error_dict:
            df['Error'] = df['Source'].map(lambda x: error_dict.get(x, 0))
        else:
            df['Error'] = 0
        
        # Sort by values, handling NaN
        df = df.sort_values('Value', na_position='first')
        
        # Plot horizontal bars with correct colors
        colors = []
        for source in df['Source']:
            if 'This study' in source:
                colors.append('#ff7f0e')  # Orange for our calculations
            elif 'CRC Handbook' in source or 'Crangle' in source:
                colors.append('#2ca02c')  # Green for experimental
            else:
                colors.append('#1f77b4')  # Blue for other computational
        
        # Filter out NaN values for plotting
        plot_df = df.dropna(subset=['Value'])
        plot_colors = [colors[i] for i in plot_df.index]
        
        # Bar height
        bars = ax.barh(plot_df['Source'], plot_df['Value'], color=plot_colors, height=0.7)
        
        # Add error bars where applicable
        for i, (idx, row) in enumerate(plot_df.iterrows()):
            if row['Error'] > 0:
                ax.errorbar(row['Value'], i, xerr=row['Error'], 
                           fmt='none', ecolor='black', capsize=4, linewidth=1.5)
        
        # Customize the plot
        ax.set_xlabel(xlabel, fontsize=14)
        ax.set_title(title, fontsize=15, pad=8)
        
        # Add subplot label
        ax.text(-0.2, 1.08, subplot_label, transform=ax.transAxes, 
                fontsize=16, fontweight='bold', va='top')
        
        # Set reasonable x-limits based on data
        values = plot_df['Value'].values
        if len(values) > 0:
            min_val, max_val = np.min(values), np.max(values)
            range_val = max_val - min_val
            # Tighter limits for better visibility
            if range_val < 0.05:  # For very small ranges
                ax.set_xlim(min_val - range_val * 0.5, max_val + range_val * 0.8)
            else:
                ax.set_xlim(min_val - range_val * 0.1, max_val + range_val * 0.3)
        
        ax.xaxis.set_minor_locator(AutoMinorLocator())
        ax.tick_params(axis='both', which='major', labelsize=12)
        ax.tick_params(axis='y', labelsize=11)
        ax.grid(True, linestyle='--', alpha=0.7, which='major', axis='x')
        ax.grid(True, linestyle=':', alpha=0.4, which='minor', axis='x')
        
        # Add value annotations
        for i, bar in enumerate(bars):
            value = plot_df['Value'].iloc[i]
            error = plot_df['Error'].iloc[i]
            if pd.notna(value):
                if 'Lattice' in xlabel:
                    text = f'{value:.3f}'
                elif 'Magnetic' in xlabel:
                    text = f'{value:.2f}'
                else:  # Energy values
                    text = f'{value:.3f}'
                
                # Calculate text position
                if error > 0:
                    text_x = value + error + range_val * 0.05
                else:
                    text_x = value + range_val * 0.02
                    
                ax.text(text_x, bar.get_y() + bar.get_height()/2, 
                       text, va='center', fontsize=11)
    
    # Use simplified, compact labels to save space
    pp1_fm_label = "PP1: kjpaw v0.2.1\n(This study)"
    pp2_fm_label = "PP2: kjpaw v1.0.0\n(This study)"
    pp3_fm_label = "PP3: rrkjus v0.2.1\n(This study)"
    
    lattice_calcs = {
        pp1_fm_label: pp_data['PP1']['FM'].get('lattice'),
        pp2_fm_label: pp_data['PP2']['FM'].get('lattice'),
        pp3_fm_label: pp_data['PP3']['FM'].get('lattice')
    }
    
    magnetic_calcs = {
        pp1_fm_label: pp_data['PP1']['FM'].get('magnetic_moment_per_atom'),
        pp2_fm_label: pp_data['PP2']['FM'].get('magnetic_moment_per_atom'),
        pp3_fm_label: pp_data['PP3']['FM'].get('magnetic_moment_per_atom')
    }
    
    energy_diff_calcs = {
        pp1_fm_label: pp_data['PP1'].get('energy_diff'),
        pp2_fm_label: pp_data['PP2'].get('energy_diff'),
        pp3_fm_label: pp_data['PP3'].get('energy_diff')
    }
    
    # Plot each subplot
    plot_on_axis(ax1, LATTICE_CONST_FM, lattice_calcs, 
                 "Lattice Constant (FM)", "Lattice Constant (Å)", 
                 ERROR_BARS, '(a)')
    
    plot_on_axis(ax2, MAGNETIC_MOMENT, magnetic_calcs,
                 "Magnetic Moment", "Magnetic Moment (μ$_B$/atom)", 
                 ERROR_BARS, '(b)')
    
    plot_on_axis(ax3, ENERGY_DIFF, energy_diff_calcs,
                 "Energy Difference (FM-NM)", "Energy Difference (eV/atom)", 
                 None, '(c)')
    
    # Add a shared legend at the bottom
    handles = [
        plt.Rectangle((0,0),1,1, color='#2ca02c'),
        plt.Rectangle((0,0),1,1, color='#1f77b4'),
        plt.Rectangle((0,0),1,1, color='#ff7f0e')
    ]
    labels = ['Experimental', 'Computational (Literature)', 'This Study']
    fig.legend(handles, labels, loc='lower center', ncol=3, fontsize=13, 
               bbox_to_anchor=(0.5, -0.07), frameon=False)
    
    # Adjust layout
    plt.tight_layout(rect=[0, 0.06, 1, 0.98])
    
    # Save the combined figure
    plt.savefig(RESULTS_DIR / 'pp_comparison_combined.png', dpi=300, bbox_inches='tight')
    plt.savefig(RESULTS_DIR / 'pp_comparison_combined.pdf', bbox_inches='tight')
    
    print("Combined figure saved as pp_comparison_combined.png/pdf")

def create_summary_table(pp_data):
    """Create and save a summary table of all pseudopotential results."""
    # Prepare data for table
    summary = {
        'Pseudopotential': [PP1_LABEL, PP2_LABEL, PP3_LABEL],
        'FM Lattice Constant (Å)': [
            pp_data['PP1']['FM'].get('lattice', np.nan) if 'PP1' in pp_data else np.nan,
            pp_data['PP2']['FM'].get('lattice', np.nan) if 'PP2' in pp_data else np.nan,
            pp_data['PP3']['FM'].get('lattice', np.nan) if 'PP3' in pp_data else np.nan
        ],
        'NM Lattice Constant (Å)': [
            pp_data['PP1']['NM'].get('lattice', np.nan) if 'PP1' in pp_data else np.nan,
            pp_data['PP2']['NM'].get('lattice', np.nan) if 'PP2' in pp_data else np.nan,
            pp_data['PP3']['NM'].get('lattice', np.nan) if 'PP3' in pp_data else np.nan
        ],
        'Magnetic Moment (μ$_B$)': [
            pp_data['PP1']['FM'].get('magnetic_moment_per_atom', np.nan) if 'PP1' in pp_data else np.nan,
            pp_data['PP2']['FM'].get('magnetic_moment_per_atom', np.nan) if 'PP2' in pp_data else np.nan,
            pp_data['PP3']['FM'].get('magnetic_moment_per_atom', np.nan) if 'PP3' in pp_data else np.nan
        ],
        'FM-NM Energy Difference (eV/atom)': [
            pp_data['PP1']['energy_diff'] if 'PP1' in pp_data and 'energy_diff' in pp_data['PP1'] else np.nan,
            pp_data['PP2']['energy_diff'] if 'PP2' in pp_data and 'energy_diff' in pp_data['PP2'] else np.nan,
            pp_data['PP3']['energy_diff'] if 'PP3' in pp_data and 'energy_diff' in pp_data['PP3'] else np.nan
        ]
    }
    
    # Create DataFrame and save to CSV
    df = pd.DataFrame(summary)
    csv_path = RESULTS_DIR / 'pseudopotential_comparison.csv'
    df.to_csv(csv_path, index=False, float_format='%.6f')
    print(f"Summary table saved to {csv_path}")
    
    return df

def main():
    """Main function to analyze and plot BCC iron data for different pseudopotentials."""
    try:
        # Extract data from output files
        pp_data = {
            'PP1': {
                'FM': extract_data_from_output(PP1_FM_OUTPUT),
                'NM': extract_data_from_output(PP1_NM_OUTPUT)
            },
            'PP2': {
                'FM': extract_data_from_output(PP2_FM_OUTPUT),
                'NM': extract_data_from_output(PP2_NM_OUTPUT)
            },
            'PP3': {
                'FM': extract_data_from_output(PP3_FM_OUTPUT),
                'NM': extract_data_from_output(PP3_NM_OUTPUT)
            }
        }
        
        # Calculate energy differences
        for pp, data in pp_data.items():
            if 'total_energy_ev' in data['FM'] and 'total_energy_ev' in data['NM'] and 'num_atoms' in data['FM']:
                data['energy_diff'] = (data['NM']['total_energy_ev'] - data['FM']['total_energy_ev']) / data['FM']['num_atoms']
        
        # Create labels for plot
        pp1_fm_label = f"{PP1_LABEL}\n(This study)"
        pp2_fm_label = f"{PP2_LABEL}\n(This study)"
        pp3_fm_label = f"{PP3_LABEL}\n(This study)"
        
        # Prepare calculated values for plots
        lattice_calcs = {
            pp1_fm_label: pp_data['PP1']['FM'].get('lattice'),
            pp2_fm_label: pp_data['PP2']['FM'].get('lattice'),
            pp3_fm_label: pp_data['PP3']['FM'].get('lattice')
        }
        
        magnetic_calcs = {
            pp1_fm_label: pp_data['PP1']['FM'].get('magnetic_moment_per_atom'),
            pp2_fm_label: pp_data['PP2']['FM'].get('magnetic_moment_per_atom'),
            pp3_fm_label: pp_data['PP3']['FM'].get('magnetic_moment_per_atom')
        }
        
        energy_diff_calcs = {
            f"{PP1_LABEL}\n(This study)": pp_data['PP1'].get('energy_diff'),
            f"{PP2_LABEL}\n(This study)": pp_data['PP2'].get('energy_diff'),
            f"{PP3_LABEL}\n(This study)": pp_data['PP3'].get('energy_diff')
        }
        
        # Create individual plots
        plot_pseudopotential_comparison(
            LATTICE_CONST_FM, 
            lattice_calcs,
            "BCC Iron Lattice Constant (FM State)", 
            "Lattice Constant (Å)", 
            "pp_comparison_lattice",
            ERROR_BARS
        )
        
        plot_pseudopotential_comparison(
            MAGNETIC_MOMENT, 
            magnetic_calcs,
            "BCC Iron Magnetic Moment", 
            "Magnetic Moment (μ$_B$/atom)", 
            "pp_comparison_magnetic",
            ERROR_BARS
        )
        
        plot_pseudopotential_comparison(
            ENERGY_DIFF, 
            energy_diff_calcs,
            "BCC Iron Energy Difference (FM-NM)", 
            "Energy Difference (eV/atom)", 
            "pp_comparison_energy_diff"
        )
        
        # Create combined figure
        plot_combined_comparison(
            {
                'lattice': LATTICE_CONST_FM,
                'magnetic': MAGNETIC_MOMENT,
                'energy': ENERGY_DIFF
            },
            pp_data
        )
        
        # Create summary table
        summary_df = create_summary_table(pp_data)
        print("\nPseudopotential Comparison Summary:")
        print(summary_df.to_string(index=False, float_format='%.4f'))
        
        print(f"\nPlots and data saved to {RESULTS_DIR}")
        
    except Exception as e:
        import traceback
        print(f"Error analyzing BCC iron data: {e}")
        print(traceback.format_exc())

if __name__ == "__main__":
    main()