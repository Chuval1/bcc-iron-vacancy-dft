# bcc-iron-vacancy-dft

## Overview

This repository contains scripts and data for density functional theory (DFT) studies of vacancies in body-centered cubic (BCC) iron. The project focuses on analyzing vacancy formation energies, structural relaxations, and magnetic properties using Quantum ESPRESSO calculations. Both PAW and USPP pseudopotentials are compared, and results are benchmarked against literature and experimental values.

## Main Features
- **Vacancy Analysis**: Calculates vacancy formation energies, analyzes atomic relaxations, and magnetic moments for BCC iron supercells.
- **Pseudopotential Comparison**: Compares results from different pseudopotentials and visualizes key properties (lattice constant, magnetic moment, energy differences) against literature.
- **Convergence Testing**: Provides scripts for k-point and supercell convergence analysis to ensure robust DFT results.

## Directory Structure
- `vacancy/` — Scripts and data for vacancy calculations and analysis.
  - `vacancy_analysis_FM.py`: Main script for vacancy formation energy and magnetic analysis.
  - Subfolders contain input/output files and results for various supercell and pseudopotential setups.
- `primitive/` — Scripts and data for bulk BCC iron calculations and pseudopotential comparisons.
  - `pseudopotential_comparison2.py`: Script for comparing different pseudopotentials.
- `convergence_tests/` — Scripts and data for k-point and supercell convergence studies.
  - `analyze_convergence.py`: Script for analyzing convergence test outputs.
- `LICENSE` — MIT License.
- `README.md` — This file.

## Getting Started

1. **Requirements**:
   - Python 3.x
   - numpy, pandas, matplotlib
2. **Usage**:
   - Run analysis scripts from their respective directories after generating Quantum ESPRESSO output files.
   - Example: `python vacancy/vacancy_analysis_FM.py`

## Acknowledgements

This project benefited from the use of AI-assisted coding and debugging tools, including:
- [GitHub Copilot](https://github.com/features/copilot)
- [ChatGPT](https://chat.openai.com/)
- [Claude](https://claude.ai/)

These tools were used for code suggestions, debugging, and documentation improvements.

## License

This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for details.