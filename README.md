# APOBEC3-Mpox-IRs

## Folders and Files

### `/input/`
Contains reference and supporting input data:
- `mpox_genome_seq.fna` – MPXV reference genome sequence  
- `palindrom_analyzer_output.txt` – Output from PalindromeAnalyzer (predicted inverted repeats / hairpins)  
- `snp_locations.xlsx` – List of single-nucleotide mutations associated with the 2022 Mpox outbreak

### `/binomial/`
This folder contains scripts for estimating statistical significance using an approximation by the binomial distribution.

- `config.yaml` – Configuration file specifying:
  - The structure type to analyze (hairpin, loop, TC, or C at the 3′ end of the loop)  
  - The method for selecting overlapping structures (e.g., single most stable, all most stable non-overlapping, maximum or minimum coverage)

Usage:
To perform the estimation, modify the desired settings in `config.yaml`, then run:
```bash
python main.py
```

### `simulations.ipynb`  
This is a Python Jupyter notebook containing Monte Carlo simulations.  
All relevant comments and explanations are provided directly within the notebook.
