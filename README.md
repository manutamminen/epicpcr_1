# epicPCR experiment December 18 2020 

This sequencing run was primarily about testing the performance differences between mineral oil and HFE7500 in emulsion PCR. Also, different spike-in concentrations of Rhodomonas and Ochromonas cells to frozen wastewater were tested.

The following samples were included:

- Rhodomonas and Ochromonas cells; emulsion PCR in HFE7500 + 5% RAN surfactant
- Frozen waterwater + Rhodo and Ochro cells; 10e3 cells / ml; emulsion PCR in HFE7500 + 5% RAN surfactant
- Frozen waterwater + Rhodo and Ochro cells; 10e4 cells / ml; emulsion PCR in HFE7500 + 5% RAN surfactant
- Frozen wastewater; emulsion PCR in HFE7500 + 5% RAN surfactant

- Rhodomonas and Ochromonas cells; emulsion PCR in mineral oil + 4% Abil EM90
- Frozen waterwater + Rhodo and Ochro cells; 10e3 cells / ml; emulsion PCR in mineral oil + 4% Abil EM90
- Frozen waterwater + Rhodo and Ochro cells; 10e4 cells / ml; emulsion PCR in mineral oil + 4% Abil EM90
- Frozen wastewater; emulsion PCR in mineral oil + 4% Abil EM90

Raw sequencing data is available at https://zenodo.org/record/4769327

Summary of the results is available at https://github.com/manutamminen/epicpcr_1/blob/main/docs/index.md

# Building

## Dependencies

- Snakemake
- VSEARCH
- SINA
- FastTree
- Tidyverse

Download the raw data into data/raw.

Start the processing pipeline by invoking `snakemake --cores all report`.



