# High-throughput-LMTO

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://github.com/balaranjan/High-throughput-LMTO/blob/main/LICENSE)
![Python 3.12](https://img.shields.io/badge/python-3.12-blue.svg)

This is a package to automate the steps involved in performing LMTO calculations. The code will take a list of .cif files as input and for each structure in the list, an LMTO calculation will be performed including optimization, band structure calculation, and density of states calculation (DOS). The outputs from band structure and DOS calculation are saved as comma-separated (.csv) files.

Outputs can be visualized using the following [plotter package](https://github.com/EmilJaffal/High-throughput-LMTO-plotter)

> The current README.md serves as a tutorial and documentation - last update January 16, 2025

## Demo

![HT-demo-gif](assets/HT_DEMO.gif)

## Getting started

The steps to perform an LMTO calculation are as follows:
1. Write the INIT file and run lminit.run to initialize CTRL and CBAK files.
2. Run lmhart.run and get the VOLSPH_by_VOL
3. If VOLSPH_by_VOL less than 100.0, then run lmes.run and lmovl.run iteratively.
       At each iteration check the output of lmovl.run for errors, get recommended
       changes (for RMAXS and OMMAX) based on the messages in the output file, and modify parameters in
       the CTRL file accordingly.
4. Run lmstr.run, check for recommended action in the output file if 
       the run failed. Make the changes (for RMAXS and OMMAX) in CTRL file and run iteratively 
       until lmstr.run finishes without errors.
5. Run lm.run until it converges.
6. Modify the number of KPOINTS (currently doubles) and optimize.
7. Perform band structure calculation.
8. Perform DOS calculation and extract total DOS and partial DOS for each element in the system.

### For Hf, Lanthanides, and Actinides
- All lanthanides are substituted with La and their 4f orbitals are downfolded.
- All actinides are substituted with Th and and their 5f orbitals are downfolded. 
- Hf is replaced with Zr.

## Installation

```bash
$ git clone https://github.com/balaranjan/High-throughput-LMTO
$ cd High-throughput-LMTO
$ pip install -r requirements.txt
$ python main.py [foldername_containing_cif_files]
```

## Contributors

- [Balaranjan Selvaratnam](https://github.com/balaranjan)
- Anton Oliynyk
- [Emil Jaffal](https://github.com/EmilJaffal)

## How to ask for help

- If you have any issues or questions, please feel free to reach out or
  [leave an issue](https://github.com/balaranjan/High-throughput-LMTO/issues).