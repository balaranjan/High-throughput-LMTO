# High-throughput-LMTO

This is a package to automate the steps involved in performing LMTO calculations.  The code will take a list of .cif files as input and for each structure in the list, an LMTO calculation will be performed including optimization, band structure calculation, and density of states calculation (DOS).  The outputs from band structure and DOS calculation are saved as comma-separated files.

The steps to perform an LMTO calculation.
1. Write the INIT file and run lminit.run to initialize CTRL and CBAK files.
2. Run lmhart.run and get the VOLSPH_by_VOL
3. If VOLSPH_by_VOL less than 100.0, then run lmes.run and lmovl.run iteratively.
       At each iteration check the output of lmovl.run for errors, get recommended
       changes based on the messages in the output file, and modify parameters in
       the CTRL file accordingly.
4. Run lmstr.run, check for recommended action in the output file if 
       the run failed. Make the changes in CTRL file and run iteratively 
       until lmstr.run finishes without errors.
5. Run lm.run until it converges.
6. Modify the number of KPOINTS (currently doubles) and optimize.
7. Perform band structure calculation.
8. Perform DOS calculation and extract total DOS and partial DOS for each element in the system.

This an ongoing development, expect frequent changes.
