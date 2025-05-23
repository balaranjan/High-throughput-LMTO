from nlmto import run_lmto, extract_data_from_cif
import multiprocessing as mp
import pandas as pd
import sys
import os

if len(sys.argv) != 2:
    print("Please provide the path to a directory containing .cif files as the only argument.")

# MP
def run_lmto_aux(*args):
    for arg in args:
        run_lmto(**arg)
    
cif_names = [c for c in os.listdir(sys.argv[1]) if c.endswith('.cif')]

total = len(cif_names)
tasks = []

for i, cname in enumerate(cif_names, 0):
        cif = cname.split(os.sep)[-1]

        if cif.endswith(".cif"):
            cif_data = extract_data_from_cif(f"{sys.argv[1]}{os.sep}{cif}")
            cif_data['calc_path'] = sys.argv[1]
            cif_data['cif_path'] = f"{sys.argv[1]}{os.sep}{cif}"
            tasks.append(cif_data)

with mp.Pool(processes=2) as p:
    p.map(func=run_lmto_aux, iterable=tasks)
    
p.close()
p.join()
