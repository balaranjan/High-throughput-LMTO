from nlmto import run_lmto, extract_data_from_cif
import multiprocessing as mp
import pandas as pd
import os


skip_cifs = [        
    "262117.cif", "1935906.cif", # VOLSPH_by_VOL=98.7
    "526719.cif", "303534.cif",  #  NiS6V3 lmstr.run Fatal  : MSTRX2: nlsqri=25gt nlsqr=16
    
    #lm.run
    "1719812.cif", # E=-7.56750D+03 DE= 9.75079D+06 after 6 runs
]

# MP
def run_lmto_aux(*args):
    for arg in args:
        run_lmto(**arg)
        
df = pd.read_csv("testfiles.csv")
total = len(df)
cif_names = df['fname'].tolist()
tasks = []

for i, cname in enumerate(cif_names, 0):
        cif = cname.split(os.sep)[-1]

        if cif.endswith(".cif"):
            if cif in skip_cifs:
                continue
            cif_data = extract_data_from_cif(f"/home/bala/research/1_LMTO/lmto_script/not_prototype_CIFs/{cif}")
            tasks.append(cif_data)

print(len(tasks))
with mp.Pool(processes=32) as p:
    p.map(func=run_lmto_aux, iterable=tasks)
    
p.close()
p.join()