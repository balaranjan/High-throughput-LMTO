import os
import sys
import time
import warnings
import pexpect
import subprocess
import pandas as pd
from typing import Any
from collections import defaultdict
from utils import *
import shutil


def print_to_console(func):
    """
    Print current step and its staus at the end to console.
    """
    def wrapper(**kwargs):
        func_name = func.__name__[4:]
        if func_name == "lm":
            print("\tRunning optimization ", end=" ")
        else:
            print(f"\tRunning {func_name}", end=" ")
        res = func(**kwargs)
        
        if func_name == "lm":
            print("error" if res[0] else " ok")
        else:
            print("error" if res[0] else "ok")
        return res
    return wrapper


def write_INIT_file(**kwargs):
    
    """
    Write INIT file for LMTO calculation using cif data obtained from utils.extract_data_from_cif().
    """
    
    atomic_numbers = {'H': 1, 'He': 2, 'Li': 3, 'Be': 4, 'B': 5, 'C': 6, 'N': 7, 'O': 8, 
                      'F': 9, 'Ne': 10, 'Na': 11, 'Mg': 12, 'Al': 13, 'Si': 14, 'P': 15, 
                      'S': 16, 'Cl': 17, 'Ar': 18, 'K': 19, 'Ca': 20, 'Sc': 21, 'Ti': 22, 
                      'V': 23, 'Cr': 24, 'Mn': 25, 'Fe': 26, 'Co': 27, 'Ni': 28, 'Cu': 29, 
                      'Zn': 30, 'Ga': 31, 'Ge': 32, 'As': 33, 'Se': 34, 'Br': 35, 'Kr': 36, 
                      'Rb': 37, 'Sr': 38, 'Y': 39, 'Zr': 40, 'Nb': 41, 'Mo': 42, 'Tc': 43, 
                      'Ru': 44, 'Rh': 45, 'Pd': 46, 'Ag': 47, 'Cd': 48, 'In': 49, 'Sn': 50, 
                      'Sb': 51, 'Te': 52, 'I': 53, 'Xe': 54, 'Cs': 55, 'Ba': 56, 'La': 57, 
                      'Ce': 58, 'Pr': 59, 'Nd': 60, 'Pm': 61, 'Sm': 62, 'Eu': 63, 'Gd': 64, 
                      'Tb': 65, 'Dy': 66, 'Ho': 67, 'Er': 68, 'Tm': 69, 'Yb': 70, 'Lu': 71, 
                      'Hf': 72, 'Ta': 73, 'W': 74, 'Re': 75, 'Os': 76, 'Ir': 77, 'Pt': 78, 
                      'Au': 79, 'Hg': 80, 'Tl': 81, 'Pb': 82, 'Bi': 83, 'Po': 84, 'At': 85, 
                      'Rn': 86, 'Fr': 87, 'Ra': 88, 'Ac': 89, 'Th': 90, 'Pa': 91, 'U': 92, 
                      'Np': 93, 'Pu': 94, 'Am': 95, 'Cm': 96, 'Bk': 97, 'Cf': 98, 'Es': 99, 
                      'Fm': 100, 'Md': 101, 'No': 102, 'Lr': 103, 'Rf': 104, 'Db': 105, 
                      'Sg': 106, 'Bh': 107, 'Hs': 108, 'Mt': 109, 'Ds': 110, 'Rg': 111, 
                      'Cn': 112, 'Nh': 113, 'Fl': 114, 'Mc': 115, 'Lv': 116, 'Ts': 117, 
                      'Og': 118}
    
    space_group_nums_params = {}
    
    # triclinic
    for i in range(1, 3):
        space_group_nums_params[i] = "A-B-C_a-b-g"
        
    # monoclinic
    for i in range(3, 16):
        space_group_nums_params[i] = "A-B-C_b"
        
    # orthorhombic
    for i in range(16, 75):
        space_group_nums_params[i] = "A-B-C"
        
    # tetragonal
    for i in range(75, 143):
        space_group_nums_params[i] = "A-C"
        
    # trigonal
    for i in range(143, 168):
        space_group_nums_params[i] = "A-C"
        
    # hexagonal
    for i in range(168, 195):
        space_group_nums_params[i] = "A-C"
        
    # cubic
    for i in range(195, 231):
        space_group_nums_params[i] = "A"
    
    param_map = {"A": "_cell_length_a", "B": "_cell_length_b", "C": "_cell_length_c",
                 "a": "_cell_angle_alpha", "b": "_cell_angle_beta", "g": "_cell_angle_gamma"}
    params = []
    sg_params = space_group_nums_params[int(kwargs['_space_group_IT_number'])]

    if '_' in sg_params:
        params, angles = sg_params.split('_')
        params = params.split('-') if '-' in params else [params]
        angles = angles.split('-') if '-' in angles else [angles]
    else:
        params = sg_params.split('-') if '-' in sg_params else [sg_params]
        angles = []
    
    params = [param_map[p] for p in params]
    angles = [param_map[a] for a in angles]
    
    lines = []
    first_line = f"          SPCGRP={kwargs['_space_group_IT_number']} IORIGIN={kwargs.get('origin', 1)} ATUNITS={kwargs.get('aomic_units', 'F')}"
    for k, param in [["A", "_cell_length_a"], ["B", "_cell_length_b"], ["C", "_cell_length_c"]]:
        if param in kwargs and param in params:
            first_line += f" {k}={kwargs[param]}"
            
    for k, param in ["ALPHA", "_cell_angle_alpha"], ["BETA", "_cell_angle_beta"], ["GAMMA", "_cell_angle_gamma"]:
        if param in kwargs and param in angles:
            first_line += f" {k}={kwargs[param]}"
            
    first_line += '\n'
    lines.append(first_line)
    
    for site in kwargs['atom_site_data']:
        e = site[0]
        for c in e:
            if c.isnumeric():
                e = e.replace(c, "")
    
        lines.append(f"          ATOM={site[0]:<3}  Z={atomic_numbers[e]:<3}   X={site[1]}\n")
            
    with open("INIT", 'w') as f:
        f.writelines(lines)
        
    return True

@print_to_console
def run_lminit(**kwargs):
    
    write_INIT_file(**kwargs)
    lminit_output_filename = 'output_lminit.txt'
    output = open(lminit_output_filename, 'a')
    subprocess.run('lminit.run', stdout=output, stderr=output)
    
    error = aborted(lminit_output_filename)

    return [error]


@print_to_console
def run_lmhart():
    
    lmhart_output_filename = 'output_lmhart.txt'
    output = open(lmhart_output_filename, 'a')
    subprocess.run('lmhart.run', stdout=output, stderr=output)
    error = aborted(lmhart_output_filename)
    
    if error:
        return [error, 0.0]
    
    VOLSPH_by_VOL = [line for line in open(lmhart_output_filename, 'r').readlines() if 'VOLSPH/VOL' in line][-1]
    try:
        VOLSPH_by_VOL = float(VOLSPH_by_VOL.split('=')[-1].replace('%', ''))
    except Exception as e:
        print("Error reading VOLSPH/VOL!")
        print(e)
    
    return [error, VOLSPH_by_VOL]


@print_to_console
def run_lmovl(iteration):
    
    lmvol_output_filename = f'output_lmovl_{iteration}.txt'
    output = open(lmvol_output_filename, 'a')
    subprocess.run('lmovl.run', stdout=output, stderr=output)
    VOLSPH_by_VOL = [line for line in open(lmvol_output_filename, 'r').readlines() if 'VOLSPH/VOL' in line][-1]
    overlap = [line for line in open(lmvol_output_filename, 'r').readlines() if 'Overlap' in line][-1]
    
    issues = find_issues(lmvol_output_filename)
    
    try:
        VOLSPH_by_VOL = float(VOLSPH_by_VOL.split('=')[-1].replace('%', ''))
    except Exception as e:
        VOLSPH_by_VOL = None
        print("Error reading VOLSPH/VOL!")
        print(e)
        
    try:
        overlap = float(overlap.split('=')[-1].replace('%', ''))
    except Exception as e:
        overlap = None
        print("Error reading overlap volume!")
        print(e)
        
    error = aborted(lmvol_output_filename)
        
    return [error, VOLSPH_by_VOL, overlap, issues]


@print_to_console
def run_lmes(iteration):
    
    lmes_output_filename = f'output_lmes_{iteration}.txt'
    output = open(lmes_output_filename, 'a')
    
    subprocess.run('lmes.run', stdout=output, stderr=output)
    error = aborted(lmes_output_filename)
    
    if error:
        VOLSPH_by_VOL = 0.
    else:
        VOLSPH_by_VOL = [line for line in open(lmes_output_filename, 'r').readlines() if 'VOLSPH/VOL' in line][-1]
        overlap = [line for line in open(lmes_output_filename, 'r').readlines() if 'Overlap' in line]
        if len(overlap):
            overlap = overlap[-1]
        
        try:
            VOLSPH_by_VOL = float(VOLSPH_by_VOL.split('=')[-1].replace('%', ''))
        except Exception as e:
            VOLSPH_by_VOL = None
            print("Error reading VOLSPH/VOL!")
            print(e)
        
    return [error, VOLSPH_by_VOL]


@print_to_console
def run_lmctl():
    
    lmctl_output_filename = 'output_lmctl.txt'
    output = open(lmctl_output_filename, 'a')
    subprocess.run('lmctl.run', stdout=output, stderr=output)
    
    error = aborted(lmctl_output_filename)
    return [error]


@print_to_console
def run_lmstr(iteration):
    lmstr_output_filename = f'output_lmstr_{iteration}.txt'
    output = open(lmstr_output_filename, 'w')
    subprocess.run("lmstr.run", stdout=output, stderr=output)
    
    error = aborted(lmstr_output_filename)
    if error:
        issues = find_issues(lmstr_output_filename)
    else:
        issues = {}
    return [error, issues]


@print_to_console
def run_lm(calc_type, num_atoms, n_try_max=5):
    
    converged = False
    n_try = 1
    etot_and_time = []
    while not converged and n_try < n_try_max:
        lm_output_filename = f'output_lm_{calc_type}_{n_try}.txt'
        t0 = time.time() 
        print(f"iteration {n_try}", end=" ")
        output = open(lm_output_filename, 'w')

        subprocess.run('lm.run', stdout=output, stderr=output)
        t = f"{time.time() - t0:.1f}"
        
        etot = extract_scf_data(lm_output_filename, n_try, natoms=num_atoms)
        with open(lm_output_filename, 'r') as f:
            converged = "Jolly good show!" in f.read()
            
        etot_and_time.append({
            'Type': 'first optimization', 'Iteration': n_try,
            'converged': converged, 'ETOT (eV)': etot, 'Eime (s)': t
        })
        n_try += 1
        
    return [not converged, etot_and_time]


@print_to_console
def run_lmdos():
    modify_CTRL_file(set_DOS_EMIN=-1.1,
                     set_DOS_EMAX=1.1,
                     set_DOS_NOPTS=1800)
    
    lmdos_output_filename = 'output_lmdos.txt'
    output = open(lmdos_output_filename, 'a')
    subprocess.call('lmdos.run', stdout=output)
    
    error = aborted(lmdos_output_filename)
    return [error]

@print_to_console
def run_lmbnd():
    lmbnd_output_filename = 'output_lmcbnd.txt'
    output = open(lmbnd_output_filename, 'a')
    subprocess.call('lmbnd.run', stdout=output)
    
    error = aborted(lmbnd_output_filename)
    return [error]


def run_lmto(**kwargs):
    os.chdir(kwargs['calc_path'])
    """
    Run the steps to perform an LMTO calculation.
     1. Write INIT file and run lminit.run to initialize CTRL and CBAK files.
     2. Run lmhart.run and get the VOLSPH_by_VOL
     3. If VOLSPH_by_VOL less than 100.0, then run lmes.run and lmovl.run iteratively.
        At each iteration check the output of lmovl.run for errors, get recommended
        changes based on the messages in the output file, and modify parameters in
        the CTRL file accordingly.
     4. Run lmstr.run, check for recommeded action in the output file if 
        the run failed. Make the changes in CTRL file and run iteratively 
        until lmstr.run finishes without errors.
     5. Run lm.run until it converges.
    """
    
    print(f"\tName: {kwargs['name']}")
    
    # create dir
    name = f"{kwargs['name']}"

    if name is not None:
        try:
            os.mkdir(name)
            os.chdir(name)
        except FileExistsError:
            os.chdir(name)
            for f in os.listdir('.'):
                os.remove(f)

    error_init = run_lminit(**kwargs)[0]
    if error_init:
        print(f"{kwargs['name']} failed")
        return True
    
    error_hart, VOLSPH_by_VOL = run_lmhart()
    if error_hart:
        print(f"{kwargs['name']} failed")
        return True
    
    modify_CTRL_file(set_IO_VERBOS=50)
    error_ctl = run_lmctl()
    error_ovl, VOLSPH_by_VOL, overlap, issues = run_lmovl(iteration=1)
    
    if VOLSPH_by_VOL >= 100:
        error_es = False
    else:
        # set VERBOSE to 50
        n_try = 1
        print(f"\t\tVOLSPH_by_VOL={VOLSPH_by_VOL}, overlap={overlap}")
        while VOLSPH_by_VOL is not None and VOLSPH_by_VOL < 100.0 and n_try < 10:
            print(f"\tRunning lmes.run and lmovl.run, iteration {n_try}")
            error_ctl = run_lmctl()
            if len(issues):
                # print(f"\tChanging params based on recomendations {issues}")
                modify_CTRL_file(**issues)
            
            error_es, VOLSPH_by_VOL = run_lmes(iteration=n_try)
            if not error_es:
                error_ovl, VOLSPH_by_VOL, overlap, issues = run_lmovl(iteration=1+n_try)
            
            if error_es or error_ovl:
                break
            
            # print(f"\t\tVOLSPH_by_VOL={VOLSPH_by_VOL}, overlap={overlap}")
            n_try += 1
        if VOLSPH_by_VOL < 100.0:
            error_ovl = True
            
    if error_ovl:
        print(f"{kwargs['name']} failed")
        return True
            
    error_str, issues = run_lmstr(iteration=1)
    n_try = 1
    if error_str and len(issues):
        while error_str and n_try <= 10:
            modify_CTRL_file(**issues)
            # print(f"\tRunning lmstr, iteration: {1+n_try}", end=" ")
            error_str, issues = run_lmstr(iteration=1+n_try)
            n_try += 1
            
    if error_str:
        print(f"{kwargs['name']} failed")
        return True
    
    etot_and_time = []
    not_converged, etot_t = run_lm(calc_type='first', num_atoms=kwargs['num_atoms'],
                                          n_try_max=5)
    etot_and_time.extend(etot_t)
    if not_converged:
        print(f"{kwargs['name']} failed")
        return True
    
    # 2x K-point, BEGMOM=F
    modify_CTRL_file(set_START_BEGMOM='T', modify_BZ_NKABC=2)
    not_converged, etot_t = run_lm(calc_type='2xKPTS', num_atoms=kwargs['num_atoms'],
                                   n_try_max=5)
    etot_and_time.extend(etot_t)
    pd.DataFrame(etot_and_time).to_csv("etot_time.csv", index=False)
    if not_converged:
        print(f"{kwargs['name']} failed")
        return True
    
    # DOS
    error_dos = run_lmdos()[0]
    if error_dos:
        print(f"{kwargs['name']} failed")
        return True
    elem_classes = process_dos_data(elements='all', name='DOS')  # TDOS
    for k, v in elem_classes.items():
        process_dos_data(elements=v, name=f"DOS-{k}")
        
    # Band structure
    error_bnd = run_lmbnd()[0]
    if error_bnd:
        print(f"{kwargs['name']} failed")
        return True
    get_band_structure(kwargs['name'].split('-')[0])
    
    if not (error_init or error_hart or error_ovl or error_es or error_str or
            not_converged or error_dos or error_bnd):
        pass
        # shutil.rmtree(os.getcwd())
        # shutil.move(os.getcwd(), "/home/bala/research/1_LMTO/lmto_script/lmto_tests/noi/")
    else:
        print(f"{kwargs['name']} failed")
        print(dict(zip(['init', 'hart', 'ovl', 'es', 'str', 'lm', 'dos', 'band'], 
                       [error_init , error_hart , error_ovl , error_es , error_str ,
            not_converged , error_dos, error_bnd])))
    return error_init or error_hart or error_ovl or error_es or error_str or not_converged or error_dos or error_bnd


def extract_scf_data(lm_output_filename, n_opt, natoms):
    # write iter, etot, dE as csv; return last etot
    with open(lm_output_filename, 'r') as f:
        lines = f.readlines()
        etot = [l for l in lines if ' ETOT' in l]

        iter_etot = [l for l in lines if 'ITER' in l and 'OUT OF' in l]
        iter_etot = [[int(l.split('OUT OF')[0][5:]), float(l.split("ETOT=")[-1][:-1])] for l in iter_etot]
        iter_etot = pd.DataFrame({'n_opt': [n_opt for _ in range(len(iter_etot))],
                                  'type': ["initial" for _ in range(len(iter_etot))], 
                                  'iteration':[i[0] for i in iter_etot], 
                                  'Energy (mRy)':[i[1] for i in iter_etot]})
        
        iter_etot['Energy (eV)'] = iter_etot['Energy (mRy)'] * 13.6057039763 / 1000.0
        iter_etot['Energy (eV/atom)'] = iter_etot['Energy (eV)'] / natoms
        
        iter_etot.to_csv(f"ETOT_values_{lm_output_filename.split('_')[2]}.csv", 
                         index=False)
        
        return float(iter_etot.iloc[-1]['Energy (eV)'])
    

if __name__ == "__main__":
    import multiprocessing as mp
    wd = "/home/bala/research/1_LMTO/lmto_script"    

    
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
                cif_data = extract_data_from_cif(f"{wd}/not_prototype_CIFs/{cif}")
                tasks.append(cif_data)
    
    print(len(tasks))
    with mp.Pool(processes=32) as p:
        p.map(func=run_lmto_aux, iterable=tasks)
        
    p.close()
    p.join()