import os
import re
import pexpect
import pandas as pd
import numpy as np
from typing import Any
from collections import defaultdict


def extract_params_from_line(line):
    """
    Extract all parameters and their values from a line from CTRL file. 
    """
    line = line.strip()
    _value = line.split()
    value = []
    if "ATOM" in line and "SIGMA" in line:
        return [f"STRATOMLINE=={line}"]
    for v in _value:
        if '=' in v:
            value.append(v)
        else:
            value[-1] += f" {v}"
    return value

def read_ctrl():
    
    """
    Read the CTRL file and return as dictionary.
    """
    ctrl = {}
    atom_counter = 1
    with open("CTRL", 'r') as f:
        lines = f.readlines()
        i = 0
        while i < len(lines):
            line = lines[i]
            if not line.startswith(" "*10):
                # new block
                block_name = line[:10].strip()

                # read block
                if block_name in ['HEADER', 'VERS']:
                    block_value = line[10:]
                    
                elif block_name in ['IO', 'OPTIONS', 'STR', 'BZ', 'DOS', 'SCALE']:
                    
                    # dictionaries maintain their insertion order as of version 3.7
                    value = extract_params_from_line(line[10:])
                    block_value = {}
                    for v in value:
                        v = v.split('=')
                        block_value[v[0]] = v[1]
                    
                    in_block = lines[i+1].startswith(" "*10)
                    j = i + 1
                    
                    while in_block:
                        value = extract_params_from_line(lines[j])
                        for v in value:
                           
                            if "STRATOMLINE" in v:
                                v = v.split('==')
                                block_value[f"{v[0]}{atom_counter:02}"] = v[1]
                                atom_counter += 1
                            else:
                                v = v.split('=')
                                block_value[v[0]] = v[1]
                        in_block = lines[j+1].startswith(" "*10) if j+1 < len(lines) else False
                        i = j
                        j += 1 
                    
                else:
                    in_block = len(lines) > i+1 and lines[i+1].startswith(" "*10)
                    j = i + 1
                    block_value = [line[10:]]
                    while in_block:
                        block_value.append(lines[j][10:])
                        in_block = len(lines) > j+1 and lines[j+1].startswith(" "*10) if j+1 < len(lines) else False
                        # lines[j+1].startswith(" "*10) if j+1 < len(lines) else False
                        i = j
                        j += 1  
                ctrl[block_name] = block_value 
                
            i += 1
    return ctrl


def write_ctrl(ctrl, **kwargs):
    
    """
    Write CTRL file.  Modifications for certain blocks in the CTRL file can be supplied as 
    keyword arguments.
    
    The keyword will have the following format:
    {modify/set}_BLOCK_PARAMETER=VALUE
    
    When the modify prefix is used, the existing value(s) will be incresed by VALUE.
    For set prefix, the existing value(s) will be replaced instead.
    """
    
    # alter params from kwargs
    for k, v in kwargs.items():
        change_type, block, param = k.split("_")
        
        # numeric can be modified, string will be set
        if change_type == "modify":
            # KPTS
            current_val = ctrl[block][param]
            if block == "SCALE":
                new_val = [float(c)+v for c in current_val.split()]
                new_val = " ".join([f"{c:.2f}"[1:] for c in new_val])
            elif block == "BZ" and param == "NKABC":
                new_val = [str(int(i)*v) for i in current_val.strip().split()]
                new_val = ' '.join(new_val)
            else:
                if current_val.strip().isdigit():
                    new_val = int(current_val) + v
                else:
                    nround = len(str(v).split('.')[-1])
                    new_val = round(float(current_val) + v, nround)
                    
            ctrl[block][param] = new_val
            # print(f"Modifying {param} in {block} block from {current_val} to {new_val}.")
        
        elif change_type == "set":
            if block == "CLASS":
                current_val = ctrl[block]
                v = float(v)
                for i, atom in enumerate(current_val):
                    if "ATOM=" not in atom:
                        continue
                    atom, radius = atom.split("R=")
                    radius = radius.split()[0].strip()
                    atom = atom.split()[0].replace("ATOM=", "").strip()
                    # print(atom, radius)
                    
                    element = atom
                    for e in element:
                        if e.isnumeric():
                            element = element.replace(e, "")
                    if element == param and float(radius) < v:
                        current_val[i] = current_val[i].replace(f"R={radius}", f"R={v}")
                    
                    # print(i, atom, radius)
                ctrl[block] = current_val
                    
            elif block == "START" and param == "BEGMOM":
                current_val = "F"
                ctrl[block][1] = ctrl[block][1].replace("BEGMOM=F", f"BEGMOM={v}")
            elif block == "COHP":
                ctrl[block] = v # "\n".join(v)
            else:
                current_val = ctrl[block][param]
                ctrl[block][param] = f"{v}"
            # print(f"Setting {param} in {block} block from {current_val} to {v}.")
    
    max_len = len("CHARGE    LMTODAT=T ELF=F ADDCOR=F SPINDENS=F CHARWIN=F EMIN=-2 EMAX=2")

    with open("CTRL", "w") as f:
        for k, v in ctrl.items():
            if k == "COHP":
                print("cohp", type(v), v)
            if isinstance(v, dict):
                block = ""
                line = f"{k:<10}"
                for bk, bv in v.items():
                    if bk[:-2] == "STRATOMLINE":
                        if line != "":
                            block += f"{line}\n"
                        line = ""
                        block += f"{' '*10}{bv}\n"
                    else:
                        kv = f"{bk}={bv} "
                        if len(line + kv) <= max_len:
                            line += kv
                        else:
                            block += f"{line}\n"
                            line = " "*10
                            line += kv
                if not block.endswith(line):
                    block += f"{line}\n"
   
                f.write(block)
            elif isinstance(v, str):
                print(f"{k:<10}{v}")
                f.write(f"{k:<10}{v}")
            else:
                f.write(f"{k:<10}{v[0]}")
                for _v in v[1:]:
                    f.write(f"          {_v.replace("\n", "")}\n")
                    
                    
def modify_CTRL_file(**kwargs):
    ctrl = read_ctrl()
    write_ctrl(ctrl, **kwargs)
    

def _parse_formula(formula: str, strict: bool = True) -> dict[str, float]:
    """
    copied from pymatgen
    
    Args:
        formula (str): A string formula, e.g. Fe2O3, Li3Fe2(PO4)3.
        strict (bool): Whether to throw an error if formula string is invalid (e.g. empty).
            Defaults to True.

    Returns:
        Composition with that formula.

    Notes:
        In the case of Metallofullerene formula (e.g. Y3N@C80),
        the @ mark will be dropped and passed to parser.
    """
    # Raise error if formula contains special characters or only spaces and/or numbers
    if "'" in formula:
        formula = formula.replace("'", "")

    if strict and re.match(r"[\s\d.*/]*$", formula):
        print(formula)
        raise ValueError(f"Invalid {formula=}")

    # For Metallofullerene like "Y3N@C80"
    formula = formula.replace("@", "")
    # Square brackets are used in formulas to denote coordination complexes (gh-3583)
    formula = formula.replace("[", "(")
    formula = formula.replace("]", ")")
    
    def get_sym_dict(form: str, factor: float) -> dict[str, float]:
        sym_dict: dict[str, float] = defaultdict(float)
        for match in re.finditer(r"([A-Z][a-z]*)\s*([-*\.e\d]*)", form):
            el = match[1]
            amt = 1.0
            if match[2].strip() != "":
                amt = float(match[2])
            sym_dict[el] += amt * factor
            form = form.replace(match.group(), "", 1)
        if form.strip():
            raise ValueError(f"{form} is an invalid formula!")
        return sym_dict

    match = re.search(r"\(([^\(\)]+)\)\s*([\.e\d]*)", formula)
    while match:
        factor = 1.0
        if match[2] != "":
            factor = float(match[2])
        unit_sym_dict = get_sym_dict(match[1], factor)
        expanded_sym = "".join(f"{el}{amt}" for el, amt in unit_sym_dict.items())
        expanded_formula = formula.replace(match.group(), expanded_sym, 1)
        formula = expanded_formula
        match = re.search(r"\(([^\(\)]+)\)\s*([\.e\d]*)", formula)
    return get_sym_dict(formula, 1)
    
    
def extract_data_from_cif(path):
    
    """
    Parse CIF and extract the information needed for LMTO calculation.
    """

    if not os.path.isfile(path):
        print(f"File {path} not found!")
        return
    
    attributes_to_read = [
        '_chemical_formula_sum',
        '_cell_length_a',
        '_cell_length_b',
        '_cell_length_c',
        '_cell_angle_alpha',
        '_cell_angle_beta',
        '_cell_angle_gamma',
        '_space_group_IT_number',
        '_space_group_name_H-M_alt',
        '#_database_code_PCD',
        '_cell_formula_units_Z'
    ]
    
    data = defaultdict(Any)
    with open(path, 'r') as f:
        lines = f.readlines()
        ln = 0
        while ln < len(lines):
            line = lines[ln].lstrip()
            if "_symmetry_space_group_name_H-M" in line:
                line = line.replace("_symmetry_space_group_name_H-M", "_space_group_name_H-M_alt")
            elif "_symmetry_Int_Tables_number" in line:
                line = line.replace("_symmetry_Int_Tables_number", "_space_group_IT_number")
            if not len(line.split()):
                ln += 1
                continue
            
            if line.split()[0] in attributes_to_read:
                next_line = lines[ln+1].lstrip()
                if next_line.startswith('_') or next_line.startswith('loop_'):
                    data[line.split()[0]] = line.split()[1:]
                else:
                    line_data = ''
                    while not next_line.startswith('_'):
                        line_data += next_line.strip().replace(';', '').replace(' ', '')
                        ln += 1
                        next_line = lines[ln+1].lstrip()

                    data[line.split()[0]] = line_data.strip().replace(';', '').replace(' ', '')
                
            if line.startswith('loop_') and lines[ln+1].lstrip().startswith('_atom_site') \
                and "aniso" not in lines[ln+1].lstrip():
                site_data = []
                keys = []
                ln += 1
                while lines[ln].lstrip().startswith('_atom_site'):
                    keys.append(lines[ln].lstrip().replace('\n', ''))
                    ln += 1

                while ln < len(lines) and not lines[ln].lstrip().startswith('_') and not lines[ln].lstrip().startswith('#'):
                    if len(lines[ln].strip()):
                        site_data.append(lines[ln].replace('\n', "").lstrip())
                    ln += 1
                data['atom_site_data'] = site_data
                
            ln += 1
    # print(data)
    data = dict(data)
    for k in ["#_database_code_PCD", "_cell_length_a", "_cell_length_b", "_cell_length_c", 
              "_cell_angle_alpha", "_cell_angle_beta", "_cell_angle_gamma",
              "_space_group_IT_number"]:
        if k in data:
            data[k] = get_float(data[k][0])
        else:
            if k == "#_database_code_PCD":
                data[k] = 000
    data["_chemical_formula_sum"] = "".join(data["_chemical_formula_sum"]).replace("'", "")
    data["atom_site_data"] = format_sites(data["atom_site_data"])
    
    if "choice" in data["_space_group_name_H-M_alt"]:
        origin = data["_space_group_name_H-M_alt"][-1]
        for c in origin:
            if not c.isnumeric():
                origin = origin.replace(c, "")
        data['origin'] = origin
        
    data['name'] = f"{data['_chemical_formula_sum']}-{data['#_database_code_PCD']}"
    data['num_atoms'] = int(data['_cell_formula_units_Z'][0]) * \
        sum(list(_parse_formula(data['_chemical_formula_sum']).values()))
    # print(data)
    return dict(data)
    
    

def format_sites(sites: str) -> list:
    """
    Receives site data from CIF and returns reformatted labels and coordinates for LMTO.
    
    Change Ce - Lu to La
    
    """
    sites = [s.split() for s in sites]
    
    labels = defaultdict(int)
    formatted_sites = []
    
    change_map = {"D": "H",
                  "Tl": "In",
                  "Hf": "Zr"}
    
    for site in sites:
        if "loop" in site[0].strip():
            continue
        # print(site)
        for ch in ["+", "-"]:
            if ch in site[1]:
                site[1] = site[1].replace(ch, "")

        sform = _parse_formula(site[1])
        site[1] = list(sform.keys())[0]
        site[1] = change_map.get(site[1], site[1])
        
        labels[site[1]] += 1
        if "Uani" in site or "Uiso" in site:
            # Ge3 Ge 4 k 0.06789 0.19306 0.5 1
            formatted_sites.append([f"{site[1]}{labels[site[1]]}", ' '.join([str(float(get_float(p))) for p in site[2:5]])])
        else:
            formatted_sites.append([f"{site[1]}{labels[site[1]]}", ' '.join([str(float(get_float(p))) for p in site[4:7]])])
        
    return formatted_sites


def get_float(s):
    if '(' in s:
        ind = s.index('(')
        return s[:ind]
    return s


def aborted(output):
    """
    Find if the calculation is terminated from the output file.
    """
    with open(output, 'r') as f:
        lines = f.readlines()
        
        for l in range(-20, 0):
            line = lines[l]

            for phrase in ["Stop in subroutine", "program aborted", "Fatal"]:
                if phrase in line:
                    return True
    return False


def find_issues(output):
    
    """
    Find common issues such as OMMAX and RMAXS and return typical solutions.
    """
    
    solutions = {}
    with open(output, 'r') as f:
        lines = f.readlines()
        total_lines = len(lines)
        for l in range(total_lines):
            line = lines[l]
            
            if "increase RMAXS" in line:
                solutions['modify_STR_RMAXS'] = 0.2
                
            if "change the oxygen muffin-tin radius to e.g. 1.85" in line:
                solutions['set_CLASS_O'] = 1.85
                
            if "increase OMMAX" in line:
                ommax_new = []
                for j in range(l+1, min(total_lines, l+15)):
                    if "%" not in lines[j] and ("OMMAX1=" in lines[j] or "OMMAX2=" in lines[j]):
                        ommax_new.append(lines[j][28:].replace("->", "").strip())
                if len(ommax_new):
                    ommax_new = list(set(ommax_new))
                    recs = ommax_new[-1].split("OMMAX2=")
                    ommax1 = recs[0][7:]
                    ommax2 = recs[1]
                    
                    solutions['set_SCALE_OMMAX1'] = ommax1.strip()
                    solutions['set_SCALE_OMMAX2'] = ommax2.strip()
                    
                # if no OMMAX increase is suggested, set manually
                if "set_SCALE_OMMAX1" not in solutions:
                    solutions["modify_SCALE_OMMAX1"] = 0.02
                    solutions["modify_SCALE_OMMAX2"] = 0.05
                    
    return solutions


def process_dos_data(elements, name):

    p = pexpect.spawn("gnudos.run")
    
    p.expect("Enter output device:")
    p.sendline("1")
    
    p.expect(" energies in Rydberg")
    p.sendline("t")
    
    p.expect("energies relative to EF")
    p.sendline("t")

    p.expect(" Examples")
    r = p.sendline(elements)
    
    classes = str(p.before.decode()).split('\n')
    classes = [c for c in classes if "classes are" in c]
    elem_classes = {}
    if len(classes):
        classes = classes[0].replace("classes are:", "").replace(" \r", "").split()
        _elements = []
        for e in classes:
            for c in e:
                if c.isnumeric():
                    e = e.replace(c, "")
            if e not in _elements:
                _elements.append(e)
        _elements = sorted(_elements, key=lambda x: len(x), reverse=True)

        for e in _elements:
            sites = []
            for c in classes:
                if e in c:
                    sites.append(c)
            elem_classes[e] = ' '.join(sites).strip()
    
    p.expect("correct total DOS")
    p.sendline("/")
    
    p.expect("if desired, enter new emin, emax")
    p.sendline("/")
    
    p.expect("if desired, enter new dosmin, dosmax")
    p.sendline("/")
    
    p.expect("Plot also IDOS")
    p.sendline("t\n")
    
    p.expect("if desired, enter new idosmin")
    p.sendline("/")
    
    p.expect("if desired, enter new title")
    p.sendline("/")
    
    p.expect(pexpect.EOF)

    # read DOS
    if not os.path.isfile('DATA.DOS'):
        return None
    
    df_dos = pd.read_csv('DATA.DOS', delimiter='\s+', names=['Energy (eV)', 'DOS', 'Intg. DOS'], engine='python')
    if elements == "all":
        name += "-total"
    df_dos.to_csv(f"{name}.csv", index=False)
    
    if elements == "all":
        return elem_classes
    
    
def get_band_structure(name):
    
    p = pexpect.spawn("gnubnd.run", timeout=2)
    
    p.expect("Enter output device:", timeout=2)
    p.sendline("1")
    
    p.expect("enter title", timeout=2)
    p.sendline(name)
    
    p.expect("energies in Rydberg", timeout=2)
    p.sendline("t")
    
    p.expect("energies relative to EF", timeout=2)
    p.sendline("t")
    
    p.expect("landscape plot", timeout=2)
    p.sendline("t")
    
    p.expect("show E_nu", timeout=2)
    p.sendline("t")
    
    p.expect("plot orbital character", timeout=2)
    p.sendline("t")
    
    p.expect("energies connected by lines", timeout=2)
    p.sendline("t")
    
    p.expect("enter emin, emax", timeout=2)
    p.sendline("/")

    p.expect(pexpect.EOF)
    
    # read BNDS.GNU
    xticks = []
    with open("BNDS.GNU", 'r') as f:
        xticks = f.readlines()
        xticks = [l for l in xticks if "xtics" in l][1]
        
        for cs in ["set xtics (", "'", ")\n"]:
            xticks = xticks.replace(cs, "")
        xticks = [l.split() for l in xticks.split(",")]
        xticks = [['|'.join(l[0:-1]), float(l[-1])] for l in xticks]
        
        points = pd.DataFrame(dict(point=[i[0] for i in xticks], values=[i[1] for i in xticks]))
        points.to_csv("band_structure_points.csv", index=False)
    
    # read BNDS.DAT
    if not os.path.isfile('BNDS.DAT'):
        return None
    
    data = []
    with open("BNDS.DAT", "r") as f:
        lines = f.readlines()
        for line in lines:
            line = line.strip().split()
            if len(line):
                data.append({'k': float(line[0]), 'Energy (eV)': float(line[1])})
            else:
                data.append({'k': np.nan, 'Energy (eV)': np.nan})
    
    df_dos = pd.DataFrame(data)
    df_dos.to_csv(f"band_structure.csv", index=False)



def process_COHP():

    p = pexpect.spawn("gnucohp.run", timeout=2)

    p.expect("Enter output device:", timeout=2)
    p.sendline("1")

    p.expect("Energy unit is Rydberg", timeout=2)
    p.sendline("t")

    p.expect("energies relative to EF", timeout=2)
    p.sendline("t")

    p.expect("emin,emax", timeout=2)
    p.sendline("/")

    p.expect("Enter class of COHP", timeout=2)
    p.sendline("/")

    p.expect("Now enter weights for each COHP", timeout=2)
    p.sendline("/")

    p.expect("If desired, enter new min, max", timeout=2)
    p.sendline("/")

    p.expect("min/max for COHP", timeout=2)
    p.sendline("/")

    p.expect(" Plot also COHP integration", timeout=2)
    p.sendline("t")

    p.expect("if desired, enter new title", timeout=2)
    p.sendline("/")

    p.expect(pexpect.EOF)
    return


# if __name__ == "__main__":
#     import os
#     f = "/home/aoliynyk/bala/lmto_script/lmto_tests/Mg2Si-1943984/output_lmovl_1.txt"
    
#     print(find_issues(f))



if __name__ == "__main__":

    import os
    f = "/home/lmto/nl/test/Ge12InRu4Y7-1537552"
    
    os.chdir(f)
    process_COHP()
    # ctrl = read_ctrl()
    # write_ctrl(ctrl, **{'modify_SCALE_OMMAX1': 0.02})
    # print(ctrl['CLASS'])
