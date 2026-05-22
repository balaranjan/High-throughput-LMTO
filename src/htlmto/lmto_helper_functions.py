from collections import defaultdict
import pandas as pd
import numpy as np
import pexpect
import shutil
import os
import re


def write_INIT_file(**kwargs):
    """Write INIT file for LMTO calculation using cif data obtained from
    utils.extract_data_from_cif()."""

    atomic_numbers = {
        "H": 1,
        "He": 2,
        "Li": 3,
        "Be": 4,
        "B": 5,
        "C": 6,
        "N": 7,
        "O": 8,
        "F": 9,
        "Ne": 10,
        "Na": 11,
        "Mg": 12,
        "Al": 13,
        "Si": 14,
        "P": 15,
        "S": 16,
        "Cl": 17,
        "Ar": 18,
        "K": 19,
        "Ca": 20,
        "Sc": 21,
        "Ti": 22,
        "V": 23,
        "Cr": 24,
        "Mn": 25,
        "Fe": 26,
        "Co": 27,
        "Ni": 28,
        "Cu": 29,
        "Zn": 30,
        "Ga": 31,
        "Ge": 32,
        "As": 33,
        "Se": 34,
        "Br": 35,
        "Kr": 36,
        "Rb": 37,
        "Sr": 38,
        "Y": 39,
        "Zr": 40,
        "Nb": 41,
        "Mo": 42,
        "Tc": 43,
        "Ru": 44,
        "Rh": 45,
        "Pd": 46,
        "Ag": 47,
        "Cd": 48,
        "In": 49,
        "Sn": 50,
        "Sb": 51,
        "Te": 52,
        "I": 53,
        "Xe": 54,
        "Cs": 55,
        "Ba": 56,
        "La": 57,
        "Ce": 58,
        "Pr": 59,
        "Nd": 60,
        "Pm": 61,
        "Sm": 62,
        "Eu": 63,
        "Gd": 64,
        "Tb": 65,
        "Dy": 66,
        "Ho": 67,
        "Er": 68,
        "Tm": 69,
        "Yb": 70,
        "Lu": 71,
        "Hf": 72,
        "Ta": 73,
        "W": 74,
        "Re": 75,
        "Os": 76,
        "Ir": 77,
        "Pt": 78,
        "Au": 79,
        "Hg": 80,
        "Tl": 81,
        "Pb": 82,
        "Bi": 83,
        "Po": 84,
        "At": 85,
        "Rn": 86,
        "Fr": 87,
        "Ra": 88,
        "Ac": 89,
        "Th": 90,
        "Pa": 91,
        "U": 92,
        # "Np": 93,
        # "Pu": 94,
        # "Am": 95,
        # "Cm": 96,
        # "Bk": 97,
        # "Cf": 98,
        # "Es": 99,
        # "Fm": 100,
        # "Md": 101,
        # "No": 102,
        # "Lr": 103,
        # "Rf": 104,
        # "Db": 105,
        # "Sg": 106,
        # "Bh": 107,
        # "Hs": 108,
        # "Mt": 109,
        # "Ds": 110,
        # "Rg": 111,
        # "Cn": 112,
        # "Nh": 113,
        # "Fl": 114,
        # "Mc": 115,
        # "Lv": 116,
        # "Ts": 117,
        # "Og": 118,
    }

    common_subs = {
        "Ce": "La",
        "Pr": "La",
        "Nd": "La",
        "Pm": "La",
        "Sm": "La",
        "Eu": "La",
        "Gd": "La",
        "Tb": "Lu",
        "Dy": "Lu",
        "Ho": "Lu",
        "Er": "Lu",
        "Tm": "Lu",
        "Yb": "Lu",
        "Lu": "Lu",
        "Hf": "Zr",
        "Th": "Ac",
        "Pa": "Ac",
        "U": "Ac",
        "Np": "Ac",
        "Pu": "Ac",
        "Am": "Ac",
        "Cm": "Ac",
        "Bk": "Ac",
        "Cf": "Ac",
        "Es": "Ac",
        "Fm": "Ac",
        "Md": "Ac",
        "No": "Ac",
        "Lr": "Ac",
    }

    element_subs = kwargs.pop("sub_elements")

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

    param_map = {
        "A": "_cell_length_a",
        "B": "_cell_length_b",
        "C": "_cell_length_c",
        "a": "_cell_angle_alpha",
        "b": "_cell_angle_beta",
        "g": "_cell_angle_gamma",
    }
    params = []
    sg_params = space_group_nums_params[int(kwargs["_space_group_IT_number"])]

    if "_" in sg_params:
        params, angles = sg_params.split("_")
        params = params.split("-") if "-" in params else [params]
        angles = angles.split("-") if "-" in angles else [angles]
    else:
        params = sg_params.split("-") if "-" in sg_params else [sg_params]
        angles = []

    params = [param_map[p] for p in params]
    angles = [param_map[a] for a in angles]

    lines = []
    sg_number = kwargs["_space_group_IT_number"]
    first_ln = f"          SPCGRP={sg_number}\
         IORIGIN={kwargs.get('origin', 1)} \
            ATUNITS={kwargs.get('aomic_units', 'F')}"
    for k, param in [
        ["A", "_cell_length_a"],
        ["B", "_cell_length_b"],
        ["C", "_cell_length_c"],
    ]:
        if param in kwargs and param in params:
            first_ln += f" {k}={kwargs[param]}"

    for k, param in (
        ["ALPHA", "_cell_angle_alpha"],
        ["BETA", "_cell_angle_beta"],
        ["GAMMA", "_cell_angle_gamma"],
    ):
        if param in kwargs and param in angles:
            first_ln += f" {k}={kwargs[param]}"

    first_ln += "\n"
    lines.append(first_ln)

    out_subs = {}
    for site in kwargs["atom_site_data"]:
        S = site["symbol"]

        if element_subs and S in element_subs:
            S = element_subs.get(S)

        elif S in common_subs:
            S = common_subs.get(S)

        if S != site["symbol"]:
            out_subs[site["symbol"]] = S
            print(f"{site['symbol']} is substituted with {S}", end=" ")

        lines.append(
            f"          ATOM={site['label']:<3}  \
                Z={atomic_numbers[S]:<3}   \
                    X={site['x']} {site['y']} {site['z']}\n"
        )

    with open("INIT", "w") as f:
        f.writelines(lines)

    if len(out_subs):
        with open("SUBS", "w") as f:
            for k, v in out_subs.items():
                f.write(f"{k}: {v}\n")

    return True


def extract_scf_data(lm_output_filename, n_opt, natoms):
    # write iter, etot, dE as csv; return last etot
    with open(lm_output_filename, "r") as f:
        lines = f.readlines()
        # etot = [l for l in lines if " ETOT" in l]

        iter_etot = [
            line for line in lines if "ITER" in line and "OUT OF" in line
        ]
        iter_etot = [
            [
                int(line.split("OUT OF")[0][5:]),
                float(line.split("ETOT=")[-1][:-1]),
            ]
            for line in iter_etot
        ]
        iter_etot = pd.DataFrame(
            {
                "n_opt": [n_opt for _ in range(len(iter_etot))],
                "type": ["initial" for _ in range(len(iter_etot))],
                "iteration": [i[0] for i in iter_etot],
                "Energy (mRy)": [i[1] for i in iter_etot],
            }
        )

        iter_etot["Energy (eV)"] = (
            iter_etot["Energy (mRy)"] * 13.6057039763 / 1000.0
        )
        iter_etot["Energy (eV/atom)"] = iter_etot["Energy (eV)"] / natoms

        iter_etot.to_csv(
            f"ETOT_values_{lm_output_filename.split('_')[2]}.csv", index=False
        )

        return float(iter_etot.iloc[-1]["Energy (eV)"])


def extract_params_from_line(line):
    """Extract all parameters and their values from a line from CTRL file."""
    line = line.strip()
    _value = line.split()
    value = []
    if "ATOM" in line and "SIGMA" in line:
        return [f"STRATOMLINE=={line}"]
    for v in _value:
        if "=" in v:
            value.append(v)
        else:
            value[-1] += f" {v}"
    return value


def read_ctrl():
    """Read the CTRL file and return as dictionary."""
    ctrl = {}
    atom_counter = 1
    with open("CTRL", "r") as f:
        lines = f.readlines()
        i = 0
        while i < len(lines):
            line = lines[i]
            if not line.startswith(" " * 10):
                # new block
                block_name = line[:10].strip()

                # read block
                if block_name in ["HEADER", "VERS"]:
                    block_value = line[10:]

                elif block_name in [
                    "IO",
                    "OPTIONS",
                    "STR",
                    "BZ",
                    "DOS",
                    "SCALE",
                ]:

                    # dictionaries maintain their insertion order
                    # as of version 3.7
                    value = extract_params_from_line(line[10:])
                    block_value = {}
                    for v in value:
                        v = v.split("=")
                        block_value[v[0]] = v[1]

                    in_block = lines[i + 1].startswith(" " * 10)
                    j = i + 1

                    while in_block:
                        value = extract_params_from_line(lines[j])
                        for v in value:

                            if "STRATOMLINE" in v:
                                v = v.split("==")
                                block_value[f"{v[0]}{atom_counter:02}"] = v[1]
                                atom_counter += 1
                            else:
                                v = v.split("=")
                                block_value[v[0]] = v[1]
                        in_block = (
                            lines[j + 1].startswith(" " * 10)
                            if j + 1 < len(lines)
                            else False
                        )
                        i = j
                        j += 1

                else:
                    in_block = len(lines) > i + 1 and lines[i + 1].startswith(
                        " " * 10
                    )
                    j = i + 1
                    block_value = [line[10:]]
                    while in_block:
                        block_value.append(lines[j][10:])
                        in_block = (
                            len(lines) > j + 1
                            and lines[j + 1].startswith(" " * 10)
                            if j + 1 < len(lines)
                            else False
                        )
                        # lines[j+1].startswith(" "*10) if j+1 <
                        # len(lines) else False
                        i = j
                        j += 1
                ctrl[block_name] = block_value

            i += 1
    return ctrl


def write_ctrl(ctrl, **kwargs):
    """Write CTRL file.  Modifications for certain blocks in the CTRL file can
    be supplied as keyword arguments.

    The keyword will have the following format:
    {modify/set}_BLOCK_PARAMETER=VALUE

    When the modify prefix is used, the existing value(s)
    will be incresed by VALUE.
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
                new_val = [float(c) + v for c in current_val.split()]
                new_val = " ".join([f"{c:.2f}"[1:] for c in new_val])
            elif block == "BZ" and param == "NKABC":
                new_val = [
                    str(int(i) * v) for i in current_val.strip().split()
                ]
                new_val = " ".join(new_val)
            else:
                if current_val.strip().isdigit():
                    new_val = int(current_val) + v
                else:
                    nround = len(str(v).split(".")[-1])
                    new_val = round(float(current_val) + v, nround)

            ctrl[block][param] = new_val

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

                    element = atom
                    for e in element:
                        if e.isnumeric():
                            element = element.replace(e, "")
                    if element == param and float(radius) < v:
                        current_val[i] = current_val[i].replace(
                            f"R={radius}", f"R={v}"
                        )

                ctrl[block] = current_val

            elif block == "START" and param == "BEGMOM":
                current_val = "F"
                ctrl[block][1] = ctrl[block][1].replace(
                    "BEGMOM=F", f"BEGMOM={v}"
                )
            elif block == "COHP":
                ctrl[block] = v  # "\n".join(v)
            else:
                current_val = ctrl[block][param]
                ctrl[block][param] = f"{v}"

    max_len = 70

    with open("CTRL", "w") as f:
        for k, v in ctrl.items():

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
                            line = " " * 10
                            line += kv
                if not block.endswith(line):
                    block += f"{line}\n"

                f.write(block)
            elif isinstance(v, str):
                # print(f"{k:<10}{v}")
                f.write(f"{k:<10}{v}")
            else:
                f.write(f"{k:<10}{v[0]}")
                for _v in v[1:]:
                    # f.write(f"          {_v.replace("\n", "")}\n")
                    f.write(f"          {_v}")


def modify_CTRL_file(**kwargs):
    ctrl = read_ctrl()
    write_ctrl(ctrl, **kwargs)


def format_sites(sites: str) -> list:
    """Receives site data from CIF and returns reformatted labels and
    coordinates for LMTO.

    Change Ce - Lu to La
    """
    sites = [s.split() for s in sites]

    labels = defaultdict(int)
    formatted_sites = []

    change_map = {"D": "H", "Tl": "In", "Hf": "Zr"}

    for site in sites:
        site[1] = change_map.get(site[1], site[1])

        labels[site[1]] += 1

        formatted_sites.append(
            [
                f"{site[1]}{labels[site[1]]}",
                " ".join([str(p) for p in site[4:7]]),
            ]
        )

    return formatted_sites


def aborted(output):
    """Find if the calculation is terminated from the output file."""
    with open(output, "r") as f:
        lines = f.readlines()

        for line_no in range(-20, 0):
            line = lines[line_no]

            for phrase in ["Stop in subroutine", "program aborted", "Fatal"]:
                if phrase in line:
                    return True
    return False


def find_issues(output):
    """Find common issues such as OMMAX and RMAXS and return typical
    solutions."""

    solutions = {}
    with open(output, "r") as f:
        lines = f.readlines()
        total_lines = len(lines)
        for line_no in range(total_lines):
            line = lines[line_no]

            if "increase RMAXS" in line:
                solutions["modify_STR_RMAXS"] = 0.2

            if "change the oxygen muffin-tin radius to e.g. 1.85" in line:
                solutions["set_CLASS_O"] = 1.85

            if "increase OMMAX" in line:
                ommax_new = []
                for j in range(line_no + 1, min(total_lines, line_no + 15)):
                    if "%" not in lines[j] and (
                        "OMMAX1=" in lines[j] or "OMMAX2=" in lines[j]
                    ):
                        ommax_new.append(
                            lines[j][28:].replace("->", "").strip()
                        )
                if len(ommax_new):
                    ommax_new = list(set(ommax_new))
                    recs = ommax_new[-1].split("OMMAX2=")
                    ommax1 = recs[0][7:]
                    ommax2 = recs[1]

                    solutions["set_SCALE_OMMAX1"] = ommax1.strip()
                    solutions["set_SCALE_OMMAX2"] = ommax2.strip()

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
    p.sendline(elements)

    classes = str(p.before.decode()).split("\n")
    classes = [c for c in classes if "classes are" in c]
    elem_classes = {}
    if len(classes):
        classes = (
            classes[0].replace("classes are:", "").replace(" \r", "").split()
        )
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
            elem_classes[e] = " ".join(sites).strip()

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
    if not os.path.isfile("DATA.DOS"):
        return None

    df_dos = pd.read_csv(
        "DATA.DOS",
        delimiter="\\s+",
        names=["Energy (eV)", "DOS", "Intg. DOS"],
        engine="python",
    )
    if elements == "all":
        name += "-total"
    df_dos.to_csv(f"{name}.csv", index=False)

    shutil.copy("DATA.DOS", f"data.dos_{name}")

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
    with open("BNDS.GNU", "r") as f:
        xticks = f.readlines()
        xticks = [line for line in xticks if "xtics" in line][1]

        for cs in ["set xtics (", "'", ")\n"]:
            xticks = xticks.replace(cs, "")
        xticks = [line.split() for line in xticks.split(",")]
        xticks = [["|".join(line[0:-1]), float(line[-1])] for line in xticks]

        points = pd.DataFrame(
            dict(point=[i[0] for i in xticks], values=[i[1] for i in xticks])
        )
        points.to_csv("band_structure_points.csv", index=False)

    # read BNDS.DAT
    if not os.path.isfile("BNDS.DAT"):
        return None

    data = []
    with open("BNDS.DAT", "r") as f:
        lines = f.readlines()
        for line in lines:
            line = line.strip().split()
            if len(line):
                data.append(
                    {"k": float(line[0]), "Energy (eV)": float(line[1])}
                )
            else:
                data.append({"k": np.nan, "Energy (eV)": np.nan})

    df_dos = pd.DataFrame(data)
    df_dos.to_csv("band_structure.csv", index=False)


def process_COHP(selected_interactions=None):

    output = []

    def capture_output(p):
        """Helper to decode and append the received text before each prompt."""
        if p.before:
            output.append(p.before.decode("utf-8", errors="ignore"))

    p = pexpect.spawn("gnucohp.run", timeout=2)
    # p.logfile_read = sys.stdout.buffer  # log all output to stdout

    p.expect("Enter output device:", timeout=2)
    p.sendline("1")

    p.expect("Energy unit is Rydberg", timeout=2)
    p.sendline("t")

    p.expect("energies relative to EF", timeout=2)
    p.sendline("t")

    p.expect(r"Examples:.*\r\n", timeout=10)  # match the examples line
    if not selected_interactions:
        capture_output(p)

    if selected_interactions:
        p.sendline(selected_interactions)
    else:
        p.sendline("/")

    p.expect("Now enter weights for each COHP", timeout=10)
    p.sendline("/")

    # Energy range
    p.expect("If desired, enter new min, max", timeout=10)
    p.sendline("/")

    # COHP range
    p.expect("If desired, enter new min, max", timeout=10)
    p.sendline("/")

    p.expect("Plot also COHP integration", timeout=2)
    p.sendline("t")  # use 'f' unless you actually need ICOHP plot

    # ICOHP range
    p.expect("If desired, enter new min, max", timeout=10)
    p.sendline("/")

    p.expect("if desired, enter new title", timeout=2)
    p.sendline("/")

    p.expect(pexpect.EOF)

    if selected_interactions:
        output = get_coph_at_ef_zero()

    return output


def parse_classes(coph_text):
    results = []
    for line in coph_text[0].split("\n")[1:]:
        line = line.replace("COHP  are:", "")
        m = re.match(
            r"[^A-Za-z]*([A-Za-z0-9]+)\s*=\s*(\w+)-\d+/\s*(\w+)-\d+\s*.+:\s*([\d.]+)\s+Ang\.",  # noqa: E501
            line,
        )
        if m:
            results.append(
                [m.group(1), f"{m.group(2)}-{m.group(3)}", float(m.group(4))]
            )

    return results


def group_cohps(bonds):
    """Group COHP bonds by site pair and distance."""
    groups = defaultdict(list)
    for label, atom_pair, distance in bonds:
        groups[f"{atom_pair}-{round(distance, 3)}"].append(label)

    return dict(groups)


def get_coph_at_ef_zero():
    if not os.path.isfile("DATA.COHP"):
        return (-2, -2)  # Return tuple indicating error

    data_text = open("DATA.COHP", "r").read()
    rows = []
    for line in data_text.split("\n"):
        parts = line.split()
        if len(parts) >= 3:
            e, col2_val, col3_val = map(float, parts[:3])
            rows.append((e, col2_val, col3_val))

    # Find crossing of 0 in energy (column 1)
    for i in range(len(rows) - 1):
        e1, c2_1, c3_1 = rows[i]
        e2, c2_2, c3_2 = rows[i + 1]

        if (e1 < 0 and e2 >= 0) or (e1 > 0 and e2 <= 0):
            # Linear interpolation for col2 and col3 at e=0
            frac = -e1 / (e2 - e1)
            interpolated_c2 = c2_1 + (c2_2 - c2_1) * frac
            interpolated_c3 = c3_1 + (c3_2 - c3_1) * frac
            return (interpolated_c2, interpolated_c3)

    return (-1, -1)


if __name__ == "__main__":
    parts = process_COHP()
    parse_classes(parts)
