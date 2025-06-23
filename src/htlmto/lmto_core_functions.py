import os
import time
import subprocess
import pandas as pd
import shutil
from collections import defaultdict
from .cif_reader.base import _parse_formula
from .utilities import print_progress_to_console
from .utilities import get_distances_from_cifkit
from .utilities import convert_cohp_files_to_csv
from .utilities import cleanup
from .lmto_helper_functions import write_INIT_file
from .lmto_helper_functions import aborted
from .lmto_helper_functions import find_issues
from .lmto_helper_functions import modify_CTRL_file
from .lmto_helper_functions import read_ctrl
from .lmto_helper_functions import process_dos_data
from .lmto_helper_functions import process_COHP
from .lmto_helper_functions import extract_scf_data
from .lmto_helper_functions import get_band_structure
from .plotting import plot_dos
from .plotting import plot_cohps
from .plotting import plot_band_structure


@print_progress_to_console
def run_lminit(**kwargs):

    write_INIT_file(**kwargs)
    lminit_output_filename = "output_lminit.log"
    output = open(lminit_output_filename, "a")
    subprocess.run("lminit.run", stdout=output, stderr=output)

    error = aborted(lminit_output_filename)

    return [error]


@print_progress_to_console
def run_lmhart():

    lmhart_output_filename = "output_lmhart.log"
    output = open(lmhart_output_filename, "a")
    subprocess.run("lmhart.run", stdout=output, stderr=output)
    error = aborted(lmhart_output_filename)

    if error:
        return [error, 0.0]

    VOLSPH_by_VOL = [
        line
        for line in open(lmhart_output_filename, "r").readlines()
        if "VOLSPH/VOL" in line
    ][-1]
    try:
        VOLSPH_by_VOL = float(VOLSPH_by_VOL.split("=")[-1].replace("%", ""))
    except Exception as e:
        print("Error reading VOLSPH/VOL!")
        print(e)

    return [error, VOLSPH_by_VOL]


@print_progress_to_console
def run_lmovl(iteration):

    lmvol_output_filename = f"output_lmovl_{iteration}.log"
    output = open(lmvol_output_filename, "a")
    subprocess.run("lmovl.run", stdout=output, stderr=output)
    VOLSPH_by_VOL = [
        line
        for line in open(lmvol_output_filename, "r").readlines()
        if "VOLSPH/VOL" in line
    ][-1]
    overlap = [
        line
        for line in open(lmvol_output_filename, "r").readlines()
        if "Overlap" in line
    ][-1]

    issues = find_issues(lmvol_output_filename)

    try:
        VOLSPH_by_VOL = float(VOLSPH_by_VOL.split("=")[-1].replace("%", ""))
    except Exception as e:
        VOLSPH_by_VOL = None
        print("Error reading VOLSPH/VOL!")
        print(e)

    try:
        overlap = float(overlap.split("=")[-1].replace("%", ""))
    except Exception as e:
        overlap = None
        print("Error reading overlap volume!")
        print(e)

    error = aborted(lmvol_output_filename)

    return [error, VOLSPH_by_VOL, overlap, issues]


@print_progress_to_console
def run_lmes(iteration):

    lmes_output_filename = f"output_lmes_{iteration}.log"
    output = open(lmes_output_filename, "a")

    subprocess.run("lmes.run", stdout=output, stderr=output)
    error = aborted(lmes_output_filename)

    if error:
        VOLSPH_by_VOL = 0.0
    else:
        VOLSPH_by_VOL = [
            line
            for line in open(lmes_output_filename, "r").readlines()
            if "VOLSPH/VOL" in line
        ][-1]
        overlap = [
            line
            for line in open(lmes_output_filename, "r").readlines()
            if "Overlap" in line
        ]
        if len(overlap):
            overlap = overlap[-1]

        try:
            VOLSPH_by_VOL = float(
                VOLSPH_by_VOL.split("=")[-1].replace("%", "")
            )
        except Exception as e:
            VOLSPH_by_VOL = None
            print("Error reading VOLSPH/VOL!")
            print(e)

    return [error, VOLSPH_by_VOL]


@print_progress_to_console
def run_lmctl():

    lmctl_output_filename = "output_lmctl.log"
    output = open(lmctl_output_filename, "a")
    subprocess.run("lmctl.run", stdout=output, stderr=output)

    error = aborted(lmctl_output_filename)
    return [error]


@print_progress_to_console
def run_lmstr(iteration):
    lmstr_output_filename = f"output_lmstr_{iteration}.log"
    output = open(lmstr_output_filename, "w")
    subprocess.run("lmstr.run", stdout=output, stderr=output)

    error = aborted(lmstr_output_filename)
    if error:
        issues = find_issues(lmstr_output_filename)
    else:
        issues = {}
    return [error, issues]


@print_progress_to_console
def run_lm(calc_type, num_atoms, n_try_max=5, get_etots=True):

    converged = False
    n_try = 1
    etot_and_time = []
    while not converged and n_try < n_try_max:
        lm_output_filename = f"output_lm_{calc_type}_{n_try}.log"
        t0 = time.time()
        print(f"iteration {n_try}", end=" ")
        output = open(lm_output_filename, "w")

        subprocess.run("lm.run", stdout=output, stderr=output)
        t = f"{time.time() - t0:.1f}"

        if get_etots:
            etot = extract_scf_data(
                lm_output_filename, n_try, natoms=num_atoms
            )
        with open(lm_output_filename, "r") as f:
            converged = "Jolly good show!" in f.read()

        if get_etots:
            etot_and_time.append(
                {
                    "Type": "first optimization",
                    "Iteration": n_try,
                    "converged": converged,
                    "ETOT (eV)": etot,
                    "Eime (s)": t,
                }
            )
        else:
            etot_and_time = []
        n_try += 1

    return [not converged, etot_and_time]


@print_progress_to_console
def run_lmdos():
    modify_CTRL_file(set_DOS_EMIN=-1.1, set_DOS_EMAX=1.1, set_DOS_NOPTS=1800)

    lmdos_output_filename = "output_lmdos.log"
    output = open(lmdos_output_filename, "a")
    subprocess.run("lmdos.run", stdout=output, stderr=output)

    error = aborted(lmdos_output_filename)
    return [error]


@print_progress_to_console
def run_lmbnd():
    lmbnd_output_filename = "output_lmcbnd.log"
    output = open(lmbnd_output_filename, "a")
    subprocess.run("lmbnd.run", stdout=output, stderr=output)
    error = aborted(lmbnd_output_filename)
    return [error]


def calc_COHPs(cifpath):
    ctrl = read_ctrl()

    class_dict = defaultdict(list)
    # sites = []
    atom_classes = [line for line in ctrl["CLASS"] if "IDMOD" not in line]
    for i, c in enumerate(atom_classes, 1):
        site = c.split("=")[1].split()[0].strip()
        el = list(_parse_formula(site).keys())[0]
        if el == "E":
            continue
        class_dict[el].append([site, i])
        # sites.append([site, i])

    max_distances = get_distances_from_cifkit(cifpath)

    error = False

    print("\tCalculated max interaction distances (A): ")

    for k, v in max_distances.items():
        print(f"\t\t{k:<3} : {v}")

    # CLASS1=1 CLASS2=1 DIMIN=.5 DIMAX=.6
    elements = list(class_dict.keys())
    for element1 in elements:
        element1_sites = class_dict[element1]
        if not element1_sites:
            continue

        ind_el1 = elements.index(element1)
        for element2 in elements[ind_el1:]:
            element2_sites = class_dict[element2]
            if not len(element2_sites):
                continue

            print(f"\n\t\tCalculating COHP for {element1:<2} - {element2:<2}")
            class_pairs = []
            added_pairs = []
            for site1, class_num1 in element1_sites:
                str_site1 = f"{site1}({class_num1})"

                for site2, class_num2 in element2_sites:
                    str_site2 = f"{site2}({class_num2})"

                    if site1 != site2:
                        s_pair = "-".join(sorted([site1, site2]))
                        if s_pair in added_pairs:
                            continue
                        else:
                            added_pairs.append(s_pair)

                    dimax = max(
                        max_distances[site1].get(site2, -1),
                        max_distances[site2].get(site1, -1),
                    )

                    if dimax == -1:
                        print(
                            f"\t\tPair: {str_site1:<6} and {str_site2:<6} \
                                - no distance found"
                        )
                        continue
                    else:
                        print(
                            f"\t\tPair: {str_site1:<6} and {str_site2:<6} \
                                - {dimax:.4f} \u212B"
                        )

                        dimax *= 1.889
                        class_pairs.append(
                            f"CLASS1={class_num1} CLASS2={class_num2} \
                                DIMIN=0.5 DIMAX={dimax:.0f} \n"
                        )

            if class_pairs:
                cohp = [ctrl["COHP"][0]]
                cohp.extend(class_pairs)

                # calculate
                modify_CTRL_file(
                    set_DOS_EMIN=-1.1,
                    set_DOS_EMAX=1.1,
                    set_DOS_NOPTS=1800,
                    set_OPTIONS_COHP="T",
                    set_COHP_ALL=cohp,
                )

                shutil.copy("CTRL", f"bak_cohp_ctrl_{element1}_{element2}")

                error, no_cohp_found = run_cohp(iteration=i)

                if not error and not no_cohp_found:
                    shutil.copy("COHP", f"cohp_{element1}_{element2}")
                    process_COHP()

                    if os.path.isfile("DATA.COHP"):
                        shutil.move(
                            "DATA.COHP", f"data.cohp_{element1}_{element2}"
                        )

                        plot_cohps(".")

    convert_cohp_files_to_csv()

    return error


@print_progress_to_console
def run_cohp(iteration):
    lmincohp_output_filename = f"output_lmincohp_{iteration}.log"
    output = open(lmincohp_output_filename, "a")
    subprocess.run("lmincohp.run", stdout=output, stderr=output)

    error = aborted(lmincohp_output_filename)

    no_cohp_found = False
    if not error:
        with open("COHP", "r") as f:
            line = " ".join(f.readlines()[:5])

            if "NUMBER OF COHPs=  0" in line:
                no_cohp_found = True

    if no_cohp_found:
        # no COHP found
        return error, True

    if not error and not no_cohp_found:
        error, _ = run_lm(
            calc_type=f"cohp{iteration}",
            num_atoms=1,
            n_try_max=5,
            get_etots=False,
        )

    error = False
    if not error:
        lmcohp_output_filename = f"output_lmcohp_{iteration}.log"
        output = open(lmcohp_output_filename, "a")
        subprocess.run("lmcohp.run", stdout=output, stderr=output)

        error = aborted(lmcohp_output_filename)
    error = False
    return [error, no_cohp_found]


def run_lmto(**kwargs):
    os.chdir(kwargs["calc_path"])
    """Run the steps to perform an LMTO calculation.

    1. Write INIT file and run lminit.run to initialize CTRL and CBAK files.
    2. Run lmhart.run and get the VOLSPH_by_VOL
    3. If VOLSPH_by_VOL less than 100.0,
       then run lmes.run and lmovl.run iteratively.
       At each iteration check the output of lmovl.run
       for errors, get recommended
       changes based on the messages in the output file,
       and modify parameters in
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
            for f in os.listdir("."):
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
    error_ctl = run_lmctl()[0]
    error_ovl, VOLSPH_by_VOL, overlap, issues = run_lmovl(iteration=1)

    if VOLSPH_by_VOL >= 100:
        error_es = False
    else:
        # set VERBOSE to 50
        n_try = 1
        print(f"\t\tVOLSPH_by_VOL={VOLSPH_by_VOL}, overlap={overlap}")
        while (
            VOLSPH_by_VOL is not None and VOLSPH_by_VOL < 100.0 and n_try < 30
        ):
            print(f"\tRunning lmes.run and lmovl.run, iteration {n_try}")
            error_ctl = run_lmctl()[0]
            if len(issues):
                # print(f"\tChanging params based on recomendations {issues}")
                modify_CTRL_file(**issues)

            error_es, VOLSPH_by_VOL = run_lmes(iteration=n_try)
            if not error_es:
                error_ovl, VOLSPH_by_VOL, overlap, issues = run_lmovl(
                    iteration=1 + n_try
                )

            if error_es or error_ovl:
                break

            print(f"\t\tVOLSPH_by_VOL={VOLSPH_by_VOL}, overlap={overlap}")
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
            error_str, issues = run_lmstr(iteration=1 + n_try)
            n_try += 1

    if error_str:
        print(f"{kwargs['name']} failed")
        return True

    etot_and_time = []
    not_converged, etot_t = run_lm(
        calc_type="first", num_atoms=kwargs["num_atoms"], n_try_max=15
    )
    etot_and_time.extend(etot_t)
    if not_converged:
        print(f"{kwargs['name']} failed")
        return True

    # 2x K-point, BEGMOM=F
    modify_CTRL_file(set_START_BEGMOM="T", modify_BZ_NKABC=2)
    not_converged, etot_t = run_lm(
        calc_type="2xKPTS", num_atoms=kwargs["num_atoms"], n_try_max=5
    )
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

    elem_classes = process_dos_data(elements="all", name="DOS")  # TDOS
    for k, v in elem_classes.items():
        process_dos_data(elements=v, name=f"DOS-{k}")
    plot_dos(".")

    # Band structure
    error_bnd = run_lmbnd()[0]
    if error_bnd:
        print(f"{kwargs['name']} failed")
        return True
    get_band_structure(kwargs["name"].split("-")[0])
    plot_band_structure(".")

    # COHP
    error_cohp = calc_COHPs(kwargs["cif_path"])

    cleanup()

    if not any(
        [
            error_init,
            error_ctl,
            error_hart,
            error_ovl,
            error_es,
            error_str,
            not_converged,
            error_dos,
            error_bnd,
            error_cohp,
        ]
    ):
        print(f"{kwargs['name']} gracefully exited!")
    else:
        print(f"{kwargs['name']} failed")
        print(
            dict(
                zip(
                    [
                        "init",
                        "ctl",
                        "hart",
                        "ovl",
                        "es",
                        "str",
                        "lm",
                        "dos",
                        "band",
                        "cohp",
                    ],
                    [
                        error_init,
                        error_ctl,
                        error_hart,
                        error_ovl,
                        error_es,
                        error_str,
                        not_converged,
                        error_dos,
                        error_bnd,
                        error_cohp,
                    ],
                )
            )
        )
    return (
        error_init
        or error_hart
        or error_ovl
        or error_es
        or error_str
        or not_converged
        or error_dos
        or error_bnd
        or error_cohp
    )
