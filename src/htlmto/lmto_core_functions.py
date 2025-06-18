import os
import time
import subprocess
import pandas as pd
from collections import defaultdict
import shutil
from .utilities import print_to_console
from .lmto_helper_functions import write_INIT_file
from .lmto_helper_functions import aborted
from .lmto_helper_functions import find_issues
from .lmto_helper_functions import modify_CTRL_file
from .lmto_helper_functions import read_ctrl
from .lmto_helper_functions import process_dos_data
from .lmto_helper_functions import process_COHP
from .lmto_helper_functions import extract_scf_data
from .lmto_helper_functions import get_band_structure
from .utilities import get_distances_from_cifkit
from .cif_reader.base import _parse_formula


@print_to_console
def run_lminit(**kwargs):

    write_INIT_file(**kwargs)
    lminit_output_filename = "output_lminit.txt"
    output = open(lminit_output_filename, "a")
    subprocess.run("lminit.run", stdout=output, stderr=output)

    error = aborted(lminit_output_filename)

    return [error]


@print_to_console
def run_lmhart():

    lmhart_output_filename = "output_lmhart.txt"
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


@print_to_console
def run_lmovl(iteration):

    lmvol_output_filename = f"output_lmovl_{iteration}.txt"
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


@print_to_console
def run_lmes(iteration):

    lmes_output_filename = f"output_lmes_{iteration}.txt"
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


@print_to_console
def run_lmctl():

    lmctl_output_filename = "output_lmctl.txt"
    output = open(lmctl_output_filename, "a")
    subprocess.run("lmctl.run", stdout=output, stderr=output)

    error = aborted(lmctl_output_filename)
    return [error]


@print_to_console
def run_lmstr(iteration):
    lmstr_output_filename = f"output_lmstr_{iteration}.txt"
    output = open(lmstr_output_filename, "w")
    subprocess.run("lmstr.run", stdout=output, stderr=output)

    error = aborted(lmstr_output_filename)
    if error:
        issues = find_issues(lmstr_output_filename)
    else:
        issues = {}
    return [error, issues]


@print_to_console
def run_lm(calc_type, num_atoms, n_try_max=5, get_etots=True):

    converged = False
    n_try = 1
    etot_and_time = []
    while not converged and n_try < n_try_max:
        lm_output_filename = f"output_lm_{calc_type}_{n_try}.txt"
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


@print_to_console
def run_lmdos():
    modify_CTRL_file(set_DOS_EMIN=-1.1, set_DOS_EMAX=1.1, set_DOS_NOPTS=1800)

    lmdos_output_filename = "output_lmdos.txt"
    output = open(lmdos_output_filename, "a")
    subprocess.run("lmdos.run", stdout=output, stderr=output)

    error = aborted(lmdos_output_filename)
    return [error]


@print_to_console
def run_lmbnd():
    lmbnd_output_filename = "output_lmcbnd.txt"
    output = open(lmbnd_output_filename, "a")
    subprocess.run("lmbnd.run", stdout=output, stderr=output)
    error = aborted(lmbnd_output_filename)
    return [error]


def calc_COHPs(cifpath):
    ctrl = read_ctrl()

    class_dict = defaultdict(list)
    sites = []
    atom_classes = [line for line in ctrl["CLASS"] if "IDMOD" not in line]
    for i, c in enumerate(atom_classes, 1):
        site = c.split("=")[1].split()[0].strip()
        el = list(_parse_formula(site).keys())[0]
        if el == "E":
            continue
        class_dict[el].append(site)
        sites.append([site, i])

    max_distances = get_distances_from_cifkit(cifpath)

    # CLASS1=1 CLASS2=1 DIMIN=.5 DIMAX=.6
    for i in range(len(sites)):
        site1, class_num1 = sites[i]
        class_pairs = []
        for j in range(i, len(sites)):
            site2, class_num2 = sites[j]

            dimax = max_distances[site1].get(site2, None)
            if dimax is None:
                print(f"\t\tNo distance data for {site1} and {site2}.")
                continue

            dimax *= 1.889
            class_pairs.append(
                f"CLASS1={class_num1} CLASS2={class_num2} \
                    DIMIN=0.5 DIMAX={dimax:.0f} \n"
            )

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

        # modify_CTRL_file(set_OPTIONS_COHP="T")
        # modify_CTRL_file(set_COHP_ALL=cohp)

        error, no_cohp_found = run_cohp(iteration=i)

        if not error and no_cohp_found:
            shutil.copy("COHP", f"COHP_{i}")
            process_COHP()

            if os.path.isfile("DATA.COHP"):
                shutil.move("DATA.COHP", f"DATA.COHP_{i}")

    return error


@print_to_console
def run_cohp(iteration):
    lmincohp_output_filename = f"output_lmincohp_{iteration}.txt"
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

    if not error and no_cohp_found:
        error, _ = run_lm(
            calc_type=f"cohp{iteration}",
            num_atoms=1,
            n_try_max=5,
            get_etots=False,
        )

    error = False
    if not error:
        lmcohp_output_filename = f"output_lmcohp_{iteration}.txt"
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

    # Band structure
    error_bnd = run_lmbnd()[0]
    if error_bnd:
        print(f"{kwargs['name']} failed")
        return True
    get_band_structure(kwargs["name"].split("-")[0])

    # COHP
    error_cohp = calc_COHPs(kwargs["cif_path"])

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


# if __name__ == "__main__":

#     wd = "/home/lmto/nl/test/Ge12InRu4Y7-1537552"
#     os.chdir(wd)
#     calc_COHPs("/home/lmto/nl/test/1537552.cif")
#     # get_distances_from_cifkit("/home/bala/Documents
# /22_lmto/High-throughput-LMTO/Cr4AlB4.cif")
#     exit(0)

#     # run 1
#     cif = sys.argv[1]
#     cif_data = extract_data_from_cif(f"{cif}")
#     cif_data["calc_path"] = "/home/bala/research/
# 1_LMTO/lmto_script"
#     run_lmto(**cif_data)

#     exit(0)

#     import multiprocessing as mp
#     import functools

#     def std_wrapper(func):
#         @functools.wraps(
#             func
#         )  # we need this to unravel the target function name
#         def caller(
#             *args, **kwargs
#         ):  # and now for the wrapper, nothing new here

#             from io import StringIO
#             import sys

#             sys.stdout, sys.stderr = (
#                 StringIO(),
#                 StringIO(),
#             )  # use our buffers instead
#             response = None  # in case a call fails
#             try:
#                 response = func(
#                     *args, **kwargs
#                 )  # call our wrapped process function
#             except Exception as e:  # too broad but good
# enough as an example
#                 print(
#                     e
#                 )  # NOTE: the exception is also printed
# to the captured STDOUT
#             # rewind our buffers:
#             sys.stdout.seek(0)
#             sys.stderr.seek(0)
#             # return everything packed as STDOUT, STDERR,
# PROCESS_RESPONSE | NONE
#             return sys.stdout.read(), sys.stderr.read(), response

#         return caller

#     wd = "/home/bala/research/1_LMTO/lmto_script"

#     skip_cifs = [
#         "262117.cif",
#         "1935906.cif",  # VOLSPH_by_VOL=98.7
#         "526719.cif",
#         "303534.cif",  #  NiS6V3 lmstr.run Fatal  :
# MSTRX2: nlsqri=25gt nlsqr=16
#         # lm.run
#         "1719812.cif",  # E=-7.56750D+03 DE= 9.75079D+06 after 6 runs
#     ]

#     # MP
#     @std_wrapper
#     def run_lmto_aux(*args):
#         for arg in args:
#             run_lmto(**arg)

#     df = pd.read_csv("testfiles.csv")
#     total = len(df)
#     cif_names = df["fname"].tolist()  # [:1]
#     tasks = []

#     for i, cname in enumerate(cif_names, 0):
#         cif = cname.split(os.sep)[-1]

#         if cif.endswith(".cif"):
#             if cif in skip_cifs:
#                 continue
#             cif_data = extract_data_from_cif
# (f"{wd}/not_prototype_CIFs/{cif}")
#             cif_data["calc_path"] = f"{wd}/lmto_tests"
#             tasks.append(cif_data)

#     print(len(tasks))

#     pool = mp.Pool(processes=32)
#     for out, err, res in pool.imap_unordered(run_lmto_aux, tasks):
#         # print(out)
#         # get name

#         lines = out.split("\n")
#         name = None
#         failed = True
#         for l in lines:
#             if "Name:" in l:
#                 name = l.split("Name:")[-1].strip()

#             if "gracefully exited!" in l:
#                 fname = l.split()[0].strip()
#                 failed = False

#         if failed:
#             print(name)
#             out_fname = f"{wd}/lmto_tests/{name}/aout.txt"
#         else:
#             out_fname = f"{wd}/lmto_tests/noi/{name}/aout.txt"

#         try:
#             with open(out_fname, "w") as f:
#                 f.write(out)
#         except:
#             print(out)
#             exit(0)

#     pool.close()

#     # with mp.Pool(processes=32) as p:
#     #     p.map(func=run_lmto_aux, iterable=tasks)

#     # p.close()
#     # p.join()
