from .lmto_core_functions import run_lmto
from .utilities import extract_data_from_cif
import multiprocessing as mp
import argparse
import os


# MP
def run_lmto_aux(*args):
    for arg in args:
        try:
            run_lmto(**arg)
        except Exception as e:
            print(f"Error while running {arg['name']}")
            print(e)


def main():

    description = """
        A High Throughput LMTO Calculator.

        This code automates the LMTO calculations using TB-LMTO program.

        """
    parser = argparse.ArgumentParser(
        description=description, formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument("input_path", help="Path to the input file / folder")
    parser.add_argument(
        "-n",
        "--num-cores",
        help="Number of CPU cores to use. \
                        Default is 2.",
        type=int,
        default=2,
        metavar="",
    )
    parser.add_argument(
        "--force-O2", help="Force origin choice 2", action="store_true"
    )

    args = parser.parse_args()

    if os.path.isfile(args.input_path):
        cif_data = extract_data_from_cif(args.input_path)

        if args.force_O2:
            cif_data["origin"] = 2

        if cif_data:
            cif_data["calc_path"] = os.getcwd()
            cif_data["cif_path"] = os.path.join(os.getcwd(), args.input_path)

            run_lmto(**cif_data)

    elif os.path.isdir(args.input_path):
        if args.force_O2:
            print(
                "Origin choice 2 will be used for all files! \
                  \nRestart the calculation without --force-O2 flag \
                  if this is not intended."
            )

        cif_names = [
            c for c in os.listdir(args.input_path) if c.endswith(".cif")
        ]

        tasks = []
        full_path = os.path.join(os.getcwd(), args.input_path)
        for i, cname in enumerate(cif_names, 0):
            cif = cname.split(os.sep)[-1]

            if cif.endswith(".cif"):
                cif_data = extract_data_from_cif(f"{full_path}{os.sep}{cif}")
                cif_data["calc_path"] = full_path
                cif_data["cif_path"] = f"{full_path}{os.sep}{cif}"
                if args.force_O2:
                    cif_data["origin"] = 2
                tasks.append(cif_data)

        with mp.Pool(processes=min(mp.cpu_count(), args.num_cores)) as p:
            p.map(func=run_lmto_aux, iterable=tasks)

        p.close()
        p.join()

    else:
        print(
            "The input path given is not recognized as a file or directory.\
               \nPlease check."
        )
