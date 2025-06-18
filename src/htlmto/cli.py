from .lmto_core_functions import run_lmto
from .utilities import extract_data_from_cif
import argparse
import os

# if len(sys.argv) != 2:
#     print(
#         "Please provide the path to a directory \
#             containing .cif files as the only argument."
#     )


# # MP
# def run_lmto_aux(*args):
#     for arg in args:
#         run_lmto(**arg)


# cif_names = [c for c in os.listdir(sys.argv[1]) if c.endswith(".cif")]

# total = len(cif_names)
# tasks = []

# for i, cname in enumerate(cif_names, 0):
#     cif = cname.split(os.sep)[-1]

#     if cif.endswith(".cif"):
#         cif_data = extract_data_from_cif(f"{sys.argv[1]}{os.sep}{cif}")
#         cif_data["calc_path"] = sys.argv[1]
#         cif_data["cif_path"] = f"{sys.argv[1]}{os.sep}{cif}"
#         tasks.append(cif_data)

# with mp.Pool(processes=2) as p:
#     p.map(func=run_lmto_aux, iterable=tasks)

# p.close()
# p.join()


def main():

    description = """
        A High Throughput LMTO Calculator.

        This code automates the LMTO calculations using TB-LMTO program.

        """
    parser = argparse.ArgumentParser(
        description=description, formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument("input_path", help="Path to the input file / folder")

    args = parser.parse_args()
    cif_data = extract_data_from_cif(args.input_path)
    cif_data["calc_path"] = os.getcwd()
    cif_data["cif_path"] = os.path.join(os.getcwd(), args.input_path)

    run_lmto(**cif_data)


if __name__ == "__main__":
    main()
