from .pcd import PCD_reader
from .icsd import ICSD_reader
from .shellxl import SXL_reader
import os


def read_cif(filename, verbose=False):

    if not os.path.isfile(filename):
        print(f"File {filename} not found!")
        return

    if verbose:
        print(filename, end="   ... ")

    contents = open(filename, "r").read()

    if "PCD" in contents:
        if verbose:
            print(" is from PCD")
        return PCD_reader(filename, verbose=verbose)
    elif "ICSD" in contents:
        if verbose:
            print(" is from ICSD")
        return ICSD_reader(filename, verbose=verbose)
    elif "shelx" in contents:
        if verbose:
            print(" is from ShellXL")
        return SXL_reader(filename, verbose=verbose)
    else:
        print(f"File {filename} form unknown source.")
        return
