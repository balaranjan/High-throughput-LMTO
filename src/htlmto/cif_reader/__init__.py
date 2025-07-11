from .pcd import PCD_reader
from .icsd import ICSD_reader
from .shellxl import SXL_reader
from .vesta import VESTA_reader
from .critic import Critic_reader
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
        return PCD_reader(filename, data_source="PCD", verbose=verbose)
    elif "ICSD" in contents:
        if verbose:
            print(" is from ICSD")
        return ICSD_reader(filename, data_source="ICSD", verbose=verbose)
    elif "shelx" in contents:
        if verbose:
            print(" is from ShellXL")
        return SXL_reader(filename, data_source="ShellXL", verbose=verbose)
    elif "VESTA" in contents:
        if verbose:
            print(" is from VESTA")
        return VESTA_reader(filename, data_source="VESTA", verbose=verbose)
    elif "critic2" in contents:
        if verbose:
            print(" is from Critic2")
        return Critic_reader(filename, data_source="Critic2", verbose=verbose)
    else:
        print(f"File {filename} form unknown source.")
        return
