# from htlmto.src.htlmto import cif_reader
from htlmto import cif_reader
import os


def check_cif_readers():
    root = "/home/bala/Documents/00_Utils/utilities/tests/test_cifs"
    cdir = "allc"
    for cif in os.listdir(f"{root}{os.sep}{cdir}"):
        if not cif.endswith("cif"):
            continue
        cif_reader.read_cif(f"{root}{os.sep}{cdir}{os.sep}{cif}", verbose=True)


if __name__ == "__main__":
    check_cif_readers()
