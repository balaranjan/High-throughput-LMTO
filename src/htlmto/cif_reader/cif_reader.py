from cif_reader import read_cif
import os


if __name__ == "__main__":

    root = "/home/bala/Documents/00_Utils/utilities/tests/test_cifs"
    cdir = "allc"
    for cif in os.listdir(f"{root}{os.sep}{cdir}"):
        if not cif.endswith("cif"):
            continue

        read_cif(f"{root}{os.sep}{cdir}{os.sep}{cif}", verbose=True)
        # if "shellxl" in cdir:
        #     print(cif)
        #     SXL_reader(f"{root}{os.sep}{cdir}{os.sep}{cif}")
