# from htlmto.src.htlmto import cif_reader
from htlmto import cif_reader
import os


def test_cif_readers():
    CIF_DIR = os.path.join(os.path.dirname(__file__), "cifs")
    for cif in os.listdir(CIF_DIR):
        if not cif.endswith("cif"):
            continue

        # cif_reader.read_cif(
        #     f"{CIF_DIR}{os.sep}{cif}", verbose=True
        # )

        assert cif_reader.read_cif(
            f"{CIF_DIR}{os.sep}{cif}", verbose=True
        ), f"Error reading {cif}"


def test_origin_choice2():
    CIF_DIR = os.path.join(os.path.dirname(__file__), "cifs")
    cif = cif_reader.read_cif(f"{CIF_DIR}{os.sep}452521.cif")
    print(cif.read())
    assert cif.has_origin_choice_2, "wrong origin choice for PCD"


# if __name__ == "__main__":
#     test_cif_readers()
