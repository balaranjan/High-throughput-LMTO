# from htlmto.src.htlmto import cif_reader
from htlmto import cif_reader
import os


def test_cif_readers():
    for cif in os.listdir("cifs"):
        print(cif)
        if not cif.endswith("cif"):
            continue
        assert cif_reader.read_cif(
            f"cifs{os.sep}{cif}"
        ), f"Error reading {cif}"


if __name__ == "__main__":
    test_cif_readers()
