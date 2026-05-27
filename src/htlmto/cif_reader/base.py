from abc import ABC, abstractmethod
from collections import defaultdict
import re
import os


class CIF_Reader(ABC):

    def __init__(self, filename, data_source, verbose=False):

        self.filename = filename
        self.read_file()

        self.headers_numeric = ["x", "y", "z", "occupancy", "multiplicity"]

        # data
        self.data_source = data_source
        self.formula_dict = dict(self.get_formula_dict())
        self.formula_dict = dict(
            sorted(
                self.formula_dict.items(),
                key=lambda item: element_data[item[0]][2],
            )
        )

        self.formula = self.get_formula()
        self.id = self.get_id()
        Z, num_atoms = self.get_no_of_atoms()
        self.Z = Z
        self.num_atoms = num_atoms
        self.space_group_number = self.get_space_group_no()
        self.structure_type = self.get_structure_type()
        self.has_origin_choice_2 = self.get_origin_choice()
        self.cell = self.get_cell()

        try:
            self.site_data = self.get_site_data()
        except Exception as e:
            print(f"Error reading site data from {self.filename}")
            print(e)
            return

        self.get_standardized_sites()

        if verbose:
            self.read()

    def read(self):
        print(f"ID       : {self.id}")
        print(f"Source   : {self.data_source}")
        print(f"Formula  : {self.formula_dict}")
        print(f"SG       : {self.space_group_number}")
        print(f"Z        : {self.Z}")
        print(f"Num atoms: {self.num_atoms}")
        print(f"Stype    : {self.structure_type}")
        print(f"has O2   : {self.has_origin_choice_2}")
        print(f"Cell     : {self.cell}")
        print(f"Sites    : {self.site_data}")
        print(f"Defect   : {self.has_defect()}")

    def get_formula(self):

        formula = self.get_formula_dict()

        value = ""
        for k, v in formula.items():
            v = abs(v)
            if v == 1:
                value += k
            elif abs(v - int(v)) > 0.01:
                value += f"{k}{v:.2f}"
            else:
                value += f"{k}{int(v)}"
        return value

    def read_file(self):
        if not os.path.isfile(self.filename):
            print(f"File {self.filename} not found!")

        with open(self.filename, "r") as f:
            self.lines = f.readlines()

    def has_defect(self):
        for site in self.site_data:
            if "occupancy" not in site:
                return None
            if site["occupancy"] != 1.0:
                return True
        return False

    def _parse_formula(self, formula):
        return _parse_formula(formula)

    def get_float(self, s):
        s = s.split("(")[0]
        return float(s)

    def get_standardized_sites(self):
        standardized_size_data = []
        for site in self.site_data:
            for k in ["x", "y", "z"]:
                v = site[k]

                if v > 1.00:
                    v -= 1.00
                elif v < 0.00:
                    v += 1.00

                assert (
                    0.00 <= v <= 1.00
                ), f"Error: Internal coordinate {k} has a value of {site[k]}."
                site[k] = v
            standardized_size_data.append(site)

        return standardized_size_data

    @abstractmethod
    def get_block(self):
        pass

    @abstractmethod
    def get_id(self):
        pass

    @abstractmethod
    def get_formula_dict(self):
        pass

    @abstractmethod
    def get_no_of_atoms(self):
        pass

    @abstractmethod
    def get_space_group_no(self):
        pass

    @abstractmethod
    def get_cell(self):
        pass

    def get_structure_type(self):
        pass

    @abstractmethod
    def get_origin_choice(self):
        pass

    @abstractmethod
    def get_site_data(self):
        pass


def _parse_formula(formula: str, strict: bool = True) -> dict[str, float]:
    """Copied from pymatgen.

    Args:
        formula (str): A string formula, e.g. Fe2O3, Li3Fe2(PO4)3.
        strict (bool): Whether to throw an error if formula string
        is invalid (e.g. empty).
            Defaults to True.

    Returns:
        Composition with that formula.

    Notes:
        In the case of Metallofullerene formula (e.g. Y3N@C80),
        the @ mark will be dropped and passed to parser.
    """
    # Raise error if formula contains special characters
    # or only spaces and/or numbers

    if "'" in formula:
        formula = formula.replace("'", "")

    if strict and re.match(r"[\s\d.*/]*$", formula):
        raise ValueError(f"Invalid formula={formula}")

    # For Metallofullerene like "Y3N@C80"
    formula = formula.replace("@", "")
    # Square brackets are used in formulas to denote coordination
    # complexes (gh-3583)

    formula = formula.replace("[", "(")
    formula = formula.replace("]", ")")

    def get_sym_dict(form: str, factor: float) -> dict[str, float]:
        sym_dict: dict[str, float] = defaultdict(float)
        for match in re.finditer(r"([A-Z][a-z]*)\s*([-*\.e\d]*)", form):
            el = match[1]
            amt = 1.0
            if match[2].strip() != "":
                amt = float(match[2])
            sym_dict[el] += amt * factor
            form = form.replace(match.group(), "", 1)
        if form.strip():
            raise ValueError(f"{form} is an invalid formula!")
        return sym_dict

    match = re.search(r"\(([^\(\)]+)\)\s*([\.e\d]*)", formula)
    while match:
        factor = 1.0
        if match[2] != "":
            factor = float(match[2])
        unit_sym_dict = get_sym_dict(match[1], factor)
        expanded_sym = "".join(
            f"{el}{amt}" for el, amt in unit_sym_dict.items()
        )
        expanded_formula = formula.replace(match.group(), expanded_sym, 1)
        formula = expanded_formula
        match = re.search(r"\(([^\(\)]+)\)\s*([\.e\d]*)", formula)
    return get_sym_dict(formula, 1)


element_data = {
    "Vac": [0, 0.0, 0],
    "H": [1, 0.25, 92],
    "He": [2, None, 98],
    "Li": [3, 1.45, 1],
    "Be": [4, 1.05, 67],
    "B": [5, 0.85, 72],
    "C": [6, 0.7, 77],
    "N": [7, 0.65, 82],
    "O": [8, 0.6, 87],
    "F": [9, 0.5, 93],
    "Ne": [10, None, 99],
    "Na": [11, 1.8, 2],
    "Mg": [12, 1.5, 68],
    "Al": [13, 1.25, 73],
    "Si": [14, 1.1, 78],
    "P": [15, 1.0, 83],
    "S": [16, 1.0, 88],
    "Cl": [17, 1.0, 94],
    "Ar": [18, 0.71, 100],
    "K": [19, 2.2, 3],
    "Ca": [20, 1.8, 7],
    "Sc": [21, 1.6, 11],
    "Ti": [22, 1.4, 43],
    "V": [23, 1.35, 46],
    "Cr": [24, 1.4, 49],
    "Mn": [25, 1.4, 52],
    "Fe": [26, 1.4, 55],
    "Co": [27, 1.35, 58],
    "Ni": [28, 1.35, 61],
    "Cu": [29, 1.35, 64],
    "Zn": [30, 1.35, 69],
    "Ga": [31, 1.3, 74],
    "Ge": [32, 1.25, 79],
    "As": [33, 1.15, 84],
    "Se": [34, 1.15, 89],
    "Br": [35, 1.15, 95],
    "Kr": [36, None, 101],
    "Rb": [37, 2.35, 4],
    "Sr": [38, 2.0, 8],
    "Y": [39, 1.8, 12],
    "Zr": [40, 1.55, 44],
    "Nb": [41, 1.45, 47],
    "Mo": [42, 1.45, 50],
    "Tc": [43, 1.35, 53],
    "Ru": [44, 1.3, 56],
    "Rh": [45, 1.35, 59],
    "Pd": [46, 1.4, 62],
    "Ag": [47, 1.6, 65],
    "Cd": [48, 1.55, 70],
    "In": [49, 1.55, 75],
    "Sn": [50, 1.45, 80],
    "Sb": [51, 1.45, 85],
    "Te": [52, 1.4, 90],
    "I": [53, 1.4, 96],
    "Xe": [54, None, 102],
    "Cs": [55, 2.6, 5],
    "Ba": [56, 2.15, 9],
    "La": [57, 1.95, 13],
    "Ce": [58, 1.85, 15],
    "Pr": [59, 1.85, 17],
    "Nd": [60, 1.85, 19],
    "Pm": [61, 1.85, 21],
    "Sm": [62, 1.85, 23],
    "Eu": [63, 1.85, 25],
    "Gd": [64, 1.8, 27],
    "Tb": [65, 1.75, 29],
    "Dy": [66, 1.75, 31],
    "Ho": [67, 1.75, 33],
    "Er": [68, 1.75, 35],
    "Tm": [69, 1.75, 37],
    "Yb": [70, 1.75, 39],
    "Lu": [71, 1.75, 41],
    "Hf": [72, 1.55, 45],
    "Ta": [73, 1.45, 48],
    "W": [74, 1.35, 51],
    "Re": [75, 1.35, 54],
    "Os": [76, 1.3, 57],
    "Ir": [77, 1.35, 60],
    "Pt": [78, 1.35, 63],
    "Au": [79, 1.35, 66],
    "Hg": [80, 1.5, 71],
    "Tl": [81, 1.9, 76],
    "Pb": [82, 1.8, 81],
    "Bi": [83, 1.6, 86],
    "Po": [84, 1.9, 91],
    "At": [85, None, 97],
    "Rn": [86, None, 103],
    "Fr": [87, None, 6],
    "Ra": [88, 2.15, 10],
    "Ac": [89, 1.95, 14],
    "Th": [90, 1.8, 16],
    "Pa": [91, 1.8, 18],
    "U": [92, 1.75, 20],
    "Np": [93, 1.75, 22],
    "Pu": [94, 1.75, 24],
    "Am": [95, 1.75, 26],
    "Cm": [96, None, 28],
    "Bk": [97, None, 30],
    "Cf": [98, None, 32],
    "Es": [99, None, 34],
}
