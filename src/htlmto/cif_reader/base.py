from abc import ABC, abstractmethod
from collections import defaultdict
import re
import os


class CIF_Reader(ABC):

    def __init__(self, filename, verbose=False):

        self.filename = filename
        self.read_file()

        self.headers_numeric = ["x", "y", "z", "occupancy", "multiplicity"]

        # data
        self.formula_dict = self.get_formula_dict()
        self.id = self.get_id()
        self.formula = self.get_formula()
        Z, num_atoms = self.get_no_of_atoms()
        self.Z = Z
        self.num_atoms = num_atoms
        self.space_group_number = self.get_space_group_no()
        self.structure_type = self.get_structure_type()
        self.has_origin_choice_2 = self.get_origin_choice()
        self.cell = self.get_cell()
        self.site_data = self.get_site_data()
        self.get_standardized_sites()

        if verbose:
            self.read()

    def read(self):
        print(f"ID       : {self.id}")
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
            if site["occupancy"] != 1.0:
                return True
        return False

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

    def _parse_formula(
        self, formula: str, strict: bool = True
    ) -> dict[str, float]:
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
