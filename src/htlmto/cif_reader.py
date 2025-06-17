from abc import ABC, abstractmethod
from collections import defaultdict
import re
import os


def read_cif(filename, verbose=False):

    if not os.path.isfile(filename):
        print(f"File {filename} not found!")

    if verbose:
        print(filename, end="   ... ")

    contents = open(filename, "r").read()

    if "PCD" in contents:
        if verbose:
            print(" is from PCD")
        return PCD_reader(filename, verbose=verbose)
    elif "ICSD" in contents:
        print(" is from ICSD")
        return ICSD_reader(filename, verbose=verbose)
    elif "shelx" in contents:
        print(" is from ShellXL")
        return SXL_reader(filename, verbose=verbose)
    else:
        print(f"File {filename} form unknown source.")
        return


class CIF_Reader(ABC):

    def __init__(self, filename, verbose=False):

        self.filename = filename
        self.read_file()

        self.headers_numeric = ["x", "y", "z", "occupancy", "multiplicity"]

        # data
        self.formula_dict = self.get_formula_dict()
        self.id = self.get_id()
        self.formula = self.get_formula()
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


class PCD_reader(CIF_Reader):

    def __init__(self, filename, verbose=False):
        super().__init__(filename, verbose=verbose)

    def get_block(self, block_name):
        i_start = None
        i_end = None
        for i, line in enumerate(self.lines):
            if block_name in line:
                i_start = i
                break

        for i, line in enumerate(self.lines[i_start + 1 :]):
            if line.startswith("_") or line.startswith("#"):
                i_end = i + i_start + 1
                break

        value = self.lines[i_start:i_end]

        if len(value) == 1:
            return " ".join(value[0].split()[1:])
        else:
            return value

    def get_id(self):
        return self.get_block("database_code_PCD")

    def get_formula_dict(self):
        value = self.get_block("chemical_formula_sum")
        return dict(self._parse_formula(value))

    def get_space_group_no(self):
        return self.get_block("space_group_IT_number")

    def get_structure_type(self):
        return self.get_block("chemical_name_structure_type")

    def get_origin_choice(self):
        hm_alt = self.get_block("space_group_name_H-M_alt")
        return "O2" in hm_alt

    def get_cell(self):

        a = self.get_block("cell_length_a")
        b = self.get_block("cell_length_b")
        c = self.get_block("cell_length_c")

        alpha = self.get_block("cell_angle_alpha")
        beta = self.get_block("cell_angle_beta")
        gamma = self.get_block("cell_angle_gamma")

        return [
            self.get_float(a),
            self.get_float(b),
            self.get_float(c),
            self.get_float(alpha),
            self.get_float(beta),
            self.get_float(gamma),
        ]

    def get_site_data(self):

        i_start, i_end = None, None

        for i, line in enumerate(self.lines):
            if "loop_" in line and "atom_site" in self.lines[i + 1]:
                i_start = i + 1
                break

        for i, line in enumerate(self.lines[i_start + 1 :]):
            if line.startswith("_") and "atom_site" not in line:
                i_end = i + i_start
                break

        values = self.lines[i_start:i_end]
        headers = [
            line.replace("_", " ")
            .strip()
            .replace("\n", "")
            .replace("atom site ", "")
            for line in values
            if "atom_site" in line
        ]

        site_labels = defaultdict(int)
        site_data = []

        header_map = {
            "label": "label",
            "type symbol": "symbol",
            "symmetry multiplicity": "multiplicity",
            "Wyckoff symbol": "Wyckoff_symbol",
            "fract x": "x",
            "fract y": "y",
            "fract z": "z",
            "occupancy": "occupancy",
        }

        for site in values[len(headers) :]:
            site = site.replace("\n", "").strip().split()
            if not len(site):
                continue
            _site_values = {}
            for header, value in zip(headers, site):
                header = header_map[header]

                if header == "label":
                    continue
                if header == "symbol":
                    label = site_labels[value] + 1
                    site_labels[value] += 1
                    _site_values["label"] = f"{value}{label}"

                    _site_values[header] = value
                elif header in self.headers_numeric:
                    _site_values[header] = float(value)
                else:
                    _site_values[header] = value

            site_data.append(_site_values)

        return site_data


class ICSD_reader(CIF_Reader):

    def __init__(self, filename, verbose=False):
        super().__init__(filename, verbose=verbose)

    def get_block(self, block_name):
        i_start = None
        i_end = None
        for i, line in enumerate(self.lines):
            if block_name in line:
                i_start = i
                break

        if i_start is None:
            return

        for i, line in enumerate(self.lines[i_start + 1 :]):
            if (
                line.startswith("_")
                or line.startswith("#")
                or line.startswith("loop_")
            ):
                i_end = i + i_start + 1
                break

        value = self.lines[i_start:i_end]

        if len(value) == 1:
            return " ".join(value[0].split()[1:])
        else:
            return value

    def get_id(self):
        return self.get_block("database_code_ICSD")

    def get_formula_dict(self):
        value = self.get_block("chemical_formula_sum")
        return dict(self._parse_formula(value))

    def get_space_group_no(self):
        value = self.get_block("space_group_IT_number")

        if not value:
            value = self.get_block("symmetry_Int_Tables_number")
            print("ELSE", value)
        return value if value else "NA"

    def get_structure_type(self):
        return self.get_block("chemical_name_structure_type")

    def get_origin_choice(self):
        hm_alt = self.get_block("space_group_name_H-M_alt")
        return "O2" in hm_alt if hm_alt else False

    def get_cell(self):

        a = self.get_block("cell_length_a")
        b = self.get_block("cell_length_b")
        c = self.get_block("cell_length_c")

        alpha = self.get_block("cell_angle_alpha")
        beta = self.get_block("cell_angle_beta")
        gamma = self.get_block("cell_angle_gamma")

        return [
            self.get_float(a),
            self.get_float(b),
            self.get_float(c),
            self.get_float(alpha),
            self.get_float(beta),
            self.get_float(gamma),
        ]

    def get_site_data(self):

        i_start, i_end = None, None

        for i, line in enumerate(self.lines):
            if "loop_" in line and "atom_site" in self.lines[i + 1]:
                i_start = i + 1
                break

        for i, line in enumerate(self.lines[i_start + 1 :]):
            if (
                line.startswith("_")
                or line.startswith("#")
                or line.startswith("loop_")
            ) and "atom_site" not in line:
                i_end = i + i_start
                break

        values = self.lines[i_start:i_end]
        headers = [
            line.replace("_", " ")
            .strip()
            .replace("\n", "")
            .replace("atom site ", "")
            for line in values
            if "atom_site" in line and "atom_site_aniso" not in line
        ]

        site_labels = defaultdict(int)
        site_data = []

        header_map = {
            "label": "label",
            "type symbol": "symbol",
            "symmetry multiplicity": "multiplicity",
            "Wyckoff symbol": "Wyckoff_symbol",
            "fract x": "x",
            "fract y": "y",
            "fract z": "z",
            "occupancy": "occupancy",
        }

        for site in values[len(headers) :]:
            site = site.replace("\n", "").strip().split()
            if not len(site):
                continue
            _site_values = {}
            for header, value in zip(headers, site):
                header = header_map.get(header)

                if not header:
                    continue

                if header == "label":
                    continue
                if header == "symbol":
                    value = value.replace("-", "").replace("+", "")
                    value = list(self._parse_formula(value).keys())[0]
                    label = site_labels[value] + 1
                    site_labels[value] += 1
                    _site_values["label"] = f"{value}{label}"

                    _site_values[header] = value
                elif header in self.headers_numeric:
                    _site_values[header] = self.get_float(value)
                else:
                    _site_values[header] = value

            if len(_site_values):
                site_data.append(_site_values)

        return site_data


class SXL_reader(CIF_Reader):

    def __init__(self, filename, verbose=False):
        super().__init__(filename, verbose=verbose)

    def get_block(self, block_name):
        i_start = None
        i_end = None
        for i, line in enumerate(self.lines):
            if block_name in line:
                i_start = i
                break

        if i_start is None:
            return None

        for i, line in enumerate(self.lines[i_start + 1 :]):
            if (
                line.startswith("_")
                or line.startswith("#")
                or line.startswith("loop_")
            ):
                i_end = i + i_start + 1
                break

        value = self.lines[i_start:i_end]

        if len(value) == 1:
            return " ".join(value[0].split()[1:])
        else:
            return value

    def get_id(self):
        return self.get_block("database_code_ID")

    def get_formula_dict(self):
        value = self.get_block("chemical_formula_sum")
        return dict(self._parse_formula(value))

    def get_space_group_no(self):
        return self.get_block("space_group_IT_number")

    def get_structure_type(self):
        value = self.get_block("chemical_name_structure_type")

        return value if value else "Not assigned"

    def get_origin_choice(self):
        hm_alt = self.get_block("space_group_name_H-M_alt")
        return "O2" in hm_alt

    def get_cell(self):

        a = self.get_block("cell_length_a")
        b = self.get_block("cell_length_b")
        c = self.get_block("cell_length_c")

        alpha = self.get_block("cell_angle_alpha")
        beta = self.get_block("cell_angle_beta")
        gamma = self.get_block("cell_angle_gamma")

        return [
            self.get_float(a),
            self.get_float(b),
            self.get_float(c),
            self.get_float(alpha),
            self.get_float(beta),
            self.get_float(gamma),
        ]

    def get_site_data(self):

        i_start, i_end = None, None

        for i, line in enumerate(self.lines):
            if "loop_" in line and "atom_site" in self.lines[i + 1]:
                i_start = i + 1
                break

        for i, line in enumerate(self.lines[i_start + 1 :]):
            if (
                line.startswith("_")
                or line.startswith("#")
                or line.startswith("loop_")
            ) and "atom_site" not in line:
                i_end = i + i_start
                break

        values = self.lines[i_start:i_end]
        headers = [
            line.replace("_", " ")
            .strip()
            .replace("\n", "")
            .replace("atom site ", "")
            for line in values
            if "atom_site" in line and "atom_site_aniso" not in line
        ]

        site_labels = defaultdict(int)
        site_data = []

        header_map = {
            "label": "label",
            "type symbol": "symbol",
            "symmetry multiplicity": "multiplicity",
            "Wyckoff symbol": "Wyckoff_symbol",
            "fract x": "x",
            "fract y": "y",
            "fract z": "z",
            "occupancy": "occupancy",
        }

        for site in values[len(headers) :]:
            site = site.replace("\n", "").strip().split()
            if not len(site):
                continue
            _site_values = {}
            for header, value in zip(headers, site):
                header = header_map.get(header)

                if not header:
                    continue

                if header == "label":
                    continue
                if header == "symbol":
                    value = value.replace("-", "").replace("+", "")
                    value = list(self._parse_formula(value).keys())[0]
                    label = site_labels[value] + 1
                    site_labels[value] += 1
                    _site_values["label"] = f"{value}{label}"

                    _site_values[header] = value
                elif header in self.headers_numeric:
                    _site_values[header] = self.get_float(value)
                else:
                    _site_values[header] = value

            if len(_site_values):
                site_data.append(_site_values)

        return site_data


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
