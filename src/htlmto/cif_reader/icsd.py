from .base import CIF_Reader
from collections import defaultdict


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

    def get_no_of_atoms(self):
        num_formula_unit = float(self.get_block("cell_formula_units_Z"))
        atoms_in_formula = sum(self.formula_dict.values())
        return num_formula_unit, num_formula_unit * atoms_in_formula

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
