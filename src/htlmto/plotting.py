from matplotlib import pyplot as plt
from matplotlib.ticker import MaxNLocator
import pandas as pd
import os
import numpy as np
import glob
import re

# Define transition metals
transition_metals = [
    "Sc",
    "Ti",
    "V",
    "Cr",
    "Mn",
    "Fe",
    "Co",
    "Ni",
    "Cu",
    "Zn",
    "Y",
    "Zr",
    "Nb",
    "Mo",
    "Tc",
    "Ru",
    "Rh",
    "Pd",
    "Ag",
    "Cd",
    "Hf",
    "Ta",
    "W",
    "Re",
    "Os",
    "Ir",
    "Pt",
    "Au",
    "Hg",
]

# Define the Mendeleev number mapping
mendeleev_numbers = {
    "H": 92,
    "He": 98,
    "Li": 1,
    "Be": 67,
    "B": 72,
    "C": 77,
    "N": 82,
    "O": 87,
    "F": 93,
    "Ne": 99,
    "Na": 2,
    "Mg": 68,
    "Al": 73,
    "Si": 78,
    "P": 83,
    "S": 88,
    "Cl": 94,
    "Ar": 100,
    "K": 3,
    "Ca": 7,
    "Sc": 11,
    "Ti": 43,
    "V": 46,
    "Cr": 49,
    "Mn": 52,
    "Fe": 55,
    "Co": 58,
    "Ni": 61,
    "Cu": 64,
    "Zn": 69,
    "Ga": 74,
    "Ge": 79,
    "As": 84,
    "Se": 89,
    "Br": 95,
    "Kr": 101,
    "Rb": 4,
    "Sr": 8,
    "Y": 12,
    "Zr": 44,
    "Nb": 47,
    "Mo": 50,
    "Tc": 53,
    "Ru": 56,
    "Rh": 59,
    "Pd": 62,
    "Ag": 65,
    "Cd": 70,
    "In": 75,
    "Sn": 80,
    "Sb": 85,
    "Te": 90,
    "I": 96,
    "Xe": 102,
    "Cs": 5,
    "Ba": 9,
    "La": 13,
    "Ce": 15,
    "Pr": 17,
    "Nd": 19,
    "Pm": 21,
    "Sm": 23,
    "Eu": 25,
    "Gd": 27,
    "Tb": 29,
    "Dy": 31,
    "Ho": 33,
    "Er": 35,
    "Tm": 37,
    "Yb": 39,
    "Lu": 41,
    "Hf": 45,
    "Ta": 48,
    "W": 51,
    "Re": 54,
    "Os": 57,
    "Ir": 60,
    "Pt": 63,
    "Au": 66,
    "Hg": 71,
    "Tl": 76,
    "Pb": 81,
    "Bi": 86,
    "Po": 91,
    "At": 97,
    "Rn": 103,
    "Fr": 6,
    "Ra": 10,
    "Ac": 14,
    "Th": 16,
    "Pa": 18,
    "U": 20,
    "Np": 22,
    "Pu": 24,
    "Am": 26,
    "Cm": 28,
    "Bk": 30,
    "Cf": 32,
    "Es": 34,
    "Fm": 36,
    "Md": 38,
    "No": 40,
    "Lr": 42,
}


def get_dos_files(folder_path):
    return glob.glob(os.path.join(folder_path, "DOS-*.csv"))


def extract_elements(label):
    return re.findall(r"[A-Z][a-z]?", label)


def get_element_color(label, elements):
    label_lower = label.lower()
    if label_lower == "total":
        return "black", "-"
    if label_lower == "e":
        return "darkgrey", "--"

    elements = [e for e in elements if e.lower() != "e"]
    element_count = len(elements)

    if element_count == 1:
        return "blue", "-"
    if element_count == 2:
        return ("blue", "-") if label == elements[0] else ("red", "-")
    if element_count == 3:
        if label == elements[0]:
            return "blue", "-"
        if label in transition_metals:
            return "grey", "-"
        if label == elements[2]:
            return "red", "-"
    if element_count == 4:
        if label == elements[0]:
            return "blue", "-"
        if label in transition_metals:
            return "grey", "-"
        if label == elements[2]:
            return "green", "-"
        if label == elements[3]:
            return "red", "-"
        return "orange", "-"
    if element_count == 5:
        if label == elements[0]:
            return "blue", "-"
        if label in transition_metals:
            return "grey", "-"
        if label == elements[2]:
            return "green", "-"
        if label == elements[3]:
            return "red", "-"
        if label == elements[4]:
            return "orange", "-"
    return "pink", "-"


def plot_dos(calc_dir):
    def plot(include_e):
        dos_files = get_dos_files(calc_dir)
        if not dos_files:
            print("No DOS files found in the directory.")
            return

        fig, ax = plt.subplots(figsize=(8, 15))
        all_y_values = []
        all_elements = []

        file_data = []
        for filename in dos_files:
            data = pd.read_csv(filename, sep=",", skiprows=1, header=None)
            data.columns = ["Energy", "DOS", "Intg_DOS"]
            data["Energy"] = pd.to_numeric(data["Energy"], errors="coerce")
            data["DOS"] = pd.to_numeric(data["DOS"], errors="coerce")
            data["Intg_DOS"] = pd.to_numeric(data["Intg_DOS"], errors="coerce")

            x = data["Energy"].values
            y = data["DOS"].values
            all_y_values.extend(y[(x >= -6) & (x <= 2)])

            label = os.path.splitext(os.path.basename(filename))[0].replace(
                "DOS-", ""
            )
            elements = extract_elements(label)
            all_elements.extend(elements)
            file_data.append((x, y, label))

        all_elements = list(set(all_elements))
        sorted_elements = sorted(
            all_elements, key=lambda e: mendeleev_numbers.get(e, float("inf"))
        )
        if all_y_values:
            max_y = max(all_y_values)
            buffer = 0.1 * max_y
        else:
            print(
                "No DOS data in the range -6 <= Energy <= 2. Skipping plot. "
                "This may be due to an empty or incorrectly formatted DOS.csv "
                "file."
            )
            return

        def sort_key(item):
            label = item[2]
            if label.lower() == "total":
                return (-1, 0)
            if label == "E" and not include_e:
                return (float("inf"), float("inf"))
            return (mendeleev_numbers.get(label, float("inf")), label)

        file_data.sort(key=sort_key)

        for x, y, label in file_data:
            if label.lower() == "e" and not include_e:
                continue
            elements = extract_elements(label)
            color, linestyle = get_element_color(label, sorted_elements)
            zorder_value = 10 if label.lower() == "total" else 0
            ax.plot(
                y,
                x,
                label=label,
                color=color,
                linewidth=5,
                linestyle=linestyle,
                zorder=zorder_value,
            )

        ax.set_ylim(-6, 2)
        ax.set_xlim(0, max_y + buffer)
        ax.axhline(0, color="black", linestyle="--", linewidth=3)
        ax.tick_params(
            axis="x", which="both", bottom=False, top=False, labelbottom=False
        )
        ax.tick_params(axis="y", labelsize=35, width=3, length=10)

        for spine in ax.spines.values():
            spine.set_linewidth(2.5)

        ax.set_ylabel("energy (eV)", fontsize=35)
        ax.yaxis.set_major_locator(MaxNLocator(integer=True))

        legend = ax.legend(
            frameon=True,
            fontsize=30,
            loc="lower right",
            handlelength=0.75,
            columnspacing=0.1,
            facecolor="white",
        )
        legend.set_zorder(99)

        folder_name = os.path.basename(calc_dir).split("-")[0]
        folder_name_cleaned = re.sub(
            r"(?<=[A-Za-z])1(?=[A-ZaZ])", "", folder_name
        )
        folder_name_cleaned = re.sub(r"(?<!\d)1$", "", folder_name_cleaned)
        ax.set_title("DOS", fontsize=35, pad=20)

        x_position = ax.get_xlim()[1]
        ax.annotate(
            r"$E_{\mathrm{F}}$",
            xy=(x_position, 0),
            xytext=(10, 0),
            textcoords="offset points",
            fontsize=35,
            va="center",
            ha="left",
            color="black",
        )

        plt.tight_layout()

        if include_e:
            output_filename = os.path.join(calc_dir, "DOS.png")
        else:
            output_filename = os.path.join(calc_dir, "DOS_without_E.png")

        plt.savefig(output_filename, dpi=300)
        print(f"Plot saved to: {output_filename}")
        plt.close(fig)

    plot(include_e=True)
    plot(include_e=False)


def plot_band_structure(calc_dir):
    try:
        points_file = os.path.join(calc_dir, "band_structure_points.csv")
        band_file = os.path.join(calc_dir, "band_structure.csv")

        if not (os.path.exists(points_file) and os.path.exists(band_file)):
            print("Band structure files not found in:", calc_dir)
            return

        points_data = pd.read_csv(points_file, delimiter=",")
        points_data.columns = points_data.columns.str.strip()

        band_data = pd.read_csv(band_file, delimiter=",")
        band_data.columns = band_data.columns.str.strip()

        ticks = points_data["values"].tolist()
        labels = points_data["point"].tolist()

        folder_name = os.path.basename(calc_dir).split("-")[0]
        folder_name_cleaned = re.sub(
            r"(?<=[A-Za-z])1(?=[A-ZaZ])", "", folder_name
        )
        folder_name_cleaned = re.sub(r"(?<!\d)1$", "", folder_name_cleaned)

        x_min = min(ticks)
        x_max = max(ticks)

        plt.rcParams.update(
            {
                "font.size": 20,
                "axes.labelsize": 20,
                "xtick.labelsize": 20,
                "ytick.labelsize": 20,
            }
        )

        fig, ax = plt.subplots(figsize=(8, 8))

        ax.plot(
            band_data["k"],
            band_data["Energy (eV)"],
            color="blue",
            linewidth=0.85,
            zorder=1,
        )

        ax.set_xticks(ticks)
        ax.set_xticklabels(labels)

        ax.set_xlabel(r"$\mathit{k}$-points")
        ax.set_ylabel("energy (eV)", labelpad=-5)

        ax.set_title("Band Structure", fontsize=20, pad=10)
        ax.set_xlim(x_min, x_max)
        ax.set_ylim(-2.5, 2.5)
        ax.axhline(0, color="black", linewidth=1.5, linestyle="--", zorder=1)

        x_position = ax.get_xlim()[1]
        ax.annotate(
            r"$E_{\mathrm{F}}$",
            xy=(x_position, 0),
            xytext=(10, 0),
            textcoords="offset points",
            fontsize=20,
            va="center",
            ha="left",
            color="black",
            zorder=1,
        )

        ax.grid(axis="x", color="black", alpha=1, zorder=2, linewidth=1.5)

        for spine in ax.spines.values():
            spine.set_linewidth(1.5)

        ax.tick_params(axis="both", width=1.5, length=8)

        plt.tight_layout()

        save_path = os.path.join(calc_dir, "Bandstructure.png")
        plt.savefig(save_path, dpi=300, bbox_inches="tight")

        print(f"Plot saved to: {save_path}")
        plt.close(fig)
    except Exception as e:
        print("Error while plotting band structure.")
        print(e)


def plot_cohps(calc_dir):
    try:
        cohp_files = [
            os.path.join(calc_dir, f)
            for f in os.listdir(calc_dir)
            if f.lower().startswith("data.cohp")
        ]
        if not cohp_files:
            print("COHP file not found in:", calc_dir)
            return

        plt.rcParams.update(
            {
                "font.size": 35,
                "axes.labelsize": 35,
                "xtick.labelsize": 35,
                "ytick.labelsize": 35,
            }
        )

        fig, ax = plt.subplots(figsize=(8, 15))
        color_cycle = ["red", "orange", "green", "blue", "purple", "pink"]

        all_x_values = []
        for idx, cohp_file in enumerate(cohp_files):
            cohp_data = pd.read_csv(
                cohp_file,
                sep=r"\s+",
                header=None,
                names=["energy", "cohp", "int_cohp"],
            )
            pair = (
                os.path.basename(cohp_file)
                .replace("data.cohp", "")
                .replace(".csv", "")
                .replace("data.cohp", "")
                .replace("_", "-")
            )
            if not pair or pair == "data.cohp":
                pair = "Total"
            if pair.startswith("-"):
                pair = pair[1:]  # Remove leading dash if present

            color = color_cycle[idx % len(color_cycle)]
            # Plot COHP (solid, with legend)
            ax.plot(
                cohp_data["cohp"],
                cohp_data["energy"],
                label=pair,
                color=color,
                linewidth=5,
            )
            # Plot ICOHP (dashed, no legend)
            ax.plot(
                cohp_data["int_cohp"],
                cohp_data["energy"],
                color=color,
                linewidth=3,
                linestyle="--",
                zorder=0,
            )
            all_x_values.extend(np.abs(cohp_data["cohp"].values))
            all_x_values.extend(np.abs(cohp_data["int_cohp"].values))

        # Set axis limits and style to match DOS plot
        ax.set_ylim(-6, 2)
        max_x = max(all_x_values) if all_x_values else 1
        buffer = 0.1 * max_x
        ax.set_xlim(-(max_x + buffer), max_x + buffer)
        ax.axhline(0, color="black", linestyle="--", linewidth=3)
        ax.axvline(0, color="black", linestyle="--", linewidth=3)
        ax.tick_params(
            axis="x", which="both", bottom=False, top=False, labelbottom=False
        )
        ax.tick_params(axis="y", labelsize=35, width=3, length=10)

        for spine in ax.spines.values():
            spine.set_linewidth(2.5)

        ax.set_ylabel("energy (eV)", fontsize=35)
        ax.set_xlabel("-COHP", fontsize=35)
        ax.yaxis.set_major_locator(MaxNLocator(integer=True))

        legend = ax.legend(
            frameon=True,
            fontsize=30,
            loc="lower left",
            handlelength=0.75,
            columnspacing=0.1,
            facecolor="white",
        )
        legend.set_zorder(99)

        folder_name = os.path.basename(calc_dir).split("-")[0]
        folder_name_cleaned = re.sub(
            r"(?<=[A-Za-z])1(?=[A-ZaZ])", "", folder_name
        )
        folder_name_cleaned = re.sub(r"(?<!\d)1$", "", folder_name_cleaned)

        # ax.set_title(folder_name_subscripted + ' COHP', fontsize=35, pad=20)
        ax.set_title("COHP", fontsize=35, pad=20)

        x_position = ax.get_xlim()[1]
        ax.annotate(
            r"$E_{\mathrm{F}}$",
            xy=(x_position, 0),
            xytext=(10, 0),
            textcoords="offset points",
            fontsize=35,
            va="center",
            ha="left",
            color="black",
        )

        plt.tight_layout()
        save_path = os.path.join(calc_dir, "COHP.png")
        plt.savefig(save_path, dpi=300, bbox_inches="tight")
        print(f"Plot saved to: {save_path}")
        plt.close(fig)
    except Exception as e:
        print("Error while plotting COHP.")
        print(e)
