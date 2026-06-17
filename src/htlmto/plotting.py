from matplotlib import pyplot as plt
from matplotlib.ticker import MaxNLocator
from collections import defaultdict
from math import ceil
import pandas as pd
import os
import numpy as np
import argparse
import glob
import sys
import re
import traceback
import matplotlib.colors as mcolors

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
        # Check for specific combination: Cu, Cd, Sn, S
        element_set = set(elements)
        if element_set == {"Cu", "Cd", "Sn", "S"}:
            if label == "Cu":
                return "green", "-"
            if label == "Cd":
                return "blue", "-"
            if label == "Sn":
                return "red", "-"
            if label == "S":
                return "yellow", "-"
        # Default behavior for element_count == 4
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


def get_nonzero_integer_ticks(x_min, x_max, n_ticks=4):
    if n_ticks < 1 or x_max <= x_min:
        return []

    nice_intervals = [
        0.1,
        0.2,
        0.25,
        1 / 3,
        0.5,
        2 / 3,
        1,
        2,
        3,
        4,
        5,
        6,
        10,
        12,
        15,
        16,
        20,
        25,
        30,
        40,
        50,
        100,
        150,
        200,
        250,
        500,
        1000,
    ]

    ideal_interval = (x_max - x_min) / n_ticks
    interval = min(nice_intervals, key=lambda x: abs(x - ideal_interval))

    ticks = []

    if x_min < 0:
        tick = -interval
        while tick >= x_min:
            if tick != 0:
                ticks.append(tick)
            tick -= interval

    tick = interval
    while tick <= x_max:
        if tick != 0:
            ticks.append(tick)
        tick += interval

    ticks = sorted(set(ticks))

    if len(ticks) <= n_ticks:
        return ticks

    step = len(ticks) // n_ticks
    selected = [ticks[i * step] for i in range(n_ticks)]

    # if ticks[-1] not in selected:
    #     selected[-1] = ticks[-1]

    return sorted(set(selected))


def rescale_dos():

    parser = argparse.ArgumentParser(
        description="Rescale and plot the DOS \
            ithin a specified energy range."
    )

    parser.add_argument(
        "calc_dir",
        type=str,
        help="Path to the calculation directory",
    )

    parser.add_argument(
        "emin",
        type=float,
        default=-6,
        help="Minimum energy limit for the y-axis (in eV). Default: -3.9",
    )

    parser.add_argument(
        "emax",
        type=float,
        default=2,
        help="Maximum energy limit for the y-axis (in eV). Default: 3.9",
    )

    args = parser.parse_args()

    calc_dir, emin, emax = args.calc_dir, args.emin, args.emax

    if calc_dir == ".":
        calc_dir = os.getcwd()

    if "plots" in os.listdir(calc_dir):
        plot_dos(f"{calc_dir}{os.sep}csv", emin, emax)
    else:
        for _calc in os.listdir(calc_dir):
            _calc = os.path.join(calc_dir, _calc)
            if not os.path.isdir(_calc):
                continue
            plot_dos(os.path.join(_calc, "csv"), emin, emax)


def plot_dos(calc_dir, emin=-6, emax=2):
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
            all_y_values.extend(y[(x >= emin) & (x <= emax)])

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
            buffer = 0.05 * max_y
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

        ax.set_ylim(emin, emax)
        ax.set_xlim(0, max_y + buffer)
        dos_ticks = get_nonzero_integer_ticks(0, max_y + buffer, n_ticks=4)
        if dos_ticks:
            ax.set_xticks(dos_ticks)
        ax.axhline(0, color="black", linestyle="--", linewidth=3)
        ax.tick_params(axis="x", labelsize=35, width=3, length=10)
        ax.tick_params(axis="y", labelsize=35, width=3, length=10)

        for spine in ax.spines.values():
            spine.set_linewidth(2.5)

        ax.set_ylabel("energy (eV)", fontsize=35)
        ax.yaxis.set_major_locator(MaxNLocator(integer=True))

        legend = ax.legend(
            frameon=False,
            fontsize=30,
            loc="best",
            handlelength=0.75,
            columnspacing=0.1,
        )
        legend.set_zorder(99)

        folder_name = os.path.basename(calc_dir).split("-")[0]
        folder_name_cleaned = re.sub(
            r"(?<=[A-Za-z])1(?=[A-ZaZ])", "", folder_name
        )
        folder_name_cleaned = re.sub(r"(?<!\d)1$", "", folder_name_cleaned)
        ax.set_title("DOS", fontsize=35, pad=20)
        ax.set_xlabel("states", fontsize=35)

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

        # save_path = None
        output_filename = f"DOS{'' if include_e else '_without_E'}.png"
        if "csv" in calc_dir.split(os.sep):
            save_path = os.sep.join(calc_dir.split(os.sep)[:-1])
            save_path = os.path.join(save_path, "plots")
            output_filename = os.path.join(save_path, output_filename)

        plt.savefig(output_filename, dpi=300)
        print(f"\t\tSaving  {output_filename}")
        plt.close(fig)

    plot(include_e=True)
    plot(include_e=False)


def rescale_band_str():

    parser = argparse.ArgumentParser(
        description="Rescale and plot the band \
            structure within a specified energy range."
    )

    parser.add_argument(
        "calc_dir",
        type=str,
        help="Path to the calculation directory containing \
            'band_structure.csv' and 'band_structure_points.csv' \
                (and optionally a 'plots' subdirectory).",
    )

    parser.add_argument(
        "emin",
        type=float,
        default=-3.9,
        help="Minimum energy limit for the y-axis (in eV). Default: -3.9",
    )

    parser.add_argument(
        "emax",
        type=float,
        default=3.9,
        help="Maximum energy limit for the y-axis (in eV). Default: 3.9",
    )

    args = parser.parse_args()

    calc_dir, emin, emax = args.calc_dir, args.emin, args.emax
    if calc_dir == ".":
        calc_dir = os.getcwd()
    if "plots" in os.listdir(calc_dir):
        plot_band_structure(f"{calc_dir}{os.sep}csv", emin, emax)
    else:
        for _calc in os.listdir(calc_dir):
            _calc = os.path.join(calc_dir, _calc)
            if not os.path.isdir(_calc):
                continue
            plot_band_structure(os.path.join(_calc, "csv"), emin, emax)


def plot_band_structure(calc_dir, emin=-3.9, emax=3.9):
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
        ax.set_ylim(emin, emax)
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

        if "csv" in calc_dir.split(os.sep):
            calc_dir = os.sep.join(calc_dir.split(os.sep)[:-1])
            save_path = os.path.join(calc_dir, "plots")
            save_path = f"{save_path}{os.sep}Bandstructure.png"

        else:
            save_path = "Bandstructure.png"
        plt.savefig(save_path, dpi=300, bbox_inches="tight")

        print(f"\t\tSaving {save_path}")
        plt.close(fig)
    except Exception as e:
        print("Error while plotting band structure.")
        print(e)


def interactive_cohp_plot():

    calc_dir = sys.argv[1]
    if os.path.isdir(os.path.join(calc_dir, "calculation")):
        calc_dir = os.path.join(calc_dir, "calculation")
        plots_dir = os.path.join(calc_dir, "plots")
    else:
        plots_dir = os.path.join(calc_dir, "../plots")
        os.makedirs(plots_dir, exist_ok=True)

    try:
        cohp_files = [
            os.path.join(calc_dir, f)
            for f in os.listdir(calc_dir)
            if f.lower().startswith("data.cohp")
        ]
        if not cohp_files:
            print("COHP file not found in:", os.getcwd())
            return

        pairs = {}
        for cohp in cohp_files:
            pair = "-".join(sorted(cohp.split("_")[-2:]))
            pairs[pair] = cohp

        keys = list(pairs.keys())
        colors = [
            "red",
            "orange",
            "green",
            "blue",
            "purple",
            "pink",
            "aqua",
            "black",
            "aquamarine",
            "beige",
            "blueviolet",
            "brown",
            "chocolate",
            "coral",
            "cyan",
            "darkblue",
            "darkcyan",
            "darkgoldenrod",
            "darkgray",
            "darkgreen",
            "darkgrey",
            "darkmagenta",
            "darkorange",
            "darkorchid",
            "darkred",
            "darksalmon",
            "darkseagreen",
            "darkviolet",
            "fuchsia",
            "gold",
            "greenyellow",
            "indigo",
            "ivory",
            "khaki",
            "lavender",
            "lime",
            "magenta",
            "maroon",
            "navy",
            "oldlace",
            "olive",
            "orchid",
            "peru",
            "pink",
            "plum",
            "silver",
            "tan",
            "teal",
            "tomato",
            "violet",
            "yellow",
            "yellowgreen",
        ]

        print("COHP for the following element pairs are available.")
        print(f"{0:<2} {'all'}")
        for i, pair in enumerate(keys, 1):
            print(f"{i:<2} {pair} {colors[(i-1) % len(colors)]}")

        selected = False
        selection_map = {}  # Stores {pair_key: color}

        while not selected:
            raw_input_str = input(
                "Select the numbers corresponding to the pairs to plot \
                    (separated by comma), or 0 for all. "
                "\nTo change color, enter 'index-color' (e.g., '1-red'): "
            )

            if raw_input_str.strip() == "0":
                # Plot all pairs with default cyclic colors
                selection_map = {
                    k: colors[(i - 1) % len(colors)]
                    for i, k in enumerate(keys, 1)
                }
                selected = True
                continue

            parts = [p.strip() for p in raw_input_str.split(",") if p.strip()]
            temp_selection = {}
            is_valid = True

            for part in parts:
                # Check if color is specified (format: "N-color")
                if "-" in part:
                    idx_str, color_name = part.rsplit("-", 1)
                else:
                    idx_str = part
                    color_name = None

                try:
                    idx = int(idx_str)
                except ValueError:
                    print(f"Error: '{idx_str}' is not a valid number.")
                    is_valid = False
                    break

                if idx < 1 or idx > len(keys):
                    print(
                        f"Error: Index {idx} is out of range (1-{len(keys)})."
                    )
                    is_valid = False
                    break

                # Validate color if provided
                if color_name is not None:
                    if color_name not in colors and not mcolors.is_color_like(
                        color_name
                    ):
                        print(
                            f"Error: '{color_name}' \
                                is not a recognized color. \
                            Use one from the list or a valid \
                                CSS/matplotlib color."
                        )
                        is_valid = False
                        break
                    temp_selection[keys[idx - 1]] = color_name
                else:
                    # Default color assignment based on index
                    temp_selection[keys[idx - 1]] = colors[
                        (idx - 1) % len(colors)
                    ]

            if is_valid and temp_selection:
                selection_map = temp_selection
                selected = True
            elif not is_valid:
                print("Invalid input. Please try again.\n")
            else:
                print("No valid pairs selected. Please try again.\n")

        # Prepare data for plotting function
        final_pairs_to_plot = []
        final_colors = []

        if 0 in [
            k.split("-")[0] for k in selection_map.keys()
        ]:  # Check if 'all' was logic handled earlier,
            # but strictly we map keys now
            pass  # Already handled above

        for pair_key, color in selection_map.items():
            final_pairs_to_plot.append(pair_key)
            final_colors.append(color)

        print(f"Selected pairs: {final_pairs_to_plot}")
        print(f"Assigned colors: {final_colors}")

        selected_files = [pairs[k] for k in final_pairs_to_plot if k in pairs]

        if not selected_files:
            print("No files selected.")
            return

        lw, lwi = 5.0, 3.0
        print(f"The linewidth for COHP is {lw} and for iCOHP it is {lwi}.")
        valid = False
        while not valid:
            try:
                value = input(
                    "Enter new widths for \
                        COHP, iCOHP (seaprated \
                            by comma, \
                                enter 0 to accept defaults): "
                )
                if value.strip() == "0":
                    valid = True
                if "," in value:
                    lw, lwi = [float(v.strip()) for v in value.split(",")]
                    valid = True
                else:
                    lw = float(value.strip())
                    lwi = 0
                    valid = True
            except Exception as e:
                print("Invalid input. Please try again.\n")
                print(e)

        emin, emax = -6.0, 2.0
        print(f"The emin, emax for COHP are {emin} and {emax}.")
        valid = False
        while not valid:
            try:
                value = input(
                    "Enter new emin, emax \
                        (seaprated by comma, enter 0 to accept defaults): "
                )
                if value.strip() == "0":
                    valid = True
                elif "," in valid:
                    emin, emax = [float(v.strip()) for v in valid.split(",")]
                    valid = True
                else:
                    print("Invalid input. Please try again.\n")
            except Exception as e:
                print("Invalid input. Please try again.\n")
                print(e)

        plot_cohps(
            calc_dir,
            plot_icohp=lwi,
            cohp_files=[pairs[pair] for pair in final_pairs_to_plot],
            colors=final_colors,
            lw=lw,
            lwi=lwi,
            emin=emin,
            emax=emax,
        )

    except Exception as e:
        print("Error in interactive selection.")
        print(e)
        print(traceback.format_exc())

    return


def plot_cohps(
    calc_dir,
    plot_icohp=True,
    cohp_files=None,
    colors=None,
    lw=5.0,
    lwi=3.0,
    emin=-6.0,
    emax=2.0,
):
    try:
        if cohp_files is None:
            cohp_files = [
                os.path.join(calc_dir, f)
                for f in os.listdir(calc_dir)
                if f.lower().startswith("data.cohp")
            ]
            if not cohp_files:
                print("COHP file not found in:", os.getcwd())
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
        color_cycle = [
            "red",
            "orange",
            "green",
            "blue",
            "purple",
            "pink",
            "aqua",
            "black",
            "aquamarine",
            "beige",
            "blueviolet",
            "brown",
            "chocolate",
            "coral",
            "cyan",
            "darkblue",
            "darkcyan",
            "darkgoldenrod",
            "darkgray",
            "darkgreen",
            "darkgrey",
            "darkmagenta",
            "darkorange",
            "darkorchid",
            "darkred",
            "darksalmon",
            "darkseagreen",
            "darkviolet",
            "fuchsia",
            "gold",
            "greenyellow",
            "indigo",
            "ivory",
            "khaki",
            "lavender",
            "lime",
            "magenta",
            "maroon",
            "navy",
            "oldlace",
            "olive",
            "orchid",
            "peru",
            "pink",
            "plum",
            "silver",
            "tan",
            "teal",
            "tomato",
            "violet",
            "yellow",
            "yellowgreen",
        ]
        if colors:
            color_cycle = colors

        all_x_values = []
        for idx, cohp_file in enumerate(cohp_files):
            cohp_data = pd.read_csv(
                cohp_file,
                sep=r"\s+",
                header=None,
                names=["energy", "cohp", "int_cohp"],
            )
            cohp_data = cohp_data[
                (cohp_data["energy"] >= emin) & (cohp_data["energy"] <= emax)
            ]

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
                linewidth=lw,
            )
            # energies = cohp_data["energy"].values
            all_x_values.append(np.abs(cohp_data["cohp"].values).max())

            # Plot ICOHP (dashed, no legend)
            if plot_icohp and lwi:
                ax.plot(
                    cohp_data["int_cohp"],
                    cohp_data["energy"],
                    color=color,
                    linewidth=lwi,
                    linestyle="--",
                    zorder=0,
                )
                all_x_values.append(np.abs(cohp_data["int_cohp"].values).max())

        # Set axis limits and style to match DOS plot
        ax.set_ylim(emin, emax)
        max_x = max(all_x_values) if all_x_values else 1
        if plot_icohp:
            buffer = 0.05 * max_x
        else:
            buffer = 0.1 * max_x
        ax.set_xlim(-(max_x + buffer), max_x + buffer)
        cohp_ticks = get_nonzero_integer_ticks(
            -(max_x + buffer), max_x + buffer, n_ticks=4
        )
        if cohp_ticks:
            ax.set_xticks(cohp_ticks)
        ax.axhline(0, color="black", linestyle="--", linewidth=lwi)
        ax.axvline(0, color="black", linestyle="--", linewidth=lwi)
        ax.tick_params(axis="x", labelsize=35, width=3, length=10)
        ax.tick_params(axis="y", labelsize=35, width=3, length=10)

        for spine in ax.spines.values():
            spine.set_linewidth(2.5)

        ax.set_ylabel("energy (eV)", fontsize=35)
        ax.set_xlabel("-COHP", fontsize=35)
        ax.yaxis.set_major_locator(MaxNLocator(integer=True))

        legend = ax.legend(
            frameon=False,
            fontsize=30,
            loc="lower left",
            handlelength=0.5,
            columnspacing=0.1,
            handletextpad=0.5,
            bbox_to_anchor=(-0.05, 0.0),
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
        if plot_icohp:
            save_path = "COHP(i).png"
        else:
            save_path = "COHP.png"
        plt.savefig(save_path, dpi=300)
        print(f"\t\tSaving  {save_path}")
        plt.close(fig)
    except Exception as e:
        print("Error while plotting COHP.")
        print(e)


def group_neighbors(cn_conns):
    """Group COHP bonds by Atom2 base name (after '-') and then by distance."""
    cn_groups = {}
    for site, neighbors in cn_conns.items():
        groups = defaultdict(lambda: defaultdict(int))
        for n_label, distance in neighbors:
            groups[distance][n_label] += 1
        cn_groups[site] = dict(groups)

    return cn_groups


def read_cohp_data(site1, site2, d, count):

    names = [
        name
        for name in os.listdir(".")
        if f"data.site_cohp_{site1}_{site2}_{count}" in name
    ]
    names = sorted(names, key=lambda name: abs(float(name.split("_")[-1]) - d))
    # name = names[0]

    if len(names) and os.path.isfile(names[0]):
        name = names[0]
        df = pd.read_csv(name, names=["energy", "cohp", "icohp"], sep=r"\s+")
        df = df[(df["energy"] >= -6.0) & (df["energy"] <= 2.0)]
        return df["energy"].values, df["cohp"].values, df["icohp"].values

    print(f"Data for {site1, site2, d, count} not found!")
    return None, None, None


def get_ticks(xmax):

    xmax = ceil(xmax)
    if xmax == 6:
        xmax = 5

    interval = int(xmax // 2)
    if interval == 0.0:
        return [-1, -0.5, 0.5, 1]

    _vals = [i for i in range(0, int(ceil(xmax)), interval) if i != 0]
    if len(_vals) == 1:
        _vals.insert(0, _vals[0] / 2)
    vals = [-v for v in _vals[::-1]]
    vals.extend(_vals)

    return vals


def plot_site_COHPs(cn_conns, label_map, plot_icohp=True):

    cn_groups = group_neighbors(cn_conns)
    color_cycle = ["red", "orange", "green", "blue", "purple", "pink"]

    for site, dist_groups in cn_groups.items():
        plt.close()
        fig, ax = plt.subplots(figsize=(8, 15))
        all_x_values = []
        idx = 0
        for d, neighbor_groups in dist_groups.items():
            for nsite, count in neighbor_groups.items():
                energy, cohp, icohp = read_cohp_data(site, nsite, d, count)
                if energy is not None:
                    color = color_cycle[idx % len(color_cycle)]
                    all_x_values.append(np.abs(cohp).max())
                    ax.plot(
                        cohp,
                        energy,
                        lw=5,
                        ls="-",
                        label=f"{nsite} (\u00D7{count})\n{d:0.3f} \u00C5",
                        c=color,
                    )
                    if plot_icohp:
                        ax.plot(icohp, energy, lw=3, ls="--", c=color)
                        all_x_values.append(icohp.max())
                    idx += 1

        # 0.05 within -6, 2
        # do for gen COHP, DOS,  0.05

        ax.set_ylim(-6, 2)
        max_x = max(all_x_values) if all_x_values else 1
        if plot_icohp:
            buffer = 0.05 * max_x
        else:
            buffer = 0.1 * max_x
        ax.set_xlim(-(max_x + buffer), max_x + buffer)
        cohp_ticks = get_ticks(max_x + buffer)
        if cohp_ticks:
            ax.set_xticks(cohp_ticks)
            tick_label = []
            for v in cohp_ticks:
                if abs(v - int(v)) == 0:
                    tick_label.append(str(int(v)))
                else:
                    tick_label.append("")
            ax.set_xticklabels(tick_label)
        ax.axhline(0, color="black", linestyle="--", linewidth=3)
        ax.axvline(0, color="black", linestyle="--", linewidth=3)
        ax.tick_params(axis="x", labelsize=35, width=3, length=10)
        ax.tick_params(axis="y", labelsize=35, width=3, length=10)

        for spine in ax.spines.values():
            spine.set_linewidth(2.5)

        ax.set_ylabel("energy (eV)", fontsize=35)
        ax.set_xlabel("-COHP", fontsize=35)
        ax.yaxis.set_major_locator(MaxNLocator(integer=True))

        legend = ax.legend(
            frameon=False,
            fontsize=30,
            loc="lower left",
            handlelength=0.5,
            columnspacing=0.1,
            handletextpad=0.5,
            bbox_to_anchor=(-0.05, 0.0),
        )
        legend.set_zorder(99)

        plt.title(site, fontsize=35)
        plt.tight_layout()

        if plot_icohp:
            fname = f"cohp(i)_{site}"
        else:
            fname = f"cohp_{site}"

        print(f"\t\tSaving site cohp for {site}: {fname}")
        plt.savefig(f"{fname}.png", dpi=300)
