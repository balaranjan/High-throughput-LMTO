from matplotlib import pyplot as plt
import pandas as pd
import os


def plot_dos(calc_dir):

    try:
        dos_files = [
            dos_file
            for dos_file in os.listdir(calc_dir)
            if dos_file.endswith("csv") and dos_file.startswith("DOS-")
        ]

        plt.close()
        for dos_file in dos_files:
            legend = dos_file[:-4].split("-")[1]
            dos = pd.read_csv(f"{calc_dir}{os.sep}{dos_file}")

            if legend == "all":
                plt.plot(dos["DOS"], dos["Energy (eV)"], label=legend, ls="--")
            else:
                plt.plot(dos["DOS"], dos["Energy (eV)"], label=legend)

            plt.xlabel("DOS (a.u.)")
            plt.ylabel("Energy (eV)")
            ax = plt.gca()
            ax.set_aspect(1.4)
            plt.legend()
            plt.tight_layout()

        plt.savefig("dos.png", bbox_inches="tight")
    except Exception as e:
        print("Error while plotting DOS.")
        print(e)


def plot_band_structure(calc_dir):

    band_csv = f"{calc_dir}{os.sep}band_structure.csv"

    plt.close()

    try:
        if os.path.isfile(band_csv):
            band = pd.read_csv(band_csv)
            bsp = pd.read_csv(f"{calc_dir}/band_structure_points.csv")
            x_ticks = bsp["values"]
            x_tick_labels = bsp["point"]

            plt.plot(band["k"], band["Energy (eV)"])
            plt.xticks(x_ticks, x_tick_labels)
            plt.xlim(min(x_ticks), max(x_ticks))
            plt.xlabel("Wave vector")
            plt.ylabel("Energy (eV)")
            plt.savefig("band.png")
    except Exception as e:
        print("Error while plotting band structure.")
        print(e)


def plot_cohps(calc_dir):

    try:
        cohps = [
            cohp
            for cohp in os.listdir(calc_dir)
            if cohp.startswith("data.cohp_")
        ]

        plt.close()
        if cohps:
            for cohp in cohps:
                name = cohp[10:].replace("_", "-")
                cohp = pd.read_csv(
                    f"{calc_dir}{os.sep}{cohp}",
                    delimiter="\\s+",
                    names=["energy", "cohp", "icohp"],
                )
                plt.plot(cohp.cohp, cohp.energy, label=name)
            plt.legend()
            plt.xlabel("DOS (a.u.)")
            plt.ylabel("Energy (eV)")
            ax = plt.gca()
            ax.set_aspect(0.5)
            plt.savefig("cohp.png", bbox_inches="tight")
    except Exception as e:
        print("Error while plotting COHP.")
        print(e)
