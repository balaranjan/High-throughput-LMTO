from matplotlib import pyplot as plt
import pandas as pd
import os


def plot_dos(calc_dir):
    dos_files = [
        dos_file
        for dos_file in os.listdir(calc_dir)
        if dos_file.endswith("csv") and dos_file.startswith("DOS-")
    ]

    plt.close()
    print(dos_files)
    for dos_file in dos_files:
        legend = dos_file[:-4].split("-")[1]
        print(dos_file)
        dos = pd.read_csv(f"{calc_dir}{os.sep}{dos_file}", delimiter="\t")
        print(dos.head())
        plt.plot(dos["Energy (eV)"], dos["DOS"], label=legend)

    plt.savefig("dos.png")


if __name__ == "__main__":
    di = "/home/bala/Documents/22_lmto/htlmto/docker_test/AlB4Cr4-38081"
    plot_dos(di)
