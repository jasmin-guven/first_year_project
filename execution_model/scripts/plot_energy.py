import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns


stages = ["ligands", "systems"]
ligands = [9, 16]

for stage in stages:
    for ligand in ligands:
        path = f"equilibration/{stage}/lig_h_{ligand}/"
        print(f"Ligand {ligand} at {stage}")
        plot_save_path = "../../plots/"
        time, energy = np.loadtxt(f"{path}/potential.xvg",comments=["#","@"],unpack=True)

        fig = plt.figure()
        sns.set(context="notebook", palette="colorblind", style="ticks", font_scale=1.8, font="Sans-serif")

        plt.plot(time, energy, color="k")
        plt.tight_layout()
        stage_savename = stage[:-1]
        plt.savefig(f"{plot_save_path}/potential_{stage_savename}_{ligand}.png", dpi=1200, transparent=True)
        
        try:
            time, temperature = np.loadtxt(f"{path}/temperature.xvg",comments=["#","@"], unpack=True)
            fig = plt.figure()
            sns.set(context="notebook", palette="colorblind", style="ticks", font_scale=1.8, font="Sans-serif")
            plt.plot(time, temperature, c="black")
            plt.tight_layout()
            stage_savename = stage[:-1]
            plt.savefig(f"{plot_save_path}/temperature_{stage_savename}_{ligand}.png", dpi=1200, transparent=True)
        except IndexError:
            print(f"error on ligand {ligand} at stage: {stage}")
        try:
            time, density = np.loadtxt(f"{path}/density.xvg",comments=["#","@"], unpack=True)
            fig = plt.figure()
            sns.set(context="notebook", palette="colorblind", style="ticks", font_scale=1.8, font="Sans-serif")
            plt.plot(time, density, c="black")
            plt.tight_layout()
            stage_savename = stage[:-1]
            plt.savefig(f"{plot_save_path}/density_{stage_savename}_{ligand}.png", dpi=1200, transparent=True)
        except FileNotFoundError:
            print(f"density file for ligand {ligand} not found")
        
        
