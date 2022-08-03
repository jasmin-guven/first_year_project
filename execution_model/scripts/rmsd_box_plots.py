import numpy as np
import MDAnalysis as mda
import MDAnalysis.analysis.rms
import matplotlib.pyplot as plt
import os
import pandas as pd
import seaborn as sns


ligands_all = np.arange(2, 17, 1)
ligands = np.delete(ligands_all, [10, 12, 13])

times = []
rmsds = []
rmsd_dictionary = {}
stage = "bound"
for i in range(len(ligands)):
    path = f"../outputs/SOMD/lig_h_1~lig_h_{ligands[i]}/{stage}/lambda_1.0000/"
    print(f"looking at {path}")
    trajectory_filename = ""
    if os.path.exists(path+"traj000000002.dcd"):
        trajectory_filename = path+"traj000000002.dcd"
    else:
        trajectory_filename = path+"traj000000001.dcd"

    frames = []
    with mda.lib.formats.libdcd.DCDFile(trajectory_filename) as trajectory:
        for frame in trajectory:
            frames.append(frame)

    first_frame = frames[0].xyz
    universe = mda.Universe(path+"somd.prm7", trajectory_filename, topology_format="PARM7")
    reference_universe = mda.Universe(path+"somd.prm7", first_frame, topology_format="PARM7")

    ligand = universe.select_atoms("resname LIG")
    reference = reference_universe.select_atoms("resname LIG")

    rmsd = mda.analysis.rms.RMSD(ligand, reference)
    rmsd.run()

    rmsd_result = rmsd.results.rmsd.T

    time = rmsd_result[1]
    rmsd_values = rmsd_result[2]

    times.append(time)
    rmsds.append(rmsd_values)

    fig = plt.figure(figsize=(8, 8))
    ax = fig.add_subplot(111)
    ax.plot(time, rmsd_values, "k-")
    ax.set_xlabel("time (ps)", fontsize = 14)
    ax.set_ylabel(r"RMSD ($\AA$)", fontsize = 14)
    plt.savefig(f"../../plots/rmsd/{stage}/rmsd_lig_{ligands[i]}.pdf")

    # rmsd_dictionary["time"] = times
    # rmsd_dictionary["rmsd"] = rmsds


fig = plt.figure(figsize=(10, 10))
sns.set(context="notebook", palette="colorblind", style="ticks", font_scale=2)

median_line_properties = dict(linestyle='-', linewidth=2.5, color="#D0006F")
xtick_positions = np.arange(1, len(ligands) + 1, 1)
plt.boxplot(rmsd_dictionary["rmsd"], medianprops=median_line_properties)
plt.xticks(ticks=xtick_positions, labels=ligands)
plt.xlabel("Ligand")
plt.ylabel(f"RMSD ($\AA$)")
plt.savefig(f"../../plots/rmsd/{stage}/box_plot_rmsd_{stage}.pdf")
