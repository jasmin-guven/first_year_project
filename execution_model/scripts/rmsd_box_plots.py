import numpy as np
import MDAnalysis as mda
import MDAnalysis.analysis.rms
import matplotlib.pyplot as plt
import os


ligands_all = np.arange(2, 17, 1)
ligands = np.delete(ligands_all, 10)

times = []
rmsds = []
rmsd_dictionary = {}

for i in range(len(ligands)):
    path = f"../outputs/SOMD/lig_h_1~lig_h_{ligands[i]}/free/lambda_1.0000/"
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

    rmsd_dictionary["time"] = times
    rmsd_dictionary["rmsd"] = rmsds





