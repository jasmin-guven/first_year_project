import BioSimSpace as bss
import glob
ligand_path = "../inputs/ligands/"
ligand_files = glob.glob(ligand_path+"*.sdf")

ligands = []
ligand_names = []
for file in ligand_files:
    ligand_name = file.split("/")[-1].replace(".sdf", "")
    ligand_names.append(ligand_name)
    print(ligand_name)
    ligand = bss.IO.readMolecules(file)[0]
    ligands.append(ligand)
    
transformations, lomap_scores = bss.Align.generateNetwork(ligands, plot_network=True, names=ligand_names)
print(transformations, lomap_scores)