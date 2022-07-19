import BioSimSpace as bss
import pandas as pd


def add_to_dictionary(dictionary, key, value) -> dict:
    """
    add values to dictionary to each key
    """
    try:
        dictionary[key].append(value)
    except KeyError:
        dictionary[key] = []
        dictionary[key].append(value)
    return dictionary


ligands_all = np.arange(1, 17, 1)
ligands = np.delete(ligands_all, 11)
print(ligands)
data_dictionary = {}
for ligand in ligands:   
    print(f"ligand: {ligand}")
    path = f"../outputs/SOMD/lig_h_1~lig_h_{ligand}/"
    free_directory = path + "free/"
    pmf_free, overlap_matrix_free = bss.FreeEnergy.Relative.analyse(free_directory)
    bound_directory = path + "bound/"
    pmf_bound, overlap_matrix_bound = bss.FreeEnergy.Relative.analyse(bound_directory)
    free_energy_difference, free_energy_error = bss.FreeEnergy.Relative.difference(pmf_bound, pmf_free)
    try:
        data_dictionary["ligand"].append(ligand)
        data_dictionary["pmf_free"].append(pmf_free)
        data_dictionary["pmf_bound"].append(pmf_bound)
        data_dictionary["om_free"].append(overlap_matrix_free)
        data_dictionary["om_bound"].append(overlap_matrix_bound)
        data_dictionary["free_energy"].append(free_energy_difference)
        data_dictionary["error"].append(free_energy_error)
    except KeyError:
        data_dictionary["ligand"] = []
        data_dictionary["pmf_free"] = []
        data_dictionary["pmf_bound"] = []
        data_dictionary["om_free"] = []
        data_dictionary["om_bound"] = []
        data_dictionary["free_energy"] = []
        data_dictionary["error"] = []
        data_dictionary["ligand"].append(ligand)
        data_dictionary["pmf_free"].append(pmf_free)        
        data_dictionary["pmf_bound"].append(pmf_bound)
        data_dictionary["om_free"].append(overlap_matrix_free)
        data_dictionary["om_bound"].append(overlap_matrix_bound)
        data_dictionary["free_energy"].append(free_energy_difference)
        data_dictionary["error"].append(free_energy_error)      

dataframe = pd.DataFrame.from_dict(data_dictionary)
savepath = "../outputs/free_energy/free_energy_molsim.csv"
dataframe.to_csv(savepath, index=False)

