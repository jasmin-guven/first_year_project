import scipy.constants
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


_K_B = scipy.constants.Boltzmann
_TEMPERATURE = 300
_N_A = scipy.constants.Avogadro


def inhibition_to_ddg(ki_a: float, ki_b: float) -> float:
    """
    convert experimental Ki values to binding free-energy difference
    :param ki_a: experimental Ki of ligand 1
    :param ki_b: experimental Ki of ligand 1
    :return: experimental RBFE value
    """
    ic50_a = 2 * ki_a
    ic50_b = 2 * ki_b

    return (_K_B * _N_A * _TEMPERATURE / 4184) * np.log(ic50_b / ic50_a)


def get_experimental_error(error_a, ki_a, error_b, ki_b):
    fraction = ki_b / ki_a
    fraction_error = fraction * np.sqrt((error_b / ki_b) ** 2 + (error_a / ki_a) ** 2)
    return (_K_B * _TEMPERATURE * fraction_error / fraction) * _N_A / 4184


experimental_dataframe = pd.read_csv("../inputs/experimental_data/exp_no_12.csv")
experimental_values = experimental_dataframe["K_i"]
experimental_errors = experimental_dataframe["K_i_err"]
compound_1 = experimental_values[0]
compound_1_error = experimental_errors[0]
ki = experimental_values
ki_errors = experimental_errors

experimental_free_energy_differences = []
experimental_free_energy_errors = []
for i in range(len(ki)):
    ddg = inhibition_to_ddg(compound_1, ki[i])
    error = get_experimental_error(compound_1_error, compound_1, ki_errors[i], ki[i])
    experimental_free_energy_differences.append(ddg)
    experimental_free_energy_errors.append(error)

experimental_free_energy_differences = np.delete(experimental_free_energy_differences, 0)
experimental_free_energy_errors = np.delete(experimental_free_energy_errors, 0)

calculation_dataframe = pd.read_csv("../outputs/free_energy/free_energy_molsim.csv")
ligand_names = calculation_dataframe["ligand"].tolist()
free_energy_differences_str = calculation_dataframe["free_energy"]
free_energy_errors_str = calculation_dataframe["error"]
free_energy_differences = []
free_energy_errors = []
for i in range(len(free_energy_differences_str)):
    stripped_value = free_energy_differences_str[i].replace(" kcal/mol", "")
    stripped_error = free_energy_errors_str[i].replace(" kcal/mol", "")
    free_energy_differences.append(float(stripped_value))
    free_energy_errors.append(float(stripped_error))

fig = plt.figure(figsize=(8, 8))
sns.set(context="notebook", palette="colorblind", style="ticks", font_scale=2)
x_label_locations = np.arange(len(experimental_free_energy_differences))
bar_width = 0.35

plt.bar(x_label_locations - bar_width / 2,
        height=experimental_free_energy_differences,
        width=bar_width,
        yerr=experimental_free_energy_errors,
        label="Experimental",
        color="#0099AB")

(_, caps, _) = plt.errorbar(x_label_locations - bar_width / 2,
                            experimental_free_energy_differences,
                            color="black",
                            yerr=experimental_free_energy_errors,
                            capsize=3,
                            linestyle="")
for cap in caps:
    cap.set_color("black")
    cap.set_markeredgewidth(1.5)

plt.bar(x_label_locations + bar_width / 2,
        height=free_energy_differences,
        width=bar_width,
        yerr=free_energy_errors,
        label="AFE calculation",
        color="#D0006F",
        linewidth=0)

(_, caps, _) = plt.errorbar(x_label_locations + bar_width / 2,
                            free_energy_differences,
                            color="black",
                            yerr=free_energy_errors,
                            capsize=3,
                            linestyle="")
for cap in caps:
    cap.set_color("black")
    cap.set_markeredgewidth(1.5)

transformations = []
for name in ligand_names:
    label = f"Ligand 1 to {str(name)}"
    transformations.append(label)

ddg_dict = dict(transformations=transformations,
                ddg_exp=experimental_free_energy_differences,
                ddg_exp_err=experimental_free_energy_errors,
                ddg_afe=free_energy_differences,
                ddg_afe_err=free_energy_errors)
ddg_dataframe = pd.DataFrame.from_dict(ddg_dict)
ddg_dataframe.to_csv("../outputs/free_energy/ddg_data.csv")

plt.axhline(color="black")
plt.xlabel("Transformations")
plt.ylabel("$\Delta\Delta$G$_{\mathrm{bind}}$ kcal$\cdot$mol$^{-1}$")
plt.xticks(x_label_locations, transformations, rotation=70, ha="right")
plt.legend(loc="lower left")
plt.ylim(-5,4)
# sns.despine()
plt.tight_layout()
plt.legend(frameon=False)
plt.savefig("../../plots/free_energy/dd_g_tr.pdf", transparent=True)
# plt.show()
