import pandas as pd

"""
Script to calculate Specific Heat (Cv) from Energy vs Time Step data (TS column is not required)

Used to calculate   1. Energy averages, Variance
                    2. Specific Heat (Cv) for particular Energy type
                    -> Cv is usually calculated for a subset of the system, like protein or nucleic acid alone (excluding water and ions)
                
NOTE: Energy vs TS data can be generated from
1. namd_energy_log.py 
    -> extracts Energy vs TS from NAMD .log file(s)   (for the entire system, may contain water or ions)
2. namd_energy_plugin.tcl  => VMD's NAMD-Energy Plugin
    -> used to calculate Energies at each TS for a subset of the system, like protein alone.
    
USAGE:
1. copy this script in working dir
2. "python specific_heat.py"
3. Enter file_name (containing energies at each timestep)
"""

print("================== Script to calculate Cv ===================")

CAL_TO_JOULE = 4.184  # 1 cal = 4.184 J
AVOGADRO_NUMBER = 6.022e23

KCAL_PER_MOL_TO_JOULE = CAL_TO_JOULE * 1000 / AVOGADRO_NUMBER      # energy 1 kcal/mol = 0.694 x 10^-20 J
Kb_J_PER_MOL_K = 8.314  # ideal gas constant in J/(mol K)
Kb_KCAL_PER_MOL_K = Kb_J_PER_MOL_K / (CAL_TO_JOULE * 1000)  # Boltzmann constant (kcal/mol/K) = 8.314 / (4.18 x 10^3)
Kb_J_PER_MOLECULE_K = Kb_J_PER_MOL_K / AVOGADRO_NUMBER   # Boltzmann Constant Kb = R/Na J/ (molecule K)

DEFAULT_TEMPERATURE = 300       # Kelvin
#DEFAULT_NUM_MOL = 1

file_name = input(" -> Specify Energy data File (with energy at each TS): ")
in_delimiter = r"\s+"  # All whitespaces
comment_token = "#"  # To skip Comments

# Data Frame
df: pd.DataFrame = pd.read_csv(file_name, sep=in_delimiter, comment=comment_token)
print("----------------------------------------")
print(f"LOG: Available Columns: [{', '.join(df.columns.values)}]")
print("----------------------------------------")

# Selecting Energy Series
energy_col_name = ''
while True:
    energy_col_name = input("\n -> Enter Energy Column (one from above): ")
    if energy_col_name in df.columns.values:
        break
    else:
        print(f"ERR: There is no column '{energy_col_name}'. Please specify one of the above listed column...")
        continue

energy_series: pd.Series = df[energy_col_name]

# Average Energy <E> (kcal/mol)
avg_energy = energy_series.mean()
avg_energy_j = avg_energy * KCAL_PER_MOL_TO_JOULE   # In Joules

# Average Squared Energy <E^2> (kcal/mol)^2
avg_sq_energy = energy_series.map(lambda x: x * x).mean()
avg_sq_energy_j2 = avg_sq_energy * KCAL_PER_MOL_TO_JOULE * KCAL_PER_MOL_TO_JOULE  # in Joules^2

# Variance in Energy (<E^2> - <E>^2) (kcal/mol)^2
var_energy = avg_sq_energy - (avg_energy * avg_energy)
var_energy_j2 = var_energy * KCAL_PER_MOL_TO_JOULE * KCAL_PER_MOL_TO_JOULE  # in Joules^2

print("\n==============================================")
print(f"Note: 1 kcal/mol = 4180/N_A J (N_A: Avogadro Number) = 0.694 x 10^-20 J")
print("----------------")
print(f"INFO: Average {energy_col_name} Energy <E>: {avg_energy} kcal/mol = {avg_energy_j} J")
print(f"INFO: Average Squared {energy_col_name} Energy <E^2>: {avg_sq_energy} (kcal/mol)^2 = {avg_sq_energy_j2} J^2")
print(f"INFO: Variance in {energy_col_name} Energy (<E^2> - <E>^2): {var_energy} (kcal/mol)^2 = {var_energy_j2} J^2")
print("==============================================")

# Specific Heat
print("\n-----------------------------------------------------")
print(f" Specific Heat (Cv) for {energy_col_name} Energy")
print("-----------------------------------------------------")
print("# NOTE")
# print(
#     "=> Energies are in kcal/mol, where mole represents Avogadro number of repetitions of the entire selected system.")
# print("   Hence, Energies must be divided by number of molecules N (of protein/dna/water) to get actual kcal/mol")
print("=> Cv = var(E) / (kb * T^2), where var(E) = <E^2> - <E>^2")


def input_temp(default_val: float = DEFAULT_TEMPERATURE):
    _temp = -1
    _temp_str = ''
    try:
        _temp_str = input(f"\n -> Enter Temperature (K) [default: {default_val} K]: ")
        if len(_temp_str.strip()) == 0:
            _temp = default_val  # default
        else:
            _temp = float(_temp_str.strip())
            if _temp <= 0:
                raise ValueError("ERR: Absolute Temperature must be > 0")
    except ValueError:
        print(f"ERR: Invalid Temperature '{_temp_str}', must be an int or float > 0", flush=True)

    return _temp


# def input_num_mol(default_val: int = DEFAULT_NUM_MOL):
#     _num_mol = -1
#     _num_mol_str = ''
#     try:
#         _num_mol_str = input("\n -> Enter Number of Molecules [default: 1]: ")
#         if len(_num_mol_str.strip()) == 0:
#             _num_mol = default_val  # default
#         else:
#             _num_mol = int(_num_mol_str.strip())
#             if _num_mol <= 0:
#                 raise ValueError("ERR: Number of molecules must be > 0")
#     except ValueError:
#         print(f"ERR: Invalid No. of Molecules: '{_num_mol_str}', must be an int > 0", flush=True)
#
#     return _num_mol


def input_mass_grams_per_mole():
    _mol_mass = -1
    _mol_mass_str = ''
    try:
        _mol_mass_str = input("\n -> Enter TOTAL MASS of selected system (in grams/mol) [OPTIONAL]: ")
        if len(_mol_mass_str.strip()) == 0:
            _mol_mass = 0  # default
        else:
            _mol_mass = float(_mol_mass_str.strip())
            if _mol_mass <= 0:
                raise ValueError("ERR: Mass must be > 0")
    except ValueError:
        print(f"ERR: Invalid Mass: '{_mol_mass_str}', must be a float > 0", flush=True)

    return _mol_mass


def _specific_heat(_temp: float, _mass_grams_per_mole: float):
    # Cv = var(E_tot) / (kb * T^2) = (<E^2> - <E>^2) / (N * kb * T^2)
    spec_heat = var_energy / (Kb_KCAL_PER_MOL_K * _temp * _temp) # in kcal/(mol K)
    spec_heat_j = spec_heat * KCAL_PER_MOL_TO_JOULE              # in J/K

    print("\n-------------------------------------------------------------------------")
    print(f" Specific Heat (Cv_n) for {energy_col_name} energy at {_temp} K = {spec_heat} kcal/(mol K) = {spec_heat_j} J/K")

    if _mass_grams_per_mole > 0:
        # Using mass (kg/molecule)
        # mass_kg_per_molecule = _mass_grams_per_mole / (AVOGADRO_NUMBER * 1000)
        mass_spec_heat = spec_heat / (_mass_grams_per_mole / 1000)   # in kcal/(Kg K)
        mass_spec_heat_j = mass_spec_heat * CAL_TO_JOULE * 1000      # in J/(kg K)
        print(f" Mass Specific Heat (Cv_m) for {energy_col_name} energy at {_temp} K = {mass_spec_heat} kcal/(Kg K) = {mass_spec_heat_j} J/(Kg K)\n")
    print("--------------------------------------------------------------------------")


## MAIN
temp = input_temp(DEFAULT_TEMPERATURE)
while temp <= 0:
    temp = input_temp(DEFAULT_TEMPERATURE)

# num_mol = input_num_mol(DEFAULT_NUM_MOL)
# while num_mol <= 0:
#     num_mol = input_num_mol(DEFAULT_NUM_MOL)

mass_grams_per_mole = input_mass_grams_per_mole()
if mass_grams_per_mole <= 0:
    print("LOG: skipping Mass Specific Heat (Cv_m) ...")

_specific_heat(_temp=temp, _mass_grams_per_mole=mass_grams_per_mole)
