import pandas as pd

"""
Script to calculate Specific Heat (Cv) from Energy vs Time Step data (TS column is not required)

Used to calculate   1. Energy averages, Variance
                    2. Specific Heat (Cv) for particular Energy type
                    -> Cv is usually calculated for a subset of the system, like protein or nucleic acid alone (excluding water and ions)
                
NOTE: Energy vs TS data can be generated from
1. namd_energy.py 
    -> extracts Energy vs TS from NAMD .log file(s)   (for the entire system, may contain water or ions)
2. from VMD "NAMD-Energy Plugin" [VMD -> Extensions -> Analysis -> NAMD Energy]
    -> used to calculate Energies at each TS for a subset of the system, like protein alone.
    
USAGE:
1. copy this script in working dir
2. "python specific_heat.py"
3. Enter file_name (containing energies at each timestep)
"""

print("================== Script to calculate Cv ===================")

CAL_TO_JOULE = 4.184  # 1 cal = 4.184 J
K_b = 8.314 / (CAL_TO_JOULE * 1000)  # ideal gas constant in kcal/(mol K)

file_name = input(" -> Specify Energy data File (with energy at each TS): ")
in_delimiter = r"\s+"   # All whitespaces
comment_token = "#"     #To skip Comments

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

# Average Squared Energy <E^2> (kcal/mol)^2
avg_sq_energy = energy_series.map(lambda x: x * x).mean()

# Variance in Energy (<E^2> - <E>^2) (kcal/mol)^2
var_energy = avg_sq_energy - (avg_energy ** 2)

print("----------------------------------------------------")
print(f"INFO: Average {energy_col_name} Energy <E>: {avg_energy} kcal/mol")
print(f"INFO: Average Squared {energy_col_name} Energy <E^2>: {avg_sq_energy} (kcal/mol)^2")
print(f"INFO: Variance in {energy_col_name} Energy (<E^2> - <E>^2): {var_energy} (kcal/mol)^2")
print("----------------------------------------------------")

# Specific Heat
print(f"\n-------------- Specific Heat (Cv) for {energy_col_name} Energy --------------")
print("NOTE: Cv = var(E) / (kb * T^2), where var(E) = <E^2> - <E>^2")


def _specific_heat(_temp: float, _molar_mass_grams: float):
    # Cv = var(E_tot) / (kb * T^2) = (<E^2> - <E>^2) / (kb * T^2)
    mol_spec_heat = var_energy / (K_b * _temp * _temp)

    print("-------------------------------------------------------------------------")
    print(f" Molar Specific Heat (Cv_n) for {energy_col_name} energy = {mol_spec_heat} kcal/(mol K)\n"
          f"                                              = {mol_spec_heat * CAL_TO_JOULE * 1000} J/(mol K)\n")

    if (_molar_mass_grams > 0):
        mass_spec_heat = mol_spec_heat / (_molar_mass_grams * 1000)
        print(f" Mass Specific Heat (Cv_m) for {energy_col_name} energy = {mass_spec_heat} kcal/(Kg K)\n"
              f"                                        = {mass_spec_heat * CAL_TO_JOULE * 1000} J/(Kg K)")
    print("--------------------------------------------------------------------------")


temp = 0
molar_mass_grams = 0
while True:
    try:
        temp_str = input("\n -> Enter Temperature (K) [0: skip]: ")
        temp = float(temp_str)
        if temp == 0:
            print("LOG: skipping Specific Heat (Cv) calculation ...")
            break
        if temp < 0:
            raise ValueError("ERR: Absolute temp must be > 0")
    except ValueError:
        print(f"ERR: Invalid Temperature '{temp_str}', must be an int or float > 0", flush=True)
        continue

    while True:
        try:
            molar_mass_str = input("\n -> Enter Molar mass of the system [or its subset] (grams) [0: skip]: ")
            molar_mass_grams = float(molar_mass_str)
            if molar_mass_grams == 0:
                print("LOG: skipping Mass Specific Heat (Cv_m) ...")
                break
            elif molar_mass_grams < 0:
                raise ValueError("ERR: Molar mass must be > 0")
            else:
                break
        except ValueError:
            print(f"ERR: Invalid Molar Mass '{molar_mass_str}', must be a number > 0 !!", flush=True)
            continue

    _specific_heat(temp, molar_mass_grams)
    break
