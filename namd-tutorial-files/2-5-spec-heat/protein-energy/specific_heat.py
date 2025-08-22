import pandas as pd

print("----------- Script to analyze NAMD-Energy Plugin output file ------------")
print("-> About Plugin: VMD -> Extensions -> Analysis -> NAMD Energy. It is generally used to determine energies at every "
      "frame for a subset of the system (like for protein only) in contrast to "
      "running simulation directly which output energies of the entire system at each frame in .log file.\n"
      "-> About this script: Used to calculate Energy averages and Specific Heat (Cv) of the protein alone\n"
      "\n-> USAGE: copy this script in working dir and use 'python3 specific_heat.py'")
print("-----------------------------------------------------------------------------------------------------")

CAL_TO_JOULE = 4.184  # 1 cal = 4.184 J
K_b = 8.314 / (CAL_TO_JOULE * 1000)  # ideal gas constant in kcal/(mol K)

file_name = input("-> Specify NAMD-Energy Plugin output file: ")

# Data Frame
df: pd.DataFrame = pd.read_csv(file_name, sep=r"\s+")
print(f"Available Columns: [{', '.join(df.columns.values)}]")

# Selecting Energy Series
energy_col_name = ''
while True:
    energy_col_name = input("-> Enter Energy Column (one from above): ")
    if energy_col_name in df.columns.values:
        break
    else:
        print(f"There is no column '{energy_col_name}'. Please specify one of the above listed column...")
        continue

energy_series: pd.Series = df[energy_col_name]

# Average Energy <E> (kcal/mol)
avg_energy = energy_series.mean()

# Average Squared Energy <E^2> (kcal/mol)^2
avg_sq_energy = energy_series.map(lambda x: x * x).mean()

# Variance in Energy (<E^2> - <E>^2) (kcal/mol)^2
var_energy = avg_sq_energy - (avg_energy ** 2)

print("----------------------------------------------------")
print(f"Average {energy_col_name} Energy <E>: {avg_energy} kcal/mol")
print(f"Average {energy_col_name} Squared Energy <E^2>: {avg_sq_energy} (kcal/mol)^2")
print(f"Variance in {energy_col_name} Energy (<E^2> - <E>^2): {var_energy} (kcal/mol)^2")
print("----------------------------------------------------")

# Specific Heat
print(f"\n-------------- Specific Heat (Cv) for {energy_col_name} Energy --------------")
print("NOTE: Cv = var(E) / (kb * T^2), where var(E) = (<E^2> - <E>^2)")


def _specific_heat(_temp: float, _molar_mass_grams: float):
    # Cv = var(E_tot) / (kb * T^2) = (<E^2> - <E>^2) / (kb * T^2)
    mol_spec_heat = var_energy / (K_b * _temp * _temp)

    print("--------------------------------------------------------")
    print(
        f"Molar Specific Heat (Cv_n) for {energy_col_name} energy = {mol_spec_heat} kcal/(mol K)"
        f" = {mol_spec_heat * CAL_TO_JOULE * 1000} J/(mol K)")

    if (_molar_mass_grams > 0):
        mass_spec_heat = mol_spec_heat / (_molar_mass_grams * 1000)
        print(f"Mass Specific Heat (Cv_m) for {energy_col_name} energy = {mass_spec_heat} kcal/(Kg K)"
              f" = {mass_spec_heat * CAL_TO_JOULE * 1000} J/(Kg K)")
    print("--------------------------------------------------------")


temp = 0
molar_mass_grams = 0
while True:
    try:
        temp = float(input("-> Enter Temperature (K) [0: skip]: "))
        if temp == 0:
            print("Skipping specific heat calculation...")
            break
        if temp < 0:
            raise ValueError("Absolute temp must be > 0")
    except ValueError:
        print(f"Invalid temperature '{temp}', must be an int or float > 0", flush=True)
        continue

    while True:
        try:
            molar_mass_grams = float(input("-> Enter Molar mass of the system [or the subset of system]"
                                           " (grams) [0: skip]: "))
            if molar_mass_grams == 0:
                print("Skipping mass specific heat...")
                break
            elif molar_mass_grams < 0:
                raise ValueError("Molar mass must be > 0")
            else:
                break
        except ValueError:
            print(f"Invalid Molar Mass '{molar_mass_grams}', must be an int or float > 0", flush=True)
            continue

    _specific_heat(temp, molar_mass_grams)
    break
