import numpy as np
import scipy as sp
import pandas as pd
import matplotlib.pyplot as plt
import seaborn
import os

CAL_TO_JOULE = 4.184  # 1 cal = 4.184 J
K_b = 8.314 / (CAL_TO_JOULE * 1000)  # ideal gas constant in kcal/(mol K)

dist_vs_frame_file = "dist_vs_frame.dat"
potential_temp_vs_ts_file = "potential-temp_vs_ts.csv"
output_sp_vs_dist_file = "sp_vs_dist1.csv"

output_potential_vs_dist_file = "potential_vs_dist.csv"       # Can be None

T = 300     # Constant Temp (K)

xf = 28     # Folded state extension (in A)
xu = 30     # Unfolded state extension (in A)

frame_dist_df = pd.read_csv(dist_vs_frame_file, sep=r"\s+", comment="#", header=None)
ts_potential_temp_df = pd.read_csv(potential_temp_vs_ts_file, sep=r"\s+", comment="#")

df = frame_dist_df[1].to_frame("DIST")
df["POTENTIAL"] = ts_potential_temp_df["POTENTIAL"]
# df["TEMP"] = ts_potential_temp_df["TEMP"]

df.dropna(axis=0)  # drops all rows with at-least one N/A
df.sort_values("DIST", inplace=True)

if output_potential_vs_dist_file and not os.path.isfile(output_potential_vs_dist_file):
    df.to_csv(output_potential_vs_dist_file, sep="\t", header=True, index=False, index_label=False)

int_df: pd.DataFrame = df.loc[df["DIST"] >= xf].loc[df["DIST"] <= xu]
# int_df["EXP"] = np.exp(int_df["POTENTIAL"].astype("float128") / (K_b * df["TEMP"]))
int_df["EXP"] = np.exp(int_df["POTENTIAL"].astype("float128") / (K_b * T))

c = sp.integrate.trapezoid(y=int_df["EXP"], x=int_df["DIST"])

dist_col_index = int_df.columns.get_loc("DIST")
sample_count = len(int_df["DIST"])

res = np.zeros((sample_count, ), dtype=np.float128)

for i in range(0, sample_count):
    _dist = int_df.iat[i, dist_col_index]

    _df = int_df.loc[int_df["DIST"] >= _dist]
    _v = sp.integrate.trapezoid(y=_df["EXP"], x=_df["DIST"])
    res[i] = _v / c

int_df["SP"] = res

int_df[["DIST", "SP"]].to_csv(output_sp_vs_dist_file, sep="\t", header=True, index=False, index_label=False)

int_df.plot(kind="line", x="DIST", y="SP")
plt.show()
