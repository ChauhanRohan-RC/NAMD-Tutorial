import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

CAL_TO_JOULE = 4.184  # 1 cal = 4.184 J
K_b = 8.314 / (CAL_TO_JOULE * 1000)  # ideal gas constant in kcal/(mol K)

split_prob_vs_dist_file = "sp_vs_dist1.csv"
output_dist_sp_pe_file = "landscape1.csv"
T = 300 # Constant Temp (K)

df = pd.read_csv(split_prob_vs_dist_file, sep=r"\s+", comment="#")

grad = -np.gradient(df["SP"], df["DIST"])
pe = -np.log(grad) * K_b * T

df["PE"] = pe

df.replace([np.inf, -np.inf], np.nan, inplace=True)
df.dropna(subset=["PE"], axis=0, how="all", inplace=True)

df.to_csv(output_dist_sp_pe_file, sep="\t", header=True, index=False, index_label=False)

df.plot(kind="line", x="DIST", y="PE")
plt.show()
