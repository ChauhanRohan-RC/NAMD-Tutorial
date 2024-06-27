import numpy as np
import scipy
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.axes

########################################################################
## Splitting Probability Sp(x) from Extension Distribution - Apparant PMF
## and reconstruct PMF from Sp(x)
########################################################################

## Usage
# 1. Copy script to working dir
# 2. INPUT: set extension vs probability file (extension dostribution)
# 3. INPUT: set Temp, LEFT and RIGHT absorbing boundares
# 4. run with "python sp_pmf.py"
# 5. Creates output file "sp_pmf.csv"

## UNITS: Energy (kcal/mol), Distance (Å), T (K)

## Constants -------------------------
CAL_TO_JOULE = 4.184  # 1 cal = 4.184 J
K_b = 8.314 / (CAL_TO_JOULE * 1000)  # ideal gas constant in kcal/(mol K)

## INPUT -----------------------
extension_pdf_file = "dist_vs_pdf.dat"    # Input Extension vs Probability file i.e Extension Probability
T = 300         # Constant Temp (K)

x_f = 0          # LEFT Absorbing Boundary - Folded state extension (in Å)
x_u = 100        # RIGHT Absorbing Boundary - Unfolded state extension (in Å)

# Whether to negate apparant PMF for Sp(x) calculation
# -> True (default) : Sp(x) will give inflection at the potential energy barrier(s)
# -> False : Sp(x) will give inflection at the potential well(s) i.e. stable states
negate_app_pmf = True  

## OUTPUT ----------------------
output_file = "sp_pmf.csv"

## ----------------------------------
# Extension Probability Distribution (PDF) DataFrame 
ext_pdf_df = pd.read_csv(extension_pdf_file, sep=r"\s+", comment="#", names=("EXT", "PDF"))

# SUbset of the dataset for integration
int_df = ext_pdf_df.loc[ext_pdf_df["EXT"] >= x_f].loc[ext_pdf_df["EXT"] <= x_u]

# Apparaent PMF - PMF from extension distribution = -K_b * T * ln(P(x)) , where P(x) is extension probability distribution
pmf_apparent = np.log(int_df["PDF"].astype("float128").values) * (-K_b * T)
int_df["PMF_APP"] = pmf_apparent

# purge infinite PMF values due to PDF = 0 (or very low) values
int_df.replace([np.inf, -np.inf], np.nan, inplace=True)
int_df.dropna(subset=["PMF_APP"], axis=0, how="all", inplace=True)

# Reciprocal of PDF
if negate_app_pmf:
    pdf_recip = (1 / int_df["PDF"].astype(np.float128).values)
    int_df["PDF_RECIP"] = pdf_recip


## Calculating Spltting Probability Sp(x) -------------------------
ext_col_index = int_df.columns.get_loc("EXT")
main_col = "PDF_RECIP" if negate_app_pmf else "PDF"   # Column which is integrated

# Integral in the denominator = Constant
c = scipy.integrate.trapezoid(y=int_df[main_col].values, x=int_df["EXT"].values)

sample_count = len(int_df["EXT"])
split_prob = np.zeros((sample_count,), dtype=np.float128)

for i in range(0, sample_count):
    _dist = int_df.iat[i, ext_col_index]

    _df = int_df.loc[int_df["EXT"] >= _dist]
    _v = scipy.integrate.trapezoid(y=_df[main_col], x=_df["EXT"])
    split_prob[i] = _v / c

int_df["SP"] = split_prob
# ------------------------------------------------------------------------

## Re-constructing PMF from Sp(x)
grad = np.gradient(int_df["SP"], int_df["EXT"])     # Gradient of Sp(fold) : ALways Negative
pmf_re = (1 if negate_app_pmf else -1) * np.log(-grad * c) * K_b * T
int_df["PMF_RE"] = pmf_re              # Reconstructed PMF from Sp(x)

# Writing Output File
int_df[["EXT", "PDF", "PMF_APP", "SP", "PMF_RE"]].to_csv(output_file, sep="\t", header=True, index=False, index_label=False)

# Plotting
fig, axes  = plt.subplots(2,2)

axes[0, 0].plot(int_df["EXT"], int_df["PDF"])
axes[0, 0].set_title("Extension DIstribution")

axes[0, 1].plot(int_df["EXT"], int_df["PMF_APP"])
axes[0, 1].set_title("Apparent PMF (from extension distribution)")

axes[1, 0].plot(int_df["EXT"], int_df["SP"])
axes[1, 0].set_title("Splitting Probability")

axes[1, 1].plot(int_df["EXT"], int_df["PMF_RE"])
axes[1, 1].set_title("Reconstructed PMF")

plt.show()
