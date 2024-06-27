import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy


########################################################################
## Splitting Probability Sp(x) from Extension Trajectory
########################################################################

## Usage
# 1. Copy script to working dir
# 2. INPUT: set extension vs frame file (ext over time)
# 3. INPUT: set Temp, LEFT-RIGHT absorbing boundaries, and BINS
# 4. run with "python sp_traj.py"
# 5. Creates output file "sp_traj.csv"

## UNITS: Energy (kcal/mol), Distance (Å), T (K)


def find_bin(val, bins_arr, last=False):
    """
    Helper method to find the bin_index of a given value
    from bin bounds array
    """

    if len(bins_arr) < 2 or val < bins_arr[0] or val > bins_arr[-1]:
        return -1

    if last:
        for _i in range(len(bins_arr) - 2, -1, -1):
            if val >= bins_arr[_i]:
                return _i

    for _i in range(1, len(bins_arr)):
        if val < bins_arr[_i]:
            return _i - 1

    return -1


## Constants -------------------------
CAL_TO_JOULE = 4.184  # 1 cal = 4.184 J
K_b = 8.314 / (CAL_TO_JOULE * 1000)  # ideal gas constant in kcal/(mol K)

## INPUT -----------------------
frame_vs_ext_file = "dist_vs_frame.dat"  # Input Frame vs Extension file
T = 300  # Constant Temp (K)

x_a = 19  # LEFT Absorbing Boundary - Folded state extension (in Å)
x_b = 40  # RIGHT Absorbing Boundary - Unfolded state extension (in Å)

## Number of Spacial and Temporal Bins
# -1 to use bin sizes (specified below) instead
ext_bin_count = -1
frame_bin_count = -1

# Bin sizes: Only used if number of bins is not given
ext_bin_size = 2.0  # Angstrom(s) in 1 bin
frame_bin_size = 20  # Number of frames in 1 bin

## OUTPUT ----------------------
output_file = "sp_traj.csv"

## ----------------------------------
# Frame vs Extension DataFrame
frame_ext_df: pd.DataFrame = pd.read_csv(frame_vs_ext_file, sep=r"\s+", comment="#", names=("FRAME", "EXT"))

# Making the dataset DISCRETE over SPACE (EXT) and TIME (FRAME)
if ext_bin_count <= 0:
    ext_bin_count = int((frame_ext_df["EXT"].max() - frame_ext_df["EXT"].min()) // ext_bin_size)

if frame_bin_count <= 0:
    frame_bin_count = int(len(frame_ext_df["FRAME"]) // frame_bin_size)

frame_ext_df["FRAME_BIN"] = pd.cut(frame_ext_df["FRAME"], bins=frame_bin_count, labels=False, right=False)

discrete_df: pd.DataFrame = frame_ext_df[["FRAME_BIN", "EXT"]].groupby(
    by=["FRAME_BIN"]).mean()  # TODO: mean of ext in each time-bin
ext_bin_series, ext_bins = pd.cut(discrete_df["EXT"], bins=ext_bin_count, labels=False, right=False, retbins=True)

discrete_df["EXT_BIN"] = ext_bin_series

x_a_bin = find_bin(x_a, ext_bins)
x_b_bin = find_bin(x_b, ext_bins)

if x_a_bin == -1:
    print(f"ERR: LEFT Boundary = {x_a} is out of bounds [{ext_bins[0]}, {ext_bins[-1]})")
    exit(-1)

if x_b_bin == -1:
    print(f"ERR: RIGHT Boundary = {x_b} is out of bounds [{ext_bins[0]}, {ext_bins[-1]})")
    exit(-1)

sample_count = x_b_bin - x_a_bin + 1
split_prob = np.zeros((sample_count,), dtype=np.float32)

_ext_bin_series = discrete_df["EXT_BIN"]
_ext_bin_range = np.arange(x_a_bin, x_b_bin + 1)

for i, x_bin in enumerate(_ext_bin_range):
    _df = _ext_bin_series.loc[_ext_bin_series == x_bin]
    _count = len(_df.index)

    if _count == 0:
        # Should not happen
        print(f"WARNING: EXT_BIN {x_bin} NOT FOUND. Skipping...")
        continue

    _v = 0
    for frame_bin in _df.index.values:
        _ahead = _ext_bin_series.iloc[frame_bin:]

        # Hitting function
        for __xbin in _ahead:
            if __xbin == x_a_bin:  # TODO: change strict equality
                _v += 1
                break

            if __xbin == x_b_bin:  # TODO: change strict equality
                break

    split_prob[i] = _v / _count

res_df = pd.DataFrame()

res_df["EXT_BIN"] = _ext_bin_range
res_df["EXT_BIN_START"] = res_df["EXT_BIN"].apply(lambda x: ext_bins[x])
res_df["EXT_BIN_END"] = res_df["EXT_BIN"].apply(lambda x: ext_bins[x + 1])
res_df["EXT_BIN_MED"] = res_df["EXT_BIN"].apply(lambda x: (ext_bins[x] + ext_bins[x + 1]) / 2)

res_df["SP"] = split_prob

# res_df.plot(kind="scatter", x="EXT_BIN", y="SP")
# plt.show()

# # ------------------------------------------------------------------------

## Re-constructing PMF from Sp(x)
grad = np.gradient(res_df["SP"], res_df["EXT_BIN"])  # Gradient of Sp(fold) : Always Negative

## TODO: smoothen gradient with Savitzky - Golay filter
# smooth_grad = scipy.signal.savgol_filter(grad, 5, 4)

pmf_re = np.log(-grad) * K_b * T
res_df["PMF_RE"] = pmf_re  # Reconstructed PMF from Sp(x)

res_df.to_csv(output_file, sep="\t", header=True, index=False, index_label=False)

## Smoothen SP and PMF
# ext_bin_smooth = np.linspace(x_a_bin, x_b_bin, num=len(res_df["EXT_BIN"]) * 3, endpoint=True)
#
# sp_interp_cubic = scipy.interpolate.interp1d(res_df["EXT_BIN"], res_df["SP"], kind="cubic")
# sp_smooth = sp_interp_cubic(ext_bin_smooth)
#
# pmf_interp_cubic = scipy.interpolate.interp1d(res_df["EXT_BIN"], res_df["PMF_RE"], kind="cubic", fill_value=(0,0))
# pmf_re_smooth = pmf_interp_cubic(ext_bin_smooth)


ext_bin_med = res_df["EXT_BIN_MED"]
ext_bin_med_len = len(ext_bin_med)

ext_smooth = np.linspace(ext_bin_med[0], ext_bin_med[ext_bin_med_len - 1], num=ext_bin_med_len * 3, endpoint=True)
sp_interp_cubic = scipy.interpolate.interp1d(ext_bin_med, res_df["SP"], kind="cubic")
sp_smooth = sp_interp_cubic(ext_smooth)

# purge infinite PMF values due to grad(Sp) = 0 (or very low) values
res_df2 = res_df.replace([np.inf, -np.inf], np.nan, inplace=False)
res_df2 = res_df2.dropna(subset=["PMF_RE"], axis=0, how="all", inplace=False)

ext_smooth2 = None
pmf_re_smooth = None

try:
    ext_bin_med2 = res_df2["EXT_BIN_MED"]
    ext_bin_med_len2 = len(ext_bin_med2)

    ext_smooth2 = np.linspace(ext_bin_med2[0], ext_bin_med2[ext_bin_med_len2 - 1], num=ext_bin_med_len2 * 3, endpoint=True)
    pmf_interp_cubic = scipy.interpolate.interp1d(ext_bin_med2, res_df2["PMF_RE"], kind="cubic")
    pmf_re_smooth = pmf_interp_cubic(ext_smooth2)
except Exception as exc:
    print("ERR: Exception in interpolating Reconstructed PMF...Ignoring")

# Plotting
fig, axes = plt.subplots(2, 2)

axes[0, 0].plot(frame_ext_df["FRAME"], frame_ext_df["EXT"])
axes[0, 0].set_title("Trajectory")

axes[0, 1].scatter(discrete_df.index, discrete_df["EXT_BIN"])
axes[0, 1].set_title("Discrete Trajectory")

axes[1, 0].scatter(res_df["EXT_BIN_MED"], res_df["SP"])
# axes[1, 0].plot(res_df["EXT_BIN"], res_df["SP"])
axes[1, 0].plot(ext_smooth, sp_smooth, '--')
axes[1, 0].set_title("SP trajectory")

## Smoothen PMF
axes[1, 1].scatter(res_df2["EXT_BIN_MED"], res_df2["PMF_RE"])
# axes[1, 1].plot(res_df["EXT_BIN"], res_df["PMF_RE"])
if ext_smooth2 is not None and pmf_re_smooth is not None:
    axes[1, 1].plot(ext_smooth2, pmf_re_smooth, '--')

axes[1, 1].set_title("Reconstructed PMF")

plt.show()
