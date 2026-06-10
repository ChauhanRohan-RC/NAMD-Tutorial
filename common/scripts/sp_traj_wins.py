import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


########################################################################
## TEST: Splitting Probability Sp(x) from Extension Trajectory in Small Windows
########################################################################

## Usage
# 1. Copy script to working dir
# 2. INPUT: set extension vs frame file (ext over time)
# 3. INPUT: set Temp, LEFT-RIGHT absorbing boundaries, and BINS
# 4. run with "python sp_traj_wins.py"
# 5. Creates output file(s) "sp_traj-win-{i}.csv"

## UNITS: Energy (kcal/mol), Distance (Å), T (K)


def find_bin(val, bins_arr, last=False):
    """
    Helper method to find the bin_index of a given value
    from bin bounds array
    """
    if last:
        for _i in range(len(bins_arr) - 1, -1, -1):
            if val > bins_arr[_i]:
                return _i

    for _i in range(len(bins_arr)):
        if val < bins_arr[_i]:
            return _i - 1

    return -1


## Constants -------------------------
CAL_TO_JOULE = 4.184  # 1 cal = 4.184 J
K_b = 8.314 / (CAL_TO_JOULE * 1000)  # ideal gas constant in kcal/(mol K)

## INPUT -----------------------
frame_vs_ext_file = "dist_vs_frame.dat"  # Input Frame vs Extension file
T = 300  # Constant Temp (K)

x_a = 10  # LEFT Absorbing Boundary - Folded state extension (in Å)
x_b = 70  # RIGHT Absorbing Boundary - Unfolded state extension (in Å)

window_count = 5    # Number of windows in Extension Range

## Number of Spacial and Temporal Bins
# -1 to use bin sizes (specified below) instead
ext_bin_count = -1
frame_bin_count = -1

# Bin sizes: Only used if number of bins is not given
ext_bin_size = 3  # Angstrom(s) in 1 bin
frame_bin_size = 20  # Number of frames in 1 bin

## OUTPUT ----------------------
output_file_name_prefix = "sp_traj"
output_file_name_suffix = ".csv"

## ----------------------------------
# Frame vs Extension DataFrame
frame_ext_df: pd.DataFrame = pd.read_csv(frame_vs_ext_file, sep=r"\s+", comment="#", names=("FRAME", "EXT"))

# Making the dataset DISCRETE over SPACE (EXT) and TIME (FRAME)
if ext_bin_count <= 0:
    ext_bin_count = int((frame_ext_df["EXT"].max() - frame_ext_df["EXT"].min()) // ext_bin_size)

if frame_bin_count <= 0:
    frame_bin_count = int(len(frame_ext_df["FRAME"]) // frame_bin_size)

frame_ext_df["FRAME_BIN"] = pd.cut(frame_ext_df["FRAME"], bins=frame_bin_count, labels=False, right=False)

discrete_df: pd.DataFrame = frame_ext_df[["FRAME_BIN", "EXT"]].groupby(by=["FRAME_BIN"]).mean()     # TODO: mean of ext in each time-bin
ext_bin_series, ext_bins = pd.cut(discrete_df["EXT"], bins=ext_bin_count, labels=False, right=False, retbins=True)
discrete_df["EXT_BIN"] = ext_bin_series

## ------------------------- MAIN Calculation -------------------------------

# x_a_bin = find_bin(x_a, ext_bins)
# x_b_bin = find_bin(x_b, ext_bins)

win_size = abs(x_b - x_a) / window_count
res_dfs = []

for i in range(0, (window_count * 2) - 1):
    ws = x_a + (i * win_size / 2)
    we = ws + win_size

    print(f"Window {i + 1} => start: {ws} | end: {we}")

    _df = discrete_df.loc[discrete_df["EXT"] >= ws].loc[discrete_df["EXT"] <= we]
    _xa_bin = _df["EXT_BIN"].min()
    _xb_bin = _df["EXT_BIN"].max()

    _sample_count = _xb_bin - _xa_bin + 1
    _split_prob = np.zeros((_sample_count, ), dtype=np.float32)

    _ext_bin_range = np.arange(_xa_bin, _xb_bin + 1)

    for i, x_bin in enumerate(_ext_bin_range):
        __df = ext_bin_series.loc[ext_bin_series == x_bin]
        __count = len(__df.index)

        __v = 0
        for frame_bin in __df.index.values:
            __ahead = ext_bin_series.iloc[frame_bin:]

            # Hitting function
            for __xbin in __ahead:
                if __xbin == _xa_bin:       # TODO: change strict equality
                    __v += 1
                    break

                if __xbin == _xb_bin:          # TODO: change strict equality
                    break

        _split_prob[i] = __v / __count

    res_df = pd.DataFrame({ "EXT_BIN": _ext_bin_range, "SP": _split_prob })

    # PMF Reconstruction
    grad = np.gradient(res_df["SP"], res_df["EXT_BIN"])

    # TODO: smoothen gradient with Savitzky - Golay filter
    # smooth_grad = scipy.signal.savgol_filter(grad, 5, 4)

    # TODO: constant c to match segments
    pmf_re = np.log(-grad * 1) * K_b * T
    res_df["PMF_RE"] = pmf_re

    ## Subsidiary Stuff
    res_df["EXT_BIN_START"] = res_df["EXT_BIN"].apply(lambda x: ext_bins[x])
    res_df["EXT_BIN_END"] = res_df["EXT_BIN"].apply(lambda x: ext_bins[x + 1])
    res_df["EXT_BIN_MED"] = res_df["EXT_BIN"].apply(lambda x: (ext_bins[x] + ext_bins[x + 1]) / 2)

    res_dfs.append(res_df)


for i, _df in enumerate(res_dfs):
    _df[["EXT_BIN", "EXT_BIN_START", "EXT_BIN_END", "EXT_BIN_MED", "SP", "PMF_RE"]].to_csv(f"{output_file_name_prefix}-win{i + 1}{output_file_name_suffix}", sep="\t", header=True, index=False, index_label=False)
    plt.plot(_df["EXT_BIN_MED"], _df["PMF_RE"])

plt.title("Reconstructed PMF in Windows")
plt.xlabel("Extension (Å)")
plt.ylabel("PMF (kcal/mol)")
plt.show()


# sample_count = x_b_bin - x_a_bin + 1
# split_prob = np.zeros((sample_count, ), dtype=np.float32)
#
# _ext_bin_series = discrete_df["EXT_BIN"]
# _ext_bin_range = np.arange(x_a_bin, x_b_bin + 1)
#
# for i, x_bin in enumerate(_ext_bin_range):
#     _df = _ext_bin_series.loc[_ext_bin_series == x_bin]
#     _count = len(_df.index)
#
#     _v = 0
#     for frame_bin in _df.index.values:
#         _ahead = _ext_bin_series.iloc[frame_bin:]
#
#         # Hitting function
#         for __xbin in _ahead:
#             if __xbin == x_a_bin:       # TODO: change strict equality
#                 _v += 1
#                 break
#
#             if __xbin == x_b_bin:          # TODO: change strict equality
#                 break
#
#     split_prob[i] = _v / _count
#
#
# res_df = pd.DataFrame()
#
# res_df ["EXT_BIN"] = _ext_bin_range
# res_df["EXT_BIN_START"] = res_df["EXT_BIN"].apply(lambda x: ext_bins[x])
# res_df["EXT_BIN_END"] = res_df["EXT_BIN"].apply(lambda x: ext_bins[x + 1])
# res_df["EXT_BIN_MED"] = res_df["EXT_BIN"].apply(lambda x: (ext_bins[x] + ext_bins[x + 1]) / 2)
#
# res_df["SP"] = split_prob
#
#
# # res_df.plot(kind="scatter", x="EXT_BIN", y="SP")
# # plt.show()
#
# # # ------------------------------------------------------------------------
#
# ## Re-constructing PMF from Sp(x)
# grad = np.gradient(res_df["SP"], res_df["EXT_BIN"])     # Gradient of Sp(fold) : ALways Negative
#
# ## TODO: smoothen gradient with Savitzky - Golay filter
# # smooth_grad = scipy.signal.savgol_filter(grad, 5, 4)
#
# pmf_re = np.log(-grad) * K_b * T
# res_df["PMF_RE"] = pmf_re              # Reconstructed PMF from Sp(x)
#
# res_df.to_csv(output_file, sep="\t", header=True, index=False, index_label=False)
#
#
# ## Smoothen SP and PMF
# # ext_bin_smooth = np.linspace(x_a_bin, x_b_bin, num=len(res_df["EXT_BIN"]) * 3, endpoint=True)
# #
# # sp_interp_cubic = scipy.interpolate.interp1d(res_df["EXT_BIN"], res_df["SP"], kind="cubic")
# # sp_smooth = sp_interp_cubic(ext_bin_smooth)
# #
# # pmf_interp_cubic = scipy.interpolate.interp1d(res_df["EXT_BIN"], res_df["PMF_RE"], kind="cubic", fill_value=(0,0))
# # pmf_re_smooth = pmf_interp_cubic(ext_bin_smooth)
#
#
# ext_bin_med = res_df["EXT_BIN_MED"]
# ext_bin_med_len = len(ext_bin_med)
#
# ext_smooth = np.linspace(ext_bin_med[0], ext_bin_med[ext_bin_med_len - 1], num=ext_bin_med_len * 3, endpoint=True)
# sp_interp_cubic = scipy.interpolate.interp1d(ext_bin_med, res_df["SP"], kind="cubic")
# sp_smooth = sp_interp_cubic(ext_smooth)
#
# pmf_interp_cubic = scipy.interpolate.interp1d(ext_bin_med, res_df["PMF_RE"], kind="cubic")
# pmf_re_smooth = pmf_interp_cubic(ext_smooth)
#
# # Plotting
# fig, axes  = plt.subplots(2,2)
#
# axes[0, 0].plot(frame_ext_df["FRAME"], frame_ext_df["EXT"])
# axes[0, 0].set_title("Trajectory")
#
# axes[0, 1].scatter(discrete_df.index, discrete_df["EXT_BIN"])
# axes[0, 1].set_title("Discrete Trajectory")
#
# axes[1, 0].scatter(res_df["EXT_BIN_MED"], res_df["SP"])
# # axes[1, 0].plot(res_df["EXT_BIN"], res_df["SP"])
# axes[1, 0].plot(ext_smooth, sp_smooth, '--')
# axes[1, 0].set_title("SP trajectory")
#
# ## Smoothen PMF
# axes[1, 1].scatter(res_df["EXT_BIN_MED"], res_df["PMF_RE"])
# # axes[1, 1].plot(res_df["EXT_BIN"], res_df["PMF_RE"])
# axes[1, 1].plot(ext_smooth, pmf_re_smooth, '--')
# axes[1, 1].set_title("Reconstructed PMF")
#
# plt.show()
