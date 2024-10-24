import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.figure import figaspect

######################################################################
## Script to plot Extension Probability Density Function (PDF)      ##
######################################################################

## Usage
# 1. Copy script to working dir
# 2. INPUT: set EXTENSION_VS_FRAME file (ext over time)
# 3. INPUT: set HISTOGRAM and ROLLING_AVG params
# 4. run with "python ext_pdf.py"
# 5. Creates output file "ext_pdf.dat", "ext_pdf-avg[bins].dat", "ext_pdf.svg"

## UNITS: Distance (Å)

COMMENT_TOKEN = "#"

## Input ---------------------------------------------------
frame_vs_ext_file = "dist_vs_frame.dat"
frame_step_fs = 2 * 100  # time (in fs) between frames = time_step (fs) * dcd_freq. -1 for NOT_DEFINED

# Frame Range (optional). NOTE: frame_index = time_fs / frame_step_fs
frame_index_start = -1  # Inclusive [-1 for no start bound]
frame_index_end = -1 # 32e6 / frame_step_fs  # Exclusive [-1 for no end bound]

# Histogram parameters
ext_start: float = 0.0
ext_end: float = 100.0
ext_bin_count: int = 1000
_ext_bin_size = (ext_end - ext_start) / ext_bin_count

# Moving average
rolling_window_bins: int = 20  # [int] Moving Average window size (in no of bins). Set to 0 to disable
rolling_window_size: float = -1  # [float] [Only used if rolling_window_bins is not set]. Moving Average window size (in Angstroms). Set to 0 to disable.


def parse_rolling_win_bins() -> tuple[int, float]:
    bins = -1

    if rolling_window_bins > 0:
        bins = rolling_window_bins
    elif rolling_window_size > 0:
        bins = int(round(rolling_window_size / ((ext_end - ext_start) / ext_bin_count)))

    return (bins, bins * _ext_bin_size) if bins > 0 else (-1, -1)


_rolling_win_bins, _rolling_win_size = parse_rolling_win_bins()

## Output --------------------------------------------------------------
output_pdf_data_file = "ext_pdf.csv"
output_pdf_avg_data_file = f"ext_pdf-avg{_rolling_win_bins}.csv" if _rolling_win_bins > 0 else ""
output_fig_file = "ext_pdf.pdf"
comment_output_header = False

# ----------------------------------------------------------------------

# Dataframe
frame_ext_df: pd.DataFrame = pd.read_csv(frame_vs_ext_file, sep=r"\s+", comment=COMMENT_TOKEN, names=("FRAME", "EXT"))
if frame_index_start >= 0:
    frame_ext_df = frame_ext_df[frame_ext_df["FRAME"] >= frame_index_start]

if frame_index_end >= 0:
    frame_ext_df = frame_ext_df[frame_ext_df["FRAME"] < frame_index_end]

# Histogram
ext_hist, ext_bin_edges = np.histogram(a=frame_ext_df["EXT"], bins=ext_bin_count, range=(ext_start, ext_end),
                                       density=True)
ext_pdf_df = pd.DataFrame()

ext_col_label = (COMMENT_TOKEN if comment_output_header else "") + "EXT"
ext_pdf_df[ext_col_label] = ext_bin_edges[0:-1]
ext_pdf_df["PDF"] = ext_hist

meta_info_str = (
    f"{COMMENT_TOKEN} INPUT Frame vs Extension file: \"{frame_vs_ext_file}\"\n"
    f"{COMMENT_TOKEN} INPUT Frame Range => frame_index_start: {frame_index_start}  |  frame_index_end: {frame_index_end}\n"
    f"{COMMENT_TOKEN} INPUT Extension Range => ext_start: {ext_start}  |  ext_end: {ext_end}\n"
    f"{COMMENT_TOKEN} INPUT Extension Bins => Count: {ext_bin_count}  | Bin Size: {_ext_bin_size}\n"
    f"{COMMENT_TOKEN} INPUT Moving Average Window Size: {_rolling_win_bins} extension bins (or {_rolling_win_size} extension units)\n"
    f"{COMMENT_TOKEN} ---------------------------------------\n"
    )

if output_pdf_data_file:
    with open(output_pdf_data_file, "w") as out_p:
        out_p.write(f"{COMMENT_TOKEN} -------------- Extension Probability Distribution (Probability Density Function - PDF) ----------------\n")
        out_p.write(meta_info_str)
        ext_pdf_df.to_csv(out_p, mode="a", sep="\t", header=True, index=False, index_label=False)

# Rolling Average

ext_pdf_avg_df = None

if _rolling_win_bins > 0:
    ext_pdf_df["PDF_AVG"] = ext_pdf_df["PDF"].rolling(_rolling_win_bins).mean()
    if output_pdf_avg_data_file:
        ext_pdf_avg_df = ext_pdf_df[[ext_col_label, "PDF_AVG"]]
        ext_pdf_avg_df = ext_pdf_avg_df.dropna()

        with open(output_pdf_avg_data_file, "w") as out_p:
            f"{COMMENT_TOKEN} -------------- Moving-Average Extension Probability Distribution (Probability Density Function - PDF) ----------------\n"
            out_p.write(meta_info_str)
            ext_pdf_avg_df.to_csv(out_p, mode="a", sep="\t", header=True, index=False, index_label=False)

# Plot
# plt.stairs(ext_hist, ext_bin_edges, fill=True)
# if pdf_avg_df is not None:
#     plt.plot(ext_pdf_avg_df[ext_col_label], ext_pdf_avg_df["PDF_AVG"])
#
# plt.xlabel("Extension (Å)")
# plt.ylabel("Relative Population")
# plt.title("Extension Distribution")

has_time = frame_step_fs > 0
if has_time:
    frame_ext_df["TIME_NS"] = frame_ext_df["FRAME"] * (frame_step_fs * 1e-6)

w, h = figaspect(9 / 23)
fig, axes = plt.subplots(1, 2, figsize=(w * 1.4, h * 1.4))
fig.tight_layout(pad=5.0)

axes[0].plot(frame_ext_df["TIME_NS" if has_time else "FRAME"], frame_ext_df["EXT"])
axes[0].set_title("Extension vs " + ("Time" if has_time else "Frame"))
axes[0].set_xlabel("Time (ns)" if has_time else "Frame")
axes[0].set_ylabel("Extension (Å)")

axes[1].stairs(ext_hist, ext_bin_edges, fill=True, label=f"PDF ({ext_bin_count} bins, {_ext_bin_size} Å/bin)")
if ext_pdf_avg_df is not None:
    axes[1].plot(ext_pdf_avg_df[ext_col_label], ext_pdf_avg_df["PDF_AVG"],
                 label=f"Moving Avg ({_rolling_win_bins} bins, {round(_rolling_win_size, 2)} Å)")
    axes[1].legend(loc="upper right")

axes[1].set_title("Extension Distribution")
axes[1].set_xlabel("Extension (Å)")
axes[1].set_ylabel("Relative Population")

if output_fig_file:
    plt.savefig(output_fig_file)
plt.show()
