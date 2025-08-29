#!/usr/bin/env python3

import os
import traceback

import numpy as np
import pandas as pd
import plotext


##########################################################################
## General script to plot Probability Density Function (PDF)            ##
##########################################################################

## Usage
# 1. Copy script to working dir
# 2. INPUT: set data_file and column names
# 3. INPUT: set HISTOGRAM and ROLLING_AVG params
# 4. run with "python pdf.py"

Kb = 1.9872036e-3   # Kb (kcal/mol)
T = 300             # Temperature (K)

# =======================================================
# INPUT
# =======================================================
data_file = "../energy-all.csv"     # Input data file
col_name_x = "TS"                   # Input column X
col_name_y = "ELECT"                # Input column Y

# Scaling
scale_x = 1
scale_y = 1 / (Kb * T)

# Output after scaling
out_col_name_x = col_name_x     # Useful after scaling
out_col_name_y = "ELECT/KbT"    # Useful after scaling
plot_label_x = "Time (fs)"      # After scaling
plot_label_y = "Elec Energy (KbT)"  # After scaling

# Data Range after scaling (OPTIONAL)
x_start: float = None  # [OPTIONAL] Inclusive ("None" for no start bound)
x_end: float = None    # [OPTIONAL] Exclusive ("None" for no end bound)

# Histogram of Y column (after scaling)
y_hist_start: float = None        # [OPTIONAL]
y_hist_end: float = None          # [OPTIONAL]
hist_bin_count: int = 1000

# Moving average
y_hist_avg_bins: int = 20    # [int] Moving Average window size (in no of bins). Set to 0 to disable
y_hist_avg_size: float = -1  # [float] [Only used if rolling_window_bins is not set]. Moving Average window size (in unit of y). Set to 0 to disable.


# =======================================================
# OUTPUT
# =======================================================
out_prefix = "erg-elec-all"

out_scaled_data_file = f"{out_prefix}.scaled.csv"       # ONLY IF scale_x or scale_y != 1. "" for None
out_pdf_file = f"{out_prefix}.pdf.csv"
out_fig_file = f"{out_prefix}.pdf.pdf"

comment_output_header = False
plot_in_terminal = True

# -------------------------------------------------------------------------

# =======================================================
# MAIN
# =======================================================

## Constants
COMMENT_TOKEN = "#"
COL_NAME_PDF = "PDF"                       # Probability Distribution Function (PDF) of Distance
COL_NAME_PDF_AVG_PREFIX = "PDF_AVG"        # Average PDF over a window

## Dataframe
df: pd.DataFrame = pd.read_csv(data_file, sep=r"\s+", comment=COMMENT_TOKEN)

# -------------------------------------------------
# Scaling
# -------------------------------------------------
if scale_x != 1 or scale_y != 1:
    if scale_x != 1:
        df[col_name_x] *= scale_x
    if scale_y != 1:
        df[col_name_y] *= scale_y

    if out_scaled_data_file:
        meta_info_str = (
            f"{COMMENT_TOKEN} -------------- Scaled Data File ----------------\n"
            f"{COMMENT_TOKEN} INPUT data file: \"{data_file}\"\n"
            f"{COMMENT_TOKEN} INPUT Columns => X: \"{col_name_x}\"  |  Y: \"{col_name_y}\"\n"
            f"{COMMENT_TOKEN} INPUT Scale => X: {scale_x}  |  Y: {scale_y}\n"
            f"{COMMENT_TOKEN} PLOT Labels (after scaling) => X: \"{plot_label_x}\"  |  Y: \"{plot_label_y}\"\n"
            f"{COMMENT_TOKEN} ---------------------------------------\n"
        )

        scaled_df = pd.DataFrame()
        scaled_df[out_col_name_x] = df[col_name_x]
        scaled_df[out_col_name_y] = df[col_name_y]

        with open(out_scaled_data_file, "w") as out_sc:
            out_sc.write(meta_info_str)
            scaled_df.to_csv(out_sc, mode="a", sep="\t", header=True, index=False, index_label=False)


# -------------------------------------------------
# Histogram
# -------------------------------------------------
if x_start is not None:
    df = df[df[col_name_x] >= x_start]

if x_end is not None:
    df = df[df[col_name_x] < x_end]

if y_hist_start is None:
    y_hist_start = df[col_name_y].min()

if y_hist_end is None:
    y_hist_end = df[col_name_y].max()

hist_bin_size = (y_hist_end - y_hist_start) / hist_bin_count

# calculate histogram density
y_hist, y_bin_edges = np.histogram(a=df[col_name_y], bins=hist_bin_count, range=(y_hist_start, y_hist_end),
                                   density=True)
pdf_df = pd.DataFrame()

out_col_y = (COMMENT_TOKEN if comment_output_header else "") + out_col_name_y
pdf_df[out_col_y] = y_bin_edges[0:-1]
pdf_df[COL_NAME_PDF] = y_hist

# -------------------------------------------------
# Histogram Moving Average
# -------------------------------------------------
def parse_hist_avg_win_bins() -> tuple: # tuple[int, float]
    bins = -1

    if y_hist_avg_bins > 0:
        bins = y_hist_avg_bins
    elif y_hist_avg_size > 0:
        bins = int(round(y_hist_avg_size / ((y_hist_end - y_hist_start) / hist_bin_count)))

    return (bins, bins * hist_bin_size) if bins > 0 else (-1, -1)


_hist_avg_win_bins, _hist_avg_win_size = parse_hist_avg_win_bins()
do_pdf_avg: bool = _hist_avg_win_bins > 0

col_name_pdf_avg: str = COL_NAME_PDF_AVG_PREFIX
if do_pdf_avg:
    col_name_pdf_avg = f"{COL_NAME_PDF_AVG_PREFIX}{_hist_avg_win_bins}"  # col_name = PDF_AVG<bins>
    pdf_df[col_name_pdf_avg] = pdf_df[COL_NAME_PDF].rolling(_hist_avg_win_bins).mean()

# -------------------------------------------------
# Output file
# -------------------------------------------------
if out_pdf_file:
    meta_info_str = (
        f"{COMMENT_TOKEN} INPUT data file: \"{data_file}\"\n"
        f"{COMMENT_TOKEN} INPUT Columns => X: \"{col_name_x}\"  | Y: \"{col_name_y}\"\n"
        f"{COMMENT_TOKEN} INPUT Scale => X: {scale_x}  | Y: {scale_y}\n"
        f"{COMMENT_TOKEN} INPUT X Range => [{x_start}, {x_end})\n"
        f"{COMMENT_TOKEN} Histogram Y Range => ({y_hist_start}, {y_hist_end})\n"
        f"{COMMENT_TOKEN} Histogram Bins => Count: {hist_bin_count}  | Bin Size: {hist_bin_size}\n"
        f"{COMMENT_TOKEN} Histogram Moving Average Window Size: {_hist_avg_win_bins} bins (or {_hist_avg_win_size} y units)\n"
        f"{COMMENT_TOKEN} ---------------------------------------\n"
    )

    with open(out_pdf_file, "w") as out_p:
        out_p.write(f"{COMMENT_TOKEN} -------------- Probability Density (PDF) ----------------\n")
        out_p.write(meta_info_str)
        pdf_df.to_csv(out_p, mode="a", sep="\t", header=True, index=False, index_label=False)


# =======================================================
# PLOTS
# =======================================================

# -------------------------------
# PLOT: Terminal
# -------------------------------
if plot_in_terminal:
    # plotext.from_matplotlib(plt.gcf())

    def print_full_line(char='-'):
        """
        Prints a line that fills the entire width of the terminal.
        """
        try:
            terminal_width = os.get_terminal_size().columns
            print(char * terminal_width)
        except OSError:
            # Handle cases where terminal size cannot be determined (e.g., not in a TTY)
            #print("Unable to determine terminal size. Printing a default line.")
            print(char * 80)  # Print a line of 80 characters as a fallback

    plotext.theme("pro")
    plotext.limit_size(True, False)
    plotext.plot_size(plotext.terminal_width(), plotext.terminal_height() * 1.6)
    plotext.subplots(2,1)

    plotext.subplot(2,1)
    plotext.plot(df[col_name_x], df[col_name_y])
    plotext.title(f"{plot_label_y} vs {plot_label_x}")
    plotext.xlabel(plot_label_x)
    plotext.ylabel(plot_label_y)

    # A basic Histogram plot on terminal
    plotext.subplot(1, 1)
    plotext.hist(df[col_name_y], bins=200, fill=True)
    plotext.title("Probability Distribution")
    plotext.xlabel(plot_label_y)
    plotext.ylabel("Relative Population")

    print(""); print_full_line("#"); print("")
    plotext.show()
    print(""); print_full_line("#"); print("")


# -------------------------------------------------
# PLOT: Interactive (Matplotlib)
# -------------------------------------------------
try:
    import matplotlib.pyplot as plt
    from matplotlib.figure import figaspect
except ImportError:
    print("ERROR: Matplotlib is not installed !!!")
    traceback.print_exc()
    exit(1)

w, h = figaspect(9 / 23)
fig, axes = plt.subplots(1, 2, figsize=(w * 1.4, h * 1.4))
fig.tight_layout(pad=5.0)

axes[0].plot(df[col_name_x], df[col_name_y])
axes[0].set_title(f"{plot_label_y} vs {plot_label_x}")
axes[0].set_xlabel(plot_label_x)
axes[0].set_ylabel(plot_label_y)

axes[1].stairs(y_hist, y_bin_edges, fill=True, label=f"PDF (bins: {hist_bin_count}, bin_size: {round(hist_bin_size, 2)})")
if do_pdf_avg:
    pdf_avg_df = pdf_df[[out_col_y, col_name_pdf_avg]]
    pdf_avg_df = pdf_avg_df.dropna()

    axes[1].plot(pdf_avg_df[out_col_y], pdf_avg_df[col_name_pdf_avg],
                 label=f"Mov-Avg (bins: {_hist_avg_win_bins}, bin_size: {round(_hist_avg_win_size, 2)})")

axes[1].set_title("Probability Distribution")
axes[1].set_xlabel(plot_label_y)
axes[1].set_ylabel("Relative Population")
axes[1].legend(loc="best")

# Output Figure file
if out_fig_file:
    plt.savefig(out_fig_file)

# Interactive plot
plt.show()
