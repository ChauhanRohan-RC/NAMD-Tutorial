#!/usr/bin/env python3

import os
import traceback

import numpy as np
import pandas as pd
import plotext


##########################################################################
## General script to plot Probability Density Function (PDF)            ##
##########################################################################
# Search for TODO

# => pouplation of i'th bin = N(i)
# => Probability Density P(i) =  c * N(i)   [c: normalization constant]
# such that
#                   integral P(i) = 1   (total area = 1)
#
# => Boltzmann Inverted PMF
#
#                   U(x)/kBT = -ln[P_eq(i)]    (i = bin index)
#

## Usage
# 1. Copy script to working dir
# 2. INPUT: set data_file and column names
# 3. INPUT: set HISTOGRAM and ROLLING_AVG params
# 4. run with "python pdf.py"

## Constants
COMMENT_TOKEN = "#"
COL_NAME_PDF = "PDF"                       # Probability Distribution Function (PDF) with total area = 1
COL_NAME_PDF_AVG_PREFIX = "PDF_AVG"        # Average PDF over a window
COL_NAME_PMF = "PMF"    # boltzmann inverted PMF

Kb = 1.9872036e-3   # Kb (kcal/mol)
T = 300             # Temperature (K) [ONLY USED in scale_y in very few cases]

# =======================================================
# INPUT
# =======================================================
data_file: str = "energy-all.csv"     # TODO Input data file
col_name_x: str = "TS"                   # TODO Input column X
col_name_y: str = "ELECT"                # TODO Input column Y
input_delimiter: str = r"\s+"            # Input delimiter

# Scaling
scale_x: float = 1
scale_y: float = 1 / (Kb * T)       # 1 / (Kb * T)   # for energies in unit of KbT

# [Optional] Labels after scaling. can be "None" or ""
out_col_name_x: str = col_name_x     # Useful after scaling
out_col_name_y: str = "ELECT/KbT"    # Useful after scaling
plot_label_x: str = "Time (fs)"      # After scaling TODO
plot_label_y: str = "Elec Energy (KbT)"  # After scaling TODO

# Data Range after scaling (OPTIONAL)
x_start: float = None  # [OPTIONAL] Inclusive ("None" for min(x))
x_end: float = None    # [OPTIONAL] Exclusive ("None" for max(x))

# Histogram of Y column (after scaling)
y_hist_start: float = None        # [OPTIONAL] "None" for min(y)
y_hist_end: float = None          # [OPTIONAL] "None" for max(y)
hist_bin_count: int = 10000

# Moving average
y_hist_avg_bins: int = 50    # [int] Moving Average window size (in no of bins). Set to 0 to disable
y_hist_avg_size: float = -1  # [float] [Only used if rolling_window_bins is not set]. Moving Average window size (in unit of y). Set to 0 to disable.

## Boltzmann Inverted PMF: U(x)/kBT = -ln[P_eq(i)]    (i = bin index)
calc_pmf: bool = True       # TODO


# =======================================================
# OUTPUT
# =======================================================
out_prefix: str = "erg-elec-all"     # TODO

out_scaled_data_file: str = f"{out_prefix}.scaled.csv"       # [OPT] ONLY IF scale_x or scale_y != 1. "" for None
out_pdf_file: str = f"{out_prefix}.pdf.csv"     # [OPT] Histogram data file (probability density function - PDF)
out_fig_file: str = f"{out_prefix}.pdf.pdf"     # [OPT] Histogram figure file

out_format: str = "%.4E"
out_delimiter: str = " "

comment_output_header: str = False
plot_in_terminal: bool= True






# =======================================================
# MAIN
# =======================================================

## Dataframe
df: pd.DataFrame = pd.read_csv(data_file,
                               usecols=[col_name_x,col_name_y],
                               sep=input_delimiter,
                               comment=COMMENT_TOKEN,
                               skipinitialspace=True)

# Checks
if not out_col_name_x:
    out_col_name_x = col_name_x

if not out_col_name_y:
    out_col_name_y = col_name_y

if not plot_label_x:
    plot_label_x = out_col_name_x

if not plot_label_y:
    plot_label_y = out_col_name_y


# -------------------------------------------------
# Scaling
# -------------------------------------------------
has_scale: bool = scale_x != 1 or scale_y != 1
if has_scale:
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
            scaled_df.to_csv(out_sc, mode="a", sep=out_delimiter, float_format=out_format, header=True, index=False, index_label=False)


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

if calc_pmf:
    pdf_df[COL_NAME_PMF] = -np.log(pdf_df[COL_NAME_PDF])    # U(x)/KbT = -ln(P_eq(x))

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
        pdf_df.to_csv(out_p, mode="a", sep=out_delimiter, float_format=out_format, header=True, index=False, index_label=False)


# =======================================================
# PLOTS
# =======================================================

# -------------------------------
# PLOT: Terminal
# -------------------------------
pdf_plot_title: str = "Probability Distribution"
pdf_plot_label_x: str = plot_label_y
pdf_plot_label_y: str = "Probability Density"

pmf_plot_title: str = "Boltzmann Inverted PMF\n $U(x)/K_{B}T = -ln(P_{eq}(i))$"
pmf_plot_label_x: str = plot_label_y
pmf_plot_label_y: str = "$U(x)$ $[K_{B}T]$"

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

    # A input data plot on terminal
    plotext.subplot(1,1)
    plotext.plot(df[col_name_x], df[col_name_y])
    plotext.title(f"{plot_label_y} vs {plot_label_x}")
    plotext.xlabel(plot_label_x)
    plotext.ylabel(plot_label_y)

    # A basic Histogram plot on terminal
    plotext.subplot(2, 1)
    plotext.hist(df[col_name_y], bins=200, fill=True)
    plotext.title(pdf_plot_title)
    plotext.xlabel(pdf_plot_label_x)
    plotext.ylabel(pdf_plot_label_y)

    # Boltzmann Inverted PMF plot on terminal
    # plotext.subplot(3, 1)
    # plotext.plot(pdf_df[col_name_y], pdf_df[COL_NAME_PMF])
    # plotext.title(pmf_plot_title)
    # plotext.xlabel(pmf_plot_label_x)
    # plotext.ylabel(pmf_plot_label_y)

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

w, h = figaspect((9/38) if calc_pmf else (9/24))
fig, axes = plt.subplots(1, 3 if calc_pmf else 2, figsize=(w * 1.4, h * 1.4))
fig.tight_layout(pad=3.0)

# x vs y plot
axes[0].plot(df[col_name_x], df[col_name_y])
axes[0].set_title(f"{plot_label_y} vs {plot_label_x}")
axes[0].set_xlabel(plot_label_x)
axes[0].set_ylabel(plot_label_y)

# hist(y) plot
axes[1].stairs(y_hist, y_bin_edges, fill=False, label=f"PDF (bins: {hist_bin_count}, bin_size: {round(hist_bin_size, 2)})")
if do_pdf_avg:
    pdf_avg_df = pdf_df[[out_col_y, col_name_pdf_avg]]
    pdf_avg_df = pdf_avg_df.dropna()

    axes[1].plot(pdf_avg_df[out_col_y], pdf_avg_df[col_name_pdf_avg],
                 label=f"Mov-Avg (bins: {_hist_avg_win_bins}, bin_size: {round(_hist_avg_win_size, 2)})")

axes[1].set_title(pdf_plot_title)
axes[1].set_xlabel(pdf_plot_label_x)
axes[1].set_ylabel(pdf_plot_label_y)
axes[1].legend(loc="best")

# PMF plot
if calc_pmf:
    axes[2].plot(pdf_df[col_name_y], pdf_df[COL_NAME_PMF])
    axes[2].set_title(pmf_plot_title)
    axes[2].set_xlabel(pmf_plot_label_x)
    axes[2].set_ylabel(pmf_plot_label_y)


# Output Figure file
if out_fig_file:
    plt.savefig(out_fig_file)

print("\n============= Probability Density ================")
if has_scale and out_scaled_data_file:
    print(f"=> Scaled input data file: {out_scaled_data_file}")

if out_pdf_file:
    print(f"=> OUTPUT data file    :  {out_pdf_file}")
if out_fig_file:
    print(f"=> OUTPUT figure file  :  {out_fig_file}")
print("=================================================\n")

# Interactive plot
plt.show()

