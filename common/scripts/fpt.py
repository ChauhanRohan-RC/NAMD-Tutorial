import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy

"""
Script to calculate First Passage Time (FPT) distribution from Simulation Trajectory

Theoretically,
    FPT(x0, t) = -d/dt { area under P(x, t|x0) in range x=x_a to x=x_b }

    where 1. P(x, t|x0) is the Conditional Probability Density function (PDF)
          2. x_a and x_b are position of potential minima's (stable states)
          
From simulation data i.e frame_vs_extension:
    FPT(frame=f) = -ve gradient w.r.t to frame at frame=f { area under ( extension_pdf upto frame=f ) in range [x_a, x_b] }
    
## Usage
# 1. Copy script to working dir
# 2. INPUT: set EXTENSION_VS_FRAME file (ext over time)
# 3. INPUT: set x_a, x_b, frame_range etc parameters
# 4. OUTPUT: set output data and plot file names
# 4. run with "python fpt.py"
# 5. Creates output file "fpt.csv", "fpt.pdf"
UNITS: Distance (Å)
"""

COMMENT_TOKEN = "#"

## INPUT -----------------------
frame_vs_ext_file = "dist_vs_frame.dat"
frame_step_fs: float = 2 * 100  # time (in fs) between frames = time_step (fs) * dcd_freq. -1 for NOT_DEFINED

# Frame Range (optional). NOTE: frame_index = time_fs / frame_step_fs
# 6.31e-9

frame__start: int = int(6.31e-9 / (frame_step_fs * 1e-15))  # Inclusive [-1 for no start bound]
frame__end: int = int(32e-9 / (frame_step_fs * 1e-15))  # Exclusive [-1 for no end bound]

# Histogram parameters
x_a: float = 28.2  # TODO: LEFT Boundary (Å)
x_b: float = 37.7  # TODO: right Boundary (Å)
ext_bin_count: int = 100
_ext_bin_size: float = (x_b - x_a) / ext_bin_count

normalize_pdf = True

## Output --------------------------------------------------------------
# TODO: Frames between successive FPT calculation
fpt_frame_step: int = int(0.5e-9 / (frame_step_fs * 1e-15))

normalize_fpt = True

output_fpt_data_file = "fpt_traj-2.csv"
output_fpt_fig_file = "fpt_traj-2.pdf"


## ----------------------------------------------------------------------------

def agg_forward(arr: np.ndarray, agg_func, out_dtype=np.float32) -> np.ndarray:
    samples = len(arr) - 1
    res = np.zeros(samples, dtype=out_dtype)

    for i in range(samples):
        res[i] = agg_func(arr[i], arr[i + 1])
    return res


def agg_forward_mean(arr: np.ndarray) -> np.ndarray:
    return agg_forward(arr, lambda a, b: (a + b) / 2)


# Dataframe
frame_ext_df: pd.DataFrame = pd.read_csv(frame_vs_ext_file, sep=r"\s+", comment=COMMENT_TOKEN, names=("FRAME", "EXT"))
if frame__start >= 0:
    frame_ext_df = frame_ext_df[frame_ext_df["FRAME"] >= frame__start]

if frame__end >= 0:
    frame_ext_df = frame_ext_df[frame_ext_df["FRAME"] < frame__end]

frame_count = len(frame_ext_df["FRAME"])
print(f"Frame Count: {frame_count}")

fpt_sample_count = int(frame_count // fpt_frame_step)
print(f"FPT Sam[les: {fpt_sample_count}")

# ext_series = frame_ext_df["EXT"].values

_start = frame__start if frame__start > 0 else 0
frame_instants_arr = _start + np.array([(i + 1) * fpt_frame_step for i in range(fpt_sample_count)])
area_arr = np.zeros((fpt_sample_count,), dtype=np.float128)

for i in range(fpt_sample_count):
    _frame_instant = frame_instants_arr[i]
    _ext = frame_ext_df[frame_ext_df["FRAME"] < _frame_instant]["EXT"].values

    # Histogram
    ext_hist, ext_bin_edges = np.histogram(a=_ext,
                                           bins=ext_bin_count,
                                           # range=(ext_start, ext_end),
                                           density=normalize_pdf)

    xa_bin = np.searchsorted(ext_bin_edges, x_a) - 1
    xb_bin = np.searchsorted(ext_bin_edges, x_b) - 1
    ext_vals = agg_forward_mean(ext_bin_edges)

    # plt.stairs(ext_hist, ext_bin_edges, label=f"{frame_end * frame_step_fs * 1e-6:g} ns")

    area = scipy.integrate.trapezoid(y=ext_hist[xa_bin: xb_bin + 1], x=ext_vals[xa_bin: xb_bin + 1])
    area_arr[i] = area

# plt.legend(loc="upper right")
# plt.savefig("fpt-1-test.pdf")
# plt.show()

# plt.plot(frame_arr, area_arr)
# plt.scatter(frame_arr, area_arr)
# plt.show()

fpt_arr = -np.gradient(area_arr, frame_instants_arr)
if normalize_fpt:
    fpt_arr -= np.min(fpt_arr)
    fpt_arr /= np.sum(fpt_arr)

# plt.plot(frame_arr, fpt_arr)
# plt.scatter(frame_arr, fpt_arr)
# plt.show()

# Output Data ----------------------------
time_arr = None
if frame_step_fs > 0:
    time_arr = frame_instants_arr * frame_step_fs * 1e-15

res_df = pd.DataFrame()
res_df["FRAME"] = frame_instants_arr
if time_arr is not None:
    res_df["TIME"] = time_arr
res_df["PDF_AREA"] = area_arr
res_df["FPTD"] = fpt_arr

if output_fpt_data_file:
    comments = [
        "-------------------- First Passage Time Distribution (FPTD) from Simulation Trajectory -------------------",
        f"INPUT frame_vs_ext file: \"{frame_vs_ext_file}\"",
        f"INPUT Frame Range: [{frame__start}, {frame__end}] | Frame Count: {frame_count} | Time b/w Frames: {f'{frame_step_fs} fs' if frame_step_fs > 0 else '<Not-Set>'}",
        f"INPUT x_a: {x_a} | x_b: {x_b} | x_bins: {ext_bin_count} | x_bin_size: {_ext_bin_size} | Normalize PDF: {normalize_pdf}",
        "-----------------------------------------------",
        f"OUTPUT -> FPT Frame Step: {fpt_frame_step} | FPT samples: {fpt_sample_count} | Normalize FPT: {normalize_fpt}",
        "-----------------------------------------------"
    ]

    with open(output_fpt_data_file, "w") as f:
        for c in comments:
            f.write(f"{COMMENT_TOKEN} {c}\n")

        res_df.to_csv(f, sep="\t", mode="a", header=True, index=False, index_label=False)

# Plot ------------------------------------
plot_x = time_arr if time_arr is not None else frame_instants_arr
plot_x_label = "$t$ (s)" if time_arr is not None else "Frame"

plot_y = fpt_arr
plot_y_label = "$P_{fpt}(t)$"

plt.scatter(plot_x, plot_y)
plt.plot(plot_x, plot_y)
plt.xlabel(plot_x_label)
plt.ylabel(plot_y_label)
plt.suptitle("First Passage Time Distribution")
plt.title(f"using Simulation Trajectory ($x_a$: {x_a:.2f}, $x_b$: {x_b:.2f})")

if output_fpt_fig_file:
    plt.savefig(output_fpt_fig_file)
plt.show()
