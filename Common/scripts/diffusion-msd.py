import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.figure import figaspect

"""
=====================================================================
Script to find self-Diffusion coefficient from MSD vs time (frame)
===================================================================== 
In general, 
                MSD ~ t^a   =>   log(MSD) ~ a log(t) 
  
In DIFFUSIVE REGIME, a = 1
                
                MSD = (1/2d) D * t          D: self-diffusivity (m^2/s), d: dimensions  

-------------------- 
PROCEDURE
--------------------

1. Plot MSD vs time (or frame) and the corresponding log-log plot
   
    =>  plot_msd(rmsd_file, time_fs_per_frame, out_fig_file=<out_fig_file>)
 
2. Divide the log-log plot into regions of straight line fits

    # Manual Bounds
    => fit_win_bounds_log = np.log( [ (frame_start1, frame_stop1), (frame_start2, frame_stop2), ... ] ) 
    OR
    # Auto bounds 
    => fit_win_count = 10   # divide log space in sequal width regions
    
    => fit_msd_loglog_linear(rmsd_file=rmsd_file,
                              fit_window_bounds_log=fit_window_bounds_log,
                              fit_window_count=fit_window_count,
                              time_fs_per_frame=time_fs_per_frame,
                              out_file_prefix=<prefix>)
 
3. Look at csv output of loglog plot linear-fit to find the slope closest to 1
4. Look at csv output MSD vs time linear fit to find Diffusion coefficient in that region

"""

# CONSTANTS

COMMENT_TOKEN = "#"

COL_FRAME = "FRAME"
COL_RMSD = "RMSD"
COL_MSD = "MSD"
COL_TIME_PS = "TIME_PS"
COL_FRAME_LOG = "FRAME_LOG"
COL_MSD_LOG = "MSD_LOG"
COL_TIME_PS_LOG = "TIME_PS_LOG"

COL_FRAME_START = "FRAME_START"
COL_FRAME_END = "FRAME_END"
COL_TIME_PS_START = "TIME_PS_START"
COL_TIME_PS_END = "TIME_PS_END"
COL_FIT_R_SQ = "R_SQ"
COL_FIT_INTERCEPT = "FIT_INTERCEPT"
COL_FIT_SLOPE = "FIT_SLOPE"
COL_FIT_SELF_DIFF_COEFF = "DIFF_COEFF"


def load_df(file, sep=r"\s+", x_column: str = None, x_start=None, x_end=None,
            positive_bounds: bool = False) -> pd.DataFrame:
    df = pd.read_csv(file, sep=sep, comment=COMMENT_TOKEN)
    if x_column:
        if x_start is not None and (not positive_bounds or x_start >= 0):
            df = df[df[x_column] >= x_start]
        if x_end is not None and (not positive_bounds or x_end >= 0):
            df = df[df[x_column] < x_end]

    return df


def plot_msd(rmsd_file,  # file or path buffer
             col_frame: str = COL_FRAME,
             col_rmsd: str = COL_RMSD,
             time_fs_per_frame: float = -1,  # time b/w frames (in fs) = dcd-frequency
             frame_start: int = None,
             frame_end: int = None,
             separator: str = r"\s+",
             plot_msd_vs_t: bool = True,
             plot_loglog: bool = True,
             uniform_logspace: bool = False,  # use plt.logspace(x,y) instead of plt.plot(np.log(x), np.log(y))
             out_fig_file: str = None):
    """
    Plot MSD vs Frame (or Time) and the corresponding Log-Log plot

    The log-log plot is useful in identifying diffusive regime
    In general,
                    MSD ~ t^a   =>   log(MSD) ~ a log(t)

    In DIFFUSIVE REGIME, a = 1, so the log-log plot will give a line with slope = 1
    in the diffusive region of time.

    """

    df = load_df(rmsd_file, sep=separator, x_column=col_frame, x_start=frame_start, x_end=frame_end,
                 positive_bounds=True)
    _frame_arr = df[col_frame]

    has_time = time_fs_per_frame is not None and time_fs_per_frame > 0
    if has_time:
        _time_ps_arr = _frame_arr * time_fs_per_frame * 1e-3  # Time in ps
        x_arr = _time_ps_arr
    else:
        x_arr = _frame_arr

    _rmsd_arr = df[col_rmsd]
    _msd_arr = np.square(_rmsd_arr)

    # Plotting

    def __do_plot(plt_or_axis, is_plt: bool,
                  x: np.ndarray, y: np.ndarray,
                  title: str, x_label: str, y_label: str,
                  plot_loglog: bool = False, series_label: str = ""):
        method = plt_or_axis.loglog if plot_loglog else plt_or_axis.plot
        method(x, y, label=series_label)

        if is_plt:
            plt_or_axis.title(title)
            plt_or_axis.xlabel(x_label)
            plt_or_axis.ylabel(y_label)
        else:
            plt_or_axis.set_title(title)
            plt_or_axis.set_xlabel(x_label)
            plt_or_axis.set_ylabel(y_label)

    title = f"MSD vs {'Time' if has_time else 'Frame'}"
    xlabel = "Time / ps" if has_time else "Frame"
    ylabel = "MSD / $Å^2$"

    if plot_loglog:
        loglog_x = x_arr if uniform_logspace else np.log(x_arr)
        loglog_y = _msd_arr if uniform_logspace else np.log(_msd_arr)

        loglog_title = "Log-Log Plot"
        _pre = 'log' if uniform_logspace else 'ln'

        loglog_xlabel = f"{_pre}({xlabel})"
        loglog_ylabel = f"{_pre}({ylabel})"

    w, h = figaspect(9 / 23)
    if plot_msd_vs_t and plot_loglog:
        fig, axes = plt.subplots(1, 2, figsize=(w * 1.4, h * 1.4))
        fig.tight_layout(pad=5.0)

        __do_plot(axes[0], is_plt=False, x=x_arr, y=_msd_arr, title=title, x_label=xlabel,
                  y_label=ylabel)

        __do_plot(axes[1], is_plt=False, x=loglog_x, y=loglog_y, title=loglog_title, x_label=loglog_xlabel,
                  y_label=loglog_ylabel, plot_loglog=uniform_logspace)
    else:
        plt.figure(figsize=(w * 1.4, h * 1.4))
        if plot_msd_vs_t:
            __do_plot(plt, is_plt=True, x=x_arr, y=_msd_arr, title=title, x_label=xlabel,
                      y_label=ylabel)
        elif plot_loglog:
            __do_plot(plt, is_plt=True, x=loglog_x, y=loglog_y, title=loglog_title, x_label=loglog_xlabel,
                      y_label=loglog_ylabel, plot_loglog=uniform_logspace)

    # if plotted somthing
    if out_fig_file and (plot_msd_vs_t or plot_loglog):
        plt.savefig(out_fig_file)
    plt.show()


def fit_msd_loglog_linear(rmsd_file,  # file or path buffer
                          # list of (frame_log_start, frame_log_end), both bounds are inclusive
                          fit_window_bounds_log: list[tuple[int, int]] | None,
                          # [Optional] num of fit windows. Used to generate frame_ranges if not specified
                          fit_window_count: int | None = 1,
                          # [Optional] num of frames in each fit window. Used to generate frame_ranges if not specified
                          fit_window_size_log: int | None = None,
                          col_frame: str = COL_FRAME,
                          col_rmsd: str = COL_RMSD,
                          frame_start: int | None = None,
                          frame_end: int | None = None,
                          time_fs_per_frame: float = -1,  # time b/w frames (in fs) = dcd-frequency
                          dimension: int = 3,
                          separator: str = r"\s+",
                          out_file_prefix: str = "msd",
                          out_fig_file_format: str = "pdf"):
    """
    Linear fit of the log(MSD) vs log(time)

    Finds the region with slope closest to 1 in loglog plot => MSD ~ t^1 (BEST DIFFUSIVE REGION)

    Finds the self-Diffusion coefficient from the MSD vs t plot in the above region
    """

    # Explicitly given frame ranges
    has_win_bounds = fit_window_bounds_log is not None and len(fit_window_bounds_log) > 0
    has_win_size = fit_window_size_log is not None and fit_window_size_log > 0
    has_win_count = fit_window_count is not None and fit_window_count > 0
    if not (has_win_bounds or has_win_size or has_win_count):
        raise ValueError("ERR: No bounds / size / count of fit window provided.")

    # Loading Data
    df = load_df(rmsd_file, sep=separator, x_column=col_frame, x_start=frame_start, x_end=frame_end,
                 positive_bounds=True)

    col_frame_log = COL_FRAME_LOG
    col_time_ps = COL_TIME_PS
    col_time_ps_log = COL_TIME_PS_LOG
    col_msd = COL_MSD
    col_msd_log = COL_MSD_LOG

    log_func = np.log  # Log function = natural log
    antilog_func = np.exp
    log_label = "ln"

    df[col_msd] = np.square(df[col_rmsd])
    df[col_frame_log] = log_func(df[col_frame])
    df[col_msd_log] = log_func(df[col_msd])

    has_time = time_fs_per_frame is not None and time_fs_per_frame > 0
    if has_time:
        df[col_time_ps] = df[col_frame] * time_fs_per_frame * 1e-3  # Time in ps
        col_x = col_time_ps
    else:
        col_x = col_frame

    # Dropping inf, -inf, NA
    df.replace([np.inf, -np.inf], np.nan, inplace=True)
    df.dropna(inplace=True)

    # Generating frame ranges if not explicitly given
    if not has_win_bounds:
        frame_log_min = np.min(df[col_frame_log])
        frame_log_max = np.max(df[col_frame_log])

        if has_win_count:
            win_size = (frame_log_max - frame_log_min) / fit_window_count
        elif has_win_size:
            win_size = fit_window_size_log
        else:
            raise AssertionError("No FIT Window Size. Should not happen!!")

        _range = np.arange(frame_log_min, frame_log_max, win_size)

        fit_window_bounds_log = [(i, i + win_size) for i in _range]
        # frame_min = np.min(df[col_frame])
        # frame_max = np.max(df[col_frame])
        #
        # if has_win_count:
        #     win_size = (frame_max - frame_min + 1) // fit_window_count
        # elif has_win_size:
        #     win_size = fit_window_size
        # else:
        #     raise AssertionError("No FIT Window Size. Should not happen!!")
        #
        # _range = list(range(frame_min, frame_max + 1, win_size - 1))
        # fit_window_bounds = [(i, min(i + win_size - 1, frame_max)) for i in _range[:-1]]

    print(f"FIT WINDOW BOUNDS: {fit_window_bounds_log}")

    # Linear Fitting of Log-Log plot ---------------------------------
    log__line_fits = []  # list of (fit_intercept, fit_slope, fit_r_sq)
    # diff_coeffs = []
    log__line_xy = []  # list of (line_x, line_y)
    for _frame_log_start, _frame_log_end in fit_window_bounds_log:
        dfi = df[df[col_frame_log].between(_frame_log_start, _frame_log_end)]

        # linear Least-Square Fit
        coeffs, _stat = np.polynomial.polynomial.polyfit(dfi[col_frame_log], dfi[col_msd_log], deg=1, full=True)
        log__line_fits.append(
            (coeffs[0], coeffs[1], _stat[0][0]))  # (intercept, slope, R^2)   R^2 is the residual square

        # # self diffusion coefficient = (1/2d) * slope of MSd vs t (in diffusive regime MSD ~ t)
        # diff_coeff = coeffs[1] / (2 * dimension)  # in Å^2/frame or Å^2/ps
        # diff_coeffs.append(diff_coeff)

        line_x = dfi[col_frame_log]
        line_y = coeffs[0] + (coeffs[1] * line_x)
        log__line_xy.append((line_x, line_y))

    log_fit_df = pd.DataFrame()
    log_fit_df[COL_FRAME_START] = [round(antilog_func(i[0])) for i in fit_window_bounds_log]
    log_fit_df[COL_FRAME_END] = [round(antilog_func(i[1])) for i in fit_window_bounds_log]
    log_fit_df[COL_FIT_R_SQ] = [i[2] for i in log__line_fits]
    log_fit_df[COL_FIT_INTERCEPT] = [i[0] for i in log__line_fits]
    log_fit_df[COL_FIT_SLOPE] = [i[1] for i in log__line_fits]
    # fit_df[COL_FIT_SELF_DIFF_COEFF] = diff_coeffs

    # Finding the slope closest to 1
    closest_idx = np.argmin(np.abs(log_fit_df[COL_FIT_SLOPE] - 1))  # Index of slope closest to 1
    print("\n==========================================================")
    print(f'BEST DIFFUSIVE BEHAVIOUR: FIT {closest_idx} (slope closest to 1 in loglog plot)')
    print("============================================================\n")

    # Info
    _info_loglog = (f"{COMMENT_TOKEN} =================================================\n"
                    f"{COMMENT_TOKEN} Linear Least-Square fit of log(MSD) vs log(Frame) \n"
                    f"{COMMENT_TOKEN} =================================================\n"
                    f"{COMMENT_TOKEN} INPUT Time between Frames: {str(time_fs_per_frame) + ' fs' if has_time else 'NOT-SPECIFIED'} \n"
                    f"{COMMENT_TOKEN} OUTPUT Slope closest to 1 (diffusive regime): FIT {closest_idx}\n"
                    f"{COMMENT_TOKEN}                => Slope: {log_fit_df[COL_FIT_SLOPE].iloc[closest_idx]:.2E} \n"
                    f"{COMMENT_TOKEN}                => Frame Range: [{log_fit_df[COL_FRAME_START].iloc[closest_idx]}, {log_fit_df[COL_FRAME_END].iloc[closest_idx]}] \n"
                    f"{COMMENT_TOKEN} -------------\n"
                    f"{COMMENT_TOKEN} NOTE: FIT columns are linear least-square fits of log(MSD) vs log(Frame) in the corresponding frame ranges\n"
                    f"{COMMENT_TOKEN}      => {COL_FIT_R_SQ}: R-Square (sum of residual square)\n"
                    f"{COMMENT_TOKEN} ------------------------------------------------------\n")

    with open(f"{out_file_prefix}.loglog-fit.csv", 'w') as _log_f:
        _log_f.write(_info_loglog)
        log_fit_df.to_csv(_log_f, sep='\t', mode='a', index=True, index_label="FIT")

    # Linear Fitting of MSD vs frame/time ---------------------------------
    line_fits = []
    line_xy = []
    diff_coeffs = []
    for _frame_start, _frame_end in zip(log_fit_df[COL_FRAME_START], log_fit_df[COL_FRAME_END]):
        dfi = df[df[col_frame].between(_frame_start, _frame_end)]

        # linear Least-Square Fit
        coeffs, _stat = np.polynomial.polynomial.polyfit(dfi[col_x], dfi[col_msd], deg=1, full=True)
        line_fits.append(
            (coeffs[0], coeffs[1], _stat[0][0]))  # (intercept, slope, R^2)   R^2 is the residual square

        # self diffusion coefficient = (1/2d) * slope of MSd vs t (in diffusive regime MSD ~ t)
        diff_coeff = coeffs[1] / (2 * dimension)  # in Å^2/frame or Å^2/ps
        diff_coeffs.append(diff_coeff)

        line_x = dfi[col_x]
        line_y = coeffs[0] + (coeffs[1] * line_x)
        line_xy.append((line_x, line_y))

    fit_df = pd.DataFrame()
    fit_df[COL_FRAME_START] = log_fit_df[COL_FRAME_START]
    fit_df[COL_FRAME_END] = log_fit_df[COL_FRAME_END]
    if has_time:
        fit_df[COL_TIME_PS_START] = fit_df[COL_FRAME_START] * time_fs_per_frame * 1e-3  # in ps
        fit_df[COL_TIME_PS_END] = fit_df[COL_FRAME_END] * time_fs_per_frame * 1e-3  # in ps
    fit_df[COL_FIT_R_SQ] = [i[2] for i in line_fits]
    fit_df[COL_FIT_INTERCEPT] = [i[0] for i in line_fits]
    fit_df[COL_FIT_SLOPE] = [i[1] for i in line_fits]
    fit_df[COL_FIT_SELF_DIFF_COEFF] = diff_coeffs

    # Info
    _info = (f"{COMMENT_TOKEN} =================================================\n"
             f"{COMMENT_TOKEN} Linear Least-Square fit of MSD vs {'Time' if has_time else 'Frame'} \n"
             f"{COMMENT_TOKEN} =================================================\n"
             f"{COMMENT_TOKEN} Note: Self-Diffusion coefficient D = (1/2d) * Slope of MSD vs t (in diffusive regime MSD ~ t^1 => slope ~ 1 in loglog plot) \n"
             f"{COMMENT_TOKEN} ------------------------\n"
             f"{COMMENT_TOKEN} INPUT Time between Frames: {str(time_fs_per_frame) + ' fs' if has_time else 'NOT-SPECIFIED'} \n"
             f"{COMMENT_TOKEN} OUTPUT Best Diffusive Behaviour (Slope closest to 1 in loglog plot): FIT {closest_idx}\n"
             f"{COMMENT_TOKEN}                => Frame Range: [{fit_df[COL_FRAME_START].iloc[closest_idx]}, {fit_df[COL_FRAME_END].iloc[closest_idx]}] \n"
             f"{COMMENT_TOKEN}                => Slope in loglog plot: {log_fit_df[COL_FIT_SLOPE].iloc[closest_idx]:.2E} \n"
             f"{COMMENT_TOKEN}                => Slope in MSD vs {'Time' if has_time else 'Frame'}: {fit_df[COL_FIT_SLOPE].iloc[closest_idx]:.2E} \n"
             f"{COMMENT_TOKEN}                => Self-Diffusion Coefficient (D): {fit_df[COL_FIT_SELF_DIFF_COEFF][closest_idx]:.4E} Å^2/{'ps' if has_time else 'frame'} \n"
             f"{COMMENT_TOKEN} ------------------------\n"
             f"{COMMENT_TOKEN} NOTE: FIT columns are linear least-square fits of MSD vs {'Time' if has_time else 'Frame'} in the corresponding frame ranges\n"
             f"{COMMENT_TOKEN}      => {COL_FIT_R_SQ}: R-Square (sum of residual square)\n"
             f"{COMMENT_TOKEN}      => {COL_FIT_SELF_DIFF_COEFF}: Self Diffusion Coefficient (in Å^2/{'ps' if has_time else 'frame'})\n"
             f"{COMMENT_TOKEN} ------------------------------------------------------\n")

    with open(f"{out_file_prefix}.fit.csv", 'w') as _f:
        _f.write(_info)
        fit_df.to_csv(_f, sep='\t', mode='a', index=True, index_label="FIT")

    # Plot
    w, h = figaspect(9 / 26)
    fig, axes = plt.subplots(1, 2, figsize=(w * 1.6, h * 1.6))
    fig.tight_layout(pad=5)
    fig.suptitle(f"Linear fit of MSD Log-Log Plot to find the Diffusive Region")

    axes[0].plot(df[col_x], df[col_msd])  # MSD vs frame/time
    for i in range(len(line_xy)):
        axes[0].plot(line_xy[i][0], line_xy[i][1],
                     label=f"Fit {i} (m: {line_fits[i][1]:.2E}, c: {line_fits[i][0]:.2E}, $R^2$: {line_fits[i][2]:.2E})")

    axes[0].set_title(f"MSD vs {'Time' if has_time else 'Frame'} Linear Fit")
    axes[0].set_xlabel(f"{'Time / fs' if has_time else 'Frame'}")
    axes[0].set_ylabel(f"MSD / $Å^2$")
    axes[0].legend(loc="best")

    axes[1].plot(df[col_frame_log], df[col_msd_log])  # log(MSD) vs log(Frame)
    for i in range(len(log__line_xy)):
        axes[1].plot(log__line_xy[i][0], log__line_xy[i][1],
                     label=f"Fit {i} (m: {log__line_fits[i][1]:.2E}, c: {log__line_fits[i][0]:.2E}, $R^2$: {log__line_fits[i][2]:.2E})")

    axes[1].set_title(f"{log_label}(MSD) vs {log_label}(Frame) Linear Fit")
    axes[1].set_xlabel(f"{log_label}(Frame)")
    axes[1].set_ylabel(f"{log_label}(MSD / $Å^2$)")
    axes[1].legend(loc="best")

    if out_file_prefix and out_fig_file_format:
        plt.savefig(f"{out_file_prefix}.fit.{out_fig_file_format}")
    plt.show()


# ----------------------------------------------------------------------------------------------------

# ===========================
# INPUT
# ===========================
rmsd_file = "../msd/rmsd.csv"
col_frame = "FRAME"
col_rmsd = "RMSD"

time_fs_per_frame = 200  # (in fs) time b/w frames = dcd-frequency
frame_start = None  # None for no bound
frame_end = None  # None for no bound

# Explicit Fit Window bounds [in Frames] (None for no explicit bounds)
fit_win_bounds = [(665, 3000)]
fit_win_count = 10      # Only if fit_win_bounds = None

# ===========================
# OUTPUT
# ===========================
do_plot_full_msd = True  # Plot MSD vs t and their log-log plots
do_loglog_linear_fit = True

out_file_prefix = "msd"
out_fig_file_format = "pdf"

# ===========================
# MAIN
# ===========================

if do_plot_full_msd:
    plot_msd(rmsd_file=rmsd_file,
             col_frame=col_frame,
             col_rmsd=col_rmsd,
             time_fs_per_frame=time_fs_per_frame,
             frame_start=frame_start,
             frame_end=frame_end,
             plot_msd_vs_t=True,
             plot_loglog=True,
             uniform_logspace=True,
             out_fig_file=f"{out_file_prefix}-full.{out_fig_file_format}")

if do_loglog_linear_fit:
    ### Specify either one [bounds / count / size] of fit window

    # Converting bounds to log(bounds)
    fit_window_bounds_log = None
    if fit_win_bounds is not None and len(fit_win_bounds) > 0:
        _bounds = np.array(fit_win_bounds)
        _bounds[_bounds == 0] = 1   # replace zeroes by 1 before taking log
        fit_window_bounds_log = np.log(_bounds)

    fit_window_count = fit_win_count
    fit_window_size_log = None  # in ln(num_frames)

    fit_msd_loglog_linear(rmsd_file=rmsd_file,
                          col_frame=col_frame,
                          col_rmsd=col_rmsd,
                          frame_start=frame_start,
                          frame_end=frame_end,
                          fit_window_bounds_log=fit_window_bounds_log,
                          fit_window_size_log=fit_window_size_log,
                          fit_window_count=fit_window_count,
                          time_fs_per_frame=time_fs_per_frame,
                          out_file_prefix=out_file_prefix,
                          out_fig_file_format=out_fig_file_format)