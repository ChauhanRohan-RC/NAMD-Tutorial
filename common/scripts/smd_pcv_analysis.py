import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.figure import figaspect

"""
Extracts and analyzes Constant-Velocity Pull (pcv) SMD information from NAMD .log file

Analyze the position and force (acceleration) acting on the center-off-mass of SMD atom(s)
at each frame due to the constant velocity pull on dummy atom connected via virtual spring

Can be used for
1. Force (or acceleration) [in direction on pull] vs time, (in pico Newton pN)
2. SMD atom(s) COM displacement [in direction on pull] vs time  (in A)
3. Fixed-SMD distance (end-to-end distance) [in direction on pull] vs time  (in A)

If you do not specify the direction of pull, every quantity will be absolute.
Otherwise, every quantity will be "dotted" by the pull direction i.e. its component along the pull direction

Plots you can make:
1. Force vs timestep
2. Force vs extension (either smd_atom_displacement or smd_to_fixed_distance)

USAGE:
1. Copy this script to your working dir
2. go to __name__ == __main__ section and set input parameters
3. run with "python smd_pcv_analysis.py"

@see analyze_smd_pcv(namd_log_files: [str, list],
                    pull_dir_vec: np.array = None,
                    fixed_atoms_com_pos: np.array = None,
                    out_file="smd_pcv_force_displacement.csv",
                    out_file_delimiter="\t")
"""

# Time Step
COLUMN_TIMESTEP = "timestep"
COLUMN_TIME_NS = "time_ns"

# Force experienced by the Center of Mass (COM) of SMD atoms
COLUMN_SMD_COM_FORCE = "smd_com_force"

# Displacement of the COM of SMD atoms (from initial position)
COLUMN_SMD_COM_DISPLACEMENT = "smd_com_displacement"

# Distance between the fixed atom and COM of SMD atoms
COLUMN_SMD_COM_FIXED_DIST = "smd_com_fixed_dist"


def vec_normalize(v):
    return v / np.linalg.norm(v)


def analyze_smd_pcv(namd_log_files: [str, list],
                    pull_dir_vec: np.array = None,
                    fixed_atoms_com_pos: np.array = None,
                    time_step_secs: float = -1,
                    out_file_path: str = "smd_pcv_force_displacement.csv",
                    out_fig_file_path: str = "smd_pcv_force_displacement.pdf",
                    out_delimiter: str = ",",
                    comment_token: str = "#",
                    comment_header: bool = False, ):
    """
    Extracts and analyzes Constant-Velocity Pull (pcv) SMD information from NAMD .log file

    At each timestep, the function calculates
        1. smd_com_force: the mag of force acting on the COM of SMD atom(s). (in pico Newton pN)
            [in the direction of applied pull, if specified]

        2. smd_com_displacement: the distance of the COM of SMD atom(s) from its initial position. (in A)
            [in the direction of applied pull, if specified]

        3. [if position of the COM of fixed atom(s) is given]
            smd_fixed_dist: the distance between COM's of SMD and fixed atom(s), (in A)
            [independent of pull direction]

    @:parameter namd_log_files: 1 or more NAMD .log files for SMD constant vel pull simulation

    @:parameter pull_dir_vec: the direction of applied pull or None. If specified, the displacement and force acting
                              on COM of SMD atom(s) will be dotted with this unit vector

    @:parameter fixed_atoms_com_pos: the position of COM of fixed atom(s), or really any position which will be used as
                                     a reference for calculating the distance of COM of SMD atom(s) at each frame

    @:parameter time_step_secs: (OPTIONAL) time step in seconds, only used for plotting. -1 for None

    @:parameter out_file: name of the output file
    @:parameter out_delimiter: delimiter for the output CSV
    @:parameter comment_token: token used for comments
    @:parameter comment_header: whether to comment the CSV header
    """

    # Handling arguments --------------------------------------
    if isinstance(namd_log_files, str):
        namd_log_files = [namd_log_files]

    if pull_dir_vec is not None and pull_dir_vec.any():
        pull_dir_unit_vec = vec_normalize(pull_dir_vec)
    else:
        pull_dir_unit_vec = None

    # columns = ["timestep", "smd_com_force", "smd_com_displacement"]
    # if fixed_atoms_com_pos:
    #     columns.append("smd_fixed_dist")
    # ---------------------------------------------------------

    has_pull_vec: bool = pull_dir_unit_vec is not None
    has_fixed_pos: bool = fixed_atoms_com_pos is not None

    ts_arr = []
    # smd_force_vec_arr = []
    # smd_pos_vec_arr = []

    smd_force_mag_arr = []
    smd_displacement_arr = []
    smd_fixed_dist_arr = []

    first_pos_vec = None
    for log_file in namd_log_files:
        with open(log_file, 'r') as f:
            for line in f:
                line = line.strip()
                if (line.startswith("SMD") and (words := line.split())) and words[0] == "SMD":
                    if len(words) != 8:
                        print(f"-> Ignoring invalid SMD entry: {line}")
                        pass  # should not happen
                    else:
                        ts = int(words[1])
                        pos_vec = np.array([float(words[2]), float(words[3]), float(words[4])])
                        force_vec = np.array([float(words[5]), float(words[6]), float(words[7])])

                        if first_pos_vec is None:
                            first_pos_vec = pos_vec

                        ts_arr.append(ts)
                        # smd_pos_vec_arr.append(pos_vec)
                        # smd_force_vec_arr.append(force_vec)

                        if has_pull_vec:
                            force_mag = np.dot(force_vec, pull_dir_unit_vec)
                            smd_com_displacement = np.dot(np.subtract(pos_vec, first_pos_vec), pull_dir_unit_vec)
                        else:
                            force_mag = np.linalg.norm(force_vec)
                            smd_com_displacement = np.linalg.norm(np.subtract(pos_vec, first_pos_vec))

                        smd_force_mag_arr.append(force_mag)
                        smd_displacement_arr.append(smd_com_displacement)

                        if has_fixed_pos:
                            smd_fix_dist = np.linalg.norm(np.subtract(pos_vec, fixed_atoms_com_pos))
                            smd_fixed_dist_arr.append(smd_fix_dist)

    # WRITING CSV -------------------------------------------------------
    df = pd.DataFrame()
    ts_col_name = f"{comment_token}{COLUMN_TIMESTEP}" if comment_header else COLUMN_TIMESTEP
    df[ts_col_name] = ts_arr

    has_time = time_step_secs is not None and time_step_secs > 0
    if has_time:
        df[COLUMN_TIME_NS] = df[ts_col_name] * (time_step_secs * 1e9)

    df[COLUMN_SMD_COM_FORCE] = smd_force_mag_arr
    df[COLUMN_SMD_COM_DISPLACEMENT] = smd_displacement_arr
    if has_fixed_pos:
        df[COLUMN_SMD_COM_FIXED_DIST] = smd_fixed_dist_arr

    comments = [
        "-------------------- Constant Vel SMD Analysis -------------------",
        f" -> INPUT Log file(s): {namd_log_files}",
        f" -> INPUT Fixed atom(s) COM position: {fixed_atoms_com_pos}",
        f" -> INPUT Pull direction unit vector: {pull_dir_unit_vec}",
        f" -> Time Step (secs): {time_step_secs}",
        "------------------------------------------------------------------",
        f" -> OUTPUT COLUMNS",
        f"    -> {COLUMN_TIMESTEP}: Time Step",
        f"    -> {COLUMN_TIME_NS}: Time (in nano-seconds)",
        f"    -> {COLUMN_SMD_COM_FORCE}: Force (pN) on COM of SMD atom(s) {'in the pulling direction' if has_pull_vec else ''}",
        f"    -> {COLUMN_SMD_COM_DISPLACEMENT}: Displacement (Å) of COM of SMD atom(s) {'in the pulling direction' if has_pull_vec else ''}",
        f"    -> {COLUMN_SMD_COM_FIXED_DIST}: Distance (Å) between COM of FIXED and SMD atom(s)",
        "-------------------------------------------------------------------"
    ]

    with open(out_file_path, 'w') as f:

        # Meta Information
        for c in comments:
            f.write(f"{comment_token}{c}\n")

        df.to_csv(f, sep=out_delimiter, header=True, index=False)

    # PLOTTING ----------------------------------------------------------------

    w, h = figaspect(9 / 23)
    fig, axes = plt.subplots(1, 2, figsize=(w * 1.4, h * 1.4))
    fig.tight_layout(pad=5.0)

    axes[0].plot(df[COLUMN_TIME_NS if has_time else COLUMN_TIMESTEP], df[COLUMN_SMD_COM_FORCE])
    axes[0].set_title("Force vs " + ("Time" if has_time else "Timestep"))
    axes[0].set_xlabel("Time (ns)" if has_time else "Timestep")
    axes[0].set_ylabel("Force (pN)")

    axes[1].plot(df[COLUMN_SMD_COM_FIXED_DIST if has_fixed_pos else COLUMN_SMD_COM_DISPLACEMENT],
                 df[COLUMN_SMD_COM_FORCE])
    axes[1].set_title("Force vs " + ("Extension" if has_fixed_pos else "Displacement"))
    axes[1].set_xlabel(("SMD-Fixed Distance" if has_fixed_pos else "SMD Displacement") + " (Å)")
    axes[1].set_ylabel("Force (pN)")

    fig.suptitle("Constant Vel SMD")

    if out_fig_file_path:
        plt.savefig(out_fig_file_path)
    plt.show()


if __name__ == "__main__":
    # TODO: set INPUT and OUTPUT params

    # INPUT ------------------
    namd_log_files = ["dna_gbis_pcv.log"]
    pull_dir_vec = np.array([1.0, 0, 0])  # (OPTIONAL) or None
    fixed_atoms_com_pos = np.array([0, 0, 0])  # (OPTIONAL) (in A) or None
    time_step_secs = 1e-15  # (OPTIONAL) Time step (in sec) or -1 for None (only used for plotting)

    # OUTPUT ----------------
    comment_token: str = "#"
    out_file_path: str = "smd_pcv_force_displacement.csv"
    out_fig_file_path: str = "smd_pcv_force_displacement.pdf"

    ## Normal CSV
    out_delimiter: str = ","
    comment_header: bool = False

    ## For XMGrace
    # out_delimiter: str = "\t"
    # comment_header: bool = True

    analyze_smd_pcv(namd_log_files=namd_log_files,
                    pull_dir_vec=pull_dir_vec,
                    fixed_atoms_com_pos=fixed_atoms_com_pos,
                    time_step_secs=time_step_secs,
                    out_file_path=out_file_path,
                    out_fig_file_path=out_fig_file_path,
                    out_delimiter=out_delimiter,
                    comment_token=comment_token,
                    comment_header=comment_header)
