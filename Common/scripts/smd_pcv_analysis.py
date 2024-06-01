"""
Extracts and analyzes Constant-Velocity Pull (pcv) SMD information from NAMD .log file

Analyze the position and force (acceleration) acting on the center-off-mass of SMD atom(s)
at each frame due to the constant velocity pull on dummy atom connected via virtual spring

Can be used for
1. Force (or accelaration) [in direction on pull] vs time, (in pico Newton pN)
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
                    pull_dir_vec: tuple = None,
                    fixed_atoms_com_pos: tuple = None,
                    out_file="smd_pcv_force_displacement.csv",
                    out_file_delimiter="\t")
"""
import math


def vec_add(v1, v2):
    return v1[0] + v2[0], v1[1] + v2[1], v1[2] + v2[2]


def vec_sub(v1, v2):
    return v1[0] - v2[0], v1[1] - v2[1], v1[2] - v2[2]


def vec_mult(v1, scalar: float):
    return v1[0] * scalar, v1[1] * scalar, v1[2] * scalar


def vec_dot(v1, v2):
    return (v1[0] * v2[0]) + (v1[1] * v2[1]) + (v1[2] * v2[2])


def vec_mag_sq(v):
    return vec_dot(v, v)


def vec_mag(v):
    return math.sqrt(vec_mag_sq(v))


def vec_normalize(v):
    return vec_mult(v, 1 / vec_mag(v))


def is_vec_nonzero(v):
    return v and (v[0] or v[1] or v[2])


def create_csv_entry(values, delimiter=","):
    return delimiter.join(values)


def calculate_smd(timestep: int,
                  force_vec: tuple,
                  smd_pos: tuple,
                  smd_pos_ref: tuple,
                  pull_dir_unit_vec: tuple = None,
                  fixed_atoms_com_pos: tuple = None):
    if pull_dir_unit_vec:
        force_mag = vec_dot(force_vec, pull_dir_unit_vec)
        smd_com_displacement = vec_dot(vec_sub(smd_pos, smd_pos_ref), pull_dir_unit_vec)
    else:
        force_mag = vec_mag(force_vec)
        smd_com_displacement = vec_mag(vec_sub(smd_pos, smd_pos_ref))

    if not fixed_atoms_com_pos:
        return timestep, force_mag, smd_com_displacement

    smd_fix_dist = vec_mag(vec_sub(smd_pos, fixed_atoms_com_pos))
    return timestep, force_mag, smd_com_displacement, smd_fix_dist


def analyze_smd_pcv(namd_log_files: [str, list],
                    pull_dir_vec: tuple = None,
                    fixed_atoms_com_pos: tuple = None,
                    out_file: str = "smd_pcv_force_displacement.csv",
                    out_file_delimiter: str = " \t "):
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
    @:parameter out_file: name of the output file
    @:parameter out_file_delimiter: delimiter for the output
    """

    # Handling arguments --------------------------------------
    if isinstance(namd_log_files, str):
        namd_log_files = [namd_log_files]

    if is_vec_nonzero(pull_dir_vec):
        pull_dir_vec = vec_normalize(pull_dir_vec)
    else:
        pull_dir_vec = None

    columns = ["timestep", "smd_com_force", "smd_com_displacement"]
    if fixed_atoms_com_pos:
        columns.append("smd_fixed_dist")
    # ---------------------------------------------------------

    with open(out_file, "w+") as out_fd:

        first_pos_vec = None
        for log_file in namd_log_files:
            with open(log_file, 'r') as f:
                for line in f:
                    line = line.strip()
                    if (line.startswith("SMD") and (words := line.split())) and words[0] == "SMD":
                        if len(words) != 8:
                            pass  # should not happen
                        else:
                            ts = int(words[1])
                            pos_vec = (float(words[2]), float(words[3]), float(words[4]))
                            force_vec = (float(words[5]), float(words[6]), float(words[7]))

                            # reference for calculating displacement of the COM of SMD atom(s)
                            if first_pos_vec is None:
                                first_pos_vec = pos_vec

                                # write csv header before writing any record
                                out_fd.write(out_file_delimiter.join(columns))

                            out = calculate_smd(timestep=ts,
                                                force_vec=force_vec,
                                                smd_pos=pos_vec,
                                                smd_pos_ref=first_pos_vec,
                                                pull_dir_unit_vec=pull_dir_vec,
                                                fixed_atoms_com_pos=fixed_atoms_com_pos)

                            out_fd.write("\n" + out_file_delimiter.join(map(str, out)))


if __name__ == "__main__":
    # TODO: set input parameters
    namd_log_files = ["ubq_wb_pcv0.log"]
    pull_dir_vec = (0.3976350723925137, 0.37333521574953327, 0.8381570055094988)  # or None
    fixed_atoms_com_pos = (24.360, 23.138, 4.649)  # or None

    out_file = "smd_pcv_force_displacement.csv"
    out_file_delimiter = " \t "

    analyze_smd_pcv(namd_log_files=namd_log_files,
                    pull_dir_vec=pull_dir_vec,
                    fixed_atoms_com_pos=fixed_atoms_com_pos,
                    out_file=out_file,
                    out_file_delimiter=out_file_delimiter)
