import datetime

"""
Script to extract ENERGIES and calculate ENERGY AVERAGES from NAMD .log file(s)

Energies output by NAMD:
TS, BOND, ANGLE, DIHED, IMPRP, ELECT, VDW, BOUNDARY, MISC, KINETIC, TOTAL, TEMP, POTENTIAL, TOTAL3, TEMPAVG

NOTE: "TS" stands for Time Step

USAGE:
1. copy this script to working dir
2. edit __name__ == "__main__" section of this script
2. run with "python namd_energy.py"
"""


def _extract_energies_internal(namd_log_files: [str, list],
                               e_titles_callback: callable,
                               energies_callback: callable,
                               start_timestep: int = -1,
                               end_timestep: int = -1):
    """
    Base method to Extract energies from NAMD .log files and recive Callbacks

    @param namd_log_files : NAMD .log file or a list of log files
    @param e_titles_callback : a function that is called with ENERGY_TITLES (as list of strings) that are present in the first log file
                            Called only ONCE at the start:
                            > e_titles_callback(energy_titles: list[str])

    @param energies_callback : a function that is called for each entry of all energies (in all log files)
                            It is guaranteed that "e_titles_callback(energy_titles)" will be called before
                            the first call to this function.
                            Energies are passed as list of strings, where each energy corresponds to the ENERGY_TITLE
                            of energy_titles list passed to e_titles_callback.
                            > energies_callback(energies: list[str])

    @param start_timestep : timestep to start reading energy_values from .log file(s).
                            -1 to start from the begining (any timestep present first)

    @param end_timestep : timestep to end reading energy_values from .log file(s).
                            -1 to read till the end of all .log files

    """

    if isinstance(namd_log_files, str):
        namd_log_files = [namd_log_files]

    e_titles = None

    for log_file_path in namd_log_files:
        with open(log_file_path, 'r') as f:
            for line in f:
                line = line.strip()
                if e_titles is None:
                    if (line.startswith("ETITLE:") and (words := line.split())) and words[0] == "ETITLE:":
                        e_titles = words[1:]  # First will be Timestep
                        e_titles_callback(e_titles)
                    continue

                # Now we have e_titles
                if (line.startswith("ENERGY:") and (words := line.split())) and words[0] == "ENERGY:":
                    energies = words[1:]  # First will be Timestep

                    if start_timestep >= 0 or end_timestep >= 0:
                        ts = int(energies[0])
                        if (start_timestep >= 0 and ts < start_timestep) or (end_timestep >= 0 and ts > end_timestep):
                            continue

                    energies_callback(energies)


def extract_energies(namd_log_files: [str, list],
                     start_timestep: int = -1,
                     end_timestep: int = -1,
                     energy_cols: [str, list] = None,
                     out_file_name: str = "energies.csv",
                     out_delimiter: str = " \t ",
                     comment_token: str = "#",
                     comment_e_titles: bool = False):
    """
    Extract specified (or all) ENERGY values from NAMD .log file(s)
    within certain timestep range or for all time steps

    Generally available ENERGY_COLUMNS:
    BOND, ANGLE, DIHED, IMPRP, ELECT, VDW, BOUNDARY, MISC, KINETIC, TOTAL, TEMP, POTENTIAL, TOTAL3, TEMPAVG

    @param namd_log_files : NAMD .log file or a list of log files
    @param start_timestep : timestep to start reading energy_values from .log file(s).
                            -1 to start from the begining (any timestep present first)
    @param end_timestep : timestep to end reading energy_values from .log file(s).
                            -1 to read till the end of all .log files
    @param energy_cols : Energy fields required.
                        can be a single string, or list of strings, or None for all energy fields present in .log file(s)
    @param out_file_name : name of the output file
    @param out_delimiter : delimiter to use for output
    @param comment_token : A token for comments. Empty string to disable comments
    @param comment_e_titles : Whether to comment Energy Titles as well. Useful for programs
                            such as Xmgrace which cannot parse headers
    """

    if isinstance(energy_cols, str):
        if energy_cols.strip():
            energy_cols = [energy_cols]
        else:
            energy_cols = None

    indices = []
    with open(out_file_name, "w") as out_fd:
        def e_titles_callback(e_titles: list):
            if energy_cols:
                if "TS" not in energy_cols and "TS" in e_titles:
                    indices.append(e_titles.index("TS"))
                indices.extend((e_titles.index(c) for c in energy_cols if c in e_titles))
            else:
                indices.clear()
                indices.extend(range(0, len(e_titles)))

            if indices:
                indices.sort()

                # Comments
                if comment_token:
                    out_fd.write(
                        f"{comment_token} ENERGIES extracted from NAMD .log file(s): [{', '.join(namd_log_files)}]\n")
                    out_fd.write(f"{comment_token} Start Timestep: {start_timestep} | End Timestep: {end_timestep}\n")
                    out_fd.write(f"{comment_token} File written on: {datetime.datetime.now()}\n")
                    out_fd.write(f"{comment_token}-----------------------------\n")

                # Header: Energy Titles
                out_fd.write((comment_token if comment_token and comment_e_titles else "")
                             + out_delimiter.join((e_titles[i] for i in indices)))

        def energies_callback(energies):
            if indices:
                # Energy Values
                out_fd.write("\n" + out_delimiter.join((energies[i] for i in indices)))

        _extract_energies_internal(namd_log_files=namd_log_files,
                                   start_timestep=start_timestep,
                                   end_timestep=end_timestep,
                                   e_titles_callback=e_titles_callback,
                                   energies_callback=energies_callback)

        print(f"INFO: Energies vs TS written to file: '{out_file_name}' | delimiter: '{out_delimiter}'"
              f" | comment_token: '{comment_token}' | comment_energy_titles: {comment_e_titles}")


def energies_average(namd_log_files: [str, list],
                     start_timestep: int = -1,
                     end_timestep: int = -1,
                     energy_cols: [str, list] = None,
                     out_file_name: str = "energies_avg.csv",
                     out_delimiter: str = " \t ",
                     comment_token: str = "#",
                     comment_e_titles: bool = False):
    """
    Compute AVERAGE ENERGIES from NAMD .log file(s)
    within certain timestep range or for all time steps

    Generally available ENERGY_COLUMNS:
    BOND, ANGLE, DIHED, IMPRP, ELECT, VDW, BOUNDARY, MISC, KINETIC, TOTAL, TEMP, POTENTIAL, TOTAL3, TEMPAVG

    @param namd_log_files : NAMD .log file or a list of log files
    @param start_timestep : timestep to start reading energy_values from .log file(s).
                            -1 to start from the begining (any timestep present first)
    @param end_timestep : timestep to end reading energy_values from .log file(s).
                            -1 to read till the end of all .log files
    @param energy_cols : Energy fields required.
                        can be a single string, or list of strings, or None for all energy fields present in .log file(s)
    @param out_file_name : name of the output file
    @param out_delimiter : delimiter to use for output
    @param comment_token : A token for comments. Empty string to disable comments
    @param comment_e_titles : Whether to comment Energy Titles as well. Useful for programs
                            such as Xmgrace which cannot parse headers
    """

    if isinstance(energy_cols, str):
        if energy_cols.strip():
            energy_cols = [energy_cols]
        else:
            energy_cols = None

    indices = []  # indices of energy cols: list[int]
    titles = []  # energy titles: list[str]
    values = []  # energy values sum: list[float]
    value_count = [0]  # hack to make value count mutable within callback

    def e_titles_callback(e_titles: list):
        if energy_cols:
            indices.extend((e_titles.index(c) for c in energy_cols if c in e_titles))
        else:
            indices.clear()
            indices.extend(range(0, len(e_titles)))

        if indices:
            indices.sort()
            titles.clear()
            titles.extend((e_titles[idx] for idx in indices))
            values.clear()
            values.extend((0 for _ in range(0, len(indices))))  # initialize values with 0

    def energies_callback(energies):
        if indices:
            for _i, idx in enumerate(indices):
                values[_i] += float(energies[idx])

            value_count[0] += 1

    _extract_energies_internal(namd_log_files=namd_log_files,
                               start_timestep=start_timestep,
                               end_timestep=end_timestep,
                               e_titles_callback=e_titles_callback,
                               energies_callback=energies_callback)

    for i in range(0, len(values)):
        values[i] /= value_count[0]

    # Console Output
    print("----------------------")
    print("AVERAGE ENERGIES:")
    print("\n".join(f"-> {titles[i]} : {values[i]}" for i in range(0, len(titles))))

    # File Output
    if out_file_name:
        if isinstance(namd_log_files, str):
            namd_log_files = [namd_log_files]

        with open(out_file_name, "w") as out_fd:
            # Comments
            if comment_token:
                out_fd.write(
                    f"{comment_token} ENERGY AVERAGES calculated from NAMD .log file(s): [{', '.join(namd_log_files)}]\n")
                out_fd.write(f"{comment_token} Start Timestep: {start_timestep} | End Timestep: {end_timestep}\n")
                out_fd.write(f"{comment_token} File written on: {datetime.datetime.now()}\n")
                out_fd.write(f"{comment_token}-----------------------------\n")

            # Header: Energy Titles
            out_fd.write((comment_token if comment_token and comment_e_titles else "")
                         + out_delimiter.join(titles) + "\n")

            # Average Energy Values
            out_fd.write(out_delimiter.join(map(str, values)))

        print("----------------------")
        print(f"INFO: Average Energies written to file: '{out_file_name}' | delimiter: '{out_delimiter}'"
              f" | comment_token: '{comment_token}' | comment_energy_titles: {comment_e_titles}")

    return dict((titles[i], values[i]) for i in range(0, len(titles)))


if __name__ == '__main__':
    #### Energies output by NAMD:
    # TS, BOND, ANGLE, DIHED, IMPRP, ELECT, VDW, BOUNDARY, MISC, KINETIC, TOTAL, TEMP, POTENTIAL, TOTAL3, TEMPAVG

    ## ====================  INPUT CONFIG  =========================
    namd_log_files = ["test.log"]  # TODO: input .log file(s)
    # energy_columns = ["BOND", "TEMP", "KINETIC"]  # TODO: Energy columns required. None or empty list for all columns
    energy_columns = None
    
    start_timestep = -1  # TODO: start timestep, or -1 to start form beginning
    end_timestep = -1  # TODO: end timestep, or -1 to read till the end of all .log file(s)

    ## =====================  OUTPUT CONFIG  ==========================
    # output_file = "energies.csv"
    # out_delimiter = " \t "
    # comment_token = ""         # Empty string to disable comments
    comment_e_titles = False  # whether to comment ENERGY TITLES in the Header, useful for Xmgrace

    ## ----------------------------------------------------------------

    ## EXTRACT ENERGIES vs TS
    extract_energies(namd_log_files=namd_log_files,
                     energy_cols=energy_columns,
                     start_timestep=start_timestep,
                     end_timestep=end_timestep,
                     # out_file_name=output_file,
                     # out_delimiter=out_delimiter,
                     # comment_token=comment_token,
                     comment_e_titles=comment_e_titles)

    ## AVERAGE ENERGIES
    energies_average(namd_log_files=namd_log_files,
                     energy_cols=energy_columns,
                     start_timestep=start_timestep,
                     end_timestep=end_timestep,
                     # out_file_name=output_file,
                     # out_delimiter=out_delimiter,
                     # comment_token=comment_token,
                     comment_e_titles=comment_e_titles)
