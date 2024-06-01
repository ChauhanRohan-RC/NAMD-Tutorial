"""
Script to extract Energies from NAMD .log file(s)

Energies output by NAMD:
BOND, ANGLE, DIHED, IMPRP, ELECT, VDW, BOUNDARY, MISC, KINETIC, TOTAL, TEMP, POTENTIAL, TOTAL3, TEMPAVG

NOTE: "TS" stands for Time Step

USAGE:
1. copy this script to working dir
2. edit __name__ == "__main__" section of this script
2. run with "python namd_energy.py"
"""


def _extract_energies(namd_log_files: [str, list],
                      e_titles_callback,
                      energies_callback,
                      start_timestep=-1,
                      end_timestep=-1):
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
                     start_timestep=-1,
                     end_timestep=-1,
                     energy_cols: [str, list] = None,
                     output_file_name="energies.csv",
                     out_delimiter=" \t "):

    if isinstance(energy_cols, str):
        if energy_cols.strip():
            energy_cols = [energy_cols]
        else:
            energy_cols = None

    indices = []
    with open(output_file_name, "w") as out_fd:

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
                out_fd.write(out_delimiter.join((e_titles[i] for i in indices)))

        def energies_callback(energies):
            if indices:
                out_fd.write("\n" + out_delimiter.join((energies[i] for i in indices)))

        _extract_energies(namd_log_files=namd_log_files,
                          start_timestep=start_timestep,
                          end_timestep=end_timestep,
                          e_titles_callback=e_titles_callback,
                          energies_callback=energies_callback)


if __name__ == '__main__':

    # Energies output by NAMD:
    # BOND, ANGLE, DIHED, IMPRP, ELECT, VDW, BOUNDARY, MISC, KINETIC, TOTAL, TEMP, POTENTIAL, TOTAL3, TEMPAVG

    namd_log_files = ["ubq_wb_pcf0.log"]
    energy_columns = ["TEMP"]      # None or empty list for all columns
    start_timestep = -1
    end_timestep = -1
    output_file = "temp.dat"
    out_delimiter = " \t "

    extract_energies(namd_log_files=namd_log_files,
                     energy_cols=energy_columns,
                     start_timestep=start_timestep,
                     end_timestep=end_timestep,
                     output_file_name=output_file,
                     out_delimiter=out_delimiter)
