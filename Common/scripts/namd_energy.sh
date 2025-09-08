#!/bin/bash

##############################################
# Energy Calculations using namd_energy.tcl
##############################################
# todo: set env vars for each run
# todo: set corresponding vars to "" in namd_energy.tcl to use environment variables

run_namd_energy() {
	vmd -dispdev text -e namd_energy.tcl
}

export NAMD_ENERGY_OUT_ENERGIES="-all"

# Run 1: Protein Self
export NAMD_ENERGY_SELECTION1="protein"
export NAMD_ENERGY_SELECTION2=""
export NAMD_ENERGY_OUT_PREFIX="prot_self"
run_namd_energy

# Run 2: Water Self
export NAMD_ENERGY_SELECTION1="water"
export NAMD_ENERGY_SELECTION2=""
export NAMD_ENERGY_OUT_PREFIX="water_self"
run_namd_energy

# Run 3: Protein-Water
export NAMD_ENERGY_SELECTION1="protein"
export NAMD_ENERGY_SELECTION2="water"
export NAMD_ENERGY_OUT_PREFIX="prot_water"
run_namd_energy
