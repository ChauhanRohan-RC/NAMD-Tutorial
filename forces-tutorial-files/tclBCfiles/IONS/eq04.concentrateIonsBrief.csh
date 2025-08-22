#!/bin/csh

#$ -q linux
#$ -e log-local-err.log
#$ -o log-local-out.log
#$ -cwd

namd2 eq04.concentrateIonsBrief.namd >! eq04.concentrateIonsBrief.log

exit 0
