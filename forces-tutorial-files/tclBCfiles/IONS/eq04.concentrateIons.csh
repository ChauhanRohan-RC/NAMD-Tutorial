#!/bin/csh

#$ -q linux
#$ -e log-local-err.log
#$ -o log-local-out.log
#$ -cwd

namd2 eq04.concentrateIons.namd >! eq04.concentrateIons.log

exit 0
