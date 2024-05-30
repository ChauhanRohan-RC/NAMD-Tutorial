#!/bin/csh

#$ -q linux
#$ -e log-local-err.log
#$ -o log-local-out.log
#$ -cwd

namd2 eq04.shear.namd >! eq04.shear.log

exit 0
