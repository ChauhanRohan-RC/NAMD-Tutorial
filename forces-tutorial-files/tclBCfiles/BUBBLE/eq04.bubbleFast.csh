#!/bin/csh

#$ -q linux
#$ -e log-local-err.log
#$ -o log-local-out.log
#$ -cwd

namd2 eq04.bubbleFast.namd >! eq04.bubbleFast.log

exit 0
