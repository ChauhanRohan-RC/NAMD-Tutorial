#!/bin/csh

#$ -q linux
#$ -e log-local-err.log
#$ -o log-local-out.log
#$ -cwd

namd2 eq04.bubble.namd >! eq04.bubble.log

exit 0
