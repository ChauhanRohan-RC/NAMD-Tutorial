## Script to calculate RMSD of a molecule at each timestep (TS) in VMD

# 1. Load the .psf and .dcd in VMD
# 2. Go to Extensions > Tk Console
# 3. enter command: source rmsd.tcl
# 4. Generates output file "rmsd.dat"

# set outfile [open rmsd.csv w]
set outfile [open rmsd.dat w]

# Number of frames
set nf [molinfo top get numframes]

# Reference Conformation (defaults to first frame) 
# selecting only backbone atoms of protein excluding Hydrogens ("noh" means no Hydrogen)
# to exculde some resudies, append "and not (resid 74 to 76)"
set refframe [atomselect top "protein and backbone and noh" frame 0]

# current frame variable, must have same atom selection as reference
set curframe [atomselect top "protein and backbone and noh"]

# RMSD loop
# puts $outfile "Frame,RMSD";	# Header for csv

for { set i 0 } { $i < $nf } { incr i } {
$curframe frame $i
$curframe move [measure fit $curframe $refframe]
puts $outfile "[measure rmsd $curframe $refframe]"
# puts $outfile "$i,[measure rmsd $curframe $refframe]"
}

close $outfile
