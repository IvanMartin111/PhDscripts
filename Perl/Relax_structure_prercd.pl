#!/usr/bin/perl
#
# RCD+ loop refinement for "native" or "modeling" scenarios
#
# by Mon (05/02/2019)
#

use strict;
use Storable;
#use Expect;
use File::Basename;


# SOME DEFINITIONS

# Path to the Rosetta's loop minimization python (e.g. ../loops_energyXX.py). Relative to the execution directory.
my $minimizer = "newmin_loops_energy13.py";
my $relax = "newmin_monrelax4.py";
my $insert = "insertpdb.pl";
my $inpref = "../";
my $col_rosettaE = 9; # Column index (for -k option of sort command)
my $col_rmsd = 12; # Column index of PyRosetta refined RMSD
my $run = 1; # =0, Disables execution and forces plain text dump of commands; otherwise run.

# PyRosetta options for minimization
my $optrelax = "--nmin 5 --initneighbor 10 --nminneigh 100 --irepack 0 --mintype dfpmin_armijo_nonmonotone --fas=0.05,0.30,0.60,1.00 --ntpsi --deleteloop ";


# INPUT PARSER
if( !($#ARGV == 2 || $#ARGV == 1) )
{
	print "USAGE:\n\t$0 <RCD_input_file> [native|modeling] [outsuff] [np] [\"parameters\"]>(Optional)\n\n";
	print "\t<RCD_input_file>  --> RCD's input text file indicating: loop PDB, start and end indices, chain, and sequence.\n";
	print "\t[native|modeling] --> Choose between \"native\" or \"modeling\" scenarios.\n";
	print "\t[outsuff]         --> Suffix for output files. E.g: \"_ref2015\"\n";
	print "\t[np]              --> Number of processes for MPI parallelization.\n";
	print "\t<parameters>      --> (Optional) Additional parameters for PyRosetta script. E.g: \"--energy ref2015\"\n";
	print "\n\tNOTE that current PyRosetta parameters are: \n";
	print "\nDESCRIPTION:\n\tRCD+ loop refinement for \"native\" or \"modeling\" scenarios.\n";
	print "\nWARNING:\n\tEnsure that the scripts: monrelax4.py loops_energy13.py monfuncs.py insertpdb.pl are on PATH.\n";
	print "\nEXAMPLE:\n";
	print "\ttime refine_loops.pl ../mini8.txt native _ref2015nat 4 \"--energy ref2015\"\n";
	print "\ttime refine_loops.pl ../mini8.txt modeling _ref2015mod 4 \"--energy ref2015\"\n";
	print "\nHELP: *Loops can be generated using, for example:\n\ttime mpirun.openmpi -np 4 rcd_mpi_gnu mini8.txt -n 10000 -t 0.99 -r --linear -d 0.5 -x dunbrack.bin -e 5 --energy_file ~/Korp6Dv1/korp6Dv1.bin --bench --loco_best 100 -o mini8_10k1H\n";
	print "\ttime mpirun.openmpi -np 4 rcd_mpi_gnu seok12.txt -n 10000 -t 0.9 -x ~/RCDprot_v103/dunbrack.bin -r --linear -d 0.5 -e 5 --energy_file ~/Korp6Dv1/korp6Dv1.bin --loco_best 100 --bench -o seok12_t90n10k1H\n\n";
	exit;
}
# Some commands:
# for i in *_ref2015nat_rosettaS.txt; do kk=`awk 'NR==3' $i`; echo -e "${i:0:15} $kk"; done > rosetta_ref2015nat.txt
# for i in *_ref2015mod_rosettaS.txt; do kk=`awk 'NR==3' $i`; echo -e "${i:0:15} $kk"; done > rosetta_ref2015mod.txt

# Loading things        
my $debug = 0;
my $prog = "refine_loops";
my $rcdfile = $ARGV[0];
my $outsuff = $ARGV[1];

# Read input file
my $data = parseTextFile($rcdfile);

my $nloop = 0;
my $base; # required outside
foreach my $loop (@{$data}) # Screen Loops
{
	printf "$prog> PROCESSING loop %3d: @{$loop}\n",$nloop+1;
	my $pdb = $loop->[0];
	my $start = $loop->[1];
	my $end = $loop->[2];
	my $chain = $loop->[3];
	my $start2 = $start - 1;
	my $end2 = $end + 1;
	$base = basename($pdb,".pdb");


	# MAIN LOOP for MODELING scenario

		# 5. Main PDB structure must be relaxed first in Rosetta potential without the native loop (--deleteloop).
		#    (This stage can be carried out in parallel with the next one.)
		my $cmd = "time $relax --base_in $base --ri $start2 --rf $end2 --chain $chain --pdb_in $inpref$pdb --base_out $base$outsuff $optrelax";

		# Using "Expect.pm" to dump program output during command execution
		# Install Expect.pm using: perl -MCPAN -e 'install Bundle::Expect'
		print "EXECUTING COMMAND: $cmd\n";
                system($cmd);
#


}



#################################################################
# FUNCTIONS
#################################################################

# Read text file line by line and separate each line into fields
sub parseTextFile
{
	my $file = shift; # Input file name
	open(IN,$file) or die "\nFailed to open $file\n";

	my @data = ();
	my $nline = 0; # Line index (counter)
	print "parseTextFile> Reading file: $file --> ";
	while(<IN>)
	{
		next if /^#/; # reading only non-# begining lines
		my @line = split(/\s+/, $_);
		# print "@line\n";
		push(@data,\@line);
		$nline++;
	}
	print " $nline lines\n";
	close IN;
	
	return \@data; # Return reference to data
}

