#!/usr/bin/perl

#Script generating the hibrid strctures with the AlfaFold model and the cristal loop.Also obtains de RMSD of the AF loop.

#Iván Martín Hernández

# Date: Jun - 2021



use strict;
use warnings;
use diagnostics;
use Data::Dumper qw(Dumper);
use feature 'say';
use Switch;
use Scalar::Util qw(reftype);
use POSIX;
use File::Basename;


#INPUT PARSER
if (!($#ARGV == 0)) {
    say "USAGE:\n\t$0  <input1> \n";
    say "\t<input1>  --> The rcd_list file with all the structures of the folder.";
    say "\t<input2>  --> The CDR region.";
    say "\nDESCRIPTION:\n\tSay something here...";
    say "\nWARNING:\n\tSomething....";
    say "\nEXAMPLE:";
    say "\t$0  AB_loop_RCD_H3.txt H3";
    say "\nHELP: some help\n";
    exit;
}

#take the inputs from the bash
my $rcdfile = $ARGV[0];
my $col_rosettaE = 9; # Column index (for -k option of sort command)
my $col_rmsd = 12; # Column index of PyRosetta refined RMSD


my %AA = (A=>'ALA',Y=>'TYR',M=>'MET',L=>'LEU',C=>'CYS',G=>'GLY',
         R=>'ARG',N=>'ASN',D=>'ASP',Q=>'GLN',E=>'GLU',H=>'HIS',W=>'TRP',
         K=>'LYS',F=>'PHE',P=>'PRO',S=>'SER',T=>'THR',I=>'ILE',V=>'VAL');


my %aa = reverse %AA;





#Read input file
my $data = parseTextFile($rcdfile);

my $nloop = 0;
my $base;#required outside

foreach my $loop(@{$data}){
    my $pdb = substr($loop -> [0],0,4);
    my $start = $loop -> [1];
    my $end = $loop -> [2];
    my $chain = $loop -> [3];
	$base=$pdb;
	
	my $RosSfiles = "$pdb*_rosettaS.txt";
	#print "R- $RosSfiles\n";
	
	#############################################
	### Block 1: not removing same loop rosettas
	#############################################
	my $outrosfile1=$pdb."_TOTAL_rosettaTS.txt";
	open(ROS1,'>',$outrosfile1);
	
	print ROS1 "#loop   E_full   E_loop  R_rcd       E_KORP  R_BB  R_HA  R_SC  E_full2  R_CA  R_BB R_BBO  R_CB R_OCB  R_HA  R_SC R_ALL \n";
	
	my $cmd = "tail -n +3 $RosSfiles |sort -k 9 -n|tail -n +10";
	print "$cmd\n";
	my $sorted=`$cmd`;
	
	#print "######\n $sorted\n";
	print ROS1 "$sorted";
	
	
	close (ROS1);
	
	
	#############################################
	### Block 2: removing same loop rosettas
	#############################################
	my @files = glob("$RosSfiles");
	
	my %loops;
	
	my $outrosfile2=$pdb."_TOTAL_rosettaRS.txt";
	open(ROS2,'>',$outrosfile2);
	print ROS2 "#loop   E_full   E_loop  R_rcd       E_KORP  R_BB  R_HA  R_SC  E_full2  R_CA  R_BB R_BBO  R_CB R_OCB  R_HA  R_SC R_ALL \n";
	
	foreach my $file(@files){
		open (PDB, $file) or die "Can't open file $file.\n";
		while(<PDB>)
        {
			next if /^#/; # reading only non-# beginning lines
			chomp($_);
			my $line=$_;
			#print "$line\n";
			my @line = split(/\s+/, $_);
			my $loop=$line[1];
			my $Renergy=$line[$col_rosettaE];
			
			if (defined $loops{$loop}{Energy} ) {
				if ($loops{$loop}{Energy} > $Renergy) {
                    $loops{$loop}{Energy}=$Renergy;
					$loops{$loop}{Line}=$line;
                }
				
                
            }
			else{
				$loops{$loop}{Energy}=$Renergy;
				$loops{$loop}{Line}=$line;
				}
            
		
		}
		
	}
	
	foreach my $loop (sort {$loops{$a}{Energy} <=> $loops{$b}{Energy}} keys %loops){
		
		#print "$loop -- $loops{$loop}{Energy} $loops{$loop}{Line}\n";
		
		print ROS2 "$loops{$loop}{Line}\n";
	}

	#print ROS2 "$sorted";
	
	
	close (ROS2);
	
	
	
	
	
	
}



#############################################
###  Resuming Block 1:
#############################################
# Get best-1 energy decoys statistics
my $kk = "head -n 1 $base"."_TOTAL_rosettaTS.txt > rosettaTS1.txt; for i in *_TOTAL_rosettaTS.txt; do var=\`awk \'NR==2\' \$i\ | sort -n -k $col_rmsd | head -n 1`; echo -e \"\${i:0:20} \$var\"; done > temp; statext.pl temp 1 18 1 all 2 >> rosettaTS1.txt";
print "EXECUTING: $kk\n";
open my $fh, '-|', qw(/bin/bash -c), $kk;
close($fh);
# Get best-2 energy decoys statistics
$kk = "head -n 1 $base"."_TOTAL_rosettaTS.txt > rosettaTS2.txt; for i in *_TOTAL_rosettaTS.txt; do var=\`awk \'NR>=2 && NR<2+2\' \$i\ | sort -n -k $col_rmsd | head -n 1`; echo -e \"\${i:0:20} \$var\"; done > temp; statext.pl temp 1 18 1 all 2  >> rosettaTS2.txt";
print "EXECUTING: $kk\n";
open $fh, '-|', qw(/bin/bash -c), $kk;
# Get best-5 energy decoys statistics
$kk = "head -n 1 $base"."_TOTAL_rosettaTS.txt > rosettaTS5.txt; for i in *_TOTAL_rosettaTS.txt; do var=\`awk \'NR>=2 && NR<2+5\' \$i\ | sort -n -k $col_rmsd | head -n 1`; echo -e \"\${i:0:20} \$var\"; done > temp; statext.pl temp 1 18 1 all 2 >> rosettaTS5.txt";
print "EXECUTING: $kk\n";
open $fh, '-|', qw(/bin/bash -c), $kk;
# Get best-10 energy decoys statistics
$kk = "head -n 1 $base"."_TOTAL_rosettaTS.txt > rosettaTS10.txt; for i in *_TOTAL_rosettaTS.txt; do var=\`awk \'NR>=2 && NR<2+10\' \$i\ | sort -n -k $col_rmsd | head -n 1`; echo -e \"\${i:0:20} \$var\"; done > temp; statext.pl temp 1 18 1 all 2 >> rosettaTS10.txt";
print "EXECUTING: $kk\n";
open $fh, '-|', qw(/bin/bash -c), $kk;
# Get best-20 energy decoys statistics
$kk = "head -n 1 $base"."_TOTAL_rosettaTS.txt > rosettaTS20.txt; for i in *_TOTAL_rosettaTS.txt; do var=\`awk \'NR>=2 && NR<2+20\' \$i\ | sort -n -k $col_rmsd | head -n 1`; echo -e \"\${i:0:20} \$var\"; done > temp; statext.pl temp 1 18 1 all 2  >> rosettaTS20.txt";
print "EXECUTING: $kk\n";
open $fh, '-|', qw(/bin/bash -c), $kk;
# Get best-50 energy decoys statistics
$kk = "head -n 1 $base"."_TOTAL_rosettaTS.txt > rosettaTS50.txt; for i in *_TOTAL_rosettaTS.txt; do var=\`awk \'NR>=2 && NR<2+50\' \$i\ | sort -n -k $col_rmsd | head -n 1`; echo -e \"\${i:0:20} \$var\"; done > temp; statext.pl temp 1 18 1 all 2  >> rosettaTS50.txt";
print "EXECUTING: $kk\n";
open $fh, '-|', qw(/bin/bash -c), $kk;
# Get best-100 energy decoys statistics
$kk = "head -n 1 $base"."_TOTAL_rosettaTS.txt > rosettaTS100.txt; for i in *_TOTAL_rosettaTS.txt; do var=\`awk \'NR>=2 && NR<2+100\' \$i\ | sort -n -k $col_rmsd | head -n 1`; echo -e \"\${i:0:20} \$var\"; done > temp; statext.pl temp 1 18 1 all 2 >> rosettaTS100.txt; rm temp";
print "EXECUTING: $kk\n";
open $fh, '-|', qw(/bin/bash -c), $kk;
close($fh);



#############################################
###  Resuming Block 2:
#############################################
# Get best-1 energy decoys statistics
$kk = "head -n 1 $base"."_TOTAL_rosettaRS.txt > rosettaRS1.txt; for i in *_TOTAL_rosettaRS.txt; do var=\`awk \'NR==2\' \$i\ | sort -n -k $col_rmsd | head -n 1`; echo -e \"\${i:0:20} \$var\"; done > temp; statext.pl temp 1 18 1 all 2 >> rosettaRS1.txt";
print "EXECUTING: $kk\n";
open $fh, '-|', qw(/bin/bash -c), $kk;
close($fh);
# Get best-2 energy decoys statistics
$kk = "head -n 1 $base"."_TOTAL_rosettaRS.txt > rosettaRS2.txt; for i in *_TOTAL_rosettaRS.txt; do var=\`awk \'NR>=2 && NR<2+2\' \$i\ | sort -n -k $col_rmsd | head -n 1`; echo -e \"\${i:0:20} \$var\"; done > temp; statext.pl temp 1 18 1 all 2  >> rosettaRS2.txt";
print "EXECUTING: $kk\n";
open $fh, '-|', qw(/bin/bash -c), $kk;
# Get best-5 energy decoys statistics
$kk = "head -n 1 $base"."_TOTAL_rosettaRS.txt > rosettaRS5.txt; for i in *_TOTAL_rosettaRS.txt; do var=\`awk \'NR>=2 && NR<2+5\' \$i\ | sort -n -k $col_rmsd | head -n 1`; echo -e \"\${i:0:20} \$var\"; done > temp; statext.pl temp 1 18 1 all 2 >> rosettaRS5.txt";
print "EXECUTING: $kk\n";
open $fh, '-|', qw(/bin/bash -c), $kk;
# Get best-10 energy decoys statistics
$kk = "head -n 1 $base"."_TOTAL_rosettaRS.txt > rosettaRS10.txt; for i in *_TOTAL_rosettaRS.txt; do var=\`awk \'NR>=2 && NR<2+10\' \$i\ | sort -n -k $col_rmsd | head -n 1`; echo -e \"\${i:0:20} \$var\"; done > temp; statext.pl temp 1 18 1 all 2 >> rosettaRS10.txt";
print "EXECUTING: $kk\n";
open $fh, '-|', qw(/bin/bash -c), $kk;
# Get best-20 energy decoys statistics
$kk = "head -n 1 $base"."_TOTAL_rosettaRS.txt > rosettaRS20.txt; for i in *_TOTAL_rosettaRS.txt; do var=\`awk \'NR>=2 && NR<2+20\' \$i\ | sort -n -k $col_rmsd | head -n 1`; echo -e \"\${i:0:20} \$var\"; done > temp; statext.pl temp 1 18 1 all 2  >> rosettaRS20.txt";
print "EXECUTING: $kk\n";
open $fh, '-|', qw(/bin/bash -c), $kk;
# Get best-50 energy decoys statistics
$kk = "head -n 1 $base"."_TOTAL_rosettaRS.txt > rosettaRS50.txt; for i in *_TOTAL_rosettaRS.txt; do var=\`awk \'NR>=2 && NR<2+50\' \$i\ | sort -n -k $col_rmsd | head -n 1`; echo -e \"\${i:0:20} \$var\"; done > temp; statext.pl temp 1 18 1 all 2  >> rosettaRS50.txt";
print "EXECUTING: $kk\n";
open $fh, '-|', qw(/bin/bash -c), $kk;
# Get best-100 energy decoys statistics
$kk = "head -n 1 $base"."_TOTAL_rosettaRS.txt > rosettaRS100.txt; for i in *_TOTAL_rosettaRS.txt; do var=\`awk \'NR>=2 && NR<2+100\' \$i\ | sort -n -k $col_rmsd | head -n 1`; echo -e \"\${i:0:20} \$var\"; done > temp; statext.pl temp 1 18 1 all 2 >> rosettaRS100.txt; rm temp";
print "EXECUTING: $kk\n";
open $fh, '-|', qw(/bin/bash -c), $kk;
close($fh);

print " That's all folks!\n";


##################################################################FUNCTIONS#################################################################

# Read text file line by line and separate each line into fields
sub parseTextFile {
    my $file = shift;#    Input file name
    open(IN, $file) or die "\nFailed to open $file\n";

    my@ data = ();
    my $nline = 0;#    Line index(counter)
    print "parseTextFile> Reading file: $file --> ";
    while ( <IN> ) {
        next
        if /^#/;#        reading only non - #begining lines
        my @line = split(/\s+/, $_);#        print "@line\n";
        push(@data, \@line);
        $nline++;
    }
    print " $nline lines\n";
    close IN;

    return\@ data;#    Return reference to data
}
