#!/usr/bin/perl

#Script for refinement of the RCD input for the txt file produced by "AB_SAbDab_sumarize.pl" and generatin the individual list (without repetitions) for each CDR

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

# INPUT PARSER
if( !($#ARGV == 0))
{
	say "USAGE:\n\t$0  <input1> \n";
	say "\t<input1>  --> The rcd_list file with all the strcutures of the folder.";
	say "\nDESCRIPTION:\n\tSay something here...";
	say "\nWARNING:\n\tSomething....";
	say "\nEXAMPLE:";
	say "\t$0  AB_loop_RCD_H3.txt ";
	say "\nHELP: some help\n"; 
	exit;
}

#take the inputs from the bash
my $file = $ARGV[0];



my %AA = (A=>'ALA',Y=>'TYR',M=>'MET',L=>'LEU',C=>'CYS',G=>'GLY',
         R=>'ARG',N=>'ASN',D=>'ASP',Q=>'GLN',E=>'GLU',H=>'HIS',W=>'TRP',
         K=>'LYS',F=>'PHE',P=>'PRO',S=>'SER',T=>'THR',I=>'ILE',V=>'VAL');


my %aa = reverse %AA;




####################################################
open (DATA,$file) or die "The file $file doesn't exit in this path.\n";

open(RKYD,'>',"H3_RK_and_D.txt");
open(RKOD,'>',"H3_RK_or_D.txt");

open(RKO,'>',"H3_RK_only.txt");
open(DO,'>',"H3_D_only.txt");

open(RKND,'>',"H3_RK_not_D.txt");






my $nline = 0; # Line index (counter)

my $nRKYD = 0;
my $nRKOD = 0;
my $nRKO = 0;
my $nDO = 0;
my $nRKND = 0;

while(<DATA>)
{
	next if /^#/; # reading only non-# beginning lines
	$_=~ s/^\t//;
	chop $_;
	my @line = split(/\s+/, $_);
	#push(@data,\@line);
	$nline++;
	my $pdb = substr($line[0],0,5)."3";
	
	my $aastart =$line[1];
	my $aaend =$line[2];
	my $chain=$line[3];
	my $seq =$line[4];
	my $length =$line[5];
	
	my $ss = $pdb.".ss";
			
	my $Nterm = substr($seq,0,4);
	my $Cterm = substr($seq,-4);
	my $inter= "";
	
	
	my $pos1 = substr($seq,1,1);
	my $pos2 = substr($Cterm,2,1);
	
	#print "$pdb $seq : $pos1 - $pos2 \n";
	
	
	if ($pos1 eq "R" || $pos1 eq "K") {
		
		if ($pos2 eq "D"){
			$nRKYD++;
			#print "   $pdb  :$pos1 - $pos2 \n";
			printf RKYD "$_\n";
		}
		else{
			$nRKOD++;
			$nRKO++;
			printf RKOD "$_\n";
			printf RKO "$_\n";
		}
    }
	else{
		if ($pos2 eq "D"){
			$nDO++;
			$nRKOD++;
			printf DO "$_\n";
			printf RKOD "$_\n";
		}
		else{
			$nRKND++;
			printf RKND "$_\n";
		}
	}
}


print "Total: $nline\n";
print "R/K y D: $nRKYD\n";
print "R/K o D: $nRKOD\n";

print "Solo R/K: $nRKO\n";
print "Solo D: $nDO\n";

print "No R/K ni D: $nRKND\n";


close(RKND);
close(DO);
close(RKO);
close(RKOD);
close(RKYD);


close(DATA);









