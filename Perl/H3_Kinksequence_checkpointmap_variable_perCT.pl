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
if( !($#ARGV == 1))
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
my $H3map = $ARGV[1];


my %AA = (A=>'ALA',Y=>'TYR',M=>'MET',L=>'LEU',C=>'CYS',G=>'GLY',
         R=>'ARG',N=>'ASN',D=>'ASP',Q=>'GLN',E=>'GLU',H=>'HIS',W=>'TRP',
         K=>'LYS',F=>'PHE',P=>'PRO',S=>'SER',T=>'THR',I=>'ILE',V=>'VAL');


my %aa = reverse %AA;

my @thr=(0.95,0.9,0.8,0.7,0.6);
my $limitsize=8;


my %dihedrals;
my %results;
my %resultscc;
my %rcontrols; 




#################### Kink list brick ###############
open (DATA,$file) or die "The file $file doesn't exit in this path.\n";

my $nline = 0; # Line index (counter)
my $nsizeline=0; # Line index with the size filter (counter)

while(<DATA>)
{
	next if /^#/; # reading only non-# beginning lines
	$_=~ s/^\t//;
	chop $_;
	my @line = split(/\s+/, $_);
	#push(@data,\@line);
	$nline++;
	my $pdb = substr($line[0],0,5)."3";
	my $pdb_name=substr($line[0],0,4);
	my $aastart =$line[1];
	my $aaend =$line[2];
	my $chain=$line[3];
	my $seq =$line[4];
	my $length =$line[5];
	
	my $ss = $pdb.".ss";
			
	my $Nterm = substr($seq,0,4);
	my $Cterm = substr($seq,-4);
	my $inter= "";
	
	#print "$seq -- $Nterm -- $inter --$Cterm\n";
	
	#if ($length => 10) {}

	
	#### Values for each loop
	
	my $a92=$aastart -1;
	my $a93=$aastart;
	my $a94=$aastart +1;
	my $a95=$aastart +2;
	
	my $a103=$aaend +1;
	my $a102=$aaend;
	my $a101=$aaend -1;
	my $a100=$aaend -2;
	
	
	
	if ($length >$limitsize) {
		#print "$pdb $length $chain $Nterm $Cterm: $a92 $a93 $a94 $a95 $a100 $a101 $a102 $a103 \n";
	
		#Take the dihedrals
		open (SS,$ss) or die "The file $ss doesn't exit in this path.\n";
		$nsizeline++;
		
		while(<SS>)
		{
			next if /^#/; # reading only non-# beginning lines
			$_=~ s/^\t//;
				chop $_;
			my @liness = split(/\s+/, $_);
			
			
			
			if($liness[1] eq $a100 && $liness[2] eq $chain){
				$dihedrals{$pdb_name}{100}{"PHI"}=$liness[3];
				$dihedrals{$pdb_name}{100}{"PSI"}=$liness[4];
			}
			elsif($liness[1] eq $a101 && $liness[2] eq $chain){
				$dihedrals{$pdb_name}{101}{"PHI"}=$liness[3];
				$dihedrals{$pdb_name}{101}{"PSI"}=$liness[4];
			}
			elsif($liness[1] eq $a102 && $liness[2] eq $chain){
				$dihedrals{$pdb_name}{102}{"PHI"}=$liness[3];
				$dihedrals{$pdb_name}{102}{"PSI"}=$liness[4];
			}
			elsif($liness[1] eq $a103 && $liness[2] eq $chain){
				$dihedrals{$pdb_name}{103}{"PHI"}=$liness[3];
				$dihedrals{$pdb_name}{103}{"PSI"}=$liness[4];
				last;
			}
			
		}
		close(SS);
	}
}
close (DATA);


print "Total sctuctures: $nline\n";
print "Total sctuctures bigger than limit size $limitsize: $nsizeline\n";


my %H3maps;
#repeat the same procedure for each thr and have data off all
foreach my $thr (@thr){

			
	########################## H3 brick ###############
	open (MAPS,$H3map) or die "The file $H3map doesn't exit in this path.\n";
	
	while(<MAPS>)
	{
		next if /^#/; # reading only non-# beginning lines
		$_=~ s/^\t//;
		chop $_;
		my @line = split(/\s+/, $_);
		#push(@data,\@line);
		my $mapname = shift @line;
		#print "$mapname - $line[0] \n";
		
		my @sortedmap=reverse sort @line;
		
		#print "$mapname - $sortedmap[0]\n";
		
		my $acum=0;
		my $limitvalue;
		foreach my $i (@sortedmap){
			if ($acum < $thr) {
				$acum =$acum + $i;
			}
			else{last;}
			$limitvalue=$i;
		}
		
		my $aaa=scalar(@sortedmap);
			
		for(my $i = 0; $i <= 71; $i++){
			for(my $j = 0; $j <= 71; $j++){
				my $value=shift @line;
				if ($value > $limitvalue) {
					$H3maps{$thr}{$mapname}{$i}{$j}=$value;
				}
				else{
					$H3maps{$thr}{$mapname}{$i}{$j}=0;
				}
			}
		}
	}
	close(MAPS);
	
	####################################################
}










foreach my $thr (@thr){
	foreach my $thr2 (@thr){
		foreach my $thr3 (@thr){
			foreach my $thr4 (@thr){	
				
				#	
				####################################################
				
				my @Nterm=(92,93,94,95);
				my @Cterm=(100,101,102,103);
				my @Tterm=(92,93,94,95,100,101,102,103);
					
				$results{$thr}{"Total"}=0;
				$results{$thr}{$thr2}{$thr3}{$thr4}{"Nterm"}=0;
				$results{$thr}{"Cterm"}=0;
					
				
				foreach my $pdb (keys %dihedrals){
					my $Nkink=1;
					my $Ckink=1;
					my $Tkink=1;

					
					#100
					my $mapi=floor(($dihedrals{$pdb}{100}{"PHI"})/5)+36;
					my $mapj=floor(($dihedrals{$pdb}{100}{"PSI"})/5)+36;
						
					if ($H3maps{$thr}{100}{$mapi}{$mapj} >0) {}
					else{$Nkink=0;}
					
					#101
					$mapi=floor(($dihedrals{$pdb}{101}{"PHI"})/5)+36;
					$mapj=floor(($dihedrals{$pdb}{101}{"PSI"})/5)+36;
						
					if ($H3maps{$thr2}{101}{$mapi}{$mapj} >0) {}
					else{$Nkink=0;}
					
					#102
					$mapi=floor(($dihedrals{$pdb}{102}{"PHI"})/5)+36;
					$mapj=floor(($dihedrals{$pdb}{102}{"PSI"})/5)+36;
						
					if ($H3maps{$thr3}{102}{$mapi}{$mapj} >0) {}
					else{$Nkink=0;}
					
					#103
					$mapi=floor(($dihedrals{$pdb}{103}{"PHI"})/5)+36;
					$mapj=floor(($dihedrals{$pdb}{103}{"PSI"})/5)+36;
						
					if ($H3maps{$thr4}{103}{$mapi}{$mapj} >0) {}
					else{$Nkink=0;}


					
					if ($Nkink) {$results{$thr}{$thr2}{$thr3}{$thr4}{"Nterm"}++}
					
					

					
				}
				
				
				
				#Printings
			
				#print "\n\tTreshold: $thr $thr2 $thr3 $thr4 ---$results{$thr}{$thr2}{$thr3}{$thr4}{Nterm} \n";

			}
		}
	}
}


#my @thr=(0.95,0.9,0.8,0.7,0.6);
#Printings
foreach my $thr (@thr){
	
	print "\n -100: $thr\n";
	
	print "\t 101: 0.95\t\t\t\t 101: 0.9\t\t\t\t 101: 0.8\t\t\t\t 101: 0.7\t\t\t\t 101: 0.6\n";
	


	printf "   %.2f %.2f %.2f %.2f %.2f \t%.2f %.2f %.2f %.2f %.2f \t%.2f %.2f %.2f %.2f %.2f \t%.2f %.2f %.2f %.2f %.2f \t%.2f %.2f %.2f %.2f %.2f \n",$results{$thr}{0.95}{0.95}{0.95}{"Nterm"}/$nsizeline*100,$results{$thr}{0.95}{0.9}{0.95}{"Nterm"}/$nsizeline*100,$results{$thr}{0.95}{0.8}{0.95}{"Nterm"}/$nsizeline*100,$results{$thr}{0.95}{0.7}{0.95}{"Nterm"}/$nsizeline*100,$results{$thr}{0.95}{0.6}{0.95}{"Nterm"}/$nsizeline*100,$results{$thr}{0.9}{0.95}{0.95}{"Nterm"}/$nsizeline*100,$results{$thr}{0.9}{0.9}{0.95}{"Nterm"}/$nsizeline*100,$results{$thr}{0.9}{0.8}{0.95}{"Nterm"}/$nsizeline*100,$results{$thr}{0.9}{0.7}{0.95}{"Nterm"}/$nsizeline*100,$results{$thr}{0.9}{0.6}{0.95}{"Nterm"}/$nsizeline*100,$results{$thr}{0.8}{0.95}{0.95}{"Nterm"}/$nsizeline*100,$results{$thr}{0.8}{0.9}{0.95}{"Nterm"}/$nsizeline*100,$results{$thr}{0.8}{0.8}{0.95}{"Nterm"}/$nsizeline*100,$results{$thr}{0.8}{0.7}{0.95}{"Nterm"}/$nsizeline*100,$results{$thr}{0.8}{0.6}{0.95}{"Nterm"}/$nsizeline*100,$results{$thr}{0.7}{0.95}{0.95}{"Nterm"}/$nsizeline*100,$results{$thr}{0.7}{0.9}{0.95}{"Nterm"}/$nsizeline*100,$results{$thr}{0.7}{0.8}{0.95}{"Nterm"}/$nsizeline*100,$results{$thr}{0.7}{0.7}{0.95}{"Nterm"}/$nsizeline*100,$results{$thr}{0.7}{0.6}{0.95}{"Nterm"}/$nsizeline*100,$results{$thr}{0.6}{0.95}{0.95}{"Nterm"}/$nsizeline*100,$results{$thr}{0.6}{0.9}{0.95}{"Nterm"}/$nsizeline*100,$results{$thr}{0.6}{0.8}{0.95}{"Nterm"}/$nsizeline*100,$results{$thr}{0.6}{0.7}{0.95}{"Nterm"}/$nsizeline*100,$results{$thr}{0.6}{0.6}{0.95}{"Nterm"}/$nsizeline*100; 

	printf "   %.2f %.2f %.2f %.2f %.2f \t%.2f %.2f %.2f %.2f %.2f \t%.2f %.2f %.2f %.2f %.2f \t%.2f %.2f %.2f %.2f %.2f \t%.2f %.2f %.2f %.2f %.2f \n",$results{$thr}{0.95}{0.95}{0.9}{"Nterm"}/$nsizeline*100,$results{$thr}{0.95}{0.9}{0.9}{"Nterm"}/$nsizeline*100,$results{$thr}{0.95}{0.8}{0.9}{"Nterm"}/$nsizeline*100,$results{$thr}{0.95}{0.7}{0.9}{"Nterm"}/$nsizeline*100,$results{$thr}{0.95}{0.6}{0.9}{"Nterm"}/$nsizeline*100,$results{$thr}{0.9}{0.95}{0.9}{"Nterm"}/$nsizeline*100,$results{$thr}{0.9}{0.9}{0.9}{"Nterm"}/$nsizeline*100,$results{$thr}{0.9}{0.8}{0.9}{"Nterm"}/$nsizeline*100,$results{$thr}{0.9}{0.7}{0.9}{"Nterm"}/$nsizeline*100,$results{$thr}{0.9}{0.6}{0.9}{"Nterm"}/$nsizeline*100,$results{$thr}{0.8}{0.95}{0.9}{"Nterm"}/$nsizeline*100,$results{$thr}{0.8}{0.9}{0.9}{"Nterm"}/$nsizeline*100,$results{$thr}{0.8}{0.8}{0.9}{"Nterm"}/$nsizeline*100,$results{$thr}{0.8}{0.7}{0.9}{"Nterm"}/$nsizeline*100,$results{$thr}{0.8}{0.6}{0.9}{"Nterm"}/$nsizeline*100,$results{$thr}{0.7}{0.95}{0.9}{"Nterm"}/$nsizeline*100,$results{$thr}{0.7}{0.9}{0.9}{"Nterm"}/$nsizeline*100,$results{$thr}{0.7}{0.8}{0.9}{"Nterm"}/$nsizeline*100,$results{$thr}{0.7}{0.7}{0.9}{"Nterm"}/$nsizeline*100,$results{$thr}{0.7}{0.6}{0.9}{"Nterm"}/$nsizeline*100,$results{$thr}{0.6}{0.95}{0.9}{"Nterm"}/$nsizeline*100,$results{$thr}{0.6}{0.9}{0.9}{"Nterm"}/$nsizeline*100,$results{$thr}{0.6}{0.8}{0.9}{"Nterm"}/$nsizeline*100,$results{$thr}{0.6}{0.7}{0.9}{"Nterm"}/$nsizeline*100,$results{$thr}{0.6}{0.6}{0.9}{"Nterm"}/$nsizeline*100; 

	printf "   %.2f %.2f %.2f %.2f %.2f \t%.2f %.2f %.2f %.2f %.2f \t%.2f %.2f %.2f %.2f %.2f \t%.2f %.2f %.2f %.2f %.2f \t%.2f %.2f %.2f %.2f %.2f \n",$results{$thr}{0.95}{0.95}{0.8}{"Nterm"}/$nsizeline*100,$results{$thr}{0.95}{0.9}{0.8}{"Nterm"}/$nsizeline*100,$results{$thr}{0.95}{0.8}{0.8}{"Nterm"}/$nsizeline*100,$results{$thr}{0.95}{0.7}{0.8}{"Nterm"}/$nsizeline*100,$results{$thr}{0.95}{0.6}{0.8}{"Nterm"}/$nsizeline*100,$results{$thr}{0.9}{0.95}{0.8}{"Nterm"}/$nsizeline*100,$results{$thr}{0.9}{0.9}{0.8}{"Nterm"}/$nsizeline*100,$results{$thr}{0.9}{0.8}{0.8}{"Nterm"}/$nsizeline*100,$results{$thr}{0.9}{0.7}{0.8}{"Nterm"}/$nsizeline*100,$results{$thr}{0.9}{0.6}{0.8}{"Nterm"}/$nsizeline*100,$results{$thr}{0.8}{0.95}{0.8}{"Nterm"}/$nsizeline*100,$results{$thr}{0.8}{0.9}{0.8}{"Nterm"}/$nsizeline*100,$results{$thr}{0.8}{0.8}{0.8}{"Nterm"}/$nsizeline*100,$results{$thr}{0.8}{0.7}{0.8}{"Nterm"}/$nsizeline*100,$results{$thr}{0.8}{0.6}{0.8}{"Nterm"}/$nsizeline*100,$results{$thr}{0.7}{0.95}{0.8}{"Nterm"}/$nsizeline*100,$results{$thr}{0.7}{0.9}{0.8}{"Nterm"}/$nsizeline*100,$results{$thr}{0.7}{0.8}{0.8}{"Nterm"}/$nsizeline*100,$results{$thr}{0.7}{0.7}{0.8}{"Nterm"}/$nsizeline*100,$results{$thr}{0.7}{0.6}{0.8}{"Nterm"}/$nsizeline*100,$results{$thr}{0.6}{0.95}{0.8}{"Nterm"}/$nsizeline*100,$results{$thr}{0.6}{0.9}{0.8}{"Nterm"}/$nsizeline*100,$results{$thr}{0.6}{0.8}{0.8}{"Nterm"}/$nsizeline*100,$results{$thr}{0.6}{0.7}{0.8}{"Nterm"}/$nsizeline*100,$results{$thr}{0.6}{0.6}{0.8}{"Nterm"}/$nsizeline*100; 

	printf "   %.2f %.2f %.2f %.2f %.2f \t%.2f %.2f %.2f %.2f %.2f \t%.2f %.2f %.2f %.2f %.2f \t%.2f %.2f %.2f %.2f %.2f \t%.2f %.2f %.2f %.2f %.2f \n",$results{$thr}{0.95}{0.95}{0.7}{"Nterm"}/$nsizeline*100,$results{$thr}{0.95}{0.9}{0.7}{"Nterm"}/$nsizeline*100,$results{$thr}{0.95}{0.8}{0.7}{"Nterm"}/$nsizeline*100,$results{$thr}{0.95}{0.7}{0.7}{"Nterm"}/$nsizeline*100,$results{$thr}{0.95}{0.6}{0.7}{"Nterm"}/$nsizeline*100,$results{$thr}{0.9}{0.95}{0.7}{"Nterm"}/$nsizeline*100,$results{$thr}{0.9}{0.9}{0.7}{"Nterm"}/$nsizeline*100,$results{$thr}{0.9}{0.8}{0.7}{"Nterm"}/$nsizeline*100,$results{$thr}{0.9}{0.7}{0.7}{"Nterm"}/$nsizeline*100,$results{$thr}{0.9}{0.6}{0.7}{"Nterm"}/$nsizeline*100,$results{$thr}{0.8}{0.95}{0.7}{"Nterm"}/$nsizeline*100,$results{$thr}{0.8}{0.9}{0.7}{"Nterm"}/$nsizeline*100,$results{$thr}{0.8}{0.8}{0.7}{"Nterm"}/$nsizeline*100,$results{$thr}{0.8}{0.7}{0.7}{"Nterm"}/$nsizeline*100,$results{$thr}{0.8}{0.6}{0.7}{"Nterm"}/$nsizeline*100,$results{$thr}{0.7}{0.95}{0.7}{"Nterm"}/$nsizeline*100,$results{$thr}{0.7}{0.9}{0.7}{"Nterm"}/$nsizeline*100,$results{$thr}{0.7}{0.8}{0.7}{"Nterm"}/$nsizeline*100,$results{$thr}{0.7}{0.7}{0.7}{"Nterm"}/$nsizeline*100,$results{$thr}{0.7}{0.6}{0.7}{"Nterm"}/$nsizeline*100,$results{$thr}{0.6}{0.95}{0.7}{"Nterm"}/$nsizeline*100,$results{$thr}{0.6}{0.9}{0.7}{"Nterm"}/$nsizeline*100,$results{$thr}{0.6}{0.8}{0.7}{"Nterm"}/$nsizeline*100,$results{$thr}{0.6}{0.7}{0.7}{"Nterm"}/$nsizeline*100,$results{$thr}{0.6}{0.6}{0.7}{"Nterm"}/$nsizeline*100; 

	printf "   %.2f %.2f %.2f %.2f %.2f \t%.2f %.2f %.2f %.2f %.2f \t%.2f %.2f %.2f %.2f %.2f \t%.2f %.2f %.2f %.2f %.2f \t%.2f %.2f %.2f %.2f %.2f \n",$results{$thr}{0.95}{0.95}{0.6}{"Nterm"}/$nsizeline*100,$results{$thr}{0.95}{0.9}{0.6}{"Nterm"}/$nsizeline*100,$results{$thr}{0.95}{0.8}{0.6}{"Nterm"}/$nsizeline*100,$results{$thr}{0.95}{0.7}{0.6}{"Nterm"}/$nsizeline*100,$results{$thr}{0.95}{0.6}{0.6}{"Nterm"}/$nsizeline*100,$results{$thr}{0.9}{0.95}{0.6}{"Nterm"}/$nsizeline*100,$results{$thr}{0.9}{0.9}{0.6}{"Nterm"}/$nsizeline*100,$results{$thr}{0.9}{0.8}{0.6}{"Nterm"}/$nsizeline*100,$results{$thr}{0.9}{0.7}{0.6}{"Nterm"}/$nsizeline*100,$results{$thr}{0.9}{0.6}{0.6}{"Nterm"}/$nsizeline*100,$results{$thr}{0.8}{0.95}{0.6}{"Nterm"}/$nsizeline*100,$results{$thr}{0.8}{0.9}{0.6}{"Nterm"}/$nsizeline*100,$results{$thr}{0.8}{0.8}{0.6}{"Nterm"}/$nsizeline*100,$results{$thr}{0.8}{0.7}{0.6}{"Nterm"}/$nsizeline*100,$results{$thr}{0.8}{0.6}{0.6}{"Nterm"}/$nsizeline*100,$results{$thr}{0.7}{0.95}{0.6}{"Nterm"}/$nsizeline*100,$results{$thr}{0.7}{0.9}{0.6}{"Nterm"}/$nsizeline*100,$results{$thr}{0.7}{0.8}{0.6}{"Nterm"}/$nsizeline*100,$results{$thr}{0.7}{0.7}{0.6}{"Nterm"}/$nsizeline*100,$results{$thr}{0.7}{0.6}{0.6}{"Nterm"}/$nsizeline*100,$results{$thr}{0.6}{0.95}{0.6}{"Nterm"}/$nsizeline*100,$results{$thr}{0.6}{0.9}{0.6}{"Nterm"}/$nsizeline*100,$results{$thr}{0.6}{0.8}{0.6}{"Nterm"}/$nsizeline*100,$results{$thr}{0.6}{0.7}{0.6}{"Nterm"}/$nsizeline*100,$results{$thr}{0.6}{0.6}{0.6}{"Nterm"}/$nsizeline*100; 




}


