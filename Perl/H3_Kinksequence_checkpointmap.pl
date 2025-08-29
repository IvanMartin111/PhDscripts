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

my @thr=(0.95,0.9,0.85,0.8,0.75,0.7,0.65,0.6,0.55,0.5);
my $limitsize=7;


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
			
			
			if ($liness[1] eq $a92 && $liness[2] eq $chain) {
				$dihedrals{$pdb_name}{92}{"PHI"}=$liness[3];
				$dihedrals{$pdb_name}{92}{"PSI"}=$liness[4];
			}
			elsif($liness[1] eq $a93 && $liness[2] eq $chain){
				$dihedrals{$pdb_name}{93}{"PHI"}=$liness[3];
				$dihedrals{$pdb_name}{93}{"PSI"}=$liness[4];
			}
			elsif($liness[1] eq $a94 && $liness[2] eq $chain){
				$dihedrals{$pdb_name}{94}{"PHI"}=$liness[3];
				$dihedrals{$pdb_name}{94}{"PSI"}=$liness[4];
			}
			elsif($liness[1] eq $a95 && $liness[2] eq $chain){
				$dihedrals{$pdb_name}{95}{"PHI"}=$liness[3];
				$dihedrals{$pdb_name}{95}{"PSI"}=$liness[4];
			}
			elsif($liness[1] eq $a100 && $liness[2] eq $chain){
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



#repeat the same procedure for each thr and have data off all
foreach my $thr (@thr){
	
	my %H3maps;

			
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
					$H3maps{$mapname}{$i}{$j}=$value;
				}
				else{
					$H3maps{$mapname}{$i}{$j}=0;
				}
			}
		}
	}
	close(MAPS);
	
	####################################################
	
	
	
	#	
	####################################################
	
	my @Nterm=(92,93,94,95);
	my @Cterm=(100,101,102,103);
	my @Tterm=(92,93,94,95,100,101,102,103);
		
	$results{$thr}{"Total"}=0;
	$results{$thr}{"Nterm"}=0;
	$results{$thr}{"Cterm"}=0;
		
	foreach my $aa (@Tterm){	
		$results{$thr}{$aa}=0;
	}
	
	foreach my $aa (@Tterm){
		foreach my $bb (@Tterm){
			$resultscc{$thr}{$aa}{$bb}=0;
		}
	}
	
	
	foreach my $pdb (keys %dihedrals){
		my $Nkink=1;
		my $Ckink=1;
		my $Tkink=1;
		my %aacontrol =(92 =>0,93 =>0,94 =>0,95 =>0);
		my %bbcontrol =(100 =>0,101 =>0,102 =>0,103 =>0);
		
		foreach my $aa (@Nterm){
			my $mapi=floor(($dihedrals{$pdb}{$aa}{"PHI"})/5)+36;
			my $mapj=floor(($dihedrals{$pdb}{$aa}{"PSI"})/5)+36;
		
			
			if ($H3maps{$aa}{$mapi}{$mapj} >0) {
				$results{$thr}{$aa}++;
				$aacontrol{$aa}=1;
            }
			else{
				#print "nope $aa: $dihedrals{$aa}{PHI}  -- $mapi , $mapj --- $H3maps{$aa}{$mapi}{$mapj}\n";
				$Nkink=0;
			}
		}
		foreach my $a(sort keys %aacontrol){
			foreach my $b(sort keys %aacontrol){
				if ($aacontrol{$a} && $aacontrol{$b}) {
                    $resultscc{$thr}{$a}{$b}++;
                }
			}
		}
		
		if ($Nkink) {$results{$thr}{"Nterm"}++}
		
		
		
		
		foreach my $aa (@Cterm){
			my $mapi=floor(($dihedrals{$pdb}{$aa}{"PHI"})/5)+36;
			my $mapj=floor(($dihedrals{$pdb}{$aa}{"PSI"})/5)+36;
		
			
			if ($H3maps{$aa}{$mapi}{$mapj} >0) {
				$results{$thr}{$aa}++;
				$bbcontrol{$aa}=1;
            }
			else{
				#print "nope $aa: $dihedrals{$aa}{PHI}  -- $mapi , $mapj --- $H3maps{$aa}{$mapi}{$mapj}\n";
				$Ckink=0;
			}
		}
		
		foreach my $a(sort keys %bbcontrol){
			foreach my $b(sort keys %bbcontrol){
				if ($bbcontrol{$a} && $bbcontrol{$b}) {
                    $resultscc{$thr}{$a}{$b}++;
                }
			}
		}
		if ($Ckink) {$results{$thr}{"Cterm"}++}
		
		
        if ($Nkink && $Ckink ) {$results{$thr}{"Total"}++}
		
	}
	
	
	
	#Printings

	print "\n\tTreshold: $thr\n";
	
	print "- Total : $results{$thr}{Total}\n";
	printf "- Nkink : $results{$thr}{Nterm}\t  -  %.2f \t\t\t\t92\t93\t94\t95\t\t\t\t\t93\t\t\t  94\t\t\t\t95\t  size:\t%.2f\n",$results{$thr}{"Nterm"}/$nline*100,($nline-$nsizeline)/$nline*100;
	
	printf "\t. 92 : $results{$thr}{92}\t  -  %.2f \t\t92\t%.2f\t%.2f\t%.2f\t%.2f\t\t92\t%.2f - (%.2f|%.2f|%.2f) %.2f - (%.2f|%.2f|%.2f) %.2f - (%.2f|%.2f|%.2f)\n",$results{$thr}{"92"}/$nline*100,$resultscc{$thr}{"92"}{"92"}/$nline*100,$resultscc{$thr}{"92"}{"93"}/$nline*100,$resultscc{$thr}{"92"}{"94"}/$nline*100,$resultscc{$thr}{"92"}{"95"}/$nline*100,($nsizeline-$resultscc{$thr}{"92"}{"93"})/$nline*100,(($nsizeline-$resultscc{$thr}{"92"}{"93"})-($nsizeline-$results{$thr}{"93"}))/($nsizeline-$resultscc{$thr}{"92"}{"93"})*100,(($nsizeline-$results{$thr}{"92"})+($nsizeline-$results{$thr}{"93"})-($nsizeline-$resultscc{$thr}{"92"}{"93"}))/($nsizeline-$resultscc{$thr}{"92"}{"93"})*100,(($nsizeline-$resultscc{$thr}{"92"}{"93"})-($nsizeline-$results{$thr}{"92"}))/($nsizeline-$resultscc{$thr}{"92"}{"93"})*100,($nsizeline-$resultscc{$thr}{"92"}{"94"})/$nline*100,(($nsizeline-$resultscc{$thr}{"92"}{"94"})-($nsizeline-$results{$thr}{"94"}))/($nsizeline-$resultscc{$thr}{"92"}{"94"})*100,(($nsizeline-$results{$thr}{"92"})+($nsizeline-$results{$thr}{"94"})-($nsizeline-$resultscc{$thr}{"92"}{"94"}))/($nsizeline-$resultscc{$thr}{"92"}{"94"})*100,(($nsizeline-$resultscc{$thr}{"92"}{"94"})-($nsizeline-$results{$thr}{"92"}))/($nsizeline-$resultscc{$thr}{"92"}{"94"})*100,($nsizeline-$resultscc{$thr}{"92"}{"95"})/$nline*100,(($nsizeline-$resultscc{$thr}{"92"}{"95"})-($nsizeline-$results{$thr}{"95"}))/($nsizeline-$resultscc{$thr}{"92"}{"95"})*100,(($nsizeline-$results{$thr}{"92"})+($nsizeline-$results{$thr}{"95"})-($nsizeline-$resultscc{$thr}{"92"}{"95"}))/($nsizeline-$resultscc{$thr}{"92"}{"95"})*100,(($nsizeline-$resultscc{$thr}{"92"}{"95"})-($nsizeline-$results{$thr}{"92"}))/($nsizeline-$resultscc{$thr}{"92"}{"95"})*100;                           
	printf "\t. 93 : $results{$thr}{93}\t  -  %.2f \t\t93\t%.2f\t%.2f\t%.2f\t%.2f\t\t93\t\t\t\t   %.2f - (%.2f|%.2f|%.2f) %.2f - (%.2f|%.2f|%.2f)\n",$results{$thr}{"93"}/$nline*100,$resultscc{$thr}{"93"}{"92"}/$nline*100,$resultscc{$thr}{"93"}{"93"}/$nline*100,$resultscc{$thr}{"93"}{"94"}/$nline*100,$resultscc{$thr}{"93"}{"95"}/$nline*100,($nsizeline-$resultscc{$thr}{"93"}{"94"})/$nline*100,(($nsizeline-$resultscc{$thr}{"93"}{"94"})-($nsizeline-$results{$thr}{"94"}))/($nsizeline-$resultscc{$thr}{"93"}{"94"})*100,(($nsizeline-$results{$thr}{"93"})+($nsizeline-$results{$thr}{"94"})-($nsizeline-$resultscc{$thr}{"93"}{"94"}))/($nsizeline-$resultscc{$thr}{"93"}{"94"})*100,(($nsizeline-$resultscc{$thr}{"93"}{"94"})-($nsizeline-$results{$thr}{"93"}))/($nsizeline-$resultscc{$thr}{"93"}{"94"})*100,($nsizeline-$resultscc{$thr}{"93"}{"95"})/$nline*100,(($nsizeline-$resultscc{$thr}{"93"}{"95"})-($nsizeline-$results{$thr}{"95"}))/($nsizeline-$resultscc{$thr}{"93"}{"95"})*100,(($nsizeline-$results{$thr}{"93"})+($nsizeline-$results{$thr}{"95"})-($nsizeline-$resultscc{$thr}{"93"}{"95"}))/($nsizeline-$resultscc{$thr}{"93"}{"95"})*100,(($nsizeline-$resultscc{$thr}{"93"}{"95"})-($nsizeline-$results{$thr}{"93"}))/($nsizeline-$resultscc{$thr}{"93"}{"95"})*100;
	printf "\t. 94 : $results{$thr}{94}\t  -  %.2f \t\t94\t%.2f\t%.2f\t%.2f\t%.2f\t\t94\t\t\t\t     \t\t\t       %.2f - (%.2f|%.2f|%.2f)\n",$results{$thr}{"94"}/$nline*100,$resultscc{$thr}{"94"}{"92"}/$nline*100,$resultscc{$thr}{"94"}{"93"}/$nline*100,$resultscc{$thr}{"94"}{"94"}/$nline*100,$resultscc{$thr}{"94"}{"95"}/$nline*100,($nsizeline-$resultscc{$thr}{"94"}{"95"})/$nline*100,(($nsizeline-$resultscc{$thr}{"94"}{"95"})-($nsizeline-$results{$thr}{"95"}))/($nsizeline-$resultscc{$thr}{"94"}{"95"})*100,(($nsizeline-$results{$thr}{"94"})+($nsizeline-$results{$thr}{"95"})-($nsizeline-$resultscc{$thr}{"94"}{"95"}))/($nsizeline-$resultscc{$thr}{"94"}{"95"})*100,(($nsizeline-$resultscc{$thr}{"94"}{"95"})-($nsizeline-$results{$thr}{"94"}))/($nsizeline-$resultscc{$thr}{"94"}{"95"})*100;
	printf "\t. 95 : $results{$thr}{95}\t  -  %.2f \t\t95\t%.2f\t%.2f\t%.2f\t%.2f\n",$results{$thr}{"95"}/$nline*100,$resultscc{$thr}{"95"}{"92"}/$nline*100,$resultscc{$thr}{"95"}{"93"}/$nline*100,$resultscc{$thr}{"95"}{"94"}/$nline*100,$resultscc{$thr}{"95"}{"95"}/$nline*100;
	
	printf "- Ckink : $results{$thr}{Cterm}\t  -  %.2f \t\t\t\t100\t101\t102\t103\t\t\t\t\t101\t\t\t  102\t\t\t\t103\t  size:\t%.2f	\n",$results{$thr}{"Cterm"}/$nline*100,($nline-$nsizeline)/$nline*100;
	
	printf "\t. 100 : $results{$thr}{100}\t  -  %.2f \t\t100\t%.2f\t%.2f\t%.2f\t%.2f\t\t100\t%.2f - (%.2f|%.2f|%.2f) %.2f - (%.2f|%.2f|%.2f) %.2f - (%.2f|%.2f|%.2f)\n",$results{$thr}{"100"}/$nline*100,$resultscc{$thr}{"100"}{"100"}/$nline*100,$resultscc{$thr}{"100"}{"101"}/$nline*100,$resultscc{$thr}{"100"}{"102"}/$nline*100,$resultscc{$thr}{"100"}{"103"}/$nline*100,($nsizeline-$resultscc{$thr}{"100"}{"101"})/$nline*100,(($nsizeline-$resultscc{$thr}{"100"}{"101"})-($nsizeline-$results{$thr}{"101"}))/($nsizeline-$resultscc{$thr}{"100"}{"101"})*100,(($nsizeline-$results{$thr}{"100"})+($nsizeline-$results{$thr}{"101"})-($nsizeline-$resultscc{$thr}{"100"}{"101"}))/($nsizeline-$resultscc{$thr}{"100"}{"101"})*100,(($nsizeline-$resultscc{$thr}{"100"}{"101"})-($nsizeline-$results{$thr}{"100"}))/($nsizeline-$resultscc{$thr}{"100"}{"101"})*100,($nsizeline-$resultscc{$thr}{"100"}{"102"})/$nline*100,(($nsizeline-$resultscc{$thr}{"100"}{"102"})-($nsizeline-$results{$thr}{"102"}))/($nsizeline-$resultscc{$thr}{"100"}{"102"})*100,(($nsizeline-$results{$thr}{"100"})+($nsizeline-$results{$thr}{"102"})-($nsizeline-$resultscc{$thr}{"100"}{"102"}))/($nsizeline-$resultscc{$thr}{"100"}{"102"})*100,(($nsizeline-$resultscc{$thr}{"100"}{"102"})-($nsizeline-$results{$thr}{"100"}))/($nsizeline-$resultscc{$thr}{"100"}{"102"})*100,($nsizeline-$resultscc{$thr}{"100"}{"103"})/$nline*100,(($nsizeline-$resultscc{$thr}{"100"}{"103"})-($nsizeline-$results{$thr}{"103"}))/($nsizeline-$resultscc{$thr}{"100"}{"103"})*100,(($nsizeline-$results{$thr}{"100"})+($nsizeline-$results{$thr}{"103"})-($nsizeline-$resultscc{$thr}{"100"}{"103"}))/($nsizeline-$resultscc{$thr}{"100"}{"103"})*100,(($nsizeline-$resultscc{$thr}{"100"}{"103"})-($nsizeline-$results{$thr}{"100"}))/($nsizeline-$resultscc{$thr}{"100"}{"103"})*100;                           
	printf "\t. 101 : $results{$thr}{101}\t  -  %.2f \t\t101\t%.2f\t%.2f\t%.2f\t%.2f\t\t101\t\t\t\t   %.2f - (%.2f|%.2f|%.2f) %.2f - (%.2f|%.2f|%.2f)\n",$results{$thr}{"101"}/$nline*100,$resultscc{$thr}{"101"}{"100"}/$nline*100,$resultscc{$thr}{"101"}{"101"}/$nline*100,$resultscc{$thr}{"101"}{"102"}/$nline*100,$resultscc{$thr}{"101"}{"103"}/$nline*100,($nsizeline-$resultscc{$thr}{"101"}{"102"})/$nline*100,(($nsizeline-$resultscc{$thr}{"101"}{"102"})-($nsizeline-$results{$thr}{"102"}))/($nsizeline-$resultscc{$thr}{"101"}{"102"})*100,(($nsizeline-$results{$thr}{"101"})+($nsizeline-$results{$thr}{"102"})-($nsizeline-$resultscc{$thr}{"101"}{"102"}))/($nsizeline-$resultscc{$thr}{"101"}{"102"})*100,(($nsizeline-$resultscc{$thr}{"101"}{"102"})-($nsizeline-$results{$thr}{"101"}))/($nsizeline-$resultscc{$thr}{"101"}{"102"})*100,($nsizeline-$resultscc{$thr}{"101"}{"103"})/$nline*100,(($nsizeline-$resultscc{$thr}{"101"}{"103"})-($nsizeline-$results{$thr}{"103"}))/($nsizeline-$resultscc{$thr}{"101"}{"103"})*100,(($nsizeline-$results{$thr}{"101"})+($nsizeline-$results{$thr}{"103"})-($nsizeline-$resultscc{$thr}{"101"}{"103"}))/($nsizeline-$resultscc{$thr}{"101"}{"103"})*100,(($nsizeline-$resultscc{$thr}{"101"}{"103"})-($nsizeline-$results{$thr}{"101"}))/($nsizeline-$resultscc{$thr}{"101"}{"103"})*100;
	printf "\t. 102 : $results{$thr}{102}\t  -  %.2f \t\t102\t%.2f\t%.2f\t%.2f\t%.2f\t\t102\t\t\t\t     \t\t\t       %.2f - (%.2f|%.2f|%.2f)\n",$results{$thr}{"102"}/$nline*100,$resultscc{$thr}{"102"}{"100"}/$nline*100,$resultscc{$thr}{"102"}{"101"}/$nline*100,$resultscc{$thr}{"102"}{"102"}/$nline*100,$resultscc{$thr}{"102"}{"103"}/$nline*100,($nsizeline-$resultscc{$thr}{"102"}{"103"})/$nline*100,(($nsizeline-$resultscc{$thr}{"102"}{"103"})-($nsizeline-$results{$thr}{"103"}))/($nsizeline-$resultscc{$thr}{"102"}{"103"})*100,(($nsizeline-$results{$thr}{"102"})+($nsizeline-$results{$thr}{"103"})-($nsizeline-$resultscc{$thr}{"102"}{"103"}))/($nsizeline-$resultscc{$thr}{"102"}{"103"})*100,(($nsizeline-$resultscc{$thr}{"102"}{"103"})-($nsizeline-$results{$thr}{"102"}))/($nsizeline-$resultscc{$thr}{"102"}{"103"})*100;
	printf "\t. 103 : $results{$thr}{103}\t  -  %.2f \t\t103\t%.2f\t%.2f\t%.2f\t%.2f\n",$results{$thr}{"103"}/$nline*100,$resultscc{$thr}{"103"}{"100"}/$nline*100,$resultscc{$thr}{"103"}{"101"}/$nline*100,$resultscc{$thr}{"103"}{"102"}/$nline*100,$resultscc{$thr}{"103"}{"103"}/$nline*100;
	
	
	
	

	
	
	

}





