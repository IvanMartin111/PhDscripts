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

my %H3maps;
my $thr=1;
my $limitsize=8;



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
	open(MAPTHR,'>',"Mapa_".$mapname."_$thr.txt");
	printf MAPTHR "PHI ; PSI ; VALUE \n";
	
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
			
			my $ibin=-180+($i*5)+2.5;
			my $jbin=-180+($j*5)+2.5;
			
			if ($value > $limitvalue) {
                $H3maps{$mapname}{$i}{$j}=$value;
				printf MAPTHR "$ibin; $jbin; ".$value*$thr." \n";
            }
			else{
				$H3maps{$mapname}{$i}{$j}=0;
				printf MAPTHR "$ibin; $jbin; 0 \n";
			}
		}
	}
	close(MAPTHR);
}

####################################################


close(MAPS);




#################### Kink list brick ###############
open (DATA,$file) or die "The file $file doesn't exit in this path.\n";

open(KINKNY,'>',"Kink_NT_yes_$thr.fa");
open(KINKNN,'>',"Kink_NT_no_$thr.fa");

open(KINKCY,'>',"Kink_CT_yes_$thr.fa");
open(KINKCN,'>',"Kink_CT_no_$thr.fa");


open(KINKTY,'>',"Kink_T_yes_$thr.fa");
open(KINKTN,'>',"Kink_T_no_$thr.fa");






open(KINCTMAP100,'>',"Kink_CT_Map100_$thr.csv");
printf KINCTMAP100 "PHI; PSI;\n";
open(KINCTMAP101,'>',"Kink_CT_Map101_$thr.csv");
printf KINCTMAP101 "PHI; PSI;\n";
open(KINCTMAP102,'>',"Kink_CT_Map102_$thr.csv");
printf KINCTMAP102 "PHI; PSI;\n";
open(KINCTMAP103,'>',"Kink_CT_Map103_$thr.csv");
printf KINCTMAP103 "PHI; PSI;\n";






open(KINNTMAP92,'>',"Kink_NT_Map92_$thr.csv");
printf KINNTMAP92 "PHI; PSI;\n";
open(KINNTMAP93,'>',"Kink_NT_Map93_$thr.csv");
printf KINNTMAP93 "PHI; PSI;\n";
open(KINNTMAP94,'>',"Kink_NT_Map94_$thr.csv");
printf KINNTMAP94 "PHI; PSI;\n";
open(KINNTMAP95,'>',"Kink_NT_Map95_$thr.csv");
printf KINNTMAP95 "PHI; PSI;\n";





#Counters

my $countCT_Y=0;
my $countCT_N=0;
my $countNT_Y=0;
my $countNT_N=0;
my $countT_Y=0;
my $countT_N=0;
my $countT_NorC=0;


my $nline = 0; # Line index (counter)

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
	
	my $a92na;
	my $a93na;
	my $a94na;
	my $a95na;
	
	my $a103na;
	my $a102na;
	my $a101na;
	my $a100na;
	
	my %dihedrals;
	
	
	if ($length >$limitsize) {
        #print "$pdb $length $chain $Nterm $Cterm: $a92 $a93 $a94 $a95 $a100 $a101 $a102 $a103 \n";
	
		#Take the dihedrals
		open (SS,$ss) or die "The file $ss doesn't exit in this path.\n";
	
		while(<SS>)
		{
			next if /^#/; # reading only non-# beginning lines
			$_=~ s/^\t//;
				chop $_;
			my @liness = split(/\s+/, $_);
			
			
			if ($liness[1] eq $a92 && $liness[2] eq $chain) {
				$dihedrals{92}{"PHI"}=$liness[3];
				$dihedrals{92}{"PSI"}=$liness[4];
				$a92na=$aa{$liness[0]};		
            }
			elsif($liness[1] eq $a93 && $liness[2] eq $chain){
				$dihedrals{93}{"PHI"}=$liness[3];
				$dihedrals{93}{"PSI"}=$liness[4];
				$a93na=$aa{$liness[0]};
			}
			elsif($liness[1] eq $a94 && $liness[2] eq $chain){
				$dihedrals{94}{"PHI"}=$liness[3];
				$dihedrals{94}{"PSI"}=$liness[4];
				$a94na=$aa{$liness[0]};
			}
			elsif($liness[1] eq $a95 && $liness[2] eq $chain){
				$dihedrals{95}{"PHI"}=$liness[3];
				$dihedrals{95}{"PSI"}=$liness[4];
				$a95na=$aa{$liness[0]};
			}
			elsif($liness[1] eq $a100 && $liness[2] eq $chain){
				$dihedrals{100}{"PHI"}=$liness[3];
				$dihedrals{100}{"PSI"}=$liness[4];
				$a100na=$aa{$liness[0]};
			}
			elsif($liness[1] eq $a101 && $liness[2] eq $chain){
				$dihedrals{101}{"PHI"}=$liness[3];
				$dihedrals{101}{"PSI"}=$liness[4];
				$a101na=$aa{$liness[0]};
			}
			elsif($liness[1] eq $a102 && $liness[2] eq $chain){
				$dihedrals{102}{"PHI"}=$liness[3];
				$dihedrals{102}{"PSI"}=$liness[4];
				$a102na=$aa{$liness[0]};
			}
			elsif($liness[1] eq $a103 && $liness[2] eq $chain){
				$dihedrals{103}{"PHI"}=$liness[3];
				$dihedrals{103}{"PSI"}=$liness[4];
				$a103na=$aa{$liness[0]};
				last;
			}
            
		}
		close(SS);
	
	
	
	
		#Identify kinks N and T terminal
	
		my $Nkink=1;
		my $Ckink=1;
		my $Tkink=1;
	
		
		my @Nterm=(92,93,94,95);
		my @Cterm=(100,101,102,103);
		my @Tterm=(92,93,94,95,100,101,102,103);
		
		foreach my $aa (@Nterm){
			my $mapi=floor(($dihedrals{$aa}{"PHI"})/5)+36;
			#my $mapj=71-(floor(($dihedrals{$aa}{"PSI"})/5)+36);
			my $mapj=floor(($dihedrals{$aa}{"PSI"})/5)+36;

			#print "$aa: $dihedrals{$aa}{PHI}  $dihedrals{$aa}{PSI} -- $mapi , $mapj || $H3maps{$aa}{$mapi}{$mapj}\n";
			
			
			if ($H3maps{$aa}{$mapi}{$mapj} >0) {
               #print "YESSS $aa: $dihedrals{$aa}{PHI} $dihedrals{$aa}{PSI} -- $mapi , $mapj  -- $H3maps{$aa}{$mapi}{$mapj}\n";
            }
			else{
				#print "nope $aa: $dihedrals{$aa}{PHI}  -- $mapi , $mapj --- $H3maps{$aa}{$mapi}{$mapj}\n";
				$Nkink=0;
			}	
		}
		
		
		if ($Nkink) {
            #print " -N- $pdb $length $chain : $a92 $a93 $a94 $a95 $a100 $a101 $a102 $a103 \n";
			
			printf KINKNY ">$pdb\n";
			printf KINKNY $a92na.$a93na.$a94na.$a95na."\n";
			$countNT_Y++;
			
			printf KINNTMAP92 $dihedrals{92}{"PHI"}."; ".$dihedrals{92}{"PSI"}.";\n";
			printf KINNTMAP93 $dihedrals{93}{"PHI"}."; ".$dihedrals{93}{"PSI"}.";\n";
			printf KINNTMAP94 $dihedrals{94}{"PHI"}."; ".$dihedrals{94}{"PSI"}.";\n";
			printf KINNTMAP95 $dihedrals{95}{"PHI"}."; ".$dihedrals{95}{"PSI"}.";\n";
			
        }
		else{
			printf KINKNN ">$pdb\n";
			printf KINKNN $a92na.$a93na.$a94na.$a95na."\n";
			$countNT_N++;
		}
		
		
        
		foreach my $aa (@Cterm){
			my $mapi=floor(($dihedrals{$aa}{"PHI"})/5)+36;
			#my $mapj=71-(floor(($dihedrals{$aa}{"PSI"})/5)+36);
			my $mapj=floor(($dihedrals{$aa}{"PSI"})/5)+36;

			#print "$aa: $dihedrals{$aa}{PHI}  -- $mapi , $mapj\n";
			
			
			if ($H3maps{$aa}{$mapi}{$mapj} >0) {
               #print "YESSS $aa: $dihedrals{$aa}{PHI}  -- $mapi , $mapj  -- $H3maps{$aa}{$mapi}{$mapj}\n";
            }
			else{
				#print "nope $aa: $dihedrals{$aa}{PHI}  -- $mapi , $mapj --- $H3maps{$aa}{$mapi}{$mapj}\n";
				$Ckink=0;
			}	
		}
		
		if ($Ckink) {
            #print " -C- $pdb $length $chain : $a92 $a93 $a94 $a95 $a100 $a101 $a102 $a103 \n";
			printf KINKCY ">$pdb\n";
			printf KINKCY $a100na.$a101na.$a102na.$a103na."\n";
			
			printf KINCTMAP100 $dihedrals{100}{"PHI"}."; ".$dihedrals{100}{"PSI"}.";\n";
			printf KINCTMAP101 $dihedrals{101}{"PHI"}."; ".$dihedrals{101}{"PSI"}.";\n";
			printf KINCTMAP102 $dihedrals{102}{"PHI"}."; ".$dihedrals{102}{"PSI"}.";\n";
			printf KINCTMAP103 $dihedrals{103}{"PHI"}."; ".$dihedrals{103}{"PSI"}.";\n";
			
			
			$countCT_Y++;
        }
		else{
			printf KINKCN ">$pdb\n";
			printf KINKCN $a100na.$a101na.$a102na.$a103na."\n";
			$countCT_N++;
		}
		
		if ($Ckink && $Nkink) {
			#print " -B- $pdb $length $chain : $a92 $a93 $a94 $a95 $a100 $a101 $a102 $a103 \n";
			printf KINKTY ">$pdb\n";
			printf KINKTY $a92na.$a93na.$a94na.$a95na."  ".$a100na.$a101na.$a102na.$a103na."\n";
			$countT_Y++;
		}
		else{
			printf KINKTN ">$pdb\n";
			printf KINKTN $a92na.$a93na.$a94na.$a95na."  ".$a100na.$a101na.$a102na.$a103na."\n";
			$countT_N++;
		}
		
		if ($Ckink || $Nkink) {
			$countT_NorC++;
		}

    }

}
close(KINNTMAP95);
close(KINNTMAP94);
close(KINNTMAP93);
close(KINNTMAP92);

close(KINCTMAP103);
close(KINCTMAP102);
close(KINCTMAP101);
close(KINCTMAP100);



close(KINKCN);
close(KINKCY);
close(KINKNN);
close(KINKNY);
close(DATA);


my $totalnum=$countNT_Y+$countNT_N;
my $percNT=($countNT_Y/$totalnum)*100;
my $percCT=($countCT_Y/$totalnum)*100;
my $percT=($countT_Y/$totalnum)*100;
my $percNorC=($countT_NorC/$totalnum)*100;

print "Threshold: $thr  -- Limit size: $limitsize\n";
print "Total structures: $totalnum - Kink N-Terminal: $countNT_Y  - Kink C-Terminal: $countCT_Y  - Kink both sides: $countT_Y  - Kink alt least 1 side: $countT_NorC\n";
printf "Percentajes            - Kink N-Terminal: %.2f - Kink C-Terminal: %.2f - Kink both sides: %.2f - Kink alt least 1 side: %.2f\n",$percNT,$percCT,$percT,$percNorC;

####################################################





