#!/usr/bin/perl

#Script for the generation of all the inputs of a CDR folder with only proteins in chothia format from the SabDab DB

#Iván Martín Hernández

# Date: Apr - 2021


use strict;
use warnings;
use diagnostics;
use Data::Dumper qw(Dumper);
use feature 'say';
use Switch;
use Scalar::Util qw(reftype);

# INPUT PARSER
if( !($#ARGV == 0))
{
	say "USAGE:\n\t$0  <input1> \n";
	say "\t<input1>  --> File with the tsv (modify) of all the proteins of SabDab DB";
	say "\nDESCRIPTION:\n\tSay something here...";
	say "\nWARNING:\n\tSomething....";
	say "\nEXAMPLE:";
	say "\tAB_SAbDab_sumarize.pl  sabdab_summary_all.tsv";
	say "\nHELP: some help\n";
	exit;
}

my %short = (
   "ALA" => "A",
	"CYS" => "C",
   "CYX" => "C",
   "ASP" => "D",
   "GLU" => "E",
   "PHE" => "F",
   "GLY" => "G",
   "HIS" => "H",
   "HSD" => "H",
   "HSE" => "H",
   "HSP" => "H",
   "HIP" => "H", # Rosseta
   "HID" => "H", # Rosseta
   "HIE" => "H", # Rosseta
   "ILE" => "I",
   "LYS" => "K",
   "LEU" => "L",
   "MET" => "M",
   "ASN" => "N",
   "PRO" => "P",
   "GLN" => "Q",
   "ARG" => "R",
   "SER" => "S",
   "THR" => "T",
   "VAL" => "V",
   "TRP" => "W",
   "TYR" => "Y",
   "  A" => "A",
   "  G" => "G",
   "  T" => "T",
   "  C" => "C",
   "  U" => "U",
   " DA" => "A",
   " DG" => "G",
   " DT" => "T",
   " DC" => "C" ) ;

#take the inputs from the bash
my $rawtable = $ARGV[0];




###########################
# Open the raw table and save all the information 
###########################



open (DATA,$rawtable) or die "The file $rawtable doesn't exit in this path.\n";

my %rawtable;
my @keys;

my $nline = 0; # Line index (counter)

while(<DATA>)
{
	next if /^#/; # reading only non-# beginning lines
	$_=~ s/^\t//;
	chop $_;
	my @line = split(/\t/, $_);
	#push(@data,\@line);
	$nline++;

   my $pdb = substr($line[0],0,4);
	
	if ($nline eq 1) {@keys=@line;}
	else{
		foreach (@keys){
			$rawtable{$pdb}{$_}=shift(@line);
			#print $_." ".$rawtable{$pdb}{$_}."\n";
		}
	}
	if ($rawtable{$pdb}{resolution} eq "NOT") {
        $rawtable{$pdb}{resolution} = 50;
    }
    
	if (!($rawtable{$pdb}{heavy_species} eq "")) {
      $rawtable{$pdb}{especie}=$rawtable{$pdb}{heavy_species};
   }
   else {
		$rawtable{$pdb}{especie}=$rawtable{$pdb}{light_species};
	}
    
}
close (DATA);

my $dimraw=scalar keys %rawtable;
print "$rawtable have information of ".$dimraw." different structures\n";




###########################
# Take the info of all the pdb in the folder in chothia format (obteined directly from de SAbDab database)
###########################

# Loading things
my $entire_name = "*.pdb";
my @files = glob("$entire_name");
my $outfile ="AB_loop_RCD.txt";
open(OUTPUT,'>',$outfile);
open(OUTPUTJSON,'>',"AB_loop.json");
open(OUTPUTCDRINCL,'>',"GAP_CDR2.txt");



printf OUTPUTJSON  "[\n";

foreach my $PDB(@files){
   #say $PDB;
	my $namepdb = substr($PDB,0,4);
	open (PDB, $PDB) or die "Can't open file $PDB.\n";
	my $newpdb=substr($PDB,0,4)."_2.pdb";
	open(OUTPUT2,'>',$newpdb);
	
	#print "Resolution: ".$rawtable{$namepdb}{resolution}."\n";
	
	#Definition of the variable of the final pdb object
	my %ABs;
	my $ABtype;
	my $numAB=0;
	my $marca1=0;
	my $marca2=0;
	my $marcaFF;
	my $haveantigen="NO";
	my $haveH="NO";
	my $haveL="NO";
	my @Hchn;
	my @Lchn;

	
	#Variables for the sequences and H3 length
	my $L1_seq="";
	my $L2_seq="";
	my $L3_seq="";
	my $H1_seq="";
	my $H2_seq="";
	my $H3_seq="";
	my $L1_len=0;
	my $L2_len=0;
	my $L3_len=0;
	my $H1_len=0;
	my $H2_len=0;
	my $H3_len=0;
	my $L1_control=0;
	my $L2_control=0;
	my $L3_control=0;
	my $H1_control=0;
	my $H2_control=0;
	my $H3_control=0;
	
	my $pre_res_isnum=1;	
	my $cur_res = -1;
	my $new_res_num=0;
	my $cur_chain ="";
	my $pre_alt =" ";
	my $InChainH=0;		#Marks if we are in H chain to procces de Hx loops
	my $InChainL=0;		#Marks if we are in L chain to procces de Lx loops
	my @pre_coords = (0,0,0);
		
		
	
	
	
	
	say OUTPUT2 "REMARK   5 REVERSE OF THE CHOTHIA NOMENCLATURE TO A CONSECUTIVE NUMERATION";

	
	while(<PDB>)
	{
		#print "$_\n";
		chop $_;
		
		
		#part that read the REMARK lines and take primeary info
		if (/^REMARK/){
			if (/PAIRED_HL/ or /SINGLE/ ){
				#print "$_";
				if (/PAIRED_HL/ ){
					$ABtype="PAIRED_HL";
					$marca1=1;
					$numAB++;
					}
				if (/SINGLE/ ){
					$ABtype="SINGLE";
					$marca2=1;
					$numAB++;
					}
				
				my @line = split(/\s+/, $_);
				foreach my $element(@line){
					my @values = split(/=/, $element);
					
					if ($values[0] eq "HCHAIN") {
						#print "$values[0]  $values[1]\n";
						$ABs{$numAB}{HCHAIN}=$values[1];
						push(@Hchn,$values[1]);
						$haveH="YES";
						#if ($values[1] eq "NA"){print "$namepdb $values[0] $values[1] BBBBBBBBBBBBBBBBB\n";}
               }
               if ($values[0] eq "LCHAIN") {
						#print "$values[0]  $values[1]\n";
						$ABs{$numAB}{LCHAIN}=$values[1];
						push(@Lchn,$values[1]);
						$haveL="YES";
						#if ($values[1] eq "NA"){print "$namepdb $values[0] $values[1] BBBBBBBBBBBBBBBBB\n";}
               }
					if ($values[0] eq "AGCHAIN") {
						#print "$values[0]  $values[1]\n";
						$ABs{$numAB}{AGCHAIN}=$values[1];
						$haveantigen="YES";
						if ($values[1] eq "NONE"){$haveantigen="NO";}
						#if ($values[1] eq "NONE"){print "$namepdb $values[0] $values[1] CCCCCCCCCCCCCCCCCCCCCCCCCCC\n";}
						#if ($values[1] eq "NA"){print "$namepdb $values[0] $values[1] NANANANAN\n";}
               }
					if ($values[0] eq "AGTYPE") {
						#print "$values[0]  $values[1]\n";
						$ABs{$numAB}{AGTYPE}=$values[1];
						#if ($values[1] eq "NONE"){print "$namepdb $values[0] $values[1] CCCCCCCCCCCCCCCCCCCCCCCCCCC\n";}
               }					
				}
			}
			say OUTPUT2 $_;
		}
		
		
			
		#part that read the ATOM lines and take the sequences info
	
		if (/^ATOM/ || /^HETATM/) {
			my $format = readPDBLine($_);
			#my @fields = split( / +/, $line); 
			#my $res_num = $fields[5];
			my $res_num=substr($_,22,5); # get residue number
			$res_num=~s/\s+//g;
			#my $chain =$fields[4];
			#my $chain =substr($line,21,1);
			#my $res =substr($line,17,3);
			#my $res =$fields[3];
			my $chain =$format->{chain};
			my $res =$format->{residueId};
			#print "$res_num - $chain - $res\n";
				
			my $alt = $format->{achar};
			my @coords = ($format->{x},$format->{y},$format->{z});
			
			#print "$res_num - $chain - $res - dis: $distance\n";
				
			#Determine if we change of chaen and the type of it to the consecuence process	
			if (!($cur_chain eq $chain)){
				#say $line;
				$cur_res = -1;
				$new_res_num=0;
				$pre_res_isnum=1;	
				$cur_chain = $chain;
					
						
				if ($chain ~~ @Hchn){
					#print "$_ \n $chain  ";
					#print "H -- $namepdb - $chain - $res - $res_num - $alt\n";
					$InChainH=1;
					$InChainL=0;
					$H1_seq="";
					$H2_seq="";
					$H3_seq="";
					$H1_len=0;
					$H2_len=0;
					$H3_len=0;

					#print "H -- $InChainH -- L --$InChainL\n";
				}
				else{
					if ($chain ~~ @Lchn){
						#print "$_ \n $chain  ";
						#print "L -- $namepdb - $chain - $res - $res_num - $alt\n";
						$InChainH=0;
						$InChainL=1;
						$L1_seq="";
						$L2_seq="";
						$L3_seq="";
						$L1_len=0;
						$L2_len=0;
						$L3_len=0;
						#print "H -- $InChainH -- L --$InChainL\n";
					}
					else{
						#print "NO CDR -- $namepdb - $chain - $res - $res_num - $alt\n";
						$InChainH=0;
						$InChainL=0;
						#print "H -- $InChainH -- L --$InChainL\n";
					}
				}
			}
					
					
					
			#Determine if the residue is new for change variables, or if we are in other atom of the same residue	
			if (!($res_num eq $cur_res)){
				#print "$cur_res - $res_num - $new_res_num - $short{$res}\n";
				
				my $distance= sqrt(($pre_coords[0]-$coords[0])**2+($pre_coords[1]-$coords[1])**2+($pre_coords[2]-$coords[2])**2);
				#print "$res_num - $new_res_num- $chain - $res - dis: $distance\n";

				
				
				#Check if the residue value is a number or if it have a letter following the chlothia sintaxis, and process the nest values of the renumeration depends of it
				if ($res_num =~ /^[+-]?\d+$/ && $pre_res_isnum && $cur_res != -1) {
					$pre_res_isnum=1;
					my $res_dif=$res_num -$cur_res;  
					$cur_res =$res_num;
					
					#print "$res_num Is a number  y la diff es $res_dif \n";

					if ($res_dif >1) {
						if ($InChainL == 1 || $InChainH == 1) {
							if ($L1_control eq 1 || $L2_control eq 1 || $L3_control eq 1 || $H1_control eq 1 || $H2_control eq 1 || $H3_control eq 1) {
								if ($distance >4.2) {
									$new_res_num =$new_res_num+$res_dif;
									printf OUTPUTCDRINCL "$namepdb".".pdb\t$chain\t$res_num\n";
                        }
                        else {
									$new_res_num =$new_res_num+1;

								}       

							}
							else{
								$new_res_num =$new_res_num+$res_dif;
							}
						}
						else{
							$new_res_num =$new_res_num+$res_dif;
	
						}
					}
					else{
						$new_res_num =$new_res_num+1;
					}

			   } else {
					if ($res_num =~ /^[+-]?\d+$/ ){
						#print "$res_num Is a number Camino 2222\n";
						$new_res_num =$new_res_num+1;
						$cur_res =$res_num;
						$pre_res_isnum=1;
						if (($L1_control eq 1 || $L2_control eq 1 || $L3_control eq 1 || $H1_control eq 1 || $H2_control eq 1 || $H3_control eq 1) && $distance >4.2) {
							printf OUTPUTCDRINCL "$namepdb".".pdb\t$chain\t$res_num\n";
						}
					}
					else {
						#print "$res_num Is not a number\n";
					
						if (($L1_control eq 1 || $L2_control eq 1 || $L3_control eq 1 || $H1_control eq 1 || $H2_control eq 1 || $H3_control eq 1) && $distance >4.2) {
							printf OUTPUTCDRINCL "$namepdb".".pdb\t$chain\t$res_num\n";
							$new_res_num =$new_res_num+1;
							$cur_res =$res_num;
							$pre_res_isnum=0;
                  }
						else{
							$new_res_num =$new_res_num+1;
							$cur_res =$res_num;
							$pre_res_isnum=0;
						}
					
					}
			   }
				
				@pre_coords=@coords;
				
            
				
				#say $res;
				$pre_alt =" ";
				
				if ($InChainL == 1){
					#print "$cur_res - $res_num - $new_res_num - $short{$res}\n";
					#say "aa -$res_num-";
					#L1
					if ($res_num eq 24){
						#say "$res - $short{$res}";
						$ABs{"Lchains"}{$chain}{"L1_start"}=$new_res_num;
						$L1_control=1;
					}
					if($L1_control){
						#say "$res - $short{$res}";
						my $new =$short{$res};
						$L1_seq = $L1_seq.$new;
						$L1_len++;
					}
					if ($L1_control eq 1 && $res_num ge 34){
						$L1_control=0;
						#print "$namepdb -- $chain -- $InChainL -- $InChainH -- $res_num -- len $L1_len -- seq $L1_seq\n";
						$ABs{"Lchains"}{$chain}{"L1seq"}=$L1_seq;
						$ABs{"Lchains"}{$chain}{"L1len"}=$L1_len;
						$ABs{"Lchains"}{$chain}{"L1_end"}=$new_res_num;
					}
					
					#L2
					if ($res_num eq 50){
						$ABs{"Lchains"}{$chain}{"L2_start"}=$new_res_num;
;
						$L2_control=1;
					}
					if($L2_control){
						#say "$res - $short{$res}";
						my $new =$short{$res};
						$L2_seq = $L2_seq.$new;
						$L2_len++;
					}
					if ($L2_control eq 1 && $res_num ge 56){
						$L2_control=0;
						#print "$namepdb -- $chain -- $InChainL -- $InChainH -- $res_num -- len $L2_len -- seq $L2_seq\n";
						$ABs{"Lchains"}{$chain}{"L2seq"}=$L2_seq;
						$ABs{"Lchains"}{$chain}{"L2len"}=$L2_len;
						$ABs{"Lchains"}{$chain}{"L2_end"}=$new_res_num;
					}
					
					#L3
					if ($res_num eq 89){
						$ABs{"Lchains"}{$chain}{"L3_start"}=$new_res_num;
						$L3_control=1;
					}
					if($L3_control){
						#say "$res - $short{$res}";
						my $new =$short{$res};
						$L3_seq = $L3_seq.$new;
						$L3_len++;
					}
					if ($L3_control eq 1 && $res_num ge 97){
						$L3_control=0;
						#print "$namepdb -- $chain -- $InChainL -- $InChainH -- $res_num -- len $L3_len -- seq $L3_seq\n";
						$ABs{"Lchains"}{$chain}{"L3seq"}=$L3_seq;
						$ABs{"Lchains"}{$chain}{"L3len"}=$L3_len;
						$ABs{"Lchains"}{$chain}{"L3_end"}=$new_res_num;
					}
					
									
				}
				elsif  ($InChainH == 1){
					#H1
					if ($res_num eq 26){
						$H1_control=1;
						$ABs{"Hchains"}{$chain}{"H1_start"}=$new_res_num;
					}
					if($H1_control){
						#say "$res - $short{$res}";
						my $new =$short{$res};
						$H1_seq = $H1_seq.$new;
						$H1_len++;
					}
					if ($H1_control eq 1 && $res_num ge 32){
						$H1_control=0;
						#print "$namepdb -- $chain -- $InChainL -- $InChainH -- $res_num -- len $H1_len -- seq $H1_seq\n";
						$ABs{"Hchains"}{$chain}{"H1seq"}=$H1_seq;
						$ABs{"Hchains"}{$chain}{"H1len"}=$H1_len;
						$ABs{"Hchains"}{$chain}{"H1_end"}=$new_res_num;
					}
					
					#H2
					if ($res_num eq 52){
						$ABs{"Hchains"}{$chain}{"H2_start"}=$new_res_num;
						$H2_control=1;
					}
					if($H2_control){
						#say "$res - $short{$res}";
						my $new =$short{$res};
						$H2_seq = $H2_seq.$new;
						$H2_len++;
					}
					if ($H2_control eq 1 && $res_num ge 56){
						$H2_control=0;
						#print "$namepdb -- $chain -- $InChainL -- $InChainH -- $res_num -- len $H2_len -- seq $H2_seq\n";
						$ABs{"Hchains"}{$chain}{"H2seq"}=$H2_seq;
						$ABs{"Hchains"}{$chain}{"H2len"}=$H2_len;
						$ABs{"Hchains"}{$chain}{"H2_end"}=$new_res_num;
					}
					
					#H3
					if ($res_num eq 93){
						$ABs{"Hchains"}{$chain}{"H3_start"}=$new_res_num;
						$H3_control=1;
					}
					if($H3_control){
						#say "$res - $short{$res}";
						my $new =$short{$res};
						$H3_seq = $H3_seq.$new;
						$H3_len++;
					}
					if ($H3_control eq 1 && $res_num eq 102){
						$H3_control=0;
						#print "$namepdb -- $chain -- $InChainL -- $InChainH -- $res_num -- len $H3_len -- seq $H3_seq\n";
						$ABs{"Hchains"}{$chain}{"H3seq"}=$H3_seq;
						$ABs{"Hchains"}{$chain}{"H3len"}=$H3_len;
						$ABs{"Hchains"}{$chain}{"H3_end"}=$new_res_num;
					}
				}
			}
					
			if (!($alt eq " ")){
				#say "NOW - $alt";
				if ($alt eq "A"){
					$pre_alt=$alt;
				}
				else {
					if ($pre_alt eq "A"){
						$pre_alt=$alt;
					}
					else{
						#$pre_alt=$alt;
						$format->{altloc} = " ";
					}
				}
				#say "ALT - $format->{altloc}";
			}
			
			$format->{resseq}=$new_res_num;
			
			printf OUTPUT2 "%6s%5s %4s%s%3s %s%4s%s   %8s%8s%8s%6s%6s\n",
	      $format->{tag},
	      $format->{seq},
	      $format->{name},
	      $format->{altloc},
	      $format->{residueId},
	      $format->{chain},
	      $format->{resseq},
	      " ",
	      $format->{x},
	      $format->{y},
	      $format->{z},
	      $format->{occ},
	      $format->{Bfact};	
		
		}
	}
	
	
	
	
	#Identify strctures with Single an Paired AB in the same file
	$marcaFF=$marca1+$marca2;
	if ($marcaFF eq 2) {
		$ABtype="DUAL";
    }

 
	 
	 
	#Part when we print all the information of the protein in the screen
	#Also printe the info of the loops for the RCO input in the OUPUT file
	print "\n\t\t\tPDB:\t $namepdb\n\n";
	print "Method: $rawtable{$namepdb}{method}\t";
	print "Organism: $rawtable{$namepdb}{organism}\n";
	print "Resolution: $rawtable{$namepdb}{resolution}\t";
	print "r_factor: $rawtable{$namepdb}{r_factor}\t";
	print "Date: $rawtable{$namepdb}{date}\t";
	print "\n\tAB type: $ABtype \tNº ABs: $numAB\t";
	print "Have antigen: $haveantigen\n";
	if ($haveH eq "YES") {
		print "\n\tChains H:\n";
      foreach my $key2 (sort keys %{$ABs{Hchains}}) {
			print "\t$key2\n";
			if ($ABs{Hchains}{$key2}{H1len}) {
            print "\t\tH1 length:  $ABs{Hchains}{$key2}{H1len}\t- sequence: $ABs{Hchains}{$key2}{H1seq}";
				print "\t\tnew start:  $ABs{Hchains}{$key2}{H1_start}\t- new end: $ABs{Hchains}{$key2}{H1_end} \n";
				printf OUTPUT  "%10s  %3d  %3d   %1s   %25s  %2d  H1\n", $newpdb , $ABs{Hchains}{$key2}{H1_start}, $ABs{Hchains}{$key2}{H1_end}, $key2,$ABs{Hchains}{$key2}{H1seq}, $ABs{Hchains}{$key2}{H1len};

         }
			
			if ($ABs{Hchains}{$key2}{H2len}) {
				print "\t\tH2 length:  $ABs{Hchains}{$key2}{H2len}\t- sequence: $ABs{Hchains}{$key2}{H2seq} ";
				print "\t\tnew start:  $ABs{Hchains}{$key2}{H2_start}\t- new end: $ABs{Hchains}{$key2}{H2_end} \n";
				printf OUTPUT  "%10s  %3d  %3d   %1s   %25s  %2d  H2\n", $newpdb , $ABs{Hchains}{$key2}{H2_start}, $ABs{Hchains}{$key2}{H2_end}, $key2,$ABs{Hchains}{$key2}{H2seq}, $ABs{Hchains}{$key2}{H2len};
			}
			if ($ABs{Hchains}{$key2}{H3len}) {
				print "\t\tH3 length:  $ABs{Hchains}{$key2}{H3len}\t- sequence: $ABs{Hchains}{$key2}{H3seq} ";
				print "\t\tnew start:  $ABs{Hchains}{$key2}{H3_start}\t- new end: $ABs{Hchains}{$key2}{H3_end} \n";
				printf OUTPUT  "%10s  %3d  %3d   %1s   %25s  %2d  H3\n", $newpdb , $ABs{Hchains}{$key2}{H3_start}, $ABs{Hchains}{$key2}{H3_end}, $key2,$ABs{Hchains}{$key2}{H3seq}, $ABs{Hchains}{$key2}{H3len};
			}
		} 
   }
	if ($haveL eq "YES") {
		print "\n\tChains L:\n";
		foreach my $key2 (sort keys %{$ABs{Lchains}}) {
			print "\t$key2\n";
			if ($ABs{Lchains}{$key2}{L1len}) {
				print "\t\tL1 length:  $ABs{Lchains}{$key2}{L1len}\t- sequence: $ABs{Lchains}{$key2}{L1seq} ";
				print "\t\tnew start:  $ABs{Lchains}{$key2}{L1_start}\t- new end: $ABs{Lchains}{$key2}{L1_end} \n";
				printf OUTPUT  "%10s  %3d  %3d   %1s   %25s  %2d  L1\n", $newpdb , $ABs{Lchains}{$key2}{L1_start}, $ABs{Lchains}{$key2}{L1_end}, $key2,$ABs{Lchains}{$key2}{L1seq}, $ABs{Lchains}{$key2}{L1len};
			}
			if ($ABs{Lchains}{$key2}{L2len}) {
				print "\t\tL2 length:  $ABs{Lchains}{$key2}{L2len}\t- sequence: $ABs{Lchains}{$key2}{L2seq} ";
				print "\t\tnew start:  $ABs{Lchains}{$key2}{L2_start}\t- new end: $ABs{Lchains}{$key2}{L2_end} \n";
				printf OUTPUT  "%10s  %3d  %3d   %1s   %25s  %2d  L2\n", $newpdb , $ABs{Lchains}{$key2}{L2_start}, $ABs{Lchains}{$key2}{L2_end}, $key2,$ABs{Lchains}{$key2}{L2seq}, $ABs{Lchains}{$key2}{L2len};
			}
			if ($ABs{Lchains}{$key2}{L3len}) {
				print "\t\tL3 length:  $ABs{Lchains}{$key2}{L3len}\t- sequence: $ABs{Lchains}{$key2}{L3seq} ";
				print "\t\tnew start:  $ABs{Lchains}{$key2}{L3_start}\t- new end: $ABs{Lchains}{$key2}{L3_end} \n";
				printf OUTPUT  "%10s  %3d  %3d   %1s   %25s  %2d  L3\n", $newpdb , $ABs{Lchains}{$key2}{L3_start}, $ABs{Lchains}{$key2}{L3_end}, $key2,$ABs{Lchains}{$key2}{L3seq}, $ABs{Lchains}{$key2}{L3len};
			}
		}
	}
	print "\n";
	
	
	
	
	# Part from printing the JSON file
	printf OUTPUTJSON "   {\n";
	printf OUTPUTJSON "      \"PDB\": \"$namepdb\",\n";
	printf OUTPUTJSON "      \"Method\": \"$rawtable{$namepdb}{method}\",\n";
	printf OUTPUTJSON "      \"Organism\": \"$rawtable{$namepdb}{especie}\",\n";
	printf OUTPUTJSON "      \"Resolution\": $rawtable{$namepdb}{resolution},\n";
	if (!(($rawtable{$namepdb}{r_factor} eq "NA") or ($rawtable{$namepdb}{r_factor} eq "unknown") or ($rawtable{$namepdb}{r_factor} eq ""))) {
      	printf OUTPUTJSON "      \"r_factor\": $rawtable{$namepdb}{r_factor},\n";
    }
	printf OUTPUTJSON "      \"Date\": \"$rawtable{$namepdb}{date}\",\n";
	printf OUTPUTJSON "      \"AB type\": \"$ABtype\", \n      \"Nº ABs\": $numAB,\n";
	printf OUTPUTJSON "      \"Have antigen\": \"$haveantigen\",\n";
	
	if ($haveH eq "YES") {
		printf OUTPUTJSON "      \"Chains H\": [\n";
		foreach my $key2 (sort keys %{$ABs{Hchains}}) {
			printf OUTPUTJSON "         {\n";
			printf OUTPUTJSON "            \"Chain\":\"$key2\",\n";
			if ($ABs{Hchains}{$key2}{H1len}) {
				printf OUTPUTJSON "            \"H1 length\": $ABs{Hchains}{$key2}{H1len},\n            \"H1 sequence\": \"$ABs{Hchains}{$key2}{H1seq}\",\n";
				printf OUTPUTJSON "            \"H1 new start\": $ABs{Hchains}{$key2}{H1_start},\n            \"H1 new end\": $ABs{Hchains}{$key2}{H1_end},\n";
			}
			if ($ABs{Hchains}{$key2}{H2len}) {
				printf OUTPUTJSON "            \"H2 length\": $ABs{Hchains}{$key2}{H2len},\n            \"H2 sequence\": \"$ABs{Hchains}{$key2}{H2seq}\",\n";
				printf OUTPUTJSON "            \"H2 new start\": $ABs{Hchains}{$key2}{H2_start},\n            \"H2 new end\": $ABs{Hchains}{$key2}{H2_end},\n";
			}
			if ($ABs{Hchains}{$key2}{H3len}) {
				printf OUTPUTJSON "            \"H3 length\": $ABs{Hchains}{$key2}{H3len},\n            \"H3 sequence\": \"$ABs{Hchains}{$key2}{H3seq}\",\n";
				printf OUTPUTJSON "            \"H3 new start\": $ABs{Hchains}{$key2}{H3_start},\n            \"H3 new end\": $ABs{Hchains}{$key2}{H3_end}\n";
			}
			printf OUTPUTJSON "         },\n";
		}
		printf OUTPUTJSON "      ],\n";
	}

	if ($haveL eq "YES") {
		printf OUTPUTJSON "      \"Chains L\": [\n";
		foreach my $key2 (sort keys %{$ABs{Lchains}}) {
			printf OUTPUTJSON "         {\n";
			printf OUTPUTJSON "            \"Chain\":\"$key2\",\n";
			if ($ABs{Lchains}{$key2}{L1len}) {
				printf OUTPUTJSON "            \"L1 length\": $ABs{Lchains}{$key2}{L1len},\n            \"L1 sequence\": \"$ABs{Lchains}{$key2}{L1seq}\",\n";
				printf OUTPUTJSON "            \"L1 new start\": $ABs{Lchains}{$key2}{L1_start},\n            \"L1 new end\": $ABs{Lchains}{$key2}{L1_end},\n";
			}
			if ($ABs{Lchains}{$key2}{L2len}) {
				printf OUTPUTJSON "            \"L2 length\": $ABs{Lchains}{$key2}{L2len},\n            \"L2 sequence\": \"$ABs{Lchains}{$key2}{L2seq}\",\n";
				printf OUTPUTJSON "            \"L2 new start\": $ABs{Lchains}{$key2}{L2_start},\n            \"L2 new end\": $ABs{Lchains}{$key2}{L2_end},\n";
			}
			if ($ABs{Lchains}{$key2}{L3len}) {
				printf OUTPUTJSON "            \"L3 length\": $ABs{Lchains}{$key2}{L3len},\n            \"L3 sequence\": \"$ABs{Lchains}{$key2}{L3seq}\",\n";
				printf OUTPUTJSON "            \"L3 new start\": $ABs{Lchains}{$key2}{L3_start},\n            \"L3 new end\": $ABs{Lchains}{$key2}{L3_end}\n";
			}
			printf OUTPUTJSON "         },\n";
		}
		printf OUTPUTJSON "      ],\n";
	}
	
	
	printf OUTPUTJSON "   }, \n";
	
	
	
	close (OUTPUT2);

}

printf OUTPUTJSON  "]";

close (OUTPUTCDRINCL);
close (OUTPUTJSON);
close (OUTPUT);









###########################
#     Refine RCD list
###########################




open (DATA,$outfile) or die "The file $outfile doesn't exit in this path.\n";

my $nline = 0; # Line index (counter)

my $newentry=0;
my $altentry=0;

my %H1;
my %H2;
my %H3;
my %L1;
my %L2;
my %L3;

my @altH1;
my @altH2;
my @altH3;
my @altL1;
my @altL2;
my @altL3;


#principal output
open(OUTPUTH1,'>',"AB_loop_RCD_H1.txt");
open(OUTPUTH2,'>',"AB_loop_RCD_H2.txt");
open(OUTPUTH3,'>',"AB_loop_RCD_H3.txt");
open(OUTPUTL1,'>',"AB_loop_RCD_L1.txt");
open(OUTPUTL2,'>',"AB_loop_RCD_L2.txt");
open(OUTPUTL3,'>',"AB_loop_RCD_L3.txt");

#Alternative output for those that have more than 1 seq in the same pdb
open(ALTOUTPUTH1,'>',"ALT_AB_loop_RCD_H1.txt");
open(ALTOUTPUTH2,'>',"ALT_AB_loop_RCD_H2.txt");
open(ALTOUTPUTH3,'>',"ALT_AB_loop_RCD_H3.txt");
open(ALTOUTPUTL1,'>',"ALT_AB_loop_RCD_L1.txt");
open(ALTOUTPUTL2,'>',"ALT_AB_loop_RCD_L2.txt");
open(ALTOUTPUTL3,'>',"ALT_AB_loop_RCD_L3.txt");



while(<DATA>)
{
	next if /^#/; # reading only non-# beginning lines
	$_=~ s/^\t//;
	chop $_;
	my @line = split(/\s+/, $_);
	#push(@data,\@line);
	$nline++;

		
	$newentry=0;
	$altentry=0;
	
	my $pdb = $line[0];
	my $loopseq = $line[4];
	my $loopCDR = $line[6];

	#print "$pdb - $loopseq - $loopCDR\n";
	
	
	
	switch($loopCDR) {
		case "H1" {
			#print "$pdb - $loopseq - $loopCDR\n";
			if (grep( /^$pdb/,keys(%H1) )) {		
				if ($H1{"$pdb"} eq $loopseq) {
                    #print "$loopseq - $H1{\"$pdb\"} -- Same sequence : Do nothing \n";
                }
				else {
					if (grep( /^$loopseq/,@altH1 )) {
						#print "ALTTTTT ------------- $loopseq - Same sequence : Do nothing \n";
					}
					else{
						push (@altH1, $loopseq);
						#print "$loopseq --- NEW ALTERNATIVE \n";
						printf ALTOUTPUTH1 "$_\n"
					}
				}
			}
			else{
				$H1{"$pdb"}=$loopseq;
				printf OUTPUTH1 "$_\n";
			}
		}
		
		case "H2" {
			#print "$pdb - $loopseq - $loopCDR\n";
			if (grep( /^$pdb/,keys(%H2) )) {		
				if ($H2{"$pdb"} eq $loopseq) {
                    #print "$loopseq - $H2{\"$pdb\"} -- Same sequence : Do nothing \n";
                }
				else {
					if (grep( /^$loopseq/,@altH2 )) {
						#print "ALTTTTT ------------- $loopseq - Same sequence : Do nothing \n";
					}
					else{
						push (@altH2, $loopseq);
						#print "$loopseq --- NEW ALTERNATIVE \n";
						printf ALTOUTPUTH2 "$_\n"
					}
				}
			}
			else{
				$H2{"$pdb"}=$loopseq;
				printf OUTPUTH2 "$_\n";
			}
		}
		
		case "H3" {
			#print "$pdb - $loopseq - $loopCDR\n";
			if (grep( /^$pdb/,keys(%H3) )) {		
				if ($H3{"$pdb"} eq $loopseq) {
                    #print "$loopseq - $H3{\"$pdb\"} -- Same sequence : Do nothing \n";
                }
				else {
					if (grep( /^$loopseq/,@altH3 )) {
						#print "ALTTTTT ------------- $loopseq - Same sequence : Do nothing \n";
					}
					else{
						push (@altH3, $loopseq);
						#print "$loopseq --- NEW ALTERNATIVE \n";
						printf ALTOUTPUTH3 "$_\n"
					}
				}
			}
			else{
				$H3{"$pdb"}=$loopseq;
				printf OUTPUTH3 "$_\n";
			}	
		}
		
		case "L1" {
			#print "$pdb - $loopseq - $loopCDR\n";
			if (grep( /^$pdb/,keys(%L1) )) {		
				if ($L1{"$pdb"} eq $loopseq) {
                    #print "$loopseq - $L1{\"$pdb\"} -- Same sequence : Do nothing \n";
                }
				else {
					if (grep( /^$loopseq/,@altL1 )) {
						#print "ALTTTTT ------------- $loopseq - Same sequence : Do nothing \n";
					}
					else{
						push (@altL1, $loopseq);
						#print "$loopseq --- NEW ALTERNATIVE \n";
						printf ALTOUTPUTL1 "$_\n"
					}
				}
			}
			else{
				$L1{"$pdb"}=$loopseq;
				printf OUTPUTL1 "$_\n";
			}
		}
		
		case "L2" {
			#print "$pdb - $loopseq - $loopCDR\n";
			if (grep( /^$pdb/,keys(%L2) )) {		
				if ($L2{"$pdb"} eq $loopseq) {
                    #print "$loopseq - $L2{\"$pdb\"} -- Same sequence : Do nothing \n";
                }
				else {
					if (grep( /^$loopseq/,@altL2 )) {
						#print "ALTTTTT ------------- $loopseq - Same sequence : Do nothing \n";
					}
					else{
						push (@altL2, $loopseq);
						#print "$loopseq --- NEW ALTERNATIVE \n";
						printf ALTOUTPUTL2 "$_\n"
					}
				}
			}
			else{
				$L2{"$pdb"}=$loopseq;
				printf OUTPUTL2 "$_\n";
			}	
		}
		
		case "L3" {
			#print "$pdb - $loopseq - $loopCDR\n";
			if (grep( /^$pdb/,keys(%L3) )) {		
				if ($L3{"$pdb"} eq $loopseq) {
                    #print "$loopseq - $L3{\"$pdb\"} -- Same sequence : Do nothing \n";
                }
				else {
					if (grep( /^$loopseq/,@altL3 )) {
						#print "ALTTTTT ------------- $loopseq - Same sequence : Do nothing \n";
					}
					else{
						push (@altL3, $loopseq);
						#print "$loopseq --- NEW ALTERNATIVE \n";
						printf ALTOUTPUTL3 "$_\n"
					}
				}
			}
			else{
				$L3{"$pdb"}=$loopseq;
				printf OUTPUTL3 "$_\n";
			}
		}
	}
}

close (DATA);

close (OUTPUTH1);
close (OUTPUTH2);
close (OUTPUTH3);
close (OUTPUTL1);
close (OUTPUTL2);
close (OUTPUTL3);
close (ALTOUTPUTH1);
close (ALTOUTPUTH2);
close (ALTOUTPUTH3);
close (OUTPUTL1);
close (OUTPUTL2);
close (OUTPUTL3);






###########################
#     FUNCTION
###########################

sub readPDBLine 
{
        my $line = shift;
        my $newAt = {};
        $newAt->{tag}=substr($line, 0,6);
        $newAt->{seq}=int(substr($line, 6,5));
        $newAt->{name}=substr($line, 12,4);
        $newAt->{altloc}=substr($line, 16,1);
        $newAt->{residueId}=substr($line,17,3);
        $newAt->{chain}=substr($line,21,1);
        $newAt->{resseq}=int(substr($line,22,4));
		  $newAt->{achar}=substr($line,26,1);
        $newAt->{x}=substr($line,30,8);
        $newAt->{y}=substr($line,38,8);
        $newAt->{z}=substr($line,46,8);
        $newAt->{occ}=substr($line, 54,6);
        $newAt->{Bfact}=substr($line, 60,6);
        $newAt->{charge}=$newAt->{occ};
        $newAt->{type}=substr($newAt->{name},1,1);
        return $newAt;
# Atomic Coordinate Entry Format Version 3.2
# http://www.wwpdb.org/documentation/format32/sect9.html
# COLUMNS        DATA  TYPE    FIELD        DEFINITION
#-------------------------------------------------------------------------------------
# 1 -  6        Record name   "ATOM  "
# 7 - 11        Integer       serial       Atom  serial number.
#13 - 16        Atom          name         Atom name.
#17             Character     altLoc       Alternate location indicator.
#18 - 20        Residue name  resName      Residue name.
#22             Character     chainID      Chain identifier.
#23 - 26        Integer       resSeq       Residue sequence number.
#27             AChar         iCode        Code for insertion of residues.
#31 - 38        Real(8.3)     x            Orthogonal coordinates for X in Angstroms.
#39 - 46        Real(8.3)     y            Orthogonal coordinates for Y in Angstroms.
#47 - 54        Real(8.3)     z            Orthogonal coordinates for Z in Angstroms.
#55 - 60        Real(6.2)     occupancy    Occupancy.
#61 - 66        Real(6.2)     tempFactor   Temperature  factor.
#77 - 78        LString(2)    element      Element symbol, right-justified.
#79 - 80        LString(2)    charge       Charge  on the atom.
}

# Writes a pdb-file form a "array of hashes" format
sub writePDB 
{
  # "%6c%5c%*c%4c%c%3c%*c%c%4d%c%*3c%8c%8c%8c%6c%6c%c%3c"
  my $data=shift;
  my $file=shift;
  my $renum=shift;
  open(PDB,">$file") or die "\nFailed to open $file\n";
  if($renum eq undef) # no renumber atom index
  {
	  foreach(@{$data})
	  {
	    printf PDB "%6s%5s %4s%s%3s %s%4s%s   %8s%8s%8s%6s%6s\n",
	      $_->{tag},
	      $_->{seq},
	      $_->{name},
	      $_->{altloc},
	      $_->{residueId},
	      $_->{chain},
	      $_->{resseq},
	      " ",
	      $_->{x},
	      $_->{y},
	      $_->{z},
	      $_->{occ},
	      $_->{Bfact};
	  }
  }
  else # renumber atom index
  {
	  my $index=1;
	  foreach(@{$data})
	  {
	    printf PDB "%6s%5s %4s%s%3s %s%4s%s   %8s%8s%8s%6s%6s\n",
	      $_->{tag},
	      $index,
	      $_->{name},
	      $_->{altloc},
	      $_->{residueId},
	      $_->{chain},
	      $_->{resseq},
	      " ",
	      $_->{x},
	      $_->{y},
	      $_->{z},
	      $_->{occ},
	      $_->{Bfact};
	    $index++;
	  }
  }
  close(PDB);
}


# Reads a pdb-file from a file
sub readPDB
{
  my $file=shift; # PDB file
  my @data=(); # The readed PDB will be stored here.
  
  my $debug=0;
  my $cont=0;
  open(PDB,"$file") or die "\nFailed to open $file\n";
  while(<PDB>)
  {
    next unless /^ATOM/; # reading only ATOM begining lines
    push(@data,readPDBLine($_));
#    printf "%5d bfact= %10.2f\n",$cont+1,$data[$cont]->{Bfact} if $debug;
    $cont++;
  }
  close PDB;
  return(\@data);
}

