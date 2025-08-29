#!/usr/bin/perl

#Script generating the hibrid strctures with the AlfaFold model and the cristal loop. Also obtains de RMSD of the AF loop.

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


# INPUT PARSER
if( !($#ARGV == 0))
{
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
my $CD =$ARGV[0];





my %AA = (A=>'ALA',Y=>'TYR',M=>'MET',L=>'LEU',C=>'CYS',G=>'GLY',
         R=>'ARG',N=>'ASN',D=>'ASP',Q=>'GLN',E=>'GLU',H=>'HIS',W=>'TRP',
         K=>'LYS',F=>'PHE',P=>'PRO',S=>'SER',T=>'THR',I=>'ILE',V=>'VAL');


my %aa = reverse %AA;
my %rmsds;
my %ancr;
my @mix=("mix-2","mix-3","mix-3m","mix-4");
my @rep1=(1,2,3);
my @repr=(1,2,3,4,5);
my @rosn=(1,5,10);
my $fout="Acnhor_RMSD_compare_".$CD.".csv";




### Main ###

foreach my $rep1(@rep1){
	foreach my $repr(@repr){
		foreach my $rosn(@rosn){
			my $rosfile="/home/ivan/disk/PhD/loops_H3/RosetaBenchark/hibAF_cr/50K_".$CD."_AF_map_rep".$rep1."/rosetta".$rosn."_modeling_50K_".$CD."_AF_map_rep".$repr.".txt";
			#print "$rosfile\n";
			my $outpdb="tempmodel".$rosn.".pdb";
			my $fistout="kktemp.pdb";
			
			
			
			
			
			open(FILE, $rosfile) or die "Can't open file $rosfile.\n";
			
			while(<FILE>)
			{
				next if /^#/; # reading only non-# beginning lines
				$_=~ s/^\t//;
				chop $_;
				my @line = split(/\s+/, $_);
            
				my $pdb = substr($line[1],0,4);
				my $model = $line[2];
				my $mname = sprintf "MODEL    %3s",$model;
				print "$mname\n";




				# Bloque guardar modelo generado por RCD para el RMSD
				my $loopfile="/home/ivan/disk/PhD/loops_H3/RosetaBenchark/hibAF_cr/50K_".$CD."_AF_map_rep".$rep1."/".$pdb."_hibAF_crloops_modeling_50K_".$CD."_AF_map_rep".$repr."_loops.pdb";
				

				open(ROS, $loopfile) or die "Can't open file $loopfile.\n";
				open(OUTT, ">", $fistout );
				
				my $control=0;
				while(<ROS>)
					{
					if (/ENDMDL/){$control=0;}
                    if ($control eq 1) {
                        print OUTT "$_";
						 
                    }
					if (/^$mname/){$control=1;}
					
				}
				close(OUTT);
				close(ROS);
				
				
				
				
				#Bloque eliminar 1 y ultimo aa desde temp.pdb
				
				my $lin = `head -1  $fistout`;
				my @lin= split(/\s+/, $lin);
				my $first=$lin[5];
				
				$lin = `tail -1  $fistout`;
				@lin= split(/\s+/, $lin);
				my $last=$lin[5];
				
				open(PDB, $fistout) or die "Can't open file $fistout.\n";
				open(OUTT, ">", $outpdb );
				
				while(<PDB>){
					my @line = split(/\s+/, $_);
					my $aa=$line[5];
					if ($aa ne $first && $aa ne $last) {
						print OUTT "$_";
                    }
                    
					
				}
				close(OUTT);				
				close (PDB);
				
				
				
				#Bloque RMSD		
				foreach my $mix(@mix){
					my $Loop_nat_mv=$mix."/".$pdb."_".$CD."_loop_native_mv.pdb";
					my $rmsdfile="rmsdCB.txt";
					my $cmd= "rmsd $outpdb $Loop_nat_mv -b --CB >$rmsdfile";
					print "$cmd\n";
					system ($cmd);
					
					open(RMSD, $rmsdfile) or die "Can't open file $rmsdfile.\n";
					my $line=<RMSD>;
					close(RMSD);
					chop $line;
					my @ele = split(/\s+/, $line);
					print "$pdb-$mix-$rep1-$repr-$rosn : $ele[2]\n";	
					$rmsds{$pdb}{$mix}{$rosn}{$rep1}{$repr}=$ele[2];
				}	
			}
			close(FILE);
		}
	}
}

#Bloque calculo totales
foreach my $pdb (sort keys %rmsds){
	foreach my $mix(@mix){
		foreach my $rosn(@rosn){
			my $sumatorio=0;
			foreach my $rep1(@rep1){
				foreach my $repr(@repr){
					$sumatorio = $sumatorio + $rmsds{$pdb}{$mix}{$rosn}{$rep1}{$repr};
				}
			}
			$rmsds{$pdb}{$mix}{$rosn}{Total}=$sumatorio/15;
		}
	}
}




    ## print CSV:
open(OUT, ">", $fout );

print OUT "PDB,";

foreach my $mix(@mix){
	foreach my $rosn(@rosn){
		print OUT "$mix-$rosn,";	
	}
}
print OUT "\n";

foreach my $pdb (sort keys %rmsds){
	print OUT "$pdb,";
	foreach my $mix(@mix){
		foreach my $rosn(@rosn){
			print OUT "$rmsds{$pdb}{$mix}{$rosn}{Total},";	
		}
	}
	print OUT "\n";
	
}


close(OUT);


#
#foreach my $loop (@{$data}) # Screen Loops
#{
#	#Variable names	
#	my $pdbCry = $loop->[0];
#	my $start = $loop->[1];
#	my $end = $loop->[2];
#	my $chain = $loop->[3];
#	my $AFchain= "B";
#	if ($chain eq "H") { $AFchain= "A";}
#	my $CDR = $loop->[6];
#	$base = basename($pdbCry,"_4.pdb");
#	
#	my $AFpdb=$base."_AFF.pdb";
#	my $AFmvpdb=$base."_".$CDR."_AFF_mv.pdb";
#	my $Loop_nat=$base."_".$CDR."_loop_native.pdb";
#	my $Loop_nat_mv=$base."_".$CDR."_loop_native_mv.pdb";
#	my $Loop_AF=$base."_".$CDR."_loop_AFF.pdb";
#	my $achAF=$base."_".$CDR."_acnAF.pdb";
#	my $achnat=$base."_".$CDR."_acncr.pdb";
#
#
#	
#	foreach my $rep1(@rep1){
#		
#		
#	}
#	
#}
#
#
#system ("rm rmsd*");
#
#
#
##printing results
#open(OUT, ">", $fout );
#print OUT "PDB,$CD R_ach_CB,$CD R_CA,$CD R_BB,$CD R_BBO,$CD R_CB,$CD R_OCB\n";
#for my $pdb (sort keys %rmsds){
#	print OUT "$pdb,";
#	print OUT"$ancr{$pdb},";
#
#	foreach my $rname (@RR){
#		print OUT "$rmsds{$pdb}{$rname},";
#	}
#	print OUT "\n";
#}
#
#
#close (OUT);







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


sub insertPDB
{
	my $loop = shift; # loop structure
	my $pdb = shift; # complete structure
	my $chain = shift; # chain id of the complete structure
	my $start = shift; # first residue index of the loop
	my $end = shift; # last residue index of the loop

	# Screen complete struture to get insertion/deletion indices...
	my $i=0;
	my ($first,$last)=(0,0);
	foreach my $full (@{$pdb}) # screen atoms
	{
		$i++;
		$first = $i if($full->{resseq} == $start-1 && $full->{chain} eq $chain);
		# $last = $i if($full->{resseq} == $end && $full->{chain} eq $chain);
	}

	for($i = $#{$pdb} ; $i >= 0 ; $i-- ) 
	{ 
		my $full = ${$pdb}[$i];
		$last = $i if($full->{resseq} == $end+1 && $full->{chain} eq $chain);
	}

	#print "$prog> start= $start  end= $end  -->  first= $first  last= $last\n";
	splice @{$pdb},$first,$last-$first,@{$loop}; # delete and insert one array into other array.
}



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






