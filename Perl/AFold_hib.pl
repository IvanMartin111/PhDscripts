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
if( !($#ARGV == 1))
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
my $rcdfile = $ARGV[0];
my $CD =$ARGV[1];


my %AA = (A=>'ALA',Y=>'TYR',M=>'MET',L=>'LEU',C=>'CYS',G=>'GLY',
         R=>'ARG',N=>'ASN',D=>'ASP',Q=>'GLN',E=>'GLU',H=>'HIS',W=>'TRP',
         K=>'LYS',F=>'PHE',P=>'PRO',S=>'SER',T=>'THR',I=>'ILE',V=>'VAL');


my %aa = reverse %AA;
my %rmsds;
my %ancr;
my @RR=("CA","BB","BBO","CB","OCB");
my $fout="AFF_thera_results_".$CD.".csv";



# Read input file
my $data = parseTextFile($rcdfile);

my $nloop = 0;
my $base; # required outside
foreach my $loop (@{$data}) # Screen Loops
{
	#Variable names	
	my $pdbCry = $loop->[0];
	my $start = $loop->[1];
	my $end = $loop->[2];
	my $chain = $loop->[3];
	my $AFchain= "B";
	if ($chain eq "H") { $AFchain= "A";}
	my $CDR = $loop->[6];
	my $startA = $start - 4;
	my $endA = $end + 4;
	$base = basename($pdbCry,"_4.pdb");
	
	my $AFpdb=$base."_AFF.pdb";
	my $AFmvpdb=$base."_".$CDR."_AFF_mv.pdb";
	my $Loop_nat=$base."_".$CDR."_loop_native.pdb";
	my $Loop_nat_mv=$base."_".$CDR."_loop_native_mv.pdb";
	my $Loop_AF=$base."_".$CDR."_loop_AFF.pdb";
	my $achAF=$base."_".$CDR."_acnAF.pdb";
	my $achnat=$base."_".$CDR."_acncr.pdb";

	my $hib=$base."_".$CDR."_hibAFF.pdb";
	
	my $hibAF_lpcr=$base."_hibAFF_crloops.pdb";
	my $hibcr_lpAF=$base."_hibcr_AFFloops.pdb";
	
	##########
	
	print "$base  $CDR- $pdbCry  - $AFpdb --$Loop_nat --$Loop_AF \n";

	
	#Change chain letter in the AF pdb
	
	my $sed="sed 's/ A / H /g' $AFpdb >kk.pdb; mv kk.pdb $AFpdb";	
	system($sed);
	$sed="sed 's/ B / L /g' $AFpdb >kk.pdb; mv kk.pdb $AFpdb";	
	system($sed);

	
	
	#anchors AF	
	open(PDB, $AFpdb) or die "Can't open file $AFpdb.\n";
	open(OUT,'>',$achAF);

	foreach my $line (<PDB>) {
		chomp($line);
		#print "$line\n";	
		
		if ($line =~ /^ATOM/ || $line =~ /^HETATM/) {
			my $format = readPDBLine($line);
			my $res_num=$format->{resseq};
			my $achain =$format->{chain};
			#print "-$chain-$achain- $res_num\n";
			if ($achain eq $chain && (($res_num >= $startA && $res_num < $start-1)||( $res_num > $end +1 &&$res_num <= $endA))) {print OUT "$line \n";}
		}
	}
	close(OUT);
	close(PDB);
	
	
	#anchors cristal
	open(PDB, $pdbCry) or die "Can't open file $pdbCry.\n";
	open(OUT,'>',$achnat);

	foreach my $line (<PDB>) {
		chomp($line);
		#print "$line\n";	
		
		if ($line =~ /^ATOM/ || $line =~ /^HETATM/) {
			my $format = readPDBLine($line);
			my $res_num=$format->{resseq};
			my $achain =$format->{chain};
#			print "$res -$chain-$achain- $res_num\n";
			if ($achain eq $chain && (($res_num >= $startA && $res_num < $start-1)||( $res_num > $end +1 &&$res_num <= $endA))) {print OUT "$line \n";}
		}
	}
	close(OUT);
	close(PDB);
	
	
	
	#Loop cristal
	open(PDB, $pdbCry) or die "Can't open file $pdbCry.\n";
	open(OUT,'>',$Loop_nat);

	foreach my $line (<PDB>) {
		chomp($line);
		#print "$line\n";	
		
		if ($line =~ /^ATOM/ || $line =~ /^HETATM/) {
			my $format = readPDBLine($line);
			my $res_num=$format->{resseq};
			my $achain =$format->{chain};
#			print "$res -$chain-$achain- $res_num\n";
			if ($achain eq $chain && $res_num >= $start && $res_num <= $end) {print OUT "$line \n";}
		}
	}
	close(OUT);
	close(PDB);
	
	
	#Align anchors and mov 
	
	my $rmsdanch="rmsdanch.txt";
	my $cmd= "rmsd $achnat $achAF  -b --CB -o temp.pdb --input2 $AFpdb >$rmsdanch; rm temp.pdb; mv temp2.pdb $AFmvpdb";
	print "$cmd\n";
	system ($cmd);
	
	open(RMSD, $rmsdanch) or die "Can't open file $rmsdanch.\n";
	my $line=<RMSD>;
	$line=<RMSD>;
	$line=<RMSD>;
	$line=<RMSD>;
	$line=<RMSD>;
	$line=<RMSD>;
	close(RMSD);
	chop $line;
	my @ele = split(/\s+/, $line);
	$ancr{$base}=$ele[2];
	
	
	$cmd= "rmsd $achAF $achnat  -b --CB -o temp.pdb --input2 $Loop_nat >$rmsdanch; rm temp.pdb; mv temp2.pdb $Loop_nat_mv";
	system ($cmd);
	
	
	#AF loop
	open(PDB, $AFmvpdb) or die "Can't open file $AFmvpdb.\n";
	open(OUT,'>',$Loop_AF);

	foreach my $line (<PDB>) {
		chomp($line);
		#print "$line\n";	
		
		if ($line =~ /^ATOM/ || $line =~ /^HETATM/) {
			my $format = readPDBLine($line);
			my $res_num=$format->{resseq};
			my $achain =$format->{chain};
#			print "$res -$chain-$achain- $res_num\n";
			if ($achain eq $chain && $res_num >= $start && $res_num <= $end){print OUT "$line \n";}
		}
	}
	close(OUT);
	close(PDB);
	
	
	

	
	
	#calculate RMSDs between AF loop and cristal moved loop
	my $rmsdfile= "rmsdbb.txt";
	$cmd= "rmsd $Loop_nat $Loop_AF -b > $rmsdfile";
	#print "\nBB: $cmd\n";
	system ($cmd);
	
	open(RMSD, $rmsdfile) or die "Can't open file $rmsdfile.\n";
	$line=<RMSD>;
	close(RMSD);
	chop $line;
	@ele = split(/\s+/, $line);
	$rmsds{$base}{BB}=$ele[2];
	
	
	$rmsdfile="rmsdca.txt";
	$cmd= "rmsd $Loop_nat $Loop_AF -c >$rmsdfile";
	#print "CA: $cmd\n";
	system ($cmd);
	
	open(RMSD, $rmsdfile) or die "Can't open file $rmsdfile.\n";
	$line=<RMSD>;
	close(RMSD);
	chop $line;
	@ele = split(/\s+/, $line);
	$rmsds{$base}{CA}=$ele[2];
	
	
	$rmsdfile="rmsdbbo.txt";
	$cmd= "rmsd $Loop_nat $Loop_AF -b --O >$rmsdfile";
	#print "CA: $cmd\n";
	system ($cmd);
	
	open(RMSD, $rmsdfile) or die "Can't open file $rmsdfile.\n";
	$line=<RMSD>;
	close(RMSD);
	chop $line;
	@ele = split(/\s+/, $line);
	$rmsds{$base}{BBO}=$ele[2];
	
	
	$rmsdfile="rmsdCB.txt";
	$cmd= "rmsd $Loop_nat $Loop_AF -b --CB >$rmsdfile";
	#print "CA: $cmd\n";
	system ($cmd);
	
	open(RMSD, $rmsdfile) or die "Can't open file $rmsdfile.\n";
	$line=<RMSD>;
	close(RMSD);
	chop $line;
	@ele = split(/\s+/, $line);
	$rmsds{$base}{CB}=$ele[2];
	
	
	
	$rmsdfile="rmsdoCB.txt";
	$cmd= "rmsd $Loop_nat $Loop_AF -b --CB --O >$rmsdfile";
	#print "CA: $cmd\n";
	system ($cmd);
	
	open(RMSD, $rmsdfile) or die "Can't open file $rmsdfile.\n";
	$line=<RMSD>;
	close(RMSD);
	chop $line;
	@ele = split(/\s+/, $line);
	$rmsds{$base}{OCB}=$ele[2];
	
	
	
	
	
	####Generate the mixed pdb
	
	my $insert="insertpdbrnb.pl $Loop_nat $AFmvpdb $hib";

	print "\t$insert\n";
	system($insert);
	
	
	### Generating the multihibridos
	if (-e $hibAF_lpcr) {
		print "the file $hibAF_lpcr exists\n";
	} else {
		print "the file $hibAF_lpcr does not exist!\n";
		
		my $cpcmd="cp $AFpdb $hibAF_lpcr";
		system($cpcmd);
	}	
	
	$insert="insertpdbrnb.pl $Loop_nat_mv $hibAF_lpcr $hibAF_lpcr";
	system($insert);
	
	
	
	if (-e $hibcr_lpAF) {
		print "the file $hibcr_lpAF exists\n";
	} else {
		print "the file $hibcr_lpAF does not exist!\n";
		
		my $cpcmd="cp $pdbCry $hibcr_lpAF";
		system($cpcmd);
	}	
	
	$insert="insertpdbrnb.pl $Loop_AF $hibcr_lpAF $hibcr_lpAF";
	system($insert);
	
	
}


system ("rm rmsd*");



#printing results
open(OUT, ">", $fout );
print OUT "PDB,$CD R_ach_CB,$CD R_CA,$CD R_BB,$CD R_BBO,$CD R_CB,$CD R_OCB\n";
for my $pdb (sort keys %rmsds){
	print OUT "$pdb,";
	print OUT"$ancr{$pdb},";

	foreach my $rname (@RR){
		print OUT "$rmsds{$pdb}{$rname},";
	}
	print OUT "\n";
}


close (OUT);







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






