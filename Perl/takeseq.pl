#!/usr/bin/perl
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

my %short = (
   "ALA" => "A",	"CYS" => "C",   "CYX" => "C",   "ASP" => "D",   "GLU" => "E",   "PHE" => "F",   "GLY" => "G",   "HIS" => "H",   "HSD" => "H",
   "HSE" => "H",  "HSP" => "H",   "HIP" => "H",   "HID" => "H",   "HIE" => "H",   "ILE" => "I",   "LYS" => "K",   "LEU" => "L",   "MET" => "M",
   "ASN" => "N",  "PRO" => "P",   "GLN" => "Q",   "ARG" => "R",   "SER" => "S",   "THR" => "T",   "VAL" => "V",   "TRP" => "W",   "TYR" => "Y",
   "  A" => "A",  "  G" => "G",  "  T" => "T",   "  C" => "C",   "  U" => "U",   " DA" => "A",   " DG" => "G",   " DT" => "T",   " DC" => "C" ) ;

my @aa = ("A","R","N","D","C","Q","E","G","H","I","L","K","M","F","P","S","T","W","Y","V");


my %pdbs;
my $file = $ARGV[0];

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
	my $pdb_name=$line[0];
  
  
  say $pdb_name;
	open (PDB, $pdb_name) or die "Can't open file $pdb_name.\n";
  my $resN=-1;
  
  foreach my $line (<PDB>) {
		chomp($line);
		#print "$line\n";	
		
		if ($line =~ /^ATOM/ || $line =~ /^HETATM/) {
			my $format = readPDBLine($line);
      my $res_num=substr($line,22,5); # get residue number
      $res_num=~ s/\s//g;
			my $chain =$format->{chain};
			my $res =$format->{residueId};
      
      #print "-$resN-$res_num-\n";
      unless ($resN eq $res_num) {
        $pdbs{$pdb_name}{$chain}=$pdbs{$pdb_name}{$chain}.$short{$res};
        $resN=$res_num;
      }
      
      
    }
  }
}
close(DATA);



foreach my $pdb(sort keys %pdbs){
  
  print "$pdb \n";
  
  foreach my $chain (sort keys %{$pdbs{$pdb}} ){
    
    print "  $chain\t  $pdbs{$pdb}{$chain}\n";
  }
}














###########################
#     FUNCTION
###########################

sub readPDBLine 
{
        my $line = shift;
        my $newAt = {};
        $newAt->{tag}=substr($line, 0,6);
        $newAt->{seq}=int(substr($line, 7,4));
        $newAt->{name}=substr($line, 12,4);
        $newAt->{altloc}=substr($line, 16,1);
        $newAt->{residueId}=substr($line,17,3);
        $newAt->{chain}=substr($line,21,1);
        $newAt->{resseq}=int(substr($line,22,4));
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

