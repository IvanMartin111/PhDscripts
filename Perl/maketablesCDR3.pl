#!/usr/bin/perl

#Script for generating the table reults of all the repicas of the CDR

#Iván Martín Hernández

# Date: Jun - 2021



use strict;
use warnings;
use diagnostics;
use Data::Dumper qw(Dumper);
use feature 'say';
use Scalar::Util qw(reftype);



#Definig variables, files and folderr

my %table;
my %rosnam;

my $fout="CDRsumary_thera_map_OCB.csv";

my $H1_rcd="20K_H1_thera_map_rep";
my $H2_rcd="20K_H2_thera_map_rep";
my $H3_rcd="50K_H3_thera_map_rep";
my $L1_rcd="20K_L1_thera_map_rep";
my $L2_rcd="20K_L2_thera_map_rep";
my $L3_rcd="20K_L3_thera_map_rep";


$rosnam{H1}="_native_20K_H1_thera_map_rep";  #rosettaX  _native_1K_H1_map_rep  X.txt
$rosnam{H2}="_native_20K_H2_thera_map_rep";  
$rosnam{H3}="_native_50K_H3_thera_map_rep";  
$rosnam{L1}="_native_20K_L1_thera_map_rep"; 
$rosnam{L2}="_native_20K_L2_thera_map_rep";  
$rosnam{L3}="_native_20K_L3_thera_map_rep"; 


my @rcds=($H1_rcd,$H2_rcd,$H3_rcd,$L1_rcd,$L2_rcd,$L3_rcd);
my @names=("H1","H2","H3","L1","L2","L3");
my @toproseta=(1,5,10);


# main
my $rcdcount=0;
foreach my $rcd (@rcds){    #Each CDR
    
    my $CDR=$names[$rcdcount];
    
    for (my $i = 1; $i < 4; $i++) {     #each replica con RDC
        my $results;
        ## reading result file
        if (($CDR eq "H3")) {
            $results="time_".$rcd."1/results.txt";
        }
        else{
        $results=$rcd.$i."/results.txt";
        }
        
        open (RES,$results) or die "The file $results doesn't exit in this path.\n";
        
        my $dummy=<RES>;
        $dummy=<RES>;
        $dummy=<RES>;
        my $count=0;
        while(<RES>)
        {
            next if /^#/; # reading only non-# beginning lines
            $_=~ s/^\t//;
            chop $_;
            my @line = split(/\s+/, $_);
            
            my $pdb = substr($line[0],0,4);
            $table{$pdb}{$CDR}{$i}{time} = $line[11];
            $table{$pdb}{$CDR}{seq}  = $line[18];
            
            $table{$pdb}{$CDR}{seq} =~ s/\A.//;
            chop $table{$pdb}{$CDR}{seq};
            
            $table{$pdb}{$CDR}{lenseq} = length($line[18])-2;
            
            $count++;
            
            if ($count == 56) {
                last;
            }
        }
        
       close (RES); 
        
        
        
        ## reading all rosettas
        
        
        for (my $j = 1; $j < 6; $j++) {     #each replica of rosetta
            foreach my $top (@toproseta){           # each best X selected of rosetta
                            
                my $roseta=$rcd.$i."/rosetta".$top.$rosnam{$CDR}.$j.".txt";
            
                open (ROS,$roseta) or die "The file $roseta doesn't exit in this path.\n";
            
                while(<ROS>)
                {
                    next if /^#/; # reading only non-# beginning lines
                    $_=~ s/^\t//;
                    chop $_;
                    my @line = split(/\s+/, $_);
                    
                    my $pdb = substr($line[1],0,4);
                    $table{$pdb}{$CDR}{$i}{$j}{$top}{E_full2} = $line[10];
                    $table{$pdb}{$CDR}{$i}{$j}{$top}{R_CA} = $line[11];
                    $table{$pdb}{$CDR}{$i}{$j}{$top}{R_BB} = $line[12];
                    $table{$pdb}{$CDR}{$i}{$j}{$top}{R_BBO} = $line[13];
                    $table{$pdb}{$CDR}{$i}{$j}{$top}{R_CB} = $line[14];
                    $table{$pdb}{$CDR}{$i}{$j}{$top}{R_OCB} = $line[15];
                }
                
               close (ROS); 
            }
            
            my $ave5=$rcd.$i."/rosettaAve5".$rosnam{$CDR}.$j.".txt";
                        
            open (AVE,$ave5) or die "The file $ave5 doesn't exit in this path.\n";
            while(<AVE>)
                {
                next if /^#/; # reading only non-# beginning lines
                $_=~ s/^\t//;
                chop $_;
                my @line = split(/\s+/, $_);
                    
                my $pdb = substr($line[1],0,4);
                $table{$pdb}{$CDR}{$i}{$j}{Ave5}{E_full2} = $line[10];
                $table{$pdb}{$CDR}{$i}{$j}{Ave5}{R_CA} = $line[11];
                $table{$pdb}{$CDR}{$i}{$j}{Ave5}{R_BB} = $line[12];
                $table{$pdb}{$CDR}{$i}{$j}{Ave5}{R_BBO} = $line[13];
                $table{$pdb}{$CDR}{$i}{$j}{Ave5}{R_CB} = $line[14];
                $table{$pdb}{$CDR}{$i}{$j}{Ave5}{R_OCB} = $line[15];
            }
            close(AVE);
        }
    }
    $rcdcount++;
}


    ## print CSV:
    



open(OUT, ">", $fout );

print OUT "PDB,H1_seq,H1_len,H1_time1,H1_time2,H1_time3,H1_Top1: 1-1,H1_Top1: 1-2,H1_Top1: 1-3,H1_Top1: 1-4,H1_Top1: 1-5,H1_Top1: 1-mean,H1_Top1: 2-1,H1_Top1: 2-2,H1_Top1: 2-3,H1_Top1: 2-4,H1_Top1: 2-5,H1_Top1: 2-mean,H1_Top1: 3-1,H1_Top1: 3-2,H1_Top1: 3-3,H1_Top1: 3-4,H1_Top1: 3-5,H1_Top1: 3-mean,H1_Top1: T-mean,H1_Top5: 1-1,H1_Top5: 1-2,H1_Top5: 1-3,H1_Top5: 1-4,H1_Top5: 1-5,H1_Top5: 1-mean,H1_Top5: 2-1,H1_Top5: 2-2,H1_Top5: 2-3,H1_Top5: 2-4,H1_Top5: 2-5,H1_Top5: 2-mean,H1_Top5: 3-1,H1_Top5: 3-2,H1_Top5: 3-3,H1_Top5: 3-4,H1_Top5: 3-5,H1_Top5: 3-mean,H1_Top5: T-mean,H1_Top10: 1-1,H1_Top10: 1-2,H1_Top10: 1-3,H1_Top10: 1-4,H1_Top10: 1-5,H1_Top10: 1-mean,H1_Top10: 2-1,H1_Top10: 2-2,H1_Top10: 2-3,H1_Top10: 2-4,H1_Top10: 2-5,H1_Top10: 2-mean,H1_Top10: 3-1,H1_Top10: 3-2,H1_Top10: 3-3,H1_Top10: 3-4,H1_Top10: 3-5,H1_Top10: 3-mean,H1_Top10: T-mean,H2_seq,H2_len,H2_time1,H2_time2,H2_time3,H2_Top1: 1-1,H2_Top1: 1-2,H2_Top1: 1-3,H2_Top1: 1-4,H2_Top1: 1-5,H2_Top1: 1-mean,H2_Top1: 2-1,H2_Top1: 2-2,H2_Top1: 2-3,H2_Top1: 2-4,H2_Top1: 2-5,H2_Top1: 2-mean,H2_Top1: 3-1,H2_Top1: 3-2,H2_Top1: 3-3,H2_Top1: 3-4,H2_Top1: 3-5,H2_Top1: 3-mean,H2_Top1: T-mean,H2_Top5: 1-1,H2_Top5: 1-2,H2_Top5: 1-3,H2_Top5: 1-4,H2_Top5: 1-5,H2_Top5: 2-1,H2_Top5: 1-mean,H2_Top5: 2-2,H2_Top5: 2-3,H2_Top5: 2-4,H2_Top5: 2-5,H2_Top5: 2-mean,H2_Top5: 3-1,H2_Top5: 3-2,H2_Top5: 3-3,H2_Top5: 3-4,H2_Top5: 3-5,H2_Top5: 3-mean,H2_Top5: T-mean,H2_Top10: 1-1,H2_Top10: 1-2,H2_Top10: 1-3,H2_Top10: 1-4,H2_Top10: 1-5,H2_Top10: 1-mean,H2_Top10: 2-1,H2_Top10: 2-2,H2_Top10: 2-3,H2_Top10: 2-4,H2_Top10: 2-5,H2_Top10: 2-mean,H2_Top10: 3-1,H2_Top10: 3-2,H2_Top10: 3-3,H2_Top10: 3-4,H2_Top10: 3-5,H2_Top10: 3-mean,H2_Top10: T-mean,H3_seq,H3_len,H3_time1,H3_time2,H3_time3,H3_Top1: 1-1,H3_Top1: 1-2,H3_Top1: 1-3,H3_Top1: 1-4,H3_Top1: 1-5,H3_Top1: 1-mean,H3_Top1: 2-1,H3_Top1: 2-2,H3_Top1: 2-3,H3_Top1: 2-4,H3_Top1: 2-5,H3_Top1: 2-mean,H3_Top1: 3-1,H3_Top1: 3-2,H3_Top1: 3-3,H3_Top1: 3-4,H3_Top1: 3-5,H3_Top1: 3-mean,H3_Top1: T-mean,H3_Top5: 1-1,H3_Top5: 1-2,H3_Top5: 1-3,H3_Top5: 1-4,H3_Top5: 1-5,H3_Top5: 1-mean,H3_Top5: 2-1,H3_Top5: 2-2,H3_Top5: 2-3,H3_Top5: 2-4,H3_Top5: 2-5,H3_Top5: 2-mean,H3_Top5: 3-1,H3_Top5: 3-2,H3_Top5: 3-3,H3_Top5: 3-4,H3_Top5: 3-5,H3_Top5: 3-mean,H3_Top5: T-mean,H3_Top10: 1-1,H3_Top10: 1-2,H3_Top10: 1-3,H3_Top10: 1-4,H3_Top10: 1-5,H3_Top10: 1-mean,H3_Top10: 2-1,H3_Top10: 2-2,H3_Top10: 2-3,H3_Top10: 2-4,H3_Top10: 2-5,H3_Top10: 2-mean,H3_Top10: 3-1,H3_Top10: 3-2,H3_Top10: 3-3,H3_Top10: 3-4,H3_Top10: 3-5,H3_Top10: 3-mean,H3_Top10: T-mean,L1_seq,L1_len,L1_time1,L1_time2,L1_time3,L1_Top1: 1-1,L1_Top1: 1-2,L1_Top1: 1-3,L1_Top1: 1-4,L1_Top1: 1-5,L1_Top1: 1-mean,L1_Top1: 2-1,L1_Top1: 2-2,L1_Top1: 2-3,L1_Top1: 2-4,L1_Top1: 2-5,L1_Top1: 2-mean,L1_Top1: 3-1,L1_Top1: 3-2,L1_Top1: 3-3,L1_Top1: 3-4,L1_Top1: 3-5,L1_Top1: 3-mean,L1_Top1: T-mean,L1_Top5: 1-1,L1_Top5: 1-2,L1_Top5: 1-3,L1_Top5: 1-4,L1_Top5: 1-5,L1_Top5: 1-mean,L1_Top5: 2-1,L1_Top5: 2-2,L1_Top5: 2-3,L1_Top5: 2-4,L1_Top5: 2-5,L1_Top5: 2-mean,L1_Top5: 3-1,L1_Top5: 3-2,L1_Top5: 3-3,L1_Top5: 3-4,L1_Top5: 3-5,L1_Top5: 3-mean,L1_Top5: T-mean,L1_Top10: 1-1,L1_Top10: 1-2,L1_Top10: 1-3,L1_Top10: 1-4,L1_Top10: 1-5,L1_Top10: 1-mean,L1_Top10: 2-1,L1_Top10: 2-2,L1_Top10: 2-3,L1_Top10: 2-4,L1_Top10: 2-5,L1_Top10: 2-mean,L1_Top10: 3-1,L1_Top10: 3-2,L1_Top10: 3-3,L1_Top10: 3-4,L1_Top10: 3-5,L1_Top10: 3-mean,L1_Top10: T-mean,L2_seq,L2_len,L2_time1,L2_time2,L2_time3,L2_Top1: 1-1,L2_Top1: 1-2,L2_Top1: 1-3,L2_Top1: 1-4,L2_Top1: 1-5,L2_Top1: 1-mean,L2_Top1: 2-1,L2_Top1: 2-2,L2_Top1: 2-3,L2_Top1: 2-4,L2_Top1: 2-5,L2_Top1: 2-mean,L2_Top1: 3-1,L2_Top1: 3-2,L2_Top1: 3-3,L2_Top1: 3-4,L2_Top1: 3-5,L2_Top1: 3-mean,L2_Top1: T-mean,L2_Top5: 1-1,L2_Top5: 1-2,L2_Top5: 1-3,L2_Top5: 1-4,L2_Top5: 1-5,L2_Top5: 1-mean,L2_Top5: 2-1,L2_Top5: 2-2,L2_Top5: 2-3,L2_Top5: 2-4,L2_Top5: 2-5,L2_Top5: 2-mean,L2_Top5: 3-1,L2_Top5: 3-2,L2_Top5: 3-3,L2_Top5: 3-4,L2_Top5: 3-5,L2_Top5: 3-mean,L2_Top5: T-mean,L2_Top10: 1-1,L2_Top10: 1-2,L2_Top10: 1-3,L2_Top10: 1-4,L2_Top10: 1-5,L2_Top10: 1-mean,L2_Top10: 2-1,L2_Top10: 2-2,L2_Top10: 2-3,L2_Top10: 2-4,L2_Top10: 2-5,L2_Top10: 2-mean,L2_Top10: 3-1,L2_Top10: 3-2,L2_Top10: 3-3,L2_Top10: 3-4,L2_Top10: 3-5,L2_Top10: 3-mean,L2_Top10: T-mean,L3_seq,L3_len,L3_time1,L3_time2,L3_time3,L3_Top1: 1-1,L3_Top1: 1-2,L3_Top1: 1-3,L3_Top1: 1-4,L3_Top1: 1-5,L3_Top1: 1-mean,L3_Top1: 2-1,L3_Top1: 2-2,L3_Top1: 2-3,L3_Top1: 2-4,L3_Top1: 2-5,L3_Top1: 2-mean,L3_Top1: 3-1,L3_Top1: 3-2,L3_Top1: 3-3,L3_Top1: 3-4,L3_Top1: 3-5,L3_Top1: 3-mean,L3_Top1: T-mean,L3_Top5: 1-1,L3_Top5: 1-2,L3_Top5: 1-3,L3_Top5: 1-4,L3_Top5: 1-5,L3_Top5: 1-mean,L3_Top5: 2-1,L3_Top5: 2-2,L3_Top5: 2-3,L3_Top5: 2-4,L3_Top5: 2-5,L3_Top5: 2-mean,L3_Top5: 3-1,L3_Top5: 3-2,L3_Top5: 3-3,L3_Top5: 3-4,L3_Top5: 3-5,L3_Top5: 3-mean,L3_Top5: T-mean,L3_Top10: 1-1,L3_Top10: 1-2,L3_Top10: 1-3,L3_Top10: 1-4,L3_Top10: 1-5,L3_Top10: 1-mean,L3_Top10: 2-1,L3_Top10: 2-2,L3_Top10: 2-3,L3_Top10: 2-4,L3_Top10: 2-5,L3_Top10: 2-mean,L3_Top10: 3-1,L3_Top10: 3-2,L3_Top10: 3-3,L3_Top10: 3-4,L3_Top10: 3-5,L3_Top10: 3-mean,L3_Top10: T-mean,";                                                                                       
foreach my $CDR (@names){ 

    for (my $i = 1; $i < 4; $i++) {
    
        for (my $j = 1; $j < 6; $j++) {  
            print OUT "$CDR-Ave5: $i-$j,";
    
        }
         print OUT "$CDR-Ave5: $i-mean,";
    }
    print OUT "$CDR-Ave5: T-mean,";
}
 print OUT "\n";




foreach my $pdb (sort keys %table){
    print OUT "$pdb,";
    
    foreach my $CDR (@names){ 
        print OUT "$table{$pdb}{$CDR}{seq},$table{$pdb}{$CDR}{lenseq},$table{$pdb}{$CDR}{1}{time},$table{$pdb}{$CDR}{2}{time},$table{$pdb}{$CDR}{3}{time},";
            
        foreach my $top (@toproseta){
            my $t=0;
            for (my $i = 1; $i < 4; $i++) {
                my $p=0;
                for (my $j = 1; $j < 6; $j++) {  
                    print OUT "$table{$pdb}{$CDR}{$i}{$j}{$top}{R_OCB},";
                    $p=$p+$table{$pdb}{$CDR}{$i}{$j}{$top}{R_OCB};
                    $t=$t+$table{$pdb}{$CDR}{$i}{$j}{$top}{R_OCB};
                }
                my $meanp=$p/5;
                print OUT "$meanp,";
            }
            my $meant=$t/15;
            print OUT "$meant,";
        }
    }
    #printing Ave5
     foreach my $CDR (@names){    
        my $t2=0;
        for (my $i = 1; $i < 4; $i++) {
            my $p2=0;
            for (my $j = 1; $j < 6; $j++) {  
                print OUT "$table{$pdb}{$CDR}{$i}{$j}{Ave5}{R_OCB},";
                $p2=$p2+$table{$pdb}{$CDR}{$i}{$j}{Ave5}{R_OCB};
                $t2=$t2+$table{$pdb}{$CDR}{$i}{$j}{Ave5}{R_OCB};
            }
            my $meanp2=$p2/5;
            print OUT "$meanp2,";
        }
        my $meant2=$t2/15;
        print OUT "$meant2,";
 
    }
    print OUT "\n";
}
close(OUT);







