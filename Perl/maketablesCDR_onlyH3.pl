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

my $fout="H3_nat_RABh_sumary_map_CB.csv";

my $H3_rcd="50K_map_ivf_RABh";

$rosnam{H3}="_native_50K_map_ivf_RABh_rep";  



my @rcds=($H3_rcd);
my @names=("H3");
my @toproseta=(1,5,10);
my $RCDrep=3;
my $rosrep=3;


# main
my $rcdcount=0;
foreach my $rcd (@rcds){    #Each CDR

    my $CDR=$names[$rcdcount];
    
    for (my $i = 1; $i < 4; $i++) {     #each replica con RDC
        
        ## reading result file
        
        my $results=$rcd.$i."/results.txt";
        
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
            
            if ($count == 45) {
                last;
            }
        }
        
       close (RES); 
        
        
        
        ## reading all rosettas
        
        
        for (my $j = 1; $j < 4; $j++) {     #each replica of rosetta
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

print OUT "PDB,H3_seq,H3_len,H3_time1,H3_time2,H3_time3,H3_Top1: 1-1,H3_Top1: 1-2,H3_Top1: 1-3,H3_Top1: 1-mean,H3_Top1: 2-1,H3_Top1: 2-2,H3_Top1: 2-3,H3_Top1: 2-mean,H3_Top1: 3-1,H3_Top1: 3-2,H3_Top1: 3-3,H3_Top1: 3-mean,H3_Top1: T-mean,H3_Top5: 1-1,H3_Top5: 1-2,H3_Top5: 1-3,H3_Top5: 1-mean,H3_Top5: 2-1,H3_Top5: 2-2,H3_Top5: 2-3,H3_Top5: 2-mean,H3_Top5: 3-1,H3_Top5: 3-2,H3_Top5: 3-3,H3_Top5: 3-mean,H3_Top5: T-mean,H3_Top10: 1-1,H3_Top10: 1-2,H3_Top10: 1-3,H3_Top10: 1-mean,H3_Top10: 2-1,H3_Top10: 2-2,H3_Top10: 2-3,H3_Top10: 2-mean,H3_Top10: 3-1,H3_Top10: 3-2,H3_Top10: 3-3,H3_Top10: 3-mean,H3_Top10: T-mean,";                                                                                       
foreach my $CDR (@names){ 

    for (my $i = 1; $i < 4; $i++) {
    
        for (my $j = 1; $j < 4; $j++) {  
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
                for (my $j = 1; $j < 4; $j++) {  
                    print OUT "$table{$pdb}{$CDR}{$i}{$j}{$top}{R_CB},";
                    $p=$p+$table{$pdb}{$CDR}{$i}{$j}{$top}{R_CB};
                    $t=$t+$table{$pdb}{$CDR}{$i}{$j}{$top}{R_CB};
                }
                my $meanp=$p/3;
                print OUT "$meanp,";
            }
            my $meant=$t/9;
            print OUT "$meant,";
        }
    }
    #printing Ave5
     foreach my $CDR (@names){    
        my $t2=0;
        for (my $i = 1; $i < 4; $i++) {
            my $p2=0;
            for (my $j = 1; $j < 4; $j++) {  
                print OUT "$table{$pdb}{$CDR}{$i}{$j}{Ave5}{R_CB},";
                $p2=$p2+$table{$pdb}{$CDR}{$i}{$j}{Ave5}{R_CB};
                $t2=$t2+$table{$pdb}{$CDR}{$i}{$j}{Ave5}{R_CB};
            }
            my $meanp2=$p2/3;
            print OUT "$meanp2,";
        }
        my $meant2=$t2/9;
        print OUT "$meant2,";
 
    }
    print OUT "\n";
}
close(OUT);







