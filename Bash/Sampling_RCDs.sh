#!/bin/bash



ncpu=20
CDR=H3

# RCD 20K
for loco_best in 250 500 1000 2000 4000 8000 16000; do
	id=$CDR'_20K_rcd_'$loco_best'_KORP'
	
	
	mpirun -np 10 rcd_mpi rcd_H3P.txt --ABmap CDRS.bin  --H3   --kink_merge 0.9 --kink_mergeA 0.9 --nokink_merge 0.5 --tkink 0.70 --tkinkA 0.70 --nterm_rama 0  -r -e 56 --energy_file ../korp3A90c18mr1f10_Allno50_6Dc14mr3r10n36nc8b4nb1yz8bf18_map.bin --bench -n 20000 --loco_best $loco_best -x dunbrack.bin -d 0.5 --linear -t 0.99 -o $id 
	
done


# RCD 40K
for loco_best in 250 500 1000 2000 4000 8000 16000; do
	id=$CDR'_40K_rcd_'$loco_best'_KORP'
	
	
	mpirun -np 10 rcd_mpi rcd_H3P.txt --ABmap CDRS.bin  --H3   --kink_merge 0.9 --kink_mergeA 0.9 --nokink_merge 0.5 --tkink 0.70 --tkinkA 0.70 --nterm_rama 0  -r -e 56 --energy_file ../korp3A90c18mr1f10_Allno50_6Dc14mr3r10n36nc8b4nb1yz8bf18_map.bin --bench -n 40000 --loco_best $loco_best -x dunbrack.bin -d 0.5 --linear -t 0.99 -o $id 
	
done



# RCD 80K
for loco_best in 250 500 1000 2000 4000 8000 16000; do
	id=$CDR'_80K_rcd_'$loco_best'_KORP'
	
	
	mpirun -np 10 rcd_mpi rcd_H3P.txt --ABmap CDRS.bin  --H3   --kink_merge 0.9 --kink_mergeA 0.9 --nokink_merge 0.5 --tkink 0.70 --tkinkA 0.70 --nterm_rama 0  -r -e 56 --energy_file ../korp3A90c18mr1f10_Allno50_6Dc14mr3r10n36nc8b4nb1yz8bf18_map.bin --bench -n 80000 --loco_best $loco_best -x dunbrack.bin -d 0.5 --linear -t 0.99 -o $id 
	
done


# RCD 160K
for loco_best in 250 500 1000 2000 4000 8000 16000; do
	id=$CDR'_160K_rcd_'$loco_best'_KORP'
	
	
	mpirun -np 10 rcd_mpi rcd_H3P.txt --ABmap CDRS.bin  --H3   --kink_merge 0.9 --kink_mergeA 0.9 --nokink_merge 0.5 --tkink 0.70 --tkinkA 0.70 --nterm_rama 0  -r -e 56 --energy_file ../korp3A90c18mr1f10_Allno50_6Dc14mr3r10n36nc8b4nb1yz8bf18_map.bin --bench -n 160000 --loco_best $loco_best -x dunbrack.bin -d 0.5 --linear -t 0.99 -o $id 
	
done
