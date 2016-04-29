 #! /bin/bash
 # original by Haoquan
 # contains some modifications to work with Amber15 MMPBSA.py
 # note that you should change the surface tension and offset 
 # accordingly for the esurf calculation (i.e. 0.0072)
 # Gil Ferreira Hoben
 echo "Processing MMGBSA data from:" $PWD
 echo _MMPBSA_complex_gb _MMPBSA_receptor_gb _MMPBSA_ligand_gb > namelist
 LIST=`cat namelist`
 mkdir Analysis
 for i in $LIST ; do
 grep VDWAALS $i.mdout.0 | awk '{print $3}' > ./Analysis/$i.vdw
 grep EGB     $i.mdout.0 | awk '{print $9}' > ./Analysis/$i.polar
 grep EGB     $i.mdout.0 | awk '{print $6}' > ./Analysis/$i.coul
 cat       ${i}_surf.dat.0 | awk 'NR>1{print $2 * 0.0072}' > ./Analysis/$i.surf
 paste -d " " ./Analysis/$i.vdw ./Analysis/$i.polar ./Analysis/$i.coul ./Analysis/$i.surf | awk '{print $1 + $2 + $3 + $4}' > Analysis/data.$i
 #rm $i.*
 done
 paste -d " " Analysis/data._MMPBSA_complex_gb Analysis/data._MMPBSA_receptor_gb Analysis/data._MMPBSA_ligand_gb | awk '{print $1 - $2 - $3}' > Analysis/data.all
 for ((j=1; j<=`wc -l Analysis/data.all | awk '{print $1}'`; j+=1)) do
 echo $j , >> Analysis/time
 done
 paste -d " " ./Analysis/time Analysis/data.all > Analysis/MMGBSA_vs_time.dat
 #rm namelist ./Analysis/time ./Analysis/data.*

