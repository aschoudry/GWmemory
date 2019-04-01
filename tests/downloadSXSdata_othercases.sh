#!/bin/bash 
#No spin different mass ratio
mkdir Spinning_Sz_antialigned
cd Spinning_Sz_antialigned

wget -c https://zenodo.org/record/1215095/files/SXS:BBH:0218/Lev5/rMPsi4_Asymptotic_GeometricUnits_CoM.h5?download=1
mv rMPsi4_Asymptotic_GeometricUnits_CoM.h5?download=1 rMPsi4_Sz1_m0p5_Sz2_0p5_q1p5data.h5

wget -c https://zenodo.org/record/1214497/files/SXS:BBH:0211/Lev5/rMPsi4_Asymptotic_GeometricUnits_CoM.h5?download=1
mv rMPsi4_Asymptotic_GeometricUnits_CoM.h5?download=1 rMPsi4_Sz1_m0p9_Sz2_0p9_q1p5data.h5

cd ..
mkdir Spinning_only_one_Sz
cd Spinning_only_one_Sz

wget -c https://zenodo.org/record/1215059/files/SXS:BBH:0210/Lev5/rMPsi4_Asymptotic_GeometricUnits_CoM.h5?download=1
mv rMPsi4_Asymptotic_GeometricUnits_CoM.h5?download=1 rMPsi4_Sz1_m0p9_q1p5data.h5

wget -c https://zenodo.org/record/1215085/files/SXS:BBH:0216/Lev5/rMPsi4_Asymptotic_GeometricUnits_CoM.h5?download=1
mv rMPsi4_Asymptotic_GeometricUnits_CoM.h5?download=1 rMPsi4_Sz1_m0p6_q1p5data.h5

wget -c https://zenodo.org/record/1214330/files/SXS:BBH:0222/Lev5/rMPsi4_Asymptotic_GeometricUnits_CoM.h5?download=1
mv rMPsi4_Asymptotic_GeometricUnits_CoM.h5?download=1 rMPsi4_Sz1_m0p3_q1p5data.h5

cd ..
mkdir Spinning_algned_with_different_mag
cd Spinning_algned_with_different_mag

wget -c https://zenodo.org/record/1214317/files/SXS:BBH:0209/Lev5/rMPsi4_Asymptotic_GeometricUnits_CoM.h5?download=1
mv rMPsi4_Asymptotic_GeometricUnits_CoM.h5?download=1 rMPsi4_Sz1_m0p9_Sz2_m0p5_q1p5data.h5

wget -c https://zenodo.org/record/1215081/files/SXS:BBH:0214/Lev5/rMPsi4_Asymptotic_GeometricUnits_CoM.h5?download=1
mv rMPsi4_Asymptotic_GeometricUnits_CoM.h5?download=1 rMPsi4_Sz1_m0p625_Sz2_m0p25_q1p5data.h5

wget -c https://zenodo.org/record/1214325/files/SXS:BBH:0221/Lev5/rMPsi4_Asymptotic_GeometricUnits_CoM.h5?download=1
mv rMPsi4_Asymptotic_GeometricUnits_CoM.h5?download=1 rMPsi4_Sz1_m0p4_Sz2_0p8_q1p5data.h5

wget -c https://zenodo.org/record/1215122/files/SXS:BBH:0226/Lev5/rMPsi4_Asymptotic_GeometricUnits_CoM.h5?download=1
mv rMPsi4_Asymptotic_GeometricUnits_CoM.h5?download=1 rMPsi4_Sz1_m0p5_Sz2_m0p9_q1p5data.h5
 











 
 

 




 




