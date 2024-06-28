# VASP-scripts
1. `BandGap` : This script computes band gap using VASP OUTCAR/EIGENVAL file. [Help: `$BandGap -h`]
2. `CreatePotcarVasp` : This scripts reads the POSCAR file and creats the POTCAR file of corresponding elements (in order) [Help: `$BandGap -h`]
3. `CheckFinishAll` : To check the convergence of a vasp run. [Help: `$CheckFinishAll -h`]
4. `Bandstructure.py` : Python scripts to plot the band structure using VASP EIGENVAL, PROCAR, POSCAR and DOSCAR file.
5. `sqs2poscarconvert.py` : Python scripts to convert sqs output from ATAT software to VASP POSCAR.
6. `MakeVacuum.py` : Python scripts to create slab+vacuum in zincblende structure.
7. `SqsGenerator` : Create SQS cell using ATAT software.
8. `subvasp` : VASP job submit scripts in cluster. [Help: `$subvasp -h`]
9. `subvaspman` : Manual of subvasp.
9. `suborca` : VASP job submit scripts in cluster. [Help: `$suborca -h`]
10. `LatticeVectors.py` : To get the lattice parameters for zincblende/cubic/strained zincblende structure from POSCAR file.
