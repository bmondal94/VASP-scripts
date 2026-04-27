# VASP-scripts

__NOTE__: To be able to use the scripts as executable, 
1. Change the mode of the files as executable using `chmod a+x <file_name>`
2. Add the folder path of all these scripts to .bashrc: e.g. `export PATH="/path/directory_scripts:$PATH"`

## Description of the scripts
1. `BandGap` : This script computes band gap using VASP OUTCAR/EIGENVAL file. [Help: `$BandGap -h`]
2. `CreatePotcarVasp` : This scripts reads the POSCAR file and creats the POTCAR file of corresponding elements (in order) [Help: `$CreatePotcarVasp -h`]\
   __NOTE__: Update the `POTDIR='/projects/p_reactivity/jalu038/VASP/potpaw_PBE.54/` line inside the file to the directory path where your Potcar files are.
4. `CheckFinishAll` : To check the convergence of a vasp run. [Help: `$CheckFinishAll -h`]
5. `Bandstructure.py` : Python scripts to plot the band structure using VASP EIGENVAL, PROCAR, POSCAR and DOSCAR file.
6. `sqs2poscarconvert.py` : Python scripts to convert sqs output from ATAT software to VASP POSCAR.
7. `MakeVacuum.py` : Python scripts to create slab+vacuum in zincblende structure.
8. `SqsGenerator` : Create SQS cell using ATAT software.
9. `subvasp` : VASP job submit scripts in cluster. [Help: `$subvasp -h`]
10. `subvaspman` : Manual of subvasp.
11. `suborca` : VASP job submit scripts in cluster. [Help: `$suborca -h`]
12. `LatticeVectors.py` : To get the lattice parameters for zincblende/cubic/strained zincblende structure from POSCAR file.
