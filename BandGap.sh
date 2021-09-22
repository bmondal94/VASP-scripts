#!/bin/bash

PAGENAME=$(basename "$0")

function usage()
{
echo " " 
echo "Description: This script is to compute band gap using VASP OUTCAR/EIGENVAL file.

Syntax:  $PAGENAME [-seh] [-f FILEPATH]

s	:	If spin nondegenerate use "s" as 2nd argument. (default: spin degenerate, double occupancy)	
e	:	If bandgap will be calculated using EIGENVAL file. (default: OUTCAR file)
f	:	Path to OUTCAR or EIGENVAL file (excluding /OUTCAR or /EIGENVAL). (default:current directory.) 
"
echo " "
}

function read_arguments()
{
    OPTIND=1
    while getopts :sehf: options ; do
        case $options in

            s) TAGSPIN="1" ;;

            e) TAGEIGENVAL="1" ;;

            h) TAGHELP="1" ;;

	    f) FILEPATH=$OPTARG ;;

           \?) echo 'ERROR: Invalid option supplied' ;exit 1;; 

        esac
    done
}


FILEPATH='.'
read_arguments "$@"

if [ -n "$TAGHELP" ] ;then
	usage
	exit 0
fi

if [ -n "$TAGEIGENVAL" ]; then
	outcar="${FILEPATH}/EIGENVAL"
	efermi=7; exline=0
	NELECT=$(awk 'NR==6 {print $1}' $outcar)
	nkpt=$(awk 'NR==6 {print $2}' $outcar)
else
	outcar="${FILEPATH}/OUTCAR"
	efermi=$(grep -n "E-fermi" $outcar | cut -d: -f1)
	NELEC=$(awk '/NELECT/ {print $3}' $outcar)
	NELECT=${NELEC%.*} # Float to integer
	nkpt=$(awk '/NKPTS/ {print $4}' $outcar)
	exline=1 # Extraline contain the  'band No.  band energies     occupation' line in OUTCAR
fi



if [ -n "$TAGSPIN" ]; then
	homo=$NELECT
else
	homo=$((NELECT/2))
fi
lumo=$((homo+1))

Ge1=`awk -v var="$efermi" 'NR>var' $outcar | grep "   $homo    " | head -1 | awk '{print $2}'`
Ge2=`awk -v var="$efermi" 'NR>var' $outcar | grep "   $lumo    " | head -1 | awk '{print $2}'`
e1=`awk -v var="$efermi" 'NR>var' $outcar | grep "   $homo    " | head -$nkpt | sort -n -k 2 | tail -1 | awk '{print $2}'`
e2=`awk -v var="$efermi" 'NR>var' $outcar | grep "   $lumo    " | head -$nkpt | sort -n -k 2 | head -1 | awk '{print $2}'`
homoline=`awk -v var="$efermi" 'NR>var' $outcar | grep -n "   $homo    " | head -$nkpt | sort -n -k 3 | tail -1 | awk '{print $1}' | cut -d: -f1`
lumoline=`awk -v var="$efermi" 'NR>var' $outcar | grep -n "   $lumo    " | head -$nkpt | sort -n -k 3 | head -1 | awk '{print $1}' | cut -d: -f1`

kpointhomo=$(sed -n "$((efermi+homoline-homo-exline)) p" $outcar)
kpointlumo=$(sed -n "$((efermi+lumoline-lumo-exline)) p" $outcar)


#diff=`bc <<< "$e2-$e1"` ;diff=`calc $e2-$e1`
diff=$(awk -v l="$e2" -v h="$e1" 'BEGIN{print l-h}')
Gdiff=$(awk -v l="$Ge2" -v h="$Ge1" 'BEGIN{print l-h}')

Efermi=$(grep "E-fermi" $outcar | awk '{print $3}')

echo ""
echo FILENAME: $outcar
echo "----------------------------------------"
echo FermiEnergy: $Efermi eV
echo "VB(G): "  " E=" $Ge1 eV
echo "CB(G): "  " E=" $Ge2 eV
echo "Band energy difference (G):" $Gdiff eV
echo "----------------------------------------"
printf "HOMO/VBM:\t$kpointhomo \tband: $homo \tE= $e1 eV \n"
printf "LUMO/CBM:\t$kpointlumo \tband: $lumo \tE= $e2 eV \n"
echo "Band gap:" $diff eV
echo ""
