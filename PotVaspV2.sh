#!/bin/sh

myPWD=$(pwd)
POTDIR='/projects/p_reactivity/jalu038/VASP/potpaw_PBE.54/'
extension=""

function usage
{
	echo " "
	echo "* Description: This scripts reads the POSCAR file and creats the POTCAR file of corresponding elements (in order)"
	echo "* Usage: $(basename "$0") [-h|--help|-gSsdGpH] [-M ELEMENTs_EXTENTION]"
	echo "  -h | --help :  For help."
	echo "  -g          :  If you want the _GW extension of the POTCAR"
	echo "  -S          :  If you want the _s extension of the POTCAR"
	echo "  -s          :  If you want the _sv extension of the POTCAR"
	echo "  -d          :  If you want the _d extension of the POTCAR"
	echo "  -G          :  If you want the _sv_GW extension of the POTCAR"
	echo "  -p          :  If you want the _pv extension of the POTCAR"
	echo "  -H          :  If you want the _h extension of the POTCAR"
	echo "  -M          :  If you want mixed extension. Pass which extension POTCAR you want by element, for all element. Order doesn't matter. e.g. -M 'Ga_d As_d N_s' 
                 Inverted commas around the passed arguments are mandatory."
	echo "                 (The default is standered version with no extension)"
	echo "* Note:  To create the POTCAR you must need the POSCAR file. You first have to move to the"
	echo "         folder that contains the POSCAR. Note, it can't create POTCAR from external path."
	echo " "
	echo "---------------- Enjoy. Have a good day. ---------------------"
	echo " "
	exit
}
function readarguments
{
	OPTIND=1
	while getopts :gSsdGpHM: options; do
		case $options in
			g) extension='_GW' ;;
			S) extension='_s' ;;
			s) extension='_sv' ;;
			d) extension='_d' ;;
			G) extension='_sv_GW' ;;
			p) extension='_pv' ;;
			H) extension='_h' ;;
			M) MIXED=1
			   ALLPOTCAR=( $OPTARG ) ;;
			\?) echo  "Invalid option: $1"; exit 1 ;;
		esac
	done
}
if [ -n "$1" ]; then
	if [ "$1" == "-h" ] || [ "$1" == "--help" ]; then
		usage
	else
		readarguments "$@"
	fi
fi
	

if [ -f POSCAR ]; then
	if [ -s POTCAR ]; then
		printf "\n\e[38;5;9;4m* Warning:\e[0m Hi $(whoami), you already have an old POTCAR. I am deleting it and creating the new one.\n\n"
		rm POTCAR
	fi
	myAGREP=$(head -6 POSCAR | tail -1)
	echo "POTCAR for these elements will be created (in order): " $myAGREP
	echo " "
	for i in $myAGREP
	do
		if [ -n "$MIXED" ];then
			for I in ${ALLPOTCAR[@]};do
				if [[ ${I%%_*} == $i ]]; then 
					extension=''
					[[ $I =~ '_' ]] && extension="_${I#*_}"
					break 
				fi
			done
		fi
		if [ ! -d "$POTDIR$i$extension" ]; then
			printf "\n\e[38;5;1m* Error: $POTDIR$i$extension doesn't exist.\e[0m \n\n "
			exit
		fi
		cat $POTDIR$i$extension/POTCAR >> $myPWD/POTCAR
	done
else
	printf "\n\e[38;5;1m* Error: No POSCAR file here! POSCAR file is mandatory. Please move to the file contains the POSCAR. \e[0m\n\n"
fi

exit

