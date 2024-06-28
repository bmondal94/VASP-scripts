#!/bin/bash

N=(00.5  00.9  01.9  02.8  03.7  04.6  05.6  07.4  09.3 11.1)
#N=(07.4)
sarray=(0.0 -0.5 -1 -1.5 -2 -2.5 -3 -3.5 -4 -4.5 -5)
bandid=1729 #Band index for which you want to get s-p contributions
unfold_point=216
BW_index=$(($((bandid-1))*unfold_point))
file='ThermExCompNearestConfig.txt'
CDDIR='/scratch/surfacechem/mondalb/project_DFT_substrate_effect/band_diagram/GaAsN/mcsqs/LatestNSQSJ302020_T1/ThermExComp'
CDDIR2='/scratch/surfacechem/mondalb/project_DFT_substrate_effect/band_diagram/GaAsN/mcsqs/LatestNSQSJ302020_T1/Eqm'
rm $CDDIR/s-p-PROCAR.dat $CDDIR/NearestConfigAfterBandgap.dat $CDDIR/BlochWeightEnergyGLDX.dat $CDDIR/BlochWeightGLDX.dat

printf "\n\e[38;5;9;4m++ Check point: Please ensure the N concentration in '$file' is in this order (1st coloumn) \n\e[0m"
echo ${N[@]}
printf "\n\e[38;5;9;4m++ Check point: Please ensure only N concentration that are switch on in '$file' are in this order (1st coloumn) \n\e[0m"
echo ${N[@]}
printf "\n\n\e[38;5;9;4m++ Check point: Please ensure the S values in '$file' is in this order (each row[2:]) \n\e[0m"
echo ${sarray[@]}
echo "** Do you want to continue? Press Yes/Y if yes, else press ENTER" 
read contin
if [[ ${contin,,} == y* ]]; then
        echo "Thanks for letting me go ahead"
else
        echo "Hi, $USER! You choose not to go further."
	echo " "
        exit
fi

CheckFinishBandgap(){
printf "++ Bandgap Convergency check\n\n"
if [ -s "OUTCAR" ] ; then
        if grep -q "METAGGA = MBJ" OUTCAR; then
		nelmoutcar=$(awk '/NELM/ {print $3}' OUTCAR | cut -d ';' -f1)
		nelmoszicar=$(tail -n 2 OSZICAR | head -n 1 | awk '{print $2}')
		nelmdiff=$(calc $nelmoutcar-$nelmoszicar)
		if [[ $nelmdiff -ne 0 ]];then
        		if [ $( grep -L 'Total CPU time used' OUTCAR) ] ; then
                		printf "\e[38;5;213m-- Warning:\e[0m OUTCAR didn't finish.\n"
				exit
			fi
		else
			printf "\e[38;5;213m-- Warning:\e[0m OUTCAR didn't finish properly. All the NELM step is consumed.\n"
			exit
		fi
		printf "++ Bandgap Convergency check successful\n\n"
	else
                printf "\e[38;5;91m++ Warning:\e[0m MBJ metagga is not used in OUTCAR.\n"
                i=$(awk '/ISIF/ {print $3}' OUTCAR)
                if [[ $i -ne 0 ]]; then
                        printf "\e[38;5;91m++ Warning:\e[0m ISIF in OUTCAR = $i (not 0).\n"
                fi
		exit
        fi
else
        printf "\e[38;5;1m** Error:\e[0m Empty/Missing OUTCAR.\n"
	exit
fi
}

st=$(date)
no_read='#'
i=0
while read -a farray
do
	if [[ ${farray[0]} == ${no_read}* ]]; then
                continue
        fi

	s=("${N[i]}|")
	BG=("${N[i]}|")
	BW=("${N[i]}|")
	BWE=("${N[i]}|")
	count=1
	for J in ${sarray[@]}
	do
		if [ -z ${farray[count]} ]; then
			echo "Configuration in $file is missing corresponding to S$J"
			exit
		fi
		if [[ $count -eq 1 ]]; then
			fname=$CDDIR2/N${N[i]}/${farray[count]}
		else
			fname=$CDDIR/S$J/N${N[i]}/${farray[count]}
		fi
		cd $fname
		echo " "
		echo $(pwd)
		#CheckFinishBandgap
		#rm conf*.out CONTCAR DOSCAR PCDAT vasprun.xml XDATCAR CHG* REPORT EIGENVAL
		#nkpts=$(awk '/NKPTS/ {print $4}' OUTCAR)
		#nbands=$(awk '/NBANDS/ {print $NF}' OUTCAR)
		### Collect Bloch weight for G-L-D_m-X
		#rm BW_*.txt
		#bwindex=$BW_index
		#for K in $(seq 1 $nkpts); do
		#        sed -n "$((bwindex+1)),$((bwindex+unfold_point))p;$((bwindex+unfold_point+1))q" WAVECAR_spinor1.f2b >> BW_CB_spinor1.txt
		#        sed -n "$((bwindex+1)),$((bwindex+unfold_point))p;$((bwindex+unfold_point+1))q" WAVECAR_spinor2.f2b >> BW_CB_spinor2.txt
		#        bwindex=$((bwindex+$((nbands*unfold_point))))
		#done
		#rm *.f2b*
		for KK in BW_*.txt; do
			#bwe1=$(head -n 1 $KK | awk '{print $4}')
			BW_G=$(grep "0.000000   0.000000   0.000000" $KK | awk '{print $NF}')
			BW_L1=$(grep "0.000000   0.000000   0.500000" $KK | awk '{print $NF}')
			BW_L2=$(grep "0.000000   0.500000   0.000000" $KK | awk '{print $NF}')
			BW_X1=$(grep "0.000000   0.500000   0.500000" $KK | awk '{print $NF}')
			BW_L3=$(grep "0.500000   0.000000   0.000000" $KK | awk '{print $NF}')
			BW_X2=$(grep "0.500000   0.000000   0.500000" $KK | awk '{print $NF}')
			BW_X3=$(grep "0.500000   0.500000   0.000000" $KK | awk '{print $NF}')
			BW_L4=$(grep "0.500000   0.500000   0.500000" $KK | awk '{print $NF}')
			BW_D=$(grep "0.000000   0.428200   0.428200" $KK | awk '{print $NF}')
			#bwe2=$(tail -n 1 $KK | awk '{print $4}')
			BW+="${BW_G},${BW_L1},${BW_L2},${BW_L3},${BW_L4},${BW_X1},${BW_X2},${BW_X3},${BW_D};"
			#BWE+="${bwe1},${bwe2};"
		done
		BW+="|"
		BWE+="|"
		## Collect Bandgap from OUTCAR
		#BG+="$(gap.sh '.' s)|"
		### Collect VASP projected s-p contribution from PROCAR
		#if [ -s PROCAR ];then
		#	lorbit=$(awk '/LORBIT/ {print $3}' OUTCAR)
		#	if [[ $lorbit -eq 12 ]] || [[ $lorbit -eq 11 ]];then
		#		s+="$(sed -n "/band  $bandid/,/band  $((bandid+1))/p" PROCAR | grep "^tot" | head -n 1 | awk '{print $2,$3,$4,$5,$11}')|"
		#	else
		#		echo "LORBIT=$lorbit found. This is not implemented. Only LORBIT=12/11 is allowed."
		#	fi
		#else
		#	echo "** Empty/Missing PROCAR"
		#fi
		((count++))
	done
	echo "$s">>$CDDIR/s-p-PROCAR.dat
	echo " ">>$CDDIR/s-p-PROCAR.dat
	echo "$BG">>$CDDIR/NearestConfigAfterBandgap.dat
	echo " ">>$CDDIR/NearestConfigAfterBandgap.dat
	echo "$BW">>$CDDIR/BlochWeightGLDX.dat
	echo " ">>$CDDIR/BlochWeightGLDX.dat
	echo "$BWE">>$CDDIR/BlochWeightEnergyGLDX.dat
	echo " ">>$CDDIR/BlochWeightEnergyGLDX.dat
((i++))
done < $CDDIR/$file

echo "Start time  : $st"
echo "Finish time : $(date)"

exit
