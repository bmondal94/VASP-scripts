#!/bin/bash
if [[ $# -eq 0 ]] ;then
echo " " 
echo "This script is to compute band gap using OUTCAR file

Syntax:  $(basename "$0") <path-to-OUTCAR> s

s			:	If spin nondegenerate use "s" as 2nd argument. (default: spin degenerate, double occupancy)	
<path-to-OUTCAR>	:	Write '.' for current directory. 
"
echo " "
exit 0
fi

outcar="$1/OUTCAR"
poscar="$1/POSCAR"

if [ "$2" == "s" ]; then
	homo=`awk '/NELECT/ {print $3+0}' $outcar`
	lumo=`awk '/NELECT/ {print $3+1}' $outcar`
else
	homo=`awk '/NELECT/ {print $3/2}' $outcar`
	lumo=`awk '/NELECT/ {print $3/2+1}' $outcar`
fi
nkpt=`awk '/NKPTS/ {print $4}' $outcar`
efermi=`grep -n "E-fermi" $outcar | cut -d: -f1`
#echo "All homo"
#awk -v var="$efermi" 'NR>var' $outcar | grep "   $homo    " | head -$nkpt | awk '{print $2}'
#echo "All lumo"
#awk -v var="$efermi" 'NR>var' $outcar | grep "   $lumo    " | head -$nkpt | awk '{print $2}'
#echo 

kall=true; multi=false
#kall=false; multi=false
if [ $kall == "true" ];then
	Ge1=`awk -v var="$efermi" 'NR>var' $outcar | grep "   $homo    " | head -1 | awk '{print $2}'`
	Ge2=`awk -v var="$efermi" 'NR>var' $outcar | grep "   $lumo    " | head -1 | awk '{print $2}'`
	e1=`awk -v var="$efermi" 'NR>var' $outcar | grep "   $homo    " | head -$nkpt | sort -n -k 2 | tail -1 | awk '{print $2}'`
	e2=`awk -v var="$efermi" 'NR>var' $outcar | grep "   $lumo    " | head -$nkpt | sort -n -k 2 | head -1 | awk '{print $2}'`
	homoline=`awk -v var="$efermi" 'NR>var' $outcar | grep -n "   $homo    " | head -$nkpt | sort -n -k 3 | tail -1 | awk '{print $1}' | cut -d: -f1`
	lumoline=`awk -v var="$efermi" 'NR>var' $outcar | grep -n "   $lumo    " | head -$nkpt | sort -n -k 3 | head -1 | awk '{print $1}' | cut -d: -f1`
	if [ $multi == "true" ]; then	
		ph1=`awk -v var="$efermi" 'NR>var' $outcar | grep "   $homo    " | head -$nkpt | sort -n -k 2 | tail -2 | head -n 1 | awk '{print $2}'`
		lh1=`awk -v var="$efermi" 'NR>var' $outcar | grep -n "   $homo    " | head -$nkpt | sort -n -k 3 | tail -2 | head -n 1 |awk '{print $1}' | cut -d: -f1`
		ph2=`awk -v var="$efermi" 'NR>var' $outcar | grep "   $homo    " | head -$nkpt | sort -n -k 2 | tail -3 | head -n 1 | awk '{print $2}'`
		lh2=`awk -v var="$efermi" 'NR>var' $outcar | grep -n "   $homo    " | head -$nkpt | sort -n -k 3 | tail -3 | head -n 1 |awk '{print $1}' | cut -d: -f1`
		ph3=`awk -v var="$efermi" 'NR>var' $outcar | grep "   $homo    " | head -$nkpt | sort -n -k 2 | tail -4 | head -n 1 | awk '{print $2}'`
		lh3=`awk -v var="$efermi" 'NR>var' $outcar | grep -n "   $homo    " | head -$nkpt | sort -n -k 3 | tail -4 | head -n 1 |awk '{print $1}' | cut -d: -f1`
		ph4=`awk -v var="$efermi" 'NR>var' $outcar | grep "   $homo    " | head -$nkpt | sort -n -k 2 | tail -5 | head -n 1 | awk '{print $2}'`
		lh4=`awk -v var="$efermi" 'NR>var' $outcar | grep -n "   $homo    " | head -$nkpt | sort -n -k 3 | tail -5 | head -n 1 |awk '{print $1}' | cut -d: -f1`
		ph5=`awk -v var="$efermi" 'NR>var' $outcar | grep "   $homo    " | head -$nkpt | sort -n -k 2 | tail -6 | head -n 1 | awk '{print $2}'`
		lh5=`awk -v var="$efermi" 'NR>var' $outcar | grep -n "   $homo    " | head -$nkpt | sort -n -k 3 | tail -6 | head -n 1 |awk '{print $1}' | cut -d: -f1`
		ph6=`awk -v var="$efermi" 'NR>var' $outcar | grep "   $homo    " | head -$nkpt | sort -n -k 2 | tail -7 | head -n 1 | awk '{print $2}'`
		lh6=`awk -v var="$efermi" 'NR>var' $outcar | grep -n "   $homo    " | head -$nkpt | sort -n -k 3 | tail -7 | head -n 1 |awk '{print $1}' | cut -d: -f1`
		ph7=`awk -v var="$efermi" 'NR>var' $outcar | grep "   $homo    " | head -$nkpt | sort -n -k 2 | tail -8 | head -n 1 | awk '{print $2}'`
		lh7=`awk -v var="$efermi" 'NR>var' $outcar | grep -n "   $homo    " | head -$nkpt | sort -n -k 3 | tail -8 | head -n 1 |awk '{print $1}' | cut -d: -f1`
		p1=`awk -v var="$efermi" 'NR>var' $outcar | grep "   $lumo    " | head -$nkpt | sort -n -k 2 | head -2 | tail -n 1 | awk '{print $2}'`
		l1=`awk -v var="$efermi" 'NR>var' $outcar | grep -n "   $lumo    " | head -$nkpt | sort -n -k 3 | head -2 | tail -n 1 |awk '{print $1}' | cut -d: -f1`
		p2=`awk -v var="$efermi" 'NR>var' $outcar | grep "   $lumo    " | head -$nkpt | sort -n -k 2 | head -3 | tail -n 1 |awk '{print $2}'`
		l2=`awk -v var="$efermi" 'NR>var' $outcar | grep -n "   $lumo    " | head -$nkpt | sort -n -k 3 | head -3 | tail -n 1 |awk '{print $1}' | cut -d: -f1`
		p3=`awk -v var="$efermi" 'NR>var' $outcar | grep "   $lumo    " | head -$nkpt | sort -n -k 2 | head -4 | tail -n 1 |awk '{print $2}'`
		l3=`awk -v var="$efermi" 'NR>var' $outcar | grep -n "   $lumo    " | head -$nkpt | sort -n -k 3 | head -4 | tail -n 1 |awk '{print $1}' | cut -d: -f1`
		p4=`awk -v var="$efermi" 'NR>var' $outcar | grep "   $lumo    " | head -$nkpt | sort -n -k 2 | head -5 | tail -n 1 | awk '{print $2}'`
		l4=`awk -v var="$efermi" 'NR>var' $outcar | grep -n "   $lumo    " | head -$nkpt | sort -n -k 3 | head -5 | tail -n 1 |awk '{print $1}' | cut -d: -f1`
		p5=`awk -v var="$efermi" 'NR>var' $outcar | grep "   $lumo    " | head -$nkpt | sort -n -k 2 | head -6 | tail -n 1 |awk '{print $2}'`
		l5=`awk -v var="$efermi" 'NR>var' $outcar | grep -n "   $lumo    " | head -$nkpt | sort -n -k 3 | head -6 | tail -n 1 |awk '{print $1}' | cut -d: -f1`
		p6=`awk -v var="$efermi" 'NR>var' $outcar | grep "   $lumo    " | head -$nkpt | sort -n -k 2 | head -7 | tail -n 1 |awk '{print $2}'`
		l6=`awk -v var="$efermi" 'NR>var' $outcar | grep -n "   $lumo    " | head -$nkpt | sort -n -k 3 | head -7 | tail -n 1 |awk '{print $1}' | cut -d: -f1`
		p7=`awk -v var="$efermi" 'NR>var' $outcar | grep "   $lumo    " | head -$nkpt | sort -n -k 2 | head -8 | tail -n 1 |awk '{print $2}'`
		l7=`awk -v var="$efermi" 'NR>var' $outcar | grep -n "   $lumo    " | head -$nkpt | sort -n -k 3 | head -8 | tail -n 1 |awk '{print $1}' | cut -d: -f1`
	fi
else
	to=`grep -n "k-point  1001 :       0.5000    0.5000    0.5000" $outcar | cut -d: -f1`
	Ge1=`awk -v var="$efermi" -v var2="$to" 'NR>var && NR<var2' $outcar | grep "   $homo    " | head -1 | awk '{print $2}'`
	Ge2=`awk -v var="$efermi" -v var2="$to" 'NR>var && NR<var2' $outcar | grep "   $lumo    " | head -1 | awk '{print $2}'`
	e1=`awk -v var="$efermi" -v var2="$to" 'NR>var && NR<var2' $outcar | grep "   $homo    " | head -$nkpt | sort -n -k 2 | tail -1 | awk '{print $2}'`
	e2=`awk -v var="$efermi" -v var2="$to" 'NR>var && NR<var2' $outcar | grep "   $lumo    " | head -$nkpt | sort -n -k 2 | head -1 | awk '{print $2}'`
	homoline=`awk -v var="$efermi" -v var2="$to" 'NR>var && NR<var2' $outcar | grep -n "   $homo    " | head -$nkpt | sort -n -k 3 | tail -1 | awk '{print $1}' | cut -d: -f1`
	lumoline=`awk -v var="$efermi" -v var2="$to" 'NR>var && NR<var2' $outcar | grep -n "   $lumo    " | head -$nkpt | sort -n -k 3 | head -1 | awk '{print $1}' | cut -d: -f1`
fi

kpointhomo=`sed -n "$((efermi+homoline-homo-1)) p" $outcar`
kpointlumo=`sed -n "$((efermi+lumoline-lumo-1)) p" $outcar`
kpointhomon=`sed -n "$((efermi+homoline-homo-1)) p" $outcar | awk '{print $2}'`
kpointlumon=`sed -n "$((efermi+lumoline-lumo-1)) p" $outcar | awk '{print $2}'`

#echo "$l1 $l2 $l3 $l4"
if [ $multi == "true" ]; then
	hp1=`sed -n "$((efermi+lh1-homo-1)) p" $outcar | awk '{print $2}'`
	hp2=`sed -n "$((efermi+lh2-homo-1)) p" $outcar | awk '{print $2}'`
	hp3=`sed -n "$((efermi+lh3-homo-1)) p" $outcar | awk '{print $2}'`
	hp4=`sed -n "$((efermi+lh4-homo-1)) p" $outcar | awk '{print $2}'`
	hp5=`sed -n "$((efermi+lh5-homo-1)) p" $outcar | awk '{print $2}'`
	hp6=`sed -n "$((efermi+lh6-homo-1)) p" $outcar | awk '{print $2}'`
	hp7=`sed -n "$((efermi+lh7-homo-1)) p" $outcar | awk '{print $2}'`
	kp1=`sed -n "$((efermi+l1-lumo-1)) p" $outcar | awk '{print $2}'`
	kp2=`sed -n "$((efermi+l2-lumo-1)) p" $outcar | awk '{print $2}'`
	kp3=`sed -n "$((efermi+l3-lumo-1)) p" $outcar | awk '{print $2}'`
	kp4=`sed -n "$((efermi+l4-lumo-1)) p" $outcar | awk '{print $2}'`
	kp5=`sed -n "$((efermi+l5-lumo-1)) p" $outcar | awk '{print $2}'`
	kp6=`sed -n "$((efermi+l6-lumo-1)) p" $outcar | awk '{print $2}'`
	kp7=`sed -n "$((efermi+l7-lumo-1)) p" $outcar | awk '{print $2}'`
fi

#diff=`bc <<< "$e2-$e1"` ;diff=`calc $e2-$e1`
diff=`awk -v l="$e2" -v h="$e1" 'BEGIN{print l-h}'`
Gdiff=`awk -v l="$Ge2" -v h="$Ge1" 'BEGIN{print l-h}'`
#echo "HOMO: " $kpointhomo " band:" $homo " E=" $e1
#echo "LUMO: " $kpointlumo " band:" $lumo " E=" $e2
#echo "Band gap:" $diff

Efermi=`grep "E-fermi" $outcar | awk '{print $3}'`

latticeparameter(){
latparam=$(python - <<PYTHON
import numpy as np
def vec(a):
        return (np.sqrt(np.sum(a**2)))
filename=open("$poscar","r")
lines=filename.readlines()
factor=float(lines[1].split()[0])
a=np.asarray(lines[2].split()).astype(np.float)
b=np.asarray(lines[3].split()).astype(np.float)
c=np.asarray(lines[4].split()).astype(np.float)
vector=np.array([vec(a),vec(b),vec(c)])
vector=vector*factor
S_factor=np.array([1., 1., 1.]) #* 6
primitive_vectors=np.divide(vector,S_factor)
nv=primitive_vectors**2
lp=  np.sqrt(2) * np.sqrt(nv[0] + nv[1] - nv[2])
#lp=  np.sqrt(2) * np.sqrt(nv[0] + nv[2] - nv[1])
#lp= np.sqrt(2) * np.sqrt(nv[1] + nv[2] - nv[0])
#lp=(lp-5.689)/5.689*100.0
print (lp)
filename.close()
PYTHON
)
echo "$latparam"
}

LP=`latticeparameter` #`head -n 2 $poscar | tail -n 1`
if [ $multi == "true" ];then
	#echo "LatticeParameter HOMOKpoint  HOMOEnergy  nLUMOpoint nLUMOEnergy ..... 7nLUMOpoint 7nLUMOEnergyLUMOKpoint  
	#LUMOEnergy nLUMOpoint nLUMOEnergy ..... 7nLUMOpoint 7nLUMOenergy  GammaHOMOEnergy  GammaLUMOEnergy"
	echo "$LP | $kpointhomon $e1 $hp1 $ph1  $hp2 $ph2  $hp3 $ph3  $hp4 $ph4  $hp5 $ph5  $hp6 $ph6 $hp7 $ph7 | $kpointlumon $e2 $kp1 $p1 $kp2 $p2 $kp3 $p3 $kp4 $p4 $kp5 $p5 $kp6 $p6 $kp7 $p7 | $Ge1 $Ge2"
else
	#echo "LatticeParameter  FermiEnergy  HOMOKpoint  HOMOEnergy  LUMOKpoint  LUMOEnergy GammaHOMOEnergy  GammaLUMOEnergy  GammaBandgap  Bandgap"
	echo $LP $Efermi $kpointhomon $e1 $kpointlumon $e2 $Ge1 $Ge2 $Gdiff $diff
fi
