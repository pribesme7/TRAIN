#!/bin/bash
# 5/02/2018 ARibes: 

if [ "$1" == "-h" ]; then
	echo " "
	echo "This script generates plots of the orbit offset and of the tune shift and saves them in a given folder."
	echo "It completes a MD cycle with the HL-LHC parameters."
	echo "The results are stored inside MDcycleimages using the stage of the cycle plus the number of the IP used"
	echo " "
	echo "    usage:  lazy.sh {IPnumber}"
	echo " "
  exit 0
fi

declare -a array=("collision_8_IP1" "flattop_B_IP1" "injection_std_IP1" "collision_8_IP1258" "flattop_B_IP1258" "injection_std_IP1258" "collision_8_IP15" "flattop_B_IP15" "injection_std_IP15" "collision_B_IP1" "flattop_std_IP1" "stable_8_IP1" "collision_B_IP1258" "flattop_std_IP1258" "stable_8_IP1258" "collision_B_IP15" "flattop_std_IP15" "stable_8_IP15" "collision_std_IP1" "injection_8_IP1" "stable_B_IP1" "collision_std_IP1258" "injection_8_IP1258" "stable_B_IP1258" "collision_std_IP15" "injection_8_IP15" "stable_B_IP15" "flattop_8_IP1" "injection_B_IP1" "stable_std_IP1" "flattop_8_IP1258" "injection_B_IP1258" "stable_std_IP1258" "flattop_8_IP15" "injection_B_IP15" "stable_std_IP15")

#array = '(*)';

for i in "${array[@]}";
do
echo "$i"
python myplots.py "$i"
python myplots2.py "$i";
done
