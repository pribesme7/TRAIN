#!/bin/bash
# 2/04/2016 AGorzawski: 
#	

if [ "$1" == "-h" ]; then
	echo " "
	echo "This script simplfies the calculations with the TRAIN code FOR AN EXISTING input file (setup.input)."
	echo "It copies the results into the folder specified in first argument (under RESULTS folder)"
	echo "Some integration in the 'setup.input' might be needed"
	echo " "
	echo "    usage:  runTrainWithInputName.sh {ResultFolder} {plot|noplot}"
	echo " "
  exit 0
fi

if [ "$#" -ne 2 ]; then
	echo " "
	echo "Illegal number of parameters - NO result folder!"
	echo "    try:  runTrainWithInputName.sh -h"
	echo " "
	exit 0
fi
#keep this name unchanged - moveResultsToSubFolder.sh uses the same subfolder
mkdir -p RESULTS
resultFolder=$1

echo 'Reminder for a current optics set:'
ls -l train.*
echo 'Running TRAIN code...'

./amtrain < setup.input

echo '...DONE WITH TRAIN!' 

./moveResultsToSubFolder.sh $resultFolder

if [ "$2" == "plot" ]; then
	python plotVerHorOffsets.py RESULTS/$resultFolder/
#	python myplot.py RESULTS/$resultFolder/
fi
