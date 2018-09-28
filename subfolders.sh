#!/bin/bash
# 16/1/18 ARibes:

if [ "$1" == "-h" ]; then
	echo " "
	echo "This script copies the results in the input subfolder. It also copies the files train.*"
	echo " "
	echo " "
  exit 0
fi

if [ "$#" -ne 1 ]; then
	echo " "
	echo "Illegal number of parameters - NO result folder!"
	echo " "
	exit 0
fi 

mkdir -p RESULTS/

cd RESULTS/

mkdir -p $1

#cd ..


if ls 25* &> /dev/null; then
   mv 25* "$1"
else
   mv 8* "$1"
fi

#path = "RESULTS/amtrain/$1"
cd ..

cp train.* RESULTS/$1/



