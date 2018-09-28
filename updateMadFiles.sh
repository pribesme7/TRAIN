#!/bin/bash

if [ "$#" -ne 1 ]; then
    echo " "
    echo "Illegal number of parameters - NO version defined!"
    echo "    try:  updateMadFiles.sh {nominal|2012|2015|2016}"
    echo " "
    exit 0
fi

version=$1
file=LHC_$version.mad

if [ ! -f MAD_PART/HLLHC/$file ]; then
    echo "Input File $file not found in the MAD_PART directory!"
    exit 0
fi

cd MAD_PART/HLLHC                                                                
#../../MAD-X/build/madx64.exe < $file
./madx64-sectormap-fixed < $file
cd ../..

#remove old files/link
rm train.manf
rm train.manb
rm train.optb
rm train.optf

#copy the newly generated files to current folder and change the name while tailoring some
cp -v MAD_PART/HLLHC/train.* .
tail -n +9 train.mapf > opt_train_$version.manf
tail -n +9 train.mapb > opt_train_$version.manb
rm train.mapf
rm train.mapb
cp train.optf opt_train_$version.optf
cp train.optb opt_train_$version.optb
cp train.surf opt_train_$version.surf
cp train.surb opt_train_$version.surb

#make proper symbolic links to the current choosen version
./updateLinkToOpticsFile.sh $version

echo 'DONE!'
