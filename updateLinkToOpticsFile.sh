#!/bin/bash

if [ "$#" -ne 1 ]; then
    echo " "
    echo "Illegal number of parameters - NO version defined!"
    echo "    try:  updateLinkToOpticsFile.sh {nominal|nominal_ip15|nominal_ip1258|2012|2015|2016|2016_ip158}"
    echo " "
    exit 0
fi
version=$1

file=opt_train_$version.optf
if [ ! -f $file ]; then
    echo "One of the input Files $file was not found in the current directory, run the updateMadFiles.sh! to re-do the optics files!"
    exit 0
fi

rm train.manf
rm train.manb
rm train.optb
rm train.optf

#echo "REMOVE not correct BB elements (MK*16) (if still persists)"
#sed '/MKIP[1258]P[LR]16/d' ./opt_train_$version.optf > file && mv file opt_train_$version.optf 
#sed '/MKIP[1258]P[LR]16/d' ./opt_train_$version.optb > file && mv file opt_train_$version.optb 
#sed '/MKIP[1258]P[LR]16/d' ./opt_train_$version.manf > file && mv file opt_train_$version.manf
#sed '/MKIP[1258]P[LR]16/d' ./opt_train_$version.manb > file && mv file opt_train_$version.manb 
#echo "done"

ln -sv opt_train_$version.optf train.optf
ln -sv opt_train_$version.optb train.optb
ln -sv opt_train_$version.manf train.manf
ln -sv opt_train_$version.manb train.manb
