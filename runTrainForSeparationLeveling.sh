#!/bin/bash
# 2/04/2016 AGorzawski: 
#######################################################
##  Here edit what to use for the all calc process
##
##  fillingSchemeFile one of the *.in files in the folder
##  version = { nominal 2015 2016 }

#
#fillingSchemeFile=train25_2040_72.in
fillingSchemeFile=train_MDLR.in.emittance
version=2016


#fillingSchemeFile=train_4440_170m.in
#version=2015

#######################################################

if [ "$#" -ne 2 ]; then
    echo " "
    echo "Illegal number of parameters no IP selected no FILE with input values!"
    echo " use: script {ip1|ip2|ip5|ip8} {filenamewithsteps}"
    echo " Actual settings used for leveling script are: filling scheme: $fillingSchemeFile and optics: $version"
	
    exit 0
fi


#######################################################
# The following part will iterate over the step (in mm) in given file and will assign the separation
# for a given IP. will run the TRAIN code
#
iptolevel=$1
filewithstepstolevel=$2
foldername=$(date +%Y%m%d)"_"$iptolevel"_optics_"$version"_sep_lev_for_$fillingSchemeFile"
mkdir -p "$foldername"
templatefile=MAD_PART/collisionConfiguration.$version.tmp
outputfile=MAD_PART/collisionConfiguration.$version

referencePlotFolder="notset"

for separation in `cat $filewithstepstolevel`;
do
	echo "################ Update coll conf $version for $separation ################" 
	echo 'GENERATING the setup.input file for '$fillingSchemeFile

	IP1SEP=0.0
	if [ $iptolevel == "ip1" ]; then
		IP1SEP=$separation
	fi
	IP2SEP=0.0
	if [ $iptolevel == "ip2" ]; then
		IP2SEP=$separation
	fi
	IP5SEP=0.200
	if [ $iptolevel == "ip5" ]; then
		IP5SEP=$separation
	fi
	IP8SEP=0.0
	if [ $iptolevel == "ip8" ]; then
		IP8SEP=$separation
	fi
	echo "$(eval "echo \"$(cat $templatefile)\"")" > $outputfile

	echo "update mad files $version"
	./updateMadFiles.sh $version
	#cat MAD_PART/collisionConfiguration.2015

	echo "update mad files $version"

	echo "Run TRAIN $fillingSchemeFile"
	./runTrainForFillingScheme.sh $fillingSchemeFile noplot
	touch RESULTS/$fillingSchemeFile/testFile.hhh
	touch RESULTS/$fillingSchemeFile/testFile.ggg
	touch RESULTS/$fillingSchemeFile/testFile

#copy the results to separate folder
	currentResultFolder=$foldername/$fillingSchemeFile.$iptolevel.$separation
	if [ $referencePlotFolder == "notset" ]; then
		referencePlotFolder=$currentResultFolder
	fi

	echo "$currentResultFolder"
	mkdir -p "$currentResultFolder"
	cp RESULTS/$fillingSchemeFile/* $currentResultFolder

#plot the results 
	python plotVerHorOffsetsComparison.py $referencePlotFolder/ $currentResultFolder/ falseForDislpay
	echo "################ DONE for $separation #################################"

done

cp $filewithstepstolevel $foldername

# create  an animated git from the all plots.
#convert -delay 250 -loop 0 *.png $foldername.gif

echo "DONE. Search the $foldername for the results"
