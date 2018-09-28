#!/bin/bash
# 2/04/2016 AGorzawski: 
#######################################################
##  Here edit what to use for the all calc process
##
##  fillingSchemeFile one of the *.in files in the folder
##  version = { nominal 2015 2016 }

#
#fillingSchemeFile=train25_2040_72.in
fillingSchemeFile=train25nom.in.emittance
version=2016


#fillingSchemeFile=train_4440_170m.in
#version=2015

#######################################################

if [ "$#" -ne 2 ]; then
    echo " "
    echo "Illegal number of parameters no IP selected no FILE with input values!"
    echo " use: script {ip1|ip5} {filenamewithsteps}"
    echo " Actual settings used for leveling script are: filling scheme: $fillingSchemeFile and optics: $version"	
    exit 0
fi


#######################################################
# The following part will iterate over the step (in mm) in given file and will assign the separation
# for a given IP. will run the TRAIN code and agreggate the code
#
iptolevel=$1
filewithstepstolevel=$2
foldername=$(date +%Y%m%d)"_"$iptolevel"_optics_"$version"_sep_scan_for_$fillingSchemeFile"
mkdir -p "$foldername"
templatefile=MAD_PART/collisionConfiguration.$version.tmp.XingPlaneSeparation
outputfile=MAD_PART/collisionConfiguration.$version

referencePlotFolder="notset"

while read -r line
do
	stringarray=($line)
	separation=${stringarray[0]}
	separationXing=${stringarray[1]}

	echo "################ UPDATE collision conf for optics $version and Separations: SEPARATION plane $separation, XING plane: $separationXing" ################" 
	echo '################ GENERATING the setup.input file for the given filling scheme: '$fillingSchemeFile

	IP1SEP=0.0
	if [ $iptolevel == "ip1" ]; then
		IP1SEP=$separation
	fi
	IP1SEPXING=0.0
	if [ $iptolevel == "ip1" ]; then
		IP1SEPXING=$separationXing
	fi
	IP2SEP=0.0
	IP8SEP=0.0
	IP5SEP=0.0
	if [ $iptolevel == "ip5" ]; then
		IP5SEP=$separation
	fi
	IP5SEPXING=0.0
	if [ $iptolevel == "ip5" ]; then
		IP5SEPXING=$separationXing
	fi

	echo "$(eval "echo \"$(cat $templatefile)\"")" > $outputfile
	echo '################ ...DONE.'; 
	echo "################ Run MAD with updated files and with optics: $version"
	./updateMadFiles.sh $version
	cat MAD_PART/collisionConfiguration.2016

	echo "################ ...DONE for MAD files with optics: $version"
	echo "################ RUN TRAIN $fillingSchemeFile"

	./runTrainForFillingScheme.sh $fillingSchemeFile noplot
	touch RESULTS/$fillingSchemeFile/testFile.hhh
	touch RESULTS/$fillingSchemeFile/testFile.ggg
	touch RESULTS/$fillingSchemeFile/testFile

##copy the results to separate folder
	currentResultFolder=$foldername/$fillingSchemeFile.$iptolevel.SEP.$separation.XING.$separationXing
	if [ $referencePlotFolder == "notset" ]; then
		referencePlotFolder=$currentResultFolder
	fi

	echo "################ MOVING TRAIN results to : $currentResultFolder"
	mkdir -p "$currentResultFolder"
	cp RESULTS/$fillingSchemeFile/* $currentResultFolder

##plot the results 
#	python plotVerHorOffsetsComparison.py $referencePlotFolder/ $currentResultFolder/ falseForDislpay
	echo "################ DONE for SEPARATION plane $separation , XING plane: $separationXing #################################"

done < "$filewithstepstolevel"

#cp $filewithstepstolevel $foldername
# create  an animated git from the all plots.
#convert -delay 250 -loop 0 *.png $foldername.gif

echo "DONE. Search the $foldername for the results"
