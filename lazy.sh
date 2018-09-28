#!/bin/bash
# 5/02/2018 ARibes: 

if [ "$1" == "-h" ]; then
	echo " "
	echo "This script simplfies the calculations with the TRAIN code."
	echo "It completes a MD cycle with the HL-LHC parameters."
	echo "The results are stored inside RESULTS/MDcycle using the stage of the cycle plus the number of the IP used"
	echo " "
	echo "    usage:  lazy.sh {IPnumber}"
	echo " "
  exit 0
fi

echo "INJECTION"

./updateMadFiles.sh HL_injection
./runTrainForFillingScheme.sh 25ns_2760b_2748_2494_2572_288bpi_13inj.in noPlot
./subfolders.sh "injection_std_IP$1"

echo " "
echo "Finished with injection standard"
echo " "

./updateMadFiles.sh HL_injection_B8
./runTrainForFillingScheme.sh 25ns_2748b_2736_2258_2374_288bpi_12inj.in noPlot
./subfolders.sh "injection_B_IP$1"

echo " "
echo "Finished with injection BCMS"
echo " "

./runTrainForFillingScheme.sh 8b4e_1972b_1967_1178_1886_224bpi_12inj.in noPlot
./subfolders.sh "injection_8_IP$1"

echo " "
echo "Finished with injection 8"
echo " "


echo "END OF RAMP"

./updateMadFiles.sh HL_flattop
./runTrainForFillingScheme.sh 25ns_2760b_2748_2494_2572_288bpi_13inj.in noPlot
./subfolders.sh "flattop_std_IP$1"

echo " "
echo "Finished with flattop standard"
echo " "

./updateMadFiles.sh HL_flattop_B8
./runTrainForFillingScheme.sh 25ns_2748b_2736_2258_2374_288bpi_12inj.in noPlot
./subfolders.sh "flattop_B_IP$1"

echo " "
echo "Finished with flattop BCMS"
echo " "

./runTrainForFillingScheme.sh 8b4e_1972b_1967_1178_1886_224bpi_12inj.in noPlot
./subfolders.sh "flattop_8_IP$1"

echo " "
echo "Finished with flattop 8"
echo " "

echo "START COLLISION (64cm)"
./updateMadFiles.sh HL_collision

./runTrainForFillingScheme.sh 25ns_2760b_2748_2494_2572_288bpi_13inj.in noPlot
./subfolders.sh "collision_std_IP$1"

echo " "
echo "Finished with collision standard"
echo " "

./runTrainForFillingScheme.sh 25ns_2748b_2736_2258_2374_288bpi_12inj.in noPlot
./subfolders.sh "collision_B_IP$1"

echo " "
echo "Finished with collision BCMS"
echo " "

./runTrainForFillingScheme.sh 8b4e_1972b_1967_1178_1886_224bpi_12inj.in noPlot
./subfolders.sh "collision_8_IP$1"

echo " "
echo "Finished with collision 8"
echo " "


echo "END OF SQUEEZE"

./updateMadFiles.sh HL_squeeze

./runTrainForFillingScheme.sh 25ns_2760b_2748_2494_2572_288bpi_13inj.in noPlot
./subfolders.sh "squeeze_std_IP$1"

echo " "
echo "Finished with squeeze standard"
echo " "

./runTrainForFillingScheme.sh 25ns_2748b_2736_2258_2374_288bpi_12inj.in noPlot
./subfolders.sh "squeeze_B_IP$1"

echo " "
echo "Finished with squeeze BCMS"
echo " "

./runTrainForFillingScheme.sh 8b4e_1972b_1967_1178_1886_224bpi_12inj.in noPlot
./subfolders.sh "squeeze_8_IP$1"

echo " "
echo "Finished with squeeze 8"
echo " "

echo "COLLISION ULTIMATE (41cm)"

./updateMadFiles.sh HL_collision_ultimate

./runTrainForFillingScheme.sh 25ns_2760b_2748_2494_2572_288bpi_13inj.in noPlot
./subfolders.sh "collision_ultimate_std_IP$1"

echo " "
echo "Finished with collision ultimate standard"
echo " "

./runTrainForFillingScheme.sh 25ns_2748b_2736_2258_2374_288bpi_12inj.in noPlot
./subfolders.sh "collision_ultimate_B_IP$1"

echo " "
echo "Finished with collision ultimate BCMS"
echo " "

./runTrainForFillingScheme.sh 8b4e_1972b_1967_1178_1886_224bpi_12inj.in noPlot
./subfolders.sh "collision_ultimate_8_IP$1"

echo " "
echo "Finished with collision ultimate 8"
echo " "


echo "COLLISION (15cm)"

./updateMadFiles.sh HL_stable

./runTrainForFillingScheme.sh 25ns_2760b_2748_2494_2572_288bpi_13inj.in noPlot
./subfolders.sh "stable_std_IP$1"

echo " "
echo "Finished with stable standard"
echo " "

./runTrainForFillingScheme.sh 25ns_2748b_2736_2258_2374_288bpi_12inj.in noPlot
./subfolders.sh "stable_B_IP$1"

echo " "
echo "Finished with stable BCMS"
echo " "

./runTrainForFillingScheme.sh 8b4e_1972b_1967_1178_1886_224bpi_12inj.in noPlot
./subfolders.sh "stable_8_IP$1"

echo " "
echo "Finished with stable 8"
echo " "
