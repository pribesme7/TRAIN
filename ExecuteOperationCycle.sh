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


echo "RAMP&SQUEEZE"

./updateMadFiles.sh HL_rampsqueeze
./runTrainForFillingScheme.sh 25ns_2760b_2748_2494_2572_288bpi_13inj.in noPlot
./subfolders.sh "rampsqueeze_std_IP$1"

echo " "
echo "Finished with rampsqueeze standard"
echo " "

./updateMadFiles.sh HL_rampsqueeze_B8
./runTrainForFillingScheme.sh 25ns_2748b_2736_2258_2374_288bpi_12inj.in noPlot
./subfolders.sh "rampsqueeze_B_IP$1"

echo " "
echo "Finished with rampsqueeze BCMS"
echo " "

./runTrainForFillingScheme.sh 8b4e_1972b_1967_1178_1886_224bpi_12inj.in noPlot
./subfolders.sh "rampsqueeze_8_IP$1"

echo " "
echo "Finished with rampsqueeze 8"
echo " "



echo "PRESQUEEZE"

./updateMadFiles.sh HL_presqueeze
./runTrainForFillingScheme.sh 25ns_2760b_2748_2494_2572_288bpi_13inj.in noPlot
./subfolders.sh "presqueeze_std_IP$1"

echo " "
echo "Finished with presqueeze standard"
echo " "

./updateMadFiles.sh HL_presqueeze_B8
./runTrainForFillingScheme.sh 25ns_2748b_2736_2258_2374_288bpi_12inj.in noPlot
./subfolders.sh "presqueeze_B_IP$1"

echo " "
echo "Finished with presqueeze BCMS"
echo " "

./runTrainForFillingScheme.sh 8b4e_1972b_1967_1178_1886_224bpi_12inj.in noPlot
./subfolders.sh "presqueeze_8_IP$1"

echo " "
echo "Finished with presqueeze 8"
echo " "




echo "START COLLISION (41cm)"
./updateMadFiles.sh HL_collision_41

./runTrainForFillingScheme.sh 25ns_2760b_2748_2494_2572_288bpi_13inj.in noPlot
./subfolders.sh "collision_41_std_IP$1"

echo " "
echo "Finished with collision 41 standard"
echo " "

./runTrainForFillingScheme.sh 25ns_2748b_2736_2258_2374_288bpi_12inj.in noPlot
./subfolders.sh "collision_41_B_IP$1"

echo " "
echo "Finished with collision 41 BCMS"
echo " "

./runTrainForFillingScheme.sh 8b4e_1972b_1967_1178_1886_224bpi_12inj.in noPlot
./subfolders.sh "collision_41_8_IP$1"

echo " "
echo "Finished with collision 41 8"
echo " "


echo "STABLE B 41"
./updateMadFiles.sh HL_stable_b_41

./runTrainForFillingScheme.sh 25ns_2760b_2748_2494_2572_288bpi_13inj.in noPlot
./subfolders.sh "stable_b_41_std_IP$1"

echo " "
echo "Finished with stable b 41 standard"
echo " "

./runTrainForFillingScheme.sh 25ns_2748b_2736_2258_2374_288bpi_12inj.in noPlot
./subfolders.sh "stable_b_41_B_IP$1"

echo " "
echo "Finished with stable b 41 BCMS"
echo " "

./runTrainForFillingScheme.sh 8b4e_1972b_1967_1178_1886_224bpi_12inj.in noPlot
./subfolders.sh "stable_b_41_8_IP$1"

echo " "
echo "Finished with stable_b 41 8"
echo " "



echo "START COLLISION (64cm)"
./updateMadFiles.sh HL_collision_64

./runTrainForFillingScheme.sh 25ns_2760b_2748_2494_2572_288bpi_13inj.in noPlot
./subfolders.sh "collision_64_std_IP$1"

echo " "
echo "Finished with collision 64 standard"
echo " "

./runTrainForFillingScheme.sh 25ns_2748b_2736_2258_2374_288bpi_12inj.in noPlot
./subfolders.sh "collision_64_B_IP$1"

echo " "
echo "Finished with collision 64 BCMS"
echo " "

./runTrainForFillingScheme.sh 8b4e_1972b_1967_1178_1886_224bpi_12inj.in noPlot
./subfolders.sh "collision_64_8_IP$1"

echo " "
echo "Finished with collision 64 8"
echo " "


echo "STABLE B 64"
./updateMadFiles.sh HL_stable_b_64

./runTrainForFillingScheme.sh 25ns_2760b_2748_2494_2572_288bpi_13inj.in noPlot
./subfolders.sh "stable_b_64_std_IP$1"

echo " "
echo "Finished with stable b 64 standard"
echo " "

./runTrainForFillingScheme.sh 25ns_2748b_2736_2258_2374_288bpi_12inj.in noPlot
./subfolders.sh "stable_b_64_B_IP$1"

echo " "
echo "Finished with stable b 64 BCMS"
echo " "

./runTrainForFillingScheme.sh 8b4e_1972b_1967_1178_1886_224bpi_12inj.in noPlot
./subfolders.sh "stable_b_64_8_IP$1"

echo " "
echo "Finished with stable_b 64 8"
echo " "



echo "STABLE E 15"
./updateMadFiles.sh HL_stable_e_15

./runTrainForFillingScheme.sh 25ns_2760b_2748_2494_2572_288bpi_13inj.in noPlot
./subfolders.sh "stable_e_15_std_IP$1"

echo " "
echo "Finished with stable e 15 standard"
echo " "

./runTrainForFillingScheme.sh 25ns_2748b_2736_2258_2374_288bpi_12inj.in noPlot
./subfolders.sh "stable_e_15_B_IP$1"

echo " "
echo "Finished with stable e 15 BCMS"
echo " "

./runTrainForFillingScheme.sh 8b4e_1972b_1967_1178_1886_224bpi_12inj.in noPlot
./subfolders.sh "stable_e_15_8_IP$1"

echo " "
echo "Finished with stable_e 15 8"
echo " "


