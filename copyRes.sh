!#/bin/bash
calcType=.emitt$2
usedFillScheme=`basename $1`
source=RESULTS/$usedFillScheme
dest=RESULTS_COMPARE_NOMINAL/$usedFillScheme$calcType

echo $source
echo $dest

mkdir $dest

cp -R $source/* $dest/