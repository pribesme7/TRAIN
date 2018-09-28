#!/bin/bash

results="RESULTS/$1"
echo "Remove old results for '$results'"
rm -r "$results"

mkdir -p "$results"
echo 'MOVE the IP result files...'
mv -t "$results" *_b.*
mv -t "$results" *_f.*
mv -t "$results" fort.*
mv -t "$results" *.IP*
echo 'MOVE general TRAIN output files...'
mv -t "$results" coll.count
mv -t "$results" train.list
mv -t "$results" train.orb
echo "!! DONE for '$results'"
