#!/bin/bash


echo 'REMOVE the IP result files...'
rm  *_b.*
rm  *_f.*
rm  fort.*
rm  *.IP*
echo 'REMOVE general TRAIN output files...'
rm  coll.count
rm  train.list
rm  train.orb
echo '!! DONE for cleaning'
