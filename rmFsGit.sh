!#/bin/bash

for file in ls FillingSchemes/*; do
 base=`basename $file`
 git rm $base
done