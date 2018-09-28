!#/bin/bash
for f in train*.in; do
  echo "File -> $f ->" $f."emittance"
  python convertFillingSchemes.py $f $f."emittance"
done
