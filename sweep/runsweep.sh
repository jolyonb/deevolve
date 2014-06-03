#!/bin/sh

# go get the executables

cp src/gini/gini .
cp src/colchi2/colchi2 .
cp src/evol.py .
cp src/main intDE

# execute in the correct order

./gini
python evol.py
./colchi2

# delete the executables (to keep everywhere clean)

rm gini
rm evol.py
rm colchi2
rm intDE