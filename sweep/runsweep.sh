#!/bin/sh

cp src/gini/gini .
cp src/colchi2/colchi2 .
cp src/evol.py .
cp src/main intDE

./gini
python evol.py
./colchi2

rm gini
rm evol.py
rm colchi2
rm intDE