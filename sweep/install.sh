#!/bin/sh

cd src/gini
make clean
make 
cd ../colchi2
make clean
make


cd ../../

cp src/gini/gini .
cp src/colchi2/colchi2 .
cp src/evol .
cp src/main intDE
cp ../intDE/SCPUnion2.1_mu_vs_z.txt .
