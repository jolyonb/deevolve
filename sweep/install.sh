#!/bin/sh

cd src/gini
make clean
make 
cd ../colchi2
make clean
make


cd ../../


cp ../intDE/SCPUnion2.1_mu_vs_z.txt .
