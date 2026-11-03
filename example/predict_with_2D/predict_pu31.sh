#!/bin/bash
BIN=../../bin
for name in pu31 
do 
  mkdir ${name}_predict && cd ${name}_predict
  ../${BIN}/Predict_RNA_structure -i ../${name}.2d -name ${name} -exclude ../${name}_exclude.dat -exe 1
done
