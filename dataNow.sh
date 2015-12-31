#!/bin/bash

./copyData.sh $1 $2

./plots/mixRatio.sh $3

cd plots
python radData.py $1 $2
cd ..

./plots/radialplots $3 DENS
mv animated.mpeg raddens.mpeg

./plots/azplots $3 DENS
mv animated.mpeg azdens.mpeg
