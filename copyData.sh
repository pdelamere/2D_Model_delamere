#!/bin/bash

cp DENS*.dat plots/.
cp MIXR*.dat plots/.
cp TEMP*.dat plots/.
#cp INTS*.dat plots/.  
cp NL2_*.dat plots/.
cp LOAD*.dat plots/.
#cp intensity*.dat plots/.  

cd plots

  python radData.py $1 $2

  cp DENS*.dat data/.
  cp MIXR*.dat data/.
  cp TEMP*.dat data/.
#  cp INTS*.dat data/.
  cp NL2_*.dat data/.
  cp LOAD*.dat data/.
#  rm intensity*.dat

  ./organize.sh

cd ..
