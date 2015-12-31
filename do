#!/bin/bash

lng=12
rad=10
npes=$(($rad * $lng))
days=200

./changeDimension.sh $rad $lng

#echo $lng $rad $npes

make all

if [ $? -eq 0 ] 
  then 

  echo "Model Compiled"
  date
  time mpirun -n $npes ./torus > runlog
   
  if [ $? -ge 0 ] 
    then 
    echo "Run Completed Successfully"

    ./moveData.sh $rad $lng

    ./transport.sh $days

    cd plots

      ./azplots $days MIXR 
      mv animated.mpeg ../azdens.mpeg

      ./3Dplots $days MIXR sp
      mv animated.mpeg ../3dspplot.mpeg

      ./3Dplots $days MIXR s3p
      mv animated.mpeg ../3ds3pplot.mpeg

      ./3Dplots $days VSUB .
      mv animated.mpeg ../3dVelPlot.mpeg

      ./3Dplots $days TEMP sp
      mv animated.mpeg ../3dTempPlot.mpeg

      ./radialplots $days DENS
      mv animated.mpeg ../raddens.mpeg

      ./radialplots $days NL2_
      mv animated.mpeg ../radnl2.mpeg

      ./radialplots $days MIXR
      mv animated.mpeg ../radmix.mpeg

      ./radialplots $days TEMP
      mv animated.mpeg ../radtemp.mpeg

#      ./azplots $days TEMP
#      mv animated.mpeg ../aztemp.mpeg

      ./azplots $days DENS
      mv animated.mpeg ../azdens.mpeg

#      ./miscPlot $days MOUT
#      mv misc.mpeg ../misc.mpeg

    cd ..
 
    python chi.py

    ./plots/mixRatio.sh $days

#    vlc dens.mpeg

  fi
  make clean
  echo "Run Complete"
fi

#gifview -a dens.gif &
#mplayer -loop -0 dens.avi &
#vlc intensity.avi &
#mplayer dens.avi &
#display peaks.jpeg & 
#display peakRatio.jpeg  
#display overlay.jpeg



