#!/bin/bash

days=$1

cd plots/

python transport.py $days

echo "set terminal jpeg" > gplot
echo "set output 'transport.jpeg'" >> gplot
echo "set key off" >> gplot
echo "set title 'Radial Transport'" >> gplot
echo "set ylabel 'Integrated Transport Time (days)'" >> gplot
echo "set xlabel 'Radial Distance (RJ)'" >> gplot
echo "set grid ytics" >> gplot
echo "set xrange [6:9]" >> gplot
echo "set grid xtics" >> gplot
echo "plot 'transport.dat' with lines" >> gplot

gnuplot gplot

#rm gplot

mv transport.jpeg ../transport.jpeg

cd ..
