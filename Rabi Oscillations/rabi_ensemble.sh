#!/bin/bash

rm *.dat
rm MolGen
rm rabi

g++ -o MolGen MolGen.cpp -L/usr/lib/x86_64-linux-gnu -lgsl -lgslcblas -lm -lconfig
g++ -o rabi rabi.cpp feynman.c -I/usr/include -L/usr/lib/x86_64-linux-gnu -lgsl -lgslcblas -lm -lpthread -lconfig

./MolGen > MolGen.dat
./rabi < MolGen.dat > out.dat

gnuplot <<- EOF

load 'script.plot'
set out 'rabi_ensemble.png'
plot 'out.dat' u 1:2 ls 2 w l t "v1"
set output

EOF
