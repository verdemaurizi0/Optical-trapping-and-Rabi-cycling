#!/bin/bash

rm *.dat
rm *.out

let vel_old=0
let vmin=1
let vmax=20

# Set the desired values of microtrap speed {vmin..vmax..dv}
for vel_new in {1..20..1}

do

  # Modify a desired parameter within Parameters.hpp
  awk -v vnew=$vel_new -v vold=$vel_old 'BEGIN {    old = "Microtrap_V "vold"."
                                                  new = "Microtrap_V "vnew"." }
  s=index($0,old) { $0 = substr($0,1,s-1) new substr($0,s+length(old)) } { print }' Parameters.hpp > temp && mv temp Parameters.hpp

  echo "Microtrap's velocity: $vel_new"

  # compile programs
  g++ -o MolDist.out MolDist.cpp Functions.cpp -O2 -lm -lgsl -lgslcblas -Wall
  g++ -o Prediction_Microtrap_Speed.out Prediction_Microtrap_Speed.cpp Functions.cpp -O2 -lm -lgsl -lgslcblas

  ./MolDist.out > MolDist.dat
  ./Prediction_Microtrap_Speed.out < MolDist.dat >> Results.dat

  vel_old=$vel_new

done

ParameterAus=Parameters.hpp
TEMP=$(awk -v line="15" -v col="3" 'NR == line { print $col }' <"$ParameterAus")
IPG_W=$(awk -v line="24" -v col="3" 'NR == line { print $col }' <"$ParameterAus")


gnuplot <<- EOF

load 'Microtrap_Speed.plot'
set out 'Microtrap_Speed.png'
set label "Microtrap temp = $TEMP mK" at $vmax/2,0.0055
set label "Optical trap w = $IPG_W µm" at $vmax/2,0.005
set xrange [$vmin:$vmax]
set xtics $vmax/4
plot 'Results.dat' u 1:2 ls 2 w l t ""
set output

EOF
