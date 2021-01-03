#!/bin/bash

rm *.dat
rm *.out

let w_old=100
let wmin=10
let wmax=100

# Set the desired values of trap waist {IPG_W_min..IPG_W_max..dIPG_W}
for w_new in {10..100..5}

do

  # Modify a desired parameter within Parameters.hpp
  awk -v wnew=$w_new -v wold=$w_old 'BEGIN {    old = "IPG_W "wold"."
                                                  new = "IPG_W "wnew"." }
  s=index($0,old) { $0 = substr($0,1,s-1) new substr($0,s+length(old)) } { print }' Parameters.hpp > temp && mv temp Parameters.hpp

  awk -v wnew=$w_new -v wold=$w_old 'BEGIN {    old = "RDL_W "wold"."
                                                  new = "RDL_W "wnew"." }
  s=index($0,old) { $0 = substr($0,1,s-1) new substr($0,s+length(old)) } { print }' Parameters.hpp > temp && mv temp Parameters.hpp
  
  echo "Optical trap waist: $w_new"

  # compile programs
  g++ -o MolDist.out MolDist.cpp Functions.cpp -O2 -lm -lgsl -lgslcblas -Wall
  g++ -o Prediction_Trap_Waist.out Prediction_Trap_Waist.cpp Functions.cpp -O2 -lm -lgsl -lgslcblas

  ./MolDist.out > MolDist.dat
  ./Prediction_Trap_Waist.out < MolDist.dat >> Results.dat

  w_old=$w_new

done

ParameterAus=Parameters.hpp
VEL=$(awk -v line="14" -v col="3" 'NR == line { print $col }' <"$ParameterAus")
TEMP=$(awk -v line="15" -v col="3" 'NR == line { print $col }' <"$ParameterAus")


gnuplot <<- EOF

load 'Trap_Waist.plot'
set out 'Trap_Waist.png'
set label "Microtrap temp = $TEMP mK" at $wmax/2,0.0055
set label "Microtrap vel = $VEL µm/µs" at $wmax/2,0.005
set xrange [$wmin:$wmax]
set xtics $wmax/10
plot 'Results.dat' u 1:2 ls 2 w l t ""
set output

EOF
