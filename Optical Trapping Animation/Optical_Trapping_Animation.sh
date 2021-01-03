#!/bin/bash

let TimeShots=100            # Number of time shots during the evolution
let TimeInterval=10          # Time interval between consecutives time shots [µs]

# compile programs
g++ -o MolDist.out MolDist.cpp Functions.cpp -O2 -lm -lgsl -lgslcblas -Wall
gcc -o Trajectories.out Trajectories.cpp Functions.cpp -O2 -lm -lgsl -lgslcblas -Wall

#MolDist.out > name t x y z vx vy vz Exc color
./MolDist.out > MolDist.dat
FileAus=MolDist.dat
ParameterAus=Parameters.hpp


# take the values from the file Parameters.hpp [line,column]
N=$(awk -v line="6" -v col="3" 'NR == line { print $col }' <"$ParameterAus")
V_IN=$(awk -v line="14" -v col="3" 'NR == line { print $col }' <"$ParameterAus")
TEMP=$(awk -v line="15" -v col="3" 'NR == line { print $col }' <"$ParameterAus")

IPG_X=$(awk -v line="19" -v col="3" 'NR == line { print $col }' <"$ParameterAus")
IPG_Y=$(awk -v line="20" -v col="3" 'NR == line { print $col }' <"$ParameterAus")
IPG_Z=$(awk -v line="21" -v col="3" 'NR == line { print $col }' <"$ParameterAus")
IPG_L=$(awk -v line="23" -v col="3" 'NR == line { print $col }' <"$ParameterAus")
IPG_W=$(awk -v line="24" -v col="3" 'NR == line { print $col }' <"$ParameterAus")

RDL_X=$(awk -v line="28" -v col="3" 'NR == line { print $col }' <"$ParameterAus")
RDL_Y=$(awk -v line="29" -v col="3" 'NR == line { print $col }' <"$ParameterAus")
RDL_Z=$(awk -v line="30" -v col="3" 'NR == line { print $col }' <"$ParameterAus")
RDL_L=$(awk -v line="31" -v col="3" 'NR == line { print $col }' <"$ParameterAus")
RDL_W=$(awk -v line="32" -v col="3" 'NR == line { print $col }' <"$ParameterAus")

for a in $(seq 1 $TimeShots);
do
echo "Cycle: $i"
let i+=1
let aus1=$a*$TimeInterval

# final time as input for Trajectories.out
./Trajectories.out $aus1 < $FileAus > Frame-$a.dat
FileAus=Frame-$a.dat


# take the values from the file Frame-$a.dat [line,column]
a3pi1_N=$(awk -v line="$((N+1))" -v col="10" 'NR == line { print $col }' <"$FileAus")
Lost_N=$(awk -v line="$((N+1))" -v col="11" 'NR == line { print $col }' <"$FileAus")
GS_N=$(awk -v line="$((N+1))" -v col="12" 'NR == line { print $col }' <"$FileAus")
GS_TRAPPED_N=$(awk -v line="$((N+1))" -v col="13" 'NR == line { print $col }' <"$FileAus")

gnuplot <<- EOF

load 'Setting_Animation.plot'
set out 'Frame_Animation_${a}.gif'

set label "Time = $aus1 µs" at -1250,2700,80
set label "  N = $((a3pi1_N+Lost_N+GS_N)) CO molecules" at -1250,2500,60
set label "  v = $V_IN µm / µs     T = $TEMP mK" at -1250,2500,40

set label "# a3pi1    = $a3pi1_N" at -100,2500,130
set label "# lost       = $Lost_N" at -100,2500,110
set label "# X1S+    = $GS_N" at -100,2500,90
set label "# Trapped = $GS_TRAPPED_N" at -100,2500,70

splot "Frame-$a.dat" u 3:4:5:14 lc rgb variable pt 7 ps 2 t "",\
$IPG_X+$IPG_W*sqrt(1+((u*$IPG_L)/(pi*$IPG_W**2))**2)*sin(v) , $IPG_Y+u , $IPG_Z+$IPG_W*sqrt(1+((u*$IPG_L)/(pi*$IPG_W**2))**2)*cos(v) t "" w l lw 1 lc rgb "#FF0000"
#$RDL_X+$RDL_W*sqrt(1+((u*$RDL_L)/(pi*$RDL_W**2))**2)*sin(v) , $RDL_Y+u , $RDL_Z+$RDL_W*sqrt(1+((u*$RDL_L)/(pi*$RDL_W**2))**2)*cos(v) t "" w l lw 1 lc rgb "#00BB00"

set output

EOF

done

convert -delay 10 Frame_Animation_?.gif Frame_Animation_??.gif Frame_Animation_???.gif -loop 0 3D_Animation.gif

rm Frame_*.gif
rm *.dat
rm *.out
