let TimeShots=100            # Number of time shots during the evolution
let TimeInterval=26          # Time interval between consecutives time shots [µs]

# compile programs
g++ -o MolDist.out MolDist.cpp Functions.cpp -O2 -lm -lgsl -lgslcblas -Wall
g++ -o Decay.out Decay.cpp Functions.cpp -O2 -lm -lgsl -lgslcblas -Wall

#MolDist.out > name t x y z vx vy vz Exc color
./MolDist.out > MolDist.dat
FileAus=MolDist.dat
ParameterAus=Parameters.hpp


# take the values from the file Parameters.h [line,column]
a3pi1_T=$(awk -v line="19" -v col="3" 'NR == line { print $col }' <"$ParameterAus")

for a in $(seq 1 $TimeShots);
do
echo "Cycle: $i"
let i+=1
let aus1=$a*$TimeInterval

# final time as input for Trajectories.out
./Decay.out $aus1 < $FileAus > Frame-$a.dat
FileAus=Frame-$a.dat

gnuplot <<- EOF

load 'Setting_Decay.plot'
set out 'Frame_Decay_${a}.gif'
set xrange [0:$TimeShots*$TimeInterval]
set xtics $TimeShots*$TimeInterval/4
plot '< tail -n 1 Frame-$a.dat' u 2:10 ls 2 t "Simulations      ",\
exp(-(x/$a3pi1_T)) w l ls 1 t "Analytical formula "
set output

EOF

done

convert -delay 10 Frame_Decay_?.gif Frame_Decay_??.gif -loop 0 Decay.gif

rm Frame_*.gif
rm *.dat
rm *.out
