# Setting

# General settings
set terminal gif font 'arial' 20 size 800,600 

# Line styles
# lt = linetype to 0/1 (dashed line/continuous line)
# lc = linecolor to a gray defined by that number
set style line 101 lt 1 lc rgb "#000000"                        # for axis
set style line 102 lt 1 lw 4 pt 7 ps 0.5         # for points #FFA500

# back means the border should be behind anything else drawn
set border 31 back ls 101

# Specify rotations about the x-axis, the z-axis, overall scaling and scaling of the z-axis 
# Normal view: 60,60 / top view: 0,90 / lateral view: 90,0
set view 60, 30, 1, 1

# The ‘ticslevel’ sets the heightfrom the bottom of the z-scale above the xy-plane
# The xtics, ytics and ztics defines the lower and the highest numbers and the space between cosecutive labels [lower, space, highest] 
set ticslevel 0
set xtics -500,250,0
set ytics -2500,2500,2500
set ztics (-50,0,50)

set xrange [-700:250]
set yrange [-3000:3000]
set zrange [-80:80]

set xlabel "Z[µm]" offset -1,-1  # label x-axis with offsets
set ylabel "Y[µm]" offset 2,-1  # label y-axis with offsets
set label "X[µm]" at -750,-3200,100 

set parametric
set isosamples 35,35 # minimum values: 2,3

set dummy u,v
set urange [-3000:3000]
set vrange [0:2*pi]


set key top
