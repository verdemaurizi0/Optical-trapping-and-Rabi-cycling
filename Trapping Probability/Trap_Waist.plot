# Setting

# General settings
set terminal png font 'arial' 20 size 800,600

# Line style for axes
# Define a line style (we're calling it 80) and set 
# lt = linetype to 0/1 (dashed line/continuous line)
# lc = linecolor to a gray defined by that number
set style line 80 lt 1 lc rgb "#808080"

# Set the border using the linestyle 80 that we defined
# 3 = 1 + 2 (1 = plot the bottom line and 2 = plot the left line)
# back means the border should be behind anything else drawn
set border 15 back ls 80

# Line style for grid
# Define a new linestyle (81)
# linetype = 0 (dashed line)
# linecolor = gray
# lw = lineweight, make it half as wide as the axes lines
set style line 81 lt 0 lc rgb "#808080" lw 0.5


# Create some linestyles for our data
# pt = point type (triangles, circles, squares, etc.)
# ps = point size
set style line 1 lt 1 lc rgb "#5060D0" lw 4 pt 9 ps 3
set style line 2 lt 1 lc rgb "#FFA500" lw 5 pt 7 ps 5

# Put X and Y labels
set xlabel "Trap waist [ µm ]"  offset 0,0
set ylabel "Trapping probability" offset 1,0


# Give the plot a title
set title ""

