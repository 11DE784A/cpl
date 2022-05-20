
set term png font "monospace, 15" linewidth 2
set pointsize 2
set grid

set view map

set xrange [-0.5:5000.5]
set xlabel "Time step"
set yrange [-0.5:20.5]
set ylabel "Position grid point"

plot 'q3_contour.txt' matrix with image
