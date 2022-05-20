
set terminal pngcairo

set xrange [0:20]
set yrange [0:20]

set view map
plot 'heatmap.txt' matrix with image
