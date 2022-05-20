
set terminal pngcairo

xmin = 0.0
xmax = 0.1
xnumtics = 40
# set xtics ()
# set for [i=0:xnumtics] xtics add (gprintf("%g", xmin + (1.0*i/xnumtics) * (xmax-xmin)) i)

ymin = 0.0
ymax = 1.0
ynumtics = 50
# set ytics ()
# set for [i=0:ynumtics] ytics add (gprintf("%g", ymin + (1.0*i/ynumtics) * (ymax-ymin)) i)

set view map
plot 'heat.txt' matrix with image
