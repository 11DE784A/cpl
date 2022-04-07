set term png font "monospace,14"

set title "Convergence Plot for Question 3"
set xlabel "Iterations"
set ylabel "Residue"
set yrange [0:0.5]
set xrange [0:4]
set xtics 1

plot 'assignment01.dat' using 1:4 title "Column 1", \
     'assignment01.dat' using 1:5 title "Column 2", \
     'assignment01.dat' using 1:6 title "Column 3"
