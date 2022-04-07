set term png font "monospace,14"

set title "Convergence Plot for Question 3"
set xlabel "Iterations"
set ylabel "Residue"
set xrange [0:50]
set xtics 5

plot 'assignment01.dat' using 1:4 title "Column 1" with line, \
     'assignment01.dat' using 1:5 title "Column 2" with line, \
     'assignment01.dat' using 1:6 title "Column 3" with line
