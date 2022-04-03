set term png font "monospace,14"

set title "Convergence Plot for Question 2"
set xlabel "Iterations"
set ylabel "Residue"
set yrange [0:5]

plot 'assignment01.dat' using 1:2 title "Jacobi", \
     'assignment01.dat' using 1:3 title "Gauss-Seidel"
