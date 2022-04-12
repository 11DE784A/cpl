set term png font "monospace,14"

set title "Data Fitting (Question 1)"
set xlabel "x"
set ylabel "y"

plot 'assign2fit.txt' using 1:2 title "Data from experiment", \
	 'assignment02.dat' using 1:2 title "Fit w/ polynomials" with line, \
	 'assignment02.dat' using 1:3 title "Fit w/ Chebyshev functions" with line
