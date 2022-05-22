#include "../code/cpl.h"
#include "../gnuplot/gnuplot_i.h"

#define gnuplot_setup(_plt, _name) ({ \
	gnuplot_cmd(_plt, "set term png font \"monospace, 15\" linewidth 2"); \
	gnuplot_cmd(_plt, "set pointsize 2"); \
	gnuplot_cmd(_plt, "set grid"); \
	gnuplot_cmd(_plt, "set output \"" _name ".png\""); \
})

double well_energy(int n) {
	return pow(n*PI, 2) / 2;
}

cpl_vector *fn_well_sch1(cpl_vector *y, double x) {
	double y1 = cpl_get(y, 1),
		   y2 = cpl_get(y, 2);

	cpl_set(y, 1, y2);
	cpl_set(y, 2, -2*well_energy(1)*y1);

	return y;
}

cpl_vector *fn_well_sch2(cpl_vector *y, double x) {
	double y1 = cpl_get(y, 1),
		   y2 = cpl_get(y, 2);

	cpl_set(y, 1, y2);
	cpl_set(y, 2, -2*well_energy(2)*y1);

	return y;
}

double rand_exp() {
	double u = cpl_rand_uniform(0, 1),
		   α = 1 - exp(-1);
	return -log(1 - u*α);
}

double integrate_gaussian(int N) {
	double Σf = 0, 
		   α = 1 - exp(-1),
		   x;
	int counter = 0;
	while (counter++ < N) {
		x = rand_exp();
		Σf += exp(-x*x) * α*exp(x);
	}

	return Σf / N;
}

double gaussian(double x) {
	return exp(-x*x);
}

int main() {
	/* Question 01 */
	printf("\nQuestion 01\n");
	int N = 1e6;

	double Σnaive = cpl_mcint1d_naive(gaussian, 0, 1, N);
	printf("Monte-Carlo integration without importance sampling ∫exp(-x^2) = % .9g\n", Σnaive);

	double Σimp = integrate_gaussian(N);
	printf("Monte-Carlo integration with importance sampling ∫exp(-x^2) = % .9g\n", Σimp);

	double Σwolf = 0.7468241328124270253994674361318530053544996868126063290276544989;
	printf("Value from Wolfram-Alpha ∫exp(-x^2) = % .9g\n", Σwolf);


	/* Question 02 */
	printf("\nQuestion 02\n");

	int Na = 51;
	double amin = 0, amax = 1;
	double da = (amax - amin) / (Na - 1);
	cpl_vector *a = cpl_vector_alloc(Na);
	for (int i = 1; i <= Na; ++i)
		cpl_set(a, i, amin + (i - 1)*da);

	cpl_matrix *ψ1 = cpl_matrix_alloc(Na, 2);
	cpl_matrix_set_all(ψ1, 0);
	cpl_set(ψ1, 1, 1, 0);
	cpl_set(ψ1, 1, 2, 1);

	cpl_diffeq_rk4(fn_well_sch1, ψ1, a, da);

	cpl_matrix *ψ2 = cpl_matrix_alloc(Na, 2);
	cpl_matrix_set_all(ψ2, 0);
	cpl_set(ψ2, 1, 1, 0);
	cpl_set(ψ2, 1, 2, 1);

	cpl_diffeq_rk4(fn_well_sch2, ψ2, a, da);

	cpl_matrix *wfplot = cpl_matrix_alloc(Na, 3);
	for (int i = 1; i <= Na; ++i) {
		cpl_set(wfplot, i, 1, cpl_get(a, i));
		cpl_set(wfplot, i, 2, cpl_get(ψ1, i, 1));
		cpl_set(wfplot, i, 3, cpl_get(ψ2, i, 1));
	}
	cpl_matrix_savetxt("q2_soln.txt", wfplot, "x \"n = 1\" \"n = 2\"");

	// Plotting
	gnuplot_ctrl *p2 = gnuplot_init();
	gnuplot_setup(p2, "q2_soln");
	gnuplot_cmd(p2, "set xlabel \"x\"");
	gnuplot_cmd(p2, "set xrange [-0.1:1.1]");
	gnuplot_cmd(p2, "set ylabel \"ψ(x)\"");
	gnuplot_cmd(p2, "set yrange [-0.2:0.4]");
	gnuplot_cmd(p2, "plot for [col=2:3] \"q2_soln.txt\" using 1:col with lines title columnheader");
	gnuplot_close(p2);

	printf("Solution for the first two wavefunctions of the infinite square well by shooting method is attached in `q2_soln.png`\n");

	/* Question 03 */
	printf("\nQuestion 03\n");

	// Number of grid points
	int Nx = 21,
		Ny = 21;

	// Region of interest
	double xmin = 0,
		   xmax = 1,
		   ymin = 0,
		   ymax = 1;

	// Grid step size
	double dx = (xmax - xmin) / (Nx - 1);
	double dy = (ymax - ymin) / (Ny - 1);

	// Setting up matrix
	cpl_matrix *ϕ = cpl_matrix_alloc(Nx, Ny);

	// Boundary conditions
	cpl_matrix_set_all(ϕ, 0);
	cpl_matrix_set_row(ϕ, 1,  1);	// ϕ(x, y = 0) = 1
	cpl_matrix_set_row(ϕ, Ny, 0);	// ϕ(x, y = 1) = 0
	cpl_matrix_set_col(ϕ, 1,  0);	// ϕ(x = 0, y) = 0
	cpl_matrix_set_col(ϕ, Nx, 0);	// ϕ(x = 1, y) = 0

	// Solution
	cpl_diffeq_laplace_dirichlet(ϕ, dx, dy);

	printf("Solution matrix:\n");
	cpl_print(ϕ);

	// Plotting
	cpl_matrix_savetxt("q3_laplace.txt", ϕ, NULL);

	gnuplot_ctrl *p3 = gnuplot_init();
	gnuplot_setup(p3, "q3_contour");
	gnuplot_cmd(p3, "set view map");
	gnuplot_cmd(p3, "set cblabel \"Potential (volts)\"");
	gnuplot_cmd(p3, "set xlabel \"x steps\"");
	gnuplot_cmd(p3, "set xrange [-0.5:20.5]");
	gnuplot_cmd(p3, "set ylabel \"y steps\"");
	gnuplot_cmd(p3, "set ylabel \"y steps\"");
	gnuplot_cmd(p3, "set yrange [-0.5:20.5]");
	gnuplot_cmd(p3, "plot \"q3_laplace.txt\" matrix with image");
	gnuplot_close(p3);

	printf("Contour plot of the solution is also attached in `q3_contour.png`\n");

	cpl_free(ϕ);
	return 0;
}

