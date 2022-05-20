#include "../code/cpl.h"
#include "../gnuplot/gnuplot_i.h"


#define gnuplot_setup(_plt, _name) ({ \
	gnuplot_cmd(_plt, "set term png font \"monospace, 15\" linewidth 2"); \
	gnuplot_cmd(_plt, "set pointsize 2"); \
	gnuplot_cmd(_plt, "set grid"); \
	gnuplot_cmd(_plt, "set output \"" _name ".png\""); \
})

double legendre(int n, double x) {
	if (n == 0) {
		return 1;
	} else if (n == 1) {
		return x;
	} else if (n == 2) {
		return (3*pow(x, 2) - 1) / 2;
	} else if (n == 3) {
		return (5*pow(x, 3) - 3*x) / 2;
	} else if (n == 4) {
		return (35*pow(x, 4) - 30*pow(x, 2) + 3) / 8;
	} else if (n == 5) {
		return (63*pow(x, 5) - 70*pow(x, 3) + 15*x) / 8;
	} else if (n == 6) {
		return (231*pow(x, 6) - 315*pow(x, 4) + 105*pow(x, 2) - 5);
	}

	return 0;
}

double q4_fn(double x, void *params) {
	return 1 / sqrt(1 + pow(x,2));
}

int main() {
	// Question 01

	printf("\n# Question 01\n\n");

	int Nsteps = 200, Nwalks = 500;
	double θ, dist_x, dist_y, d_rms, Σ = 0;
	for (int j = 1; j <= Nwalks; ++j) {
		dist_x = 0;
		dist_y = 0;
		for (int i = 1; i <= Nsteps; ++i) {
			θ = cpl_rand_uniform(0, 2*PI);
			dist_x += cos(θ);
			dist_y += sin(θ);
		}

		Σ += pow(dist_x, 2) + pow(dist_y, 2);
	}

	d_rms = sqrt(Σ / Nwalks);

	printf("RMS distance = % .4f\n", d_rms);
	printf("√N = % .4f\n", sqrt(Nsteps));

	// Question 02
	printf("\n# Question 02\n\n");

	cpl_vector *xdata = cpl_vector_loadtxt("esem4fit.txt", 1, 26, 1);
	cpl_vector *ydata = cpl_vector_loadtxt("esem4fit.txt", 1, 26, 2);

	int params = 5;

	cpl_vector *a_lg = cpl_vector_alloc(params);
	cpl_vector_set_all(a_lg, 1.0);

	cpl_matrix *Cov_lg = cpl_matrix_alloc(params, params);
	cpl_matrix_set_all(Cov_lg, 1.0);

	cpl_stats_linreg(params, legendre, xdata, ydata, NULL, a_lg, Cov_lg);

	printf("Parameters for fitting with polynomials:\n");
	cpl_print(a_lg);

	// Plotting
	int Npts = 51;
	double dx = 2.0 / (Npts - 1), x, y;
	cpl_matrix *q2plot = cpl_matrix_alloc(Npts, 2);
	for (int i = 1; i <= Npts; ++i) {
		x = -1 + dx*(i - 1);
		y = 0;
		for (int n = 0; n < params; ++n)
			y += cpl_get(a_lg, n + 1) * legendre(n, x);

		cpl_set(q2plot, i, 1, x);
		cpl_set(q2plot, i, 2, y);
	}

	cpl_matrix_savetxt("q2_fit.txt", q2plot, "x y");

	gnuplot_ctrl *p2 = gnuplot_init();
	gnuplot_setup(p2, "q2_fit");
	gnuplot_cmd(p2, "set xlabel \"x\"");
	gnuplot_cmd(p2, "set ylabel \"y\"");
	gnuplot_cmd(p2, "plot \"esem4fit.txt\" using 1:2 with points title \"Data\", \"q2_fit.txt\" using 1:2 with lines title \"Fit\"");
	gnuplot_close(p2);

	// Question 03
	printf("\n# Question 03\n\n");
	double κ = 1;
	int Nx = 21, Nt = 5001;
	double xmin = 0,
		   xmax = 2;
	dx = (xmax - xmin) / (Nx - 1);
	double tmin = 0,
		   tmax = 4,
		   dt = (tmax - tmin) / (Nt - 1);

	cpl_matrix *U = cpl_matrix_alloc(Nx, Nt);

	// Boundary conditions
	cpl_matrix_set_all(U, 0);
	for (int k = 1; k <= Nx; ++k) {
		cpl_set(U, k, 1, 20*fabs(sin(PI*(k - 1)*dx)));
	}

	// Solution by explicit (forward) method
	cpl_diffeq_heat_forward(U, κ, dx, dt);

	// Plotting shenanigans
	int Tstep[7] = { 0, 10, 20, 50, 100, 200, 500 };
	cpl_matrix *Uplot = cpl_matrix_alloc(Nx, 8);
	cpl_vector *v = cpl_vector_alloc(Nx);
	for (int i = 1; i <= Nx; ++i) {
		cpl_set(Uplot, i, 1, (i - 1)*dx);
		for (int j = 0; j < 7; ++j) {
			cpl_matrix_get_col(U, Tstep[j] + 1, v);
			cpl_matrix_set_col(Uplot, j + 2, v);
		}
	}

	cpl_matrix_savetxt("q3_profiles.txt", Uplot, "x \"T = 0C\" \"T = 10C\" \"T = 20C\" \"T = 50C\" \"T = 100C\" \"T = 200C\" \"T = 500C\"");

	gnuplot_ctrl *p3 = gnuplot_init();
	gnuplot_setup(p3, "q3_profiles");
	gnuplot_cmd(p3, "set xlabel \"Position, x\"");
	gnuplot_cmd(p3, "set xrange [-0.1:2.1]");
	gnuplot_cmd(p3, "set ylabel \"Temperature (°C)\"");
	gnuplot_cmd(p3, "set yrange [-1:25]");
	gnuplot_cmd(p3, "plot for [col=2:8] \"q3_profiles.txt\" using 1:col with lines title columnheader");
	gnuplot_close(p3);

	// Question 04

	printf("\n# Question 04\n\n");
	
	printf("For this problem we have to do the integral: ∫dx/√(1 + x²) from -1 to +1\n");
	printf("The exact answer is ln(1 + √2) - ln(-1 + √2) = % .10g\n", log(1+sqrt(2)) - log(-1+sqrt(2)));

	cpl_function F = { q4_fn, NULL };
	printf("Gaussian quadrature with N = 4: % .10g\n", cpl_integrate_quadgl(F, -1, 1, 4));
	printf("Gaussian quadrature with N = 5: % .10g\n", cpl_integrate_quadgl(F, -1, 1, 5));
	printf("Gaussian quadrature with N = 6: % .10g\n", cpl_integrate_quadgl(F, -1, 1, 6));

	return 0;
}
