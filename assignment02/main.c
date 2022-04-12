#include "../code/cpl.h"
#include "../code/cpl_bench.h"

double modified_chebyshev(int n, double x) {
	if (n == 0) {
		return 1;
	} else if (n == 1) {
		return 2*x - 1;
	} else if (n == 2) {
		return 8*pow(x, 2) - 8*x + 1;
	} else if (n == 3) {
		return 32*pow(x, 3) - 48*pow(x, 2) + 18*x -1;
	}

	return 0;
}

double condition_number(cpl_matrix *M) {
	int dim = cpl_matrix_rows(M);

	cpl_matrix *id = cpl_matrix_id(dim);

	cpl_matrix *Minv = cpl_matrix_alloc(dim, dim);
	cpl_matrix_set_all(Minv, 1.0);

	cpl_linalg_seidel(M, id, Minv, NULL);
	double cond = cpl_matrix_frobnorm(M) * cpl_matrix_frobnorm(Minv);

	cpl_free(id);
	cpl_free(Minv);

	return cond;
}

double estimate_pi(int max_iters) {
	double x, y;
	int points_inside_circle = 0;

	for (int i = 1; i <= max_iters; ++i) {
		x = cpl_rand_uniform(0, 1);
		y = cpl_rand_uniform(0, 1);

		if (x*x + y*y < 1) 
			++points_inside_circle;
	}

	return 4 * (double) points_inside_circle / max_iters;
}

double f(double x) { 
	return sqrt(1 - x*x); 
}

void check_lcg_pi(int a, int m, int s, int max_iters) {
	cpl_lcg_modulus(m);
	cpl_lcg_multiplier(a);
	cpl_lcg_seed(s);

	printf("With a = %3d, m = %5d, and num_points = %g: π = %.8f\n", 
			a, m, (float) max_iters, estimate_pi(max_iters));
}

void check_lcg_int(int a, int m, int s, int max_iters) {
	cpl_lcg_modulus(m);
	cpl_lcg_multiplier(a);
	cpl_lcg_seed(s);

	printf("With a = %3d, m = %5d, and num_points = %g: 4⋅∫√(1 - x²)dx = %.8f\n", 
			a, m, (float) max_iters, 4.0 * cpl_mcint1d_naive(f, 0, 1, max_iters));
}

double steinmetz_volume(double r, int N) {
	/* The Steinmetz solid is defined by
	 * -r ≤ x ≤ r
	 *  x² + y² ≤ r²
	 *  x² + z² ≤ r²
	 * We can estimate the volume by throwing points */

	double x, y, z;
	int points_inside_solid = 0;
	for (int i = 0; i < N; ++i) {
		x = cpl_rand_uniform(-1, 1);
		y = cpl_rand_uniform(-1, 1);
		z = cpl_rand_uniform(-1, 1);

		if (x*x + y*y <= 1.0 && x*x + z*z <= 1.0)
			++points_inside_solid;
	}

	return ((double) points_inside_solid / N) * 8 * pow(r, 3);
}

int main() {
	/*** Question 1 ***/
	cpl_vector *xdata = cpl_vector_loadtxt("assign2fit.txt", 2, 22, 1);
	cpl_vector *ydata = cpl_vector_loadtxt("assign2fit.txt", 2, 22, 2);

	int params = 4;

	printf("\n# Question 1\n");

	/* For cubic polynomial */
	cpl_vector *a_poly = cpl_vector_alloc(params);
	cpl_vector_set_all(a_poly, 1.0);

	cpl_matrix *Cov_poly = cpl_matrix_alloc(params, params);
	cpl_matrix_set_all(Cov_poly, 1.0);

	cpl_stats_linreg(params, polynomial, xdata, ydata, NULL, a_poly, Cov_poly);

	printf("\nParameters for fitting with polynomials:\n");
	cpl_print(a_poly);

	printf("Condition number: %.4g\n", condition_number(Cov_poly));

	/* For modified Chebyshev functions */
	cpl_vector *a_cheby = cpl_vector_alloc(params);
	cpl_vector_set_all(a_cheby, 1.0);

	cpl_matrix *Cov_cheby = cpl_matrix_alloc(params, params);
	cpl_matrix_set_all(Cov_cheby, 1.0);

	cpl_stats_linreg(params, modified_chebyshev, xdata, ydata, NULL, a_cheby, Cov_cheby);

	printf("\nParameters for fitting with modified Chebyshev functions:\n");
	cpl_print(a_cheby);

	printf("Condition number: %.4g\n", condition_number(Cov_cheby));

	printf("\nLower condition number for Chebyshev functions implies greater numerical stability and roughly upto 4 digits of additional numerical accuracy that is lost when fitting with regular polynomials.\n");

	cpl_free(a_cheby);
	cpl_free(Cov_cheby);
	cpl_free(a_poly);
	cpl_free(Cov_poly);
	cpl_free(xdata);
	cpl_free(ydata);

	/*** Question 2 ***/
	int max_iters = 1e6;

	printf("\n# Question 2\n");

	printf("\n## Estimating π by throwing points\n");
	check_lcg_pi(65,  1021,  1, max_iters);
	check_lcg_pi(572, 16381, 1, max_iters);

	printf("\n## Estimating π by Monte-Carlo evaluation of ∫√(1 - x²)dx\n");
	check_lcg_int(65,  1021,  1, max_iters);
	check_lcg_int(572, 16381, 1, max_iters);

	printf("\nWe see that the choice a = 527 and m = 16381 gives better estimate of π in both cases.\n");

	/*** Question 3 ***/
	printf("\n# Question 3\n");

	cpl_lcg_modulus(16381); cpl_lcg_multiplier(572);
	printf("\nVolume of Steinmetz solid computed by Monte-Carlo: %.8f\n", 
		   steinmetz_volume(1, max_iters));
	
	return 0;
}
