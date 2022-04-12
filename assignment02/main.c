#include "../code/cpl.h"
#include "../code/cpl_bench.h"

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
	int max_iters = 1e6;

	/*** Question 2 ***/
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
