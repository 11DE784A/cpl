#include "cpl_defines.h"
#include "cpl_includes.h"

#include "cpl_arrays.h"
#include "cpl_random.h"

static int SEED = 1;
static int MODULUS = 16381;
static int MULTIPLIER = 572;

void cpl_lcg_seed(int s) {
	SEED = s;
}

void cpl_lcg_modulus(int m) {
	MODULUS = m;
}

void cpl_lcg_multiplier(int a) {
	MULTIPLIER = a;
}

int cpl_rand() {
	SEED = (SEED * MULTIPLIER) % MODULUS;
	return SEED;
}

double cpl_rand_uniform(double a, double b) {
	return ((double) cpl_rand() / MODULUS) * (b - a) + a;
}

double cpl_mcint1d_naive(double (*f)(double), double a, double b, int N) {
	double x, Σf = 0;
	for (int i = 0; i < N; ++i) {
		x = cpl_rand_uniform(a, b);
		Σf += f(x);
	}

	return Σf * (b - a) / N;
}

double cpl_mcint_box_naive(double (*f)(cpl_vector *), int dim, double a, double b, int N) {
	double Σf = 0;
	cpl_vector *x = cpl_vector_alloc(dim);
	for (int i = 0; i < N; ++i) {
		for (int j = 1; j <= dim; ++j)
			cpl_set(x, j, cpl_rand_uniform(a, b));
		Σf += f(x);
	}

	return Σf * pow(b - a, dim) / N;
}
