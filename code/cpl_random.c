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
