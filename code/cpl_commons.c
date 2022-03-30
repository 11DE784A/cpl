#include "cpl_defines.h"
#include "cpl_includes.h"

#include "cpl_commons.h"

/* Assertions */
void cpl_check(int cond, char* message) {
	if (!cond) {
		fprintf(stderr, "Assertion failed: ");
		fprintf(stderr, message);
		fprintf(stderr, "\nProcess terminated.\n");
		exit(EXIT_FAILURE);
	}
}

/* Random numbers */
double cpl_rand() {
	return (double) rand() / RAND_MAX;
}

double cpl_randn() {
	return sqrt(-2*log(cpl_rand())) * cos(2*PI*cpl_rand());
}

double cpl_rand_uniform(double a, double b) {
	return (b - a) * cpl_rand() + a;
}

