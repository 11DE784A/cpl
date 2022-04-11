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

/* Delta function */
int cpl_delta(int i, int j) {
	return (i == j ? 1 : 0);
}

