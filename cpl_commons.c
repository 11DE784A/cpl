#include "cpl_includes.h"

void cpl_check(int cond, char* message) {
	if (!cond) {
		fprintf(stderr, "Assertion failed: ");
		fprintf(stderr, message);
		fprintf(stderr, "\nProcess terminated.\n");
		exit(EXIT_FAILURE);
	}
}

