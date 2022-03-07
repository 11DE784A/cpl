#include "cpl.h"

scalar cpl_get2(void *X, ...) {
	va_list ap;
	va_start(ap, X);
	int i = va_arg(ap, int);
	if (X->type == VECTOR)
		return cpl_vector_get(X, i);

	int j = va_arg(ap, int);
	if (X->type == MATRIX)
		return cpl_matrix_get(X, i, j);
}

scalar generic_solve(void *A) {
	return cpl_get2(A, 1) + 1;
}

scalar func(int i, int j) {
	return i + j;
}

int main() {
	cpl_matrix *M = cpl_matrix_calloc(3, 3);
	printf("%g\n", generic_solve(M));
	cpl_free(M);
	return 0;
}
