#include "cpl.h"

#include "cpl_bench.h"

int main() {
/**	int dim = 6;

	cpl_matrix *A = cpl_matrix_alloc(dim, dim);
	cpl_matrix_build(A, 1.0, -1.0, 4.0, 0.0, 2.0, 9.0,
						0.0, 5.0, -2.0, 7.0, 8.0, 4.0,
						1.0, 0.0, 5.0, 7.0, 3.0, -2.0,
						6.0, -1.0, 2.0, 3.0, 0.0, 8.0,
						-4.0, 2.0, 0.0, 5.0, -5.0, 3.0,
						0.0, 7.0, -1.0, 5.0, 4.0, -2.0);
	cpl_print(A);

	cpl_vector *b = cpl_vector_alloc(dim);
	cpl_vector_build(b, -5.0/2, 2.0/3, 3.0, -4.0/3, -1.0/3, 5.0/3);
	cpl_print(b);

	cpl_matrix *B = cpl_matrix_alloc(dim, 1);
	cpl_matrix_build(B, -5.0/2, 2.0/3, 3.0, -4.0/3, -1.0/3, 5.0/3);
	cpl_print(B);

	time(cpl_linalg_gaussjordan(A, B));
	cpl_print(B);

	time(cpl_linalg_lusolve(A, b));
	cpl_print(b);

	scalar x = cpl_linalg_ludet(A);
	printf("det: %g\n", x); **/

	return 0;
}
