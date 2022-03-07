#include "cpl.h"
#include "cpl_bench.h"

int main() {
	int dim = 6;

	cpl_matrix *A = cpl_matrix_alloc(dim, dim);
	cpl_matrix_build(A,  2.0, -3.0,  0.0,  0.0,  0.0,  0.0,
						-1.0,  4.0, -1.0,  0.0, -1.0,  0.0,
						 0.0, -1.0,  4.0,  0.0,  0.0, -1.0,
						 0.0,  0.0,  0.0,  2.0, -3.0,  0.0,
						 0.0, -1.0,  0.0, -1.0,  4.0, -1.0,
						 0.0,  0.0, -1.0,  0.0, -1.0,  4.0);

	cpl_vector *x = cpl_vector_calloc(dim);
	cpl_vector_set(x, 1, 1.0);

	scalar λ = cpl_eigen_power(A, x);
	cpl_print(A);
	cpl_print(x);
	printf("Eigenvalue, λ = %.4g\n", λ);
	return 0;
}
