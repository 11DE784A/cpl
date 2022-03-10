#include "../../cpl.h"

int main() {
	int dim = 6;
	cpl_matrix *A = cpl_matrix_alloc(dim, dim);
	cpl_matrix_build(A, -2.0,  0.0,  0.0, -1.0,  0.0,  0.5,
						 0.0,  4.0,  0.5,  0.0,  1.0,  0.0,
						 0.0,  0.5,  1.5,  0.0,  0.0,  0.0,
						-1.0,  0.0,  0.0, -2.0,  0.0,  1.0,
						 0.0,  1.0,  0.0,  0.0, -2.5,  0.0,
						 0.5,  0.0,  0.0,  1.0,  0.0, -3.75);

	cpl_matrix *B = cpl_matrix_alloc(dim, 1);
	cpl_matrix_build(B, -1.0, 0.0, 2.75, 2.5, -3.0, 2.0);

	cpl_matrix *X_jac = cpl_matrix_alloc(dim, 1);
	cpl_linalg_jacobi(A, B, X_jac, NULL);
	printf("Solution by Jacobi method:\n");
	cpl_print(X_jac);

	cpl_matrix *X_sei = cpl_matrix_alloc(dim, 1);
	cpl_linalg_seidel(A, B, X_sei);
	printf("Solution by Gauss-Seidel method:\n");
	cpl_print(X_sei);

	cpl_free(A);
	cpl_free(B);
	cpl_free(X_jac);
	cpl_free(X_sei);
}
