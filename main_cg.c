#include "cpl.h"
#include "cpl_bench.h"


scalar Afunc(int x, int y) {
	float m = 0.2;
	scalar Axy = ((x + 1 == y ? 1 : 0) 
					 - (x - 1 == y ? 1 : 0) 
					 + 2 * (x == y ? 1 : 0)) / 2.0 + m*m * (x == y ? 1.0 : 0.0);
	return Axy;
}

int main() {
	int dim = 6;

	cpl_matrix *A = cpl_matrix_alloc(dim, dim);
	cpl_matrix_build(A, 2.0, -3.0, 0.0, 0.0, 0.0, 0.0,
						-1.0, 4.0, -1.0, 0.0, -1.0, 0.0,
						0.0, -1.0, 4.0, 0.0, 0.0, -1.0,
						0.0, 0.0, 0.0, 2.0, -3.0, 0.0,
						0.0, -1.0, 0.0, -1.0, 4.0, -1.0,
						0.0, 0.0, -1.0, 0.0, -1.0, 4.0);
	// cpl_print(A);

	cpl_matrix *B = cpl_matrix_alloc(dim, 1);
	cpl_matrix_build(B, 5.0/2, 2.0/3, 3.0, -4.0/3, -1.0/3, 5.0/3);

	cpl_matrix *X = cpl_matrix_copy(B);
	cpl_linalg_gaussjordan(A, X);

	cpl_matrix *Ainv = cpl_matrix_id(dim);
	cpl_print(Ainv);
	cpl_linalg_gaussjordan(A, Ainv);

	cpl_mult(A, Ainv);
	cpl_print(Ainv);

	// Testing Conjugate Gradient
	dim = 4;
	cpl_matrix *C = cpl_matrix_alloc(dim, dim);
	cpl_matrix_build(C, 4.0, -2.0, 2.0, 1.0,
						-2.0, 2.0, -4.0, 2.0,
						2.0, -4.0, 11.0, -3.0,
						1.0, 2.0, -3.0,  13.0);
	cpl_print(C);

	cpl_vector *b = cpl_vector_alloc(dim);
	cpl_vector_build(b, 1.0, 2.0, 2.0, 1.0);

	cpl_print(b);

	cpl_vector *x = cpl_vector_calloc(dim);
	cpl_linalg_conjgrad_solve(C, b, x);
	cpl_print(x);

	cpl_mult(C, x);
	cpl_print(x);

	// Testing Jacobi
	// cpl_print(A);
	// cpl_vector *res = cpl_vector_alloc(0);
	// cpl_matrix *Y = cpl_matrix_id(dim);
	// cpl_matrix *id = cpl_matrix_id(dim);
	// time(cpl_linalg_jacobi(A, id, Y, res));
	// cpl_print(Y);

	// cpl_print(res);

	// cpl_free(Y);
	cpl_free(x);
	cpl_free(b);
	cpl_free(C);
	cpl_free(Ainv);
	cpl_free(X);
	cpl_free(B);
	cpl_free(A);

	return 0;
}
