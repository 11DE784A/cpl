#include "../code/cpl.h"
#include "../code/cpl_bench.h"

int main() {
	/*** Question 1 ***/
	printf("# Question 1\n\n");

	int dim = 6;

	cpl_matrix *A = cpl_matrix_alloc(dim, dim);
	cpl_matrix_build(A,  1.0, -1.0,  4.0,  0.0,  2.0,  9.0,
						 0.0,  5.0, -2.0,  7.0,  8.0,  4.0,
						 1.0,  0.0,  5.0,  7.0,  3.0, -2.0,
						 6.0, -1.0,  2.0,  3.0,  0.0,  8.0,
						-4.0,  2.0,  0.0,  5.0, -5.0,  3.0,
						 0.0,  7.0, -1.0,  5.0,  4.0, -2.0);

	printf("## Matrix A\n");
	cpl_print(A);

	cpl_matrix *B = cpl_matrix_alloc(dim, 1);
	cpl_matrix_build(B, 19.0, 2.0, 13.0, -7.0, -9.0, 2.0);

	printf("## Matrix B\n");
	cpl_print(B);

	cpl_matrix *X = cpl_matrix_copy(B);
	cpl_time(cpl_linalg_gaussjordan(A, X));

	printf("## Solution by Gauss-Jordan\n");
	cpl_print(X);

	cpl_vector *b = cpl_matrix_flatten(B);
	cpl_time(cpl_linalg_lusolve(A, b));

	printf("## Solution by LU decomposition\n");
	cpl_print(b);


	/*** Question 2 ***/
	printf("# Question 2\n\n");

	dim = 6;

	cpl_matrix_build(A,  2.0, -3.0,  0.0,  0.0,  0.0,  0.0,
						-1.0,  4.0, -1.0,  0.0, -1.0,  0.0,
						 0.0, -1.0,  4.0,  0.0,  0.0, -1.0,
						 0.0,  0.0,  0.0,  2.0, -3.0,  0.0,
						 0.0, -1.0,  0.0, -1.0,  4.0, -1.0,
						 0.0,  0.0, -1.0,  0.0, -1.0,  4.0);

	printf("## Matrix A\n");
	cpl_print(A);

	cpl_matrix_build(B, -5./3, 2./3, 3., -4./3, -1./3, 5./3);
	printf("## Matrix B\n");
	cpl_print(B);

	/* Solution by LU */
	cpl_vector *x = cpl_matrix_flatten(B);
	cpl_time(cpl_linalg_lusolve(A, x));

	printf("## Solution by LU decomposition\n");
	cpl_print(x);

	/* Solution by Jacobi */
	cpl_matrix_set_all(X, 0.0);
	cpl_time(cpl_linalg_jacobi(A, B, X, NULL));

	printf("## Solution by Jacobi iteration\n");
	cpl_print(X);

	/* Solution by Gauss-Seidel */
	cpl_matrix_set_all(X, 0.0);
	cpl_time(cpl_linalg_seidel(A, B, X, NULL));

	printf("## Solution by Gauss-Seidel\n");
	cpl_print(X);

	/** Inversion **/
	cpl_matrix *id = cpl_matrix_id(dim);

	/* Inverse by Jacobi */
	cpl_vector *jac_res = cpl_vector_alloc(0);
	cpl_matrix *A_inv = cpl_matrix_calloc(dim, dim);
	cpl_time(cpl_linalg_jacobi(A, id, A_inv, jac_res));

	printf("## Inverse by Jacobi iteration\n");
	cpl_print(A_inv);

	printf("Residue versus iterations\n");
	cpl_print(jac_res);

	/* Inverse by Gauss-Seidel */
	cpl_vector *sed_res = cpl_vector_alloc(0);
	cpl_matrix_set_all(A_inv, 0.0);
	cpl_time(cpl_linalg_seidel(A, id, A_inv, sed_res));

	printf("## Inverse by Gauss-Seidel method\n");
	cpl_print(A_inv);

	printf("Residue versus iterations\n");
	cpl_print(sed_res);

	/*** Question 3 ***/
	printf("# Question 3\n\n");
	dim = 20;

	int μ = 1;
	scalar m = 0.2;
	scalar A_fn (int x, int y) {
		return (cpl_delta(x + μ, y) + cpl_delta(x - μ, y)) / 2.0 
				+ (m*m - 1) * cpl_delta(x, y); }

	int N = 3;

	cpl_vector *e = cpl_vector_alloc(dim);
	cpl_vector *f = cpl_vector_alloc(dim);
	cpl_vector *res3 = cpl_vector_alloc(0);

	printf("We will find the inverse column-by-column for %d columns\n\n", N);
	for (int i = 1; i <= N; ++i) {
		cpl_vector_set_all(f, 1.0);

		cpl_vector_set_all(e, 0.0);
		cpl_vector_set(e, i, 1.0);

		cpl_time(cpl_linalg_conjgrad_fly(A_fn, e, f, res3));

		printf("Column %d of the inverse:\n", i);
		cpl_print(f);

	}

	printf("Residue versus iterations\n");
	cpl_print(res3);

	cpl_free(res3);
	cpl_free(f);
	cpl_free(e);
	cpl_free(sed_res);
	cpl_free(jac_res);
	cpl_free(x);
	cpl_free(b);
	cpl_free(id);
	cpl_free(A_inv);
	cpl_free(X);
	cpl_free(B);
	cpl_free(A);

	return 0;
}
