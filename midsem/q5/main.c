#include "../../cpl.h"

const scalar a = -1;
const scalar c = a;
const scalar b = 2;
const scalar n = 5;

scalar calculated_eigvals(int k) {
	scalar λk = b + 2 * sqrt(a * c) * cos(k*PI / (n + 1));
	return λk;
}

cpl_vector *calculated_eigvecs(int k, cpl_vector *v) {
	for (int i = 1; i <= cpl_vector_dim(v); ++i)
		cpl_set(v, i, 2 * pow(sqrt(c/a), k) * sin(i*k*PI / (n + 1)));

	return v;
}

int main() {
	cpl_matrix *A = cpl_matrix_alloc(n, n);
	cpl_matrix_build(A, 2.0,  -1.0,   0.0,   0.0,   0.0,
					   -1.0,   2.0,  -1.0,   0.0,   0.0,
						0.0,  -1.0,   2.0,  -1.0,   0.0,
						0.0,   0.0,  -1.0,   2.0,  -1.0,
						0.0,   0.0,   0.0,  -1.0,   2.0);

	cpl_vector *v1 = cpl_vector_1hot(n, 1);
	scalar λ1 = cpl_eigen_power(A, v1);
	printf("Largest eigenvalue: %g\n", λ1);
	printf("Corresponding eigenvector:\n");
	cpl_print(v1);

	cpl_matrix *A_def = cpl_matrix_alloc(n, n);
	for (int i = 1; i <= n; ++i) {
		for (int j = 1; j <= n; ++j) {
			cpl_set(A_def, i, j, cpl_get(v1, i) * cpl_get(v1, j));
		}
	}

	cpl_matrix_add(A, cpl_matrix_scale(A_def, -λ1));

	cpl_vector *v2 = cpl_vector_1hot(n, 1);
	scalar λ2 = cpl_eigen_power(A_def, v2);
	printf("Second largest eigenvalue: %g\n", λ2);
	printf("Corresponding eigenvector:\n");
	cpl_print(v2);


	printf("Calculated eigenvalues and eigenvectors\n");
	scalar λ1_calc = calculated_eigvals(1);
	printf("Largest eigenvalue: %g\n", λ1_calc);
	cpl_vector *v1_calc = cpl_vector_alloc(n);
	calculated_eigvecs(1, v1_calc);
	cpl_vector_normalize(v1_calc);
	cpl_print(v1_calc);

	scalar λ2_calc = calculated_eigvals(2);
	printf("Second largest eigenvalue: %g\n", λ2_calc);
	cpl_vector *v2_calc = cpl_vector_alloc(n);
	calculated_eigvecs(2, v2_calc);
	cpl_vector_normalize(v2_calc);
	cpl_print(v2_calc);
}
