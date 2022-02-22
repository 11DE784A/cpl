#include "cpl_includes.h"
#include "cpl_defines.h"

#include "cpl.h"

int main() {
	int dim = 4;
	cpl_tensor *A = cpl_matrix_alloc(dim, dim);
	for (int i = 0; i < cpl_tensor_length(A); ++i)
		A->array[i] = i+100;

	cpl_matrix_print(A);

	cpl_tensor *LU = cpl_linalg_ludecomp(A);
	cpl_matrix_print(LU);

	cpl_tensor *L = cpl_matrix_id(dim);

	cpl_tensor *U = cpl_matrix_alloc(dim, dim);
	cpl_tensor_set_all(U, 0);

	for (int i = 1; i <= dim; ++i) {
		for (int j = 1; j <= dim; ++j) {
			if (i <= j) {
				cpl_tensor_set(U, i, j, cpl_tensor_get(LU, i, j));
			} else {
				cpl_tensor_set(L, i, j, cpl_tensor_get(LU, i, j));
			}
		}
	}

	cpl_matrix_print(L);
	cpl_matrix_print(U);

	cpl_matrix_print(cpl_matrix_mult(L, U));

	cpl_tensor_free(A);
	cpl_tensor_free(LU);

	return 0;
}
