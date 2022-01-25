#include "cpl_defines.h"
#include "cpl_includes.h"

#include "cpl.h"

cpl_tensor *cpl_linalg_gaussjordan(cpl_tensor* U, cpl_tensor *B) {
	/* Returns NULL if U is singular */

	cpl_check(cpl_matrix_issquare(U),
			  "Matrix passed to cpl_linalg_gaussjordan should be square");
	cpl_tensor *V = cpl_tensor_copy(U);

	cpl_tensor *X = cpl_tensor_copy(B);
	if (cpl_tensor_rank(X) == 1)
		cpl_vector_to_matrix(X);
	cpl_check(cpl_matrix_rows(U) == cpl_matrix_rows(B),
			  "UX = B is not a well defined system");

	int swaps = 0;

	for (int j = 1; j <= cpl_matrix_cols(V); ++j) {
		for (int i = j; i <= cpl_matrix_rows(V); ++i) {
			scalar Vij = cpl_tensor_get(V, i, j);
			if (Vij) {
				swaps += cpl_matrix_swap_rows(V, i, j);
				cpl_matrix_scale_row(V, j, 1/Vij);

				cpl_matrix_swap_rows(X, i, j);
				cpl_matrix_scale_row(X, j, 1/Vij);

				break;

			} else if (i == cpl_matrix_rows(V)) {
				cpl_tensor_free(X);
				cpl_tensor_free(V);
				return NULL;
			}
		}

		cpl_tensor *v = cpl_matrix_get_row(V, j);
		cpl_tensor *x = cpl_matrix_get_row(X, j);
		for (int i = 1; i <= cpl_matrix_rows(V); ++i) {
			scalar Vij = cpl_tensor_get(V, i, j);
			if (Vij == 0.0 || i == j) continue;

			cpl_tensor *v_copy = cpl_tensor_copy(v);
			cpl_tensor_scale(v_copy, -Vij);
			cpl_matrix_add_to_row(V, i, v_copy);
			cpl_tensor_free(v_copy);

			cpl_tensor *x_copy = cpl_tensor_copy(x);
			cpl_tensor_scale(x_copy, -Vij);
			cpl_matrix_add_to_row(X, i, x_copy);
			cpl_tensor_free(x_copy);
		}

		cpl_tensor_free(v);
		cpl_tensor_free(x);
	}

	cpl_tensor_free(V);

	return X;
}

