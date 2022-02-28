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
	cpl_check(cpl_matrix_rows(U) == cpl_matrix_rows(X),
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

cpl_tensor *cpl_linalg_ludecomp(cpl_tensor *A) {
	cpl_check(cpl_matrix_issquare(A),
			  "Matrix passed to cpl_linalg_ludecomp should be square");

	cpl_tensor *LU = cpl_tensor_alloc_shape(cpl_tuple_copy(cpl_tensor_shape(A)));
	int dim = cpl_tuple_get(cpl_tensor_shape(A), 1);
	for (int j = 1; j <= dim; ++j) {
		for (int i = 1; i <= j; ++i) {
			scalar Lik_Ukj = 0;
			for (int k = 1; k <= i - 1; ++k)
				Lik_Ukj += cpl_tensor_get(LU, i, k) * cpl_tensor_get(LU, k, j);

			cpl_tensor_set(LU, i, j, cpl_tensor_get(A, i, j) - Lik_Ukj);
		}

		for (int i = j + 1; i <= dim; ++i) {
			scalar Lik_Ukj = 0;
			for (int k = 1; k <= j - 1; ++k)
				Lik_Ukj += cpl_tensor_get(LU, i, k) * cpl_tensor_get(LU, k, j);

			cpl_tensor_set(LU, i, j, (cpl_tensor_get(A, i, j) - Lik_Ukj) 
										/ cpl_tensor_get(LU, j, j));
		}
	}

	return LU;
}

cpl_tensor *cpl_linalg_cholesky(cpl_tensor *A) {
	cpl_check(cpl_matrix_issquare(A),
			  "Matrix passed to cpl_linalg_cholesky should be square");

	int dim = cpl_tuple_get(cpl_tensor_shape(A), 1);
	cpl_tensor *L = cpl_tensor_alloc_shape(cpl_tuple_copy(cpl_tensor_shape(A)));
	cpl_tensor_set_all(L, 0.0);

	for (int i = 1; i <= dim; ++i) {
		scalar Lik2 = 0;
		for (int k = 1; k <= i - 1; ++k)
			Lik2 += pow(cpl_tensor_get(L, i, k), 2);

		cpl_tensor_set(L, i, i, sqrt(cpl_tensor_get(A, i, i) - Lik2));

		for (int j = i + 1; j <= dim; ++j) {
			scalar Lik_Ljk = 0;
			for (int k = 1; k <= i - 1; ++k)
				Lik_Ljk += cpl_tensor_get(L, i, k) * cpl_tensor_get(L, j, k);

			cpl_tensor_set(L, j, i, (cpl_tensor_get(A, i, j) - Lik_Ljk)
									/ cpl_tensor_get(L, i, i));
		}
	}
	return L;
}

void cpl_linalg_jacobi(cpl_tensor *U, cpl_tensor *B, cpl_tensor *x0) {
	cpl_check(cpl_matrix_issquare(U),
			  "Matrix passed to cpl_linalg_jacobi should be square");
	cpl_check(cpl_matrix_rows(U) == cpl_vector_dim(B)
				&& cpl_matrix_cols(U) == cpl_vector_dim(x0),
			  "UX = B is not a well defined system");

	cpl_tensor *x1 = cpl_tensor_copy(x0);

	int iters = 0;
	scalar residual = 1.0 / TOL;
	while (residual  > TOL) {
		if (++iters > MAX_ITERS) {
			cpl_tensor_free(x1);
			fprintf(stderr, "Failed to converge after %d iterations\n", MAX_ITERS);
			return;
		}

		for (int i = 1; i <= cpl_vector_dim(x1); ++i) {
			scalar Uij_xj = 0;
			for (int j = 1; j <= cpl_vector_dim(x0); ++j) {
				if (i == j) continue;
				Uij_xj += cpl_tensor_get(U, i, j) * cpl_tensor_get(x0, j);
			}

			cpl_tensor_set(x1, i,
						   (cpl_tensor_get(B, i) - Uij_xj) / cpl_tensor_get(U, i, i));

			residual = cpl_tensor_distance(x0, x1);
			cpl_tensor_overwrite(x0, x1);
		}
	}

	cpl_tensor_free(x1);
}

void cpl_linalg_gaussseidel(cpl_tensor *U, cpl_tensor *B, cpl_tensor *x0) {
	cpl_check(cpl_matrix_issquare(U),
			  "Matrix passed to cpl_linalg_gaussseidel should be square");
	cpl_check(cpl_matrix_rows(U) == cpl_vector_dim(B)
				&& cpl_matrix_cols(U) == cpl_vector_dim(x0),
			  "UX = B is not a well defined system");

	cpl_tensor *x1 = cpl_tensor_copy(x0);

	int iters = 0;
	scalar residual = 1.0 / TOL;
	while (residual > TOL) {
		if (++iters > MAX_ITERS) {
			cpl_tensor_free(x1);
			fprintf(stderr, "Failed to converge after %d iterations\n", MAX_ITERS);
			return;
		}

		for (int i = 1; i <= cpl_vector_dim(x1); ++i) {
			scalar Uij_xj_1 = 0, Uij_xj_2 = 0;
			for (int j = 1; j <= i - 1; ++j)
				Uij_xj_1 += cpl_tensor_get(U, i, j) * cpl_tensor_get(x1, j);

			for (int j = i+1; j <= cpl_vector_dim(x0); ++j)
				Uij_xj_2 += cpl_tensor_get(U, i, j) * cpl_tensor_get(x0, j);

			cpl_tensor_set(x1, i,
						   (cpl_tensor_get(B, i) - Uij_xj_1 - Uij_xj_2) 
								/ cpl_tensor_get(U, i, i));

			residual = cpl_tensor_distance(x0, x1);
			cpl_tensor_overwrite(x0, x1);
		}
	}

	cpl_tensor_free(x1);
}

