#include "cpl_defines.h"
#include "cpl_includes.h"

#include "cpl_commons.h"

#include "cpl_arrays.h"
#include "cpl_linalg.h"

cpl_matrix *cpl_linalg_gaussjordan(cpl_matrix *U, cpl_matrix *B) {
	/* Returns NULL if U is singular */

	cpl_check(cpl_matrix_issquare(U),
			  "Matrix passed to cpl_linalg_gaussjordan should be square");
	cpl_check(cpl_matrix_rows(U) == cpl_matrix_rows(B),
			  "UX = B is not a well defined system");

	cpl_matrix *V = cpl_matrix_copy(U);

	cpl_vector *v = cpl_vector_alloc(cpl_matrix_cols(V));
	cpl_vector *v_copy = cpl_vector_like(v);

	cpl_vector *x = cpl_vector_alloc(cpl_matrix_cols(B));
	cpl_vector *x_copy = cpl_vector_like(x);

	int swaps = 0;
	for (int j = 1; j <= cpl_matrix_cols(V); ++j) {
		for (int i = j; i <= cpl_matrix_rows(V); ++i) {
			scalar Vij = cpl_get(V, i, j);
			if (Vij) {
				swaps += cpl_matrix_swap_rows(V, i, j);
				cpl_matrix_scale_row(V, j, 1/Vij);

				cpl_matrix_swap_rows(B, i, j);
				cpl_matrix_scale_row(B, j, 1/Vij);

				break;

			} else if (i == cpl_matrix_rows(V)) {
				cpl_free(x_copy);
				cpl_free(x);
				cpl_free(v_copy);
				cpl_free(v);
				cpl_free(V);
				cpl_free(B);
				return NULL;
			}
		}

		cpl_matrix_get_row(V, j, v);
		cpl_matrix_get_row(B, j, x);

		for (int i = 1; i <= cpl_matrix_rows(V); ++i) {
			scalar Vij = cpl_get(V, i, j);
			if (Vij == 0.0 || i == j) continue;

			cpl_vector_overwrite(v_copy, v);
			cpl_vector_scale(v_copy, -Vij);
			cpl_matrix_add_to_row(V, i, v_copy);

			cpl_vector_overwrite(x_copy, x);
			cpl_vector_scale(x_copy, -Vij);
			cpl_matrix_add_to_row(B, i, x_copy);
		}
	}

	cpl_free(x_copy);
	cpl_free(x);
	cpl_free(v_copy);
	cpl_free(v);
	cpl_free(V);

	return B;
}


cpl_matrix *cpl_linalg_ludecomp(cpl_matrix *A) {
	cpl_check(cpl_matrix_issquare(A),
			  "Matrix passed to cpl_linalg_ludecomp should be square");

	int dim = cpl_matrix_rows(A);
	cpl_matrix *LU = cpl_matrix_alloc(dim, dim);
	for (int j = 1; j <= dim; ++j) {
		for (int i = 1; i <= j; ++i) {
			scalar Lik_Ukj = 0;
			for (int k = 1; k <= i - 1; ++k)
				Lik_Ukj += cpl_get(LU, i, k) * cpl_get(LU, k, j);

			cpl_set(LU, i, j, cpl_get(A, i, j) - Lik_Ukj);
		}

		for (int i = j + 1; i <= dim; ++i) {
			scalar Lik_Ukj = 0;
			for (int k = 1; k <= j - 1; ++k)
				Lik_Ukj += cpl_get(LU, i, k) * cpl_get(LU, k, j);

			cpl_set(LU, i, j, (cpl_get(A, i, j) - Lik_Ukj) 
										/ cpl_get(LU, j, j));
		}
	}

	return LU;
}

cpl_matrix *cpl_linalg_cholesky(cpl_matrix *A) {
	cpl_check(cpl_matrix_issquare(A),
			  "Matrix passed to cpl_linalg_cholesky should be square");

	int dim = cpl_matrix_rows(A);
	cpl_matrix *L = cpl_matrix_calloc(dim, dim);

	for (int i = 1; i <= dim; ++i) {
		scalar Lik2 = 0;
		for (int k = 1; k <= i - 1; ++k)
			Lik2 += pow(cpl_get(L, i, k), 2);

		cpl_set(L, i, i, sqrt(cpl_get(A, i, i) - Lik2));

		for (int j = i + 1; j <= dim; ++j) {
			scalar Lik_Ljk = 0;
			for (int k = 1; k <= i - 1; ++k)
				Lik_Ljk += cpl_get(L, i, k) * cpl_get(L, j, k);

			cpl_set(L, j, i, (cpl_get(A, i, j) - Lik_Ljk) / cpl_get(L, i, i));
		}
	}
	return L;
}

int cpl_linalg_jacobi(cpl_matrix *U, cpl_matrix *B, cpl_matrix *X0) {
	cpl_matrix *X1 = cpl_matrix_copy(X0);

	int iter = 0;
	scalar residual = 1 / TOL;
	while (residual > TOL) {
		if (++iter > MAX_ITERS) {
			fprintf(stderr, "Failed to converge after %d iterations", MAX_ITERS);
			break;
		}

		residual = 0;
		for (int i = 1; i <= cpl_matrix_rows(X1); ++i) {
			for (int j = 1; j <= cpl_matrix_cols(X1); ++j) {
				scalar Uik_Xkj = 0;
				for (int k = 1; k <= cpl_matrix_rows(X0); ++k) {
					if (k == i) continue;
					Uik_Xkj += cpl_get(U, i, k) * cpl_get(X0, k, j);
				}

				cpl_set(X1, i, j, (cpl_get(B, i, j) - Uik_Xkj) / cpl_get(U, i, i));
				residual += sabs(cpl_get(X0, i, j) - cpl_get(X1, i, j));
				cpl_set(X0, i, j, cpl_get(X1, i, j));
			}
		}
	}

	cpl_free(X1);

	return iter;
}

int cpl_linalg_seidel(cpl_matrix *U, cpl_matrix *B, cpl_matrix *X) {

	int iter = 0;
	scalar residual = 1 / TOL;
	while (residual > TOL) {
		if (++iter > MAX_ITERS) {
			fprintf(stderr, "Failed to converge after %d iterations", MAX_ITERS);
			break;
		}

		residual = 0;
		for (int i = 1; i <= cpl_matrix_rows(X); ++i) {
			for (int j = 1; j <= cpl_matrix_cols(X); ++j) {
				scalar Uik_Xkj_1 = 0;
				for (int k = 1; k <= i - 1; ++k)
					Uik_Xkj_1 += cpl_get(U, i, k) * cpl_get(X, k, j);

				scalar Uik_Xkj_2 = 0;
				for (int k = i + 1; k <= cpl_matrix_rows(X); ++k)
					Uik_Xkj_2 += cpl_get(U, i, k) * cpl_get(X, k, j);

				scalar Xij = (cpl_get(B, i, j) - Uik_Xkj_1 - Uik_Xkj_2) 
															/ cpl_get(U, i, i);
				residual += sabs(Xij - cpl_get(X, i, j));
				cpl_set(X, i, j, Xij);
			}
		}
	}

	return iter;
}

void cpl_linalg_conjgrad(cpl_matrix *U, cpl_vector *b, cpl_vector *x) {
	int dim = cpl_vector_dim(x);

	cpl_vector *r = cpl_add(b, cpl_scale(cpl_mult_alloc(U, x), -1));
	cpl_vector *d = cpl_vector_copy(r);

	scalar rold = cpl_vector_inner(r, r);
	scalar alpha, rnew;

	for (int i = 1; i <= dim; ++i) {
		cpl_vector *Ud = cpl_mult_alloc(U, d);
		alpha = rold / cpl_vector_inner(d, Ud);
		for (int j = 1; j <= dim; ++j) {
			cpl_set(x, j, cpl_get(x, j) + alpha * cpl_get(d, j));
			cpl_set(r, j, cpl_get(r, j) - alpha * cpl_get(Ud, j));
		}

		cpl_free(Ud);

		rnew = cpl_vector_inner(r, r);
		if (rnew < TOL) break;


		for (int j = 1; j <= dim; ++j)
			cpl_set(d, j, cpl_get(r, j) + (rnew / rold) * cpl_get(d, j));

		rold = rnew;
	}

	cpl_free(r);
	cpl_free(d);
}

