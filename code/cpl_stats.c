#include "cpl_defines.h"
#include "cpl_includes.h"

#include "cpl_commons.h"

#include "cpl_arrays.h"
#include "cpl_linalg.h"
#include "cpl_stats.h"

scalar cpl_stats_mean(cpl_vector *v) {
	scalar Σ = 0;
	int N = cpl_vector_dim(v);
	for (int i = 1; i <= N; ++i)
		Σ += cpl_get(v, i);

	return Σ / N;
}

scalar cpl_stats_var(cpl_vector *v) {
	int N = cpl_vector_dim(v);
	scalar μ = cpl_stats_mean(v);

	scalar Σ = 0;
	for (int i = 1; i <= N; ++i)
		Σ += pow(cpl_get(v, i) - μ, 2);

	return Σ / (N - 1);
}

scalar cpl_stats_varm(cpl_vector *v, scalar μ) {
	int N = cpl_vector_dim(v);

	scalar Σ = 0;
	for (int i = 1; i <= N; ++i)
		Σ += pow(cpl_get(v, i) - μ, 2);

	return Σ / N;
}

scalar cpl_stats_std(cpl_vector *v) {
	return sqrt(cpl_stats_var(v));
}

scalar cpl_stats_stdm(cpl_vector *v, scalar μ) {
	return sqrt(cpl_stats_varm(v, μ));
}

cpl_tuple cpl_stats_jackknife(scalar (*f)(scalar), cpl_vector *x) {
	int N = cpl_vector_dim(x);

	scalar xi_JK;
	cpl_vector *x_JK = cpl_vector_calloc(N);
	for (int i = 1; i <= N; ++i) {
		xi_JK = 0;
		for (int j = 1; j <= N; ++j) {
			if (j == i) continue;
			xi_JK += cpl_get(x, j);
		}

		cpl_set(x_JK, i, xi_JK / (N - 1));
	}

	cpl_vector *f_JK = cpl_vector_calloc(N);
	for (int i = 1; i <= N; ++i)
		cpl_set(f_JK, i, f(cpl_get(x_JK, i)));

	scalar mean_JK = cpl_stats_mean(f_JK);
	scalar var_JK = 0;
	for (int i = 1; i <= N; ++i)
		var_JK += pow(mean_JK - cpl_get(f_JK, i), 2);

	var_JK *=  (N - 1.0) / N;

	cpl_free(f_JK);
	cpl_free(x_JK);

	return (cpl_tuple) {.first = mean_JK, .second = var_JK};
}

scalar cpl_stats_linfit(cpl_vector *x, cpl_vector *y, cpl_vector *σ, cpl_vector *params, cpl_matrix *Cov) {
	return cpl_stats_polyfit(1, x, y, σ, params, Cov);
}

scalar cpl_stats_polyfit(int deg, cpl_vector *x, cpl_vector *y, cpl_vector *σ, cpl_vector *params, cpl_matrix *Cov) {

	int N = cpl_vector_dim(y);

	cpl_vector *b = cpl_vector_calloc(deg + 1);
	cpl_matrix *A = cpl_matrix_calloc(deg + 1, deg + 1);

	scalar bk, Akl;
	for (int k = 1; k <= deg + 1; ++k) {
		bk = 0;
		for (int i = 1; i <= N; ++i)
			bk += pow(cpl_get(x, i), k - 1) * cpl_get(y, i) 
					/ (σ ? pow(cpl_get(σ, i), 2) : 1);

		cpl_set(b, k, bk);

		for (int l = 1; l <= deg + 1; ++l) {
			Akl = 0;
			for (int i = 1; i <= N; ++i)
				Akl += pow(cpl_get(x, i), k + l - 2) 
						/ (σ ? pow(cpl_get(σ, i), 2) : 1);

			cpl_set(A, k, l, Akl);
		}

	}

	/* Fitting parameters */
	cpl_linalg_conjgrad_solve(A, b, params, NULL);

	/* Covariance matrix */
	if (Cov) {
		cpl_matrix *id = cpl_matrix_id(deg + 1);
		cpl_linalg_seidel(A, id, Cov, NULL);
		cpl_free(id);
	}

	cpl_free(A);
	cpl_free(b);

	/* Goodness of fit */
	scalar χ2 = 0, yi;
	for (int i = 1; i <= N; ++i) {
		yi = 0;
		for (int k = 0; k <= deg; ++k)
			yi += cpl_get(params, k + 1) * pow(cpl_get(x, i), k);
		χ2 += pow((cpl_get(y, i) - yi) / (σ ? cpl_get(σ, i) : 1), 2);
	}

	χ2 /= (N - deg - 1);

	return χ2;
}
