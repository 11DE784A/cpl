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

scalar polynomial(int n, scalar x) {
	return pow(x, n);
}

scalar hermite(int n, scalar x) {
	if (n == 0) {
		return 1;
	} else if (n == 1) {
		return 2*x;
	} else {
		return 2*x*hermite(n - 1, x) - 2*n*hermite(n - 2, x);
	}
}

scalar cpl_stats_linfit(cpl_vector *x, cpl_vector *y, cpl_vector *σ, cpl_vector *params, cpl_matrix *Cov) {
	return cpl_stats_linreg(2, polynomial, x, y, σ, params, Cov);
}

scalar cpl_stats_linreg(int n, scalar (*F)(int, scalar),
		cpl_vector *x, cpl_vector *y, cpl_vector *σ, 
		cpl_vector *params, cpl_matrix *Cov) {

	int N = cpl_vector_dim(y);

	cpl_vector *b = cpl_vector_calloc(n);
	cpl_matrix *A = cpl_matrix_calloc(n, n);

	scalar bk, Akl;
	for (int k = 1; k <= n; ++k) {
		bk = 0;
		for (int i = 1; i <= N; ++i)
			bk += F(k - 1, cpl_get(x, i)) * cpl_get(y, i)
					/ (σ ? pow(cpl_get(σ, i), 2) : 1);

		cpl_set(b, k, bk);

		for (int l = 1; l <= n; ++l) {
			Akl = 0;
			for (int i = 1; i <= N; ++i)
				Akl += F(k - 1, cpl_get(x, i)) * F(l - 1, cpl_get(x, i))
						/ (σ ? pow(cpl_get(σ, i), 2) : 1);

			cpl_set(A, k, l, Akl);
		}

	}

	/* Fitting parameters */
	cpl_linalg_conjgrad_solve(A, b, params, NULL);

	/* Covariance matrix */
	if (Cov) {
		cpl_matrix *id = cpl_matrix_id(n);
		cpl_linalg_seidel(A, id, Cov, NULL);
		cpl_free(id);
	}

	cpl_free(A);
	cpl_free(b);

	/* Goodness of fit */
	scalar χ2 = 0, yi;
	for (int i = 1; i <= N; ++i) {
		yi = 0;
		for (int k = 0; k < n; ++k)
			yi += cpl_get(params, k + 1) * F(k, cpl_get(x, i));
		χ2 += pow((cpl_get(y, i) - yi) / (σ ? cpl_get(σ, i) : 1), 2);
	}

	χ2 /= (N - n);

	return χ2;
}
