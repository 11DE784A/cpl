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
	int N = cpl_vector_dim(y);

	scalar U00 = 0, U10 = 0, U11 = 0, v0 = 0, v1 = 0, σi2;
	for (int i = 1; i <= N; ++i) {
		σi2 = pow(cpl_get(σ, i), 2);

		U00 += 1.0 / σi2;
		U10 += cpl_get(x, i) / σi2;
		U11 += pow(cpl_get(x, i), 2) / σi2;

		v0 += cpl_get(y, i) / σi2;
		v1 += cpl_get(y, i) * cpl_get(x, i) / σi2;
	}

	scalar Δ = U00*U11 - U10*U10;

	cpl_set(params, 1, (U11*v0 - U10*v1) / Δ);
	cpl_set(params, 2, (-U10*v0 + U00*v1) / Δ);

	cpl_set(Cov, 1, 1, U11);
	cpl_set(Cov, 1, 2, -U10);
	cpl_set(Cov, 2, 1, -U10);
	cpl_set(Cov, 2, 2, U00);
	cpl_matrix_scale(Cov, 1.0 / Δ);

	scalar χ2 = 0;
	for (int i = 1; i <= N; ++i) {
		χ2 += pow((cpl_get(y, i) - cpl_get(params, 1) - cpl_get(params, 2) * cpl_get(x, i)) / cpl_get(σ, i), 2);
	}

	return χ2 / (N - 2.0);
}

scalar cpl_stats_polyfit(int degree, cpl_vector *x, cpl_vector *y, cpl_vector *σ, 
					   cpl_vector *params, cpl_matrix *Cov) {

	int N = cpl_vector_dim(y);

	cpl_vector *b = cpl_vector_calloc(degree + 1);
	cpl_matrix *A = cpl_matrix_calloc(degree + 1, degree + 1);

	for (int i = 0; i <= degree; ++i) {
		scalar bi = 0;
		for (int j = 0; j <= degree; ++j) {
			scalar Aij = 0;
			for (int k = 1; k <= N; ++k)
				Aij += pow(cpl_get(x, k), i + j - 2) / pow(cpl_get(σ, k), 2);
			cpl_set(A, i, j, Aij);

			bi += (cpl_get(y, j) * pow(cpl_get(x, j), i)) / pow(cpl_get(σ, j), 2);
		}
		cpl_set(b, i, bi);
	}

	// Parameters
	cpl_linalg_conjgrad_solve(A, b, params, NULL);

	// Covariance matrix
	cpl_matrix *id = cpl_matrix_id(degree + 1);
	cpl_linalg_conjgrad(A, id, Cov);

	// Goodness of fit
	scalar χ2 = 0;
	for (int i = 1; i <= N; ++i) {
		scalar yi = 0;
		for (int k = 0; k <= degree; ++k)
			yi += cpl_get(params, k) * pow(cpl_get(x, i), k);
		χ2 += pow((cpl_get(y, i) - yi) / cpl_get(σ, i), 2);
	}

	int Ndof = N - degree - 1;
	χ2 /= Ndof;

	cpl_free(id);
	cpl_free(b);
	cpl_free(A);

	return χ2;
}
