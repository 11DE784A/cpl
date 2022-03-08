#include "cpl_defines.h"
#include "cpl_includes.h"

#include "cpl_commons.h"

#include "cpl_arrays.h"
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
