#include "../../cpl.h"

int main() {
	int dim = 10;
	cpl_vector *time = cpl_vector_alloc(dim);
	cpl_vector_build(time, 1.0, 15.0, 30.0, 45.0, 60.0, 75.0, 90.0, 105.0, 120.0, 135.0);

	cpl_vector *counts = cpl_vector_alloc(dim);
	cpl_vector_build(counts, 106.0, 80.0, 98.0, 75.0, 74.0, 73.0, 49.0, 38.0, 37.0, 22.0);

	cpl_vector *σ = cpl_vector_alloc(dim);
	cpl_vector_build(σ, 10.0, 9.0, 10.0, 9.0, 8.0, 8.0, 7.0, 6.0, 6.0, 5.0);

	cpl_vector *time_log = cpl_vector_alloc(dim);
	cpl_vector *counts_log = cpl_vector_alloc(dim);
	cpl_vector *σ_log = cpl_vector_alloc(dim);
	for (int i = 1; i <= dim; ++i) {
		cpl_set(time_log, i, log(cpl_get(time, i)));
		cpl_set(counts_log, i, log(cpl_get(counts, i)));
		cpl_set(σ_log, i, 1.0 / sqrt(cpl_get(σ, i)));
	}

	cpl_vector *params = cpl_vector_calloc(2);
	cpl_vector_build(params, 1.0, 1.0);
	cpl_matrix *Cov = cpl_matrix_calloc(2, 2);
	scalar χ2 = cpl_stats_linfit(time_log, counts_log, σ_log, params, Cov);

	printf("Fit equation: log(A) = %.4g%#.4gt\n", cpl_get(params, 1), cpl_get(params, 2));
	printf("χ2 value for the fitting: %g\n", χ2);
	printf("Covariance matrix for the fit:\n");
	cpl_print(Cov);

	scalar τ = sabs(1.0 / cpl_get(params, 2));
	scalar στ = τ * sqrt(cpl_get(Cov, 2, 2));
	printf("Average lifetime: %.4g ± %.4g seconds \n", τ, στ);
}
