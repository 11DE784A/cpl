#include "cpl_defines.h"
#include "cpl_includes.h"

#include "cpl_commons.h"

#include "cpl_arrays.h"
#include "cpl_eigen.h"

scalar cpl_eigen_power(cpl_matrix *A, cpl_vector *x0) {
	scalar λ_old, λ_new;
	cpl_vector *y = cpl_vector_copy(x0);
	cpl_vector *x1 = cpl_vector_copy(x0);

	scalar residue = 1 / TOL;
	int iters = 0;
	while (residue > TOL) {
		if (++iters > MAX_ITERS) {
			fprintf(stderr, "Failed to converge after %d iterations", MAX_ITERS);
			break;
		}

		cpl_mult(A, x1);
		λ_new = cpl_vector_inner(x1, y) / cpl_vector_inner(x0, y);
		residue = sabs(λ_new - λ_old);

		cpl_vector_normalize(x1);
		cpl_overwrite(x0, x1);
		λ_old = λ_new;
	}

	cpl_free(x1);
	cpl_free(y);

	return λ_new;
}
