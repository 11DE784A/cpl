#include "cpl_defines.h"
#include "cpl_includes.h"

#include "cpl_commons.h"

#include "cpl_arrays.h"
#include "cpl_transform.h"

cpl_vector *cpl_transform_dft(cpl_vector *v) {
	int N = cpl_vector_dim(v);
	cpl_vector *x = cpl_vector_alloc(N);

	for (int k = 1; k <= N; ++k) {
		for
	}
}
