#include "cpl_defines.h"
#include "cpl_includes.h"

#include "cpl.h"

int main() {
	int dim = 10;
	cpl_tensor *A = cpl_matrix_alloc(dim, dim);
	for (int i = 1; i <= dim; ++i) {
		for (int j = 1; j <= dim; ++j) {
			cpl_tensor_set(A, i, j, 1.0 / (i+j-1));
		}
	}

	cpl_tensor *B = cpl_matrix_alloc(dim-1, dim);
	cpl_tensor_set_all(B, 1.0);

	cpl_linalg_gaussjordan(A, B);

	cpl_tensor_free(B);
	cpl_tensor_free(A);

	return 0;
}
