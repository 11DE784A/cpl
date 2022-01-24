#include "cpl_defines.h"
#include "cpl_includes.h"

#include <time.h>

#include "cpl.h"

int main() {
	int dim = 100;
	cpl_tensor *A = cpl_matrix_alloc(dim, dim);
	for (int i = 1; i <= dim; ++i) {
		for (int j = 1; j <= dim; ++j) {
			cpl_tensor_set(A, i, j, 1.0 / (i+j-1));
		}
	}

	cpl_tensor *B = cpl_matrix_alloc(dim, dim);
	cpl_tensor_set_all(B, 1.0);

	double start = clock();
	cpl_tensor *AB = cpl_matrix_mult(A, B);
	double end = clock();

	double time_elapsed = ((double) (end - start)) / CLOCKS_PER_SEC;
	printf("Time to multiply two %dx%d matrices: %f sec", dim, dim, time_elapsed);

	cpl_tensor_free(AB);
	cpl_tensor_free(B);
	cpl_tensor_free(A);

	return 0;
}
