#include <stdio.h>
#include "cpl.h"

int main() {
	cpl_tensor *T = cpl_matrix_alloc(4, 5);
	for (int i = 0; i < 20; ++i) {
		T->array[i] = i;
	}
	cpl_matrix_set(T, 1, 1, -1e7);
	cpl_matrix_set(T, 4, 4, 1e-7);
	cpl_matrix_print(T);

	cpl_tensor *v = cpl_matrix_zalloc(1, 4);
	cpl_matrix_set(v, 1, 2, 1.);

	cpl_tensor *w = cpl_vector_alloc(5);
	cpl_vector_set(w, 2, 1);

	cpl_matrix_print(v);

	cpl_tensor *S = cpl_matrix_mult(v, T);

	cpl_matrix_print(S);

	cpl_tensor_free(S);
	cpl_tensor_free(T);
	cpl_tensor_free(v);
	cpl_tensor_free(w);

	return 0;
}
