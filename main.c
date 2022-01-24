#include "cpl_defines.h"
#include "cpl_includes.h"

#include "cpl.h"

int main() {
	cpl_tensor *A = cpl_matrix_alloc(4, 4);
	for (int i = 0; i < cpl_tensor_length(A); ++i)
		A->array[i] = i;
	cpl_matrix_print(A);

	cpl_tensor *id = cpl_matrix_id(4);
	cpl_matrix_print(id);

	cpl_tensor *B = cpl_matrix_mult(A, id);
	cpl_matrix_print(B);

	cpl_tensor_free(A);
	cpl_tensor_free(B);
	cpl_tensor_free(id);

	return 0;
}
