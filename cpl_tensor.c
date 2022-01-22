#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>

#include "cpl.h"

#define scalar double

cpl_tuple *cpl_tensor_shape(cpl_tensor *T) {
	return T->shape;
}

int cpl_tensor_rank(cpl_tensor *T) {
	return cpl_tuple_length(cpl_tensor_shape(T));
}

int cpl_tensor_length(cpl_tensor *T) {
	return cpl_tuple_mult(cpl_tensor_shape(T), 1);
}

void cpl_tensor_free(cpl_tensor* T) {
	cpl_tuple_free(T->shape);
	free(T->array);
	free(T);
}

cpl_tensor *cpl_tensor_talloc(cpl_tuple *shape) {
	cpl_tensor *T = malloc(sizeof(cpl_tensor));
	T->shape = shape;

	scalar *array = malloc(cpl_tensor_length(T) * sizeof(scalar));
	T->array = array;

	return T;
}

cpl_tensor *cpl_tensor_alloc(int rank, ...) {
	cpl_tuple *shape = cpl_tuple_alloc(rank);

	va_list ap;
	va_start(ap, rank);
	for (int i = 1; i <= rank; ++i)
		cpl_tuple_set(shape, i, va_arg(ap, int));
	va_end(ap);

	return cpl_tensor_talloc(shape);
}

scalar cpl_tensor_get(cpl_tensor *T, ...) {
	va_list ap;
	va_start(ap, T);
	int pos = 0;
	for (int i = 1; i <= cpl_tensor_rank(T); ++i) {
		int steps = va_arg(ap, int) - 1;
		int stride = cpl_tuple_mult(cpl_tensor_shape(T), i + 1);
		pos += steps * stride;
	}

	return T->array[pos];
}

scalar cpl_tensor_set(cpl_tensor *T, ...) {
	va_list ap;
	va_start(ap, T);
	int pos = 0;
	for (int i = 1; i <= cpl_tensor_rank(T); ++i) {
		int steps = va_arg(ap, int) - 1;
		int stride = cpl_tuple_mult(cpl_tensor_shape(T), i + 1);
		pos += steps * stride;
	}

	T->array[pos] = va_arg(ap, scalar);

	return T->array[pos];
}

cpl_tensor *cpl_tensor_add(cpl_tensor *R, cpl_tensor *S) {
	cpl_tensor *T;
	if (cpl_tuple_equal(cpl_tensor_shape(R), cpl_tensor_shape(S))) {
		T = cpl_tensor_talloc(cpl_tuple_copy(cpl_tensor_shape(R)));
		for (int i = 0; i < cpl_tensor_length(R); ++i)
			T->array[i] = R->array[i] + S->array[i];
	}

	return T;
}

scalar cpl_tensor_hadamard(cpl_tensor *S, cpl_tensor *T) {
	scalar had = 0;
	if (cpl_tuple_equal(cpl_tensor_shape(S), cpl_tensor_shape(T))) {
		for (int i = 0; i < cpl_tensor_length(S); ++i) {
			had += S->array[i] + T->array[i];
		}
	}

	return had;
}

cpl_tensor *cpl_vector_alloc(int dim) {
	return cpl_tensor_alloc(1, dim);
}

scalar cpl_vector_get(cpl_tensor *v, int i) {
	return cpl_tensor_get(v, i);
}

scalar cpl_vector_set(cpl_tensor *v, int i, scalar vi) {
	return cpl_tensor_set(v, i, vi);
}

int cpl_vector_dim(cpl_tensor *v) {
	return cpl_tensor_length(v);
}

cpl_tensor *cpl_vector_add(cpl_tensor *v, cpl_tensor *w) {
	return cpl_tensor_add(v, w);
}

scalar cpl_vector_hadamard(cpl_tensor *v, cpl_tensor *w) {
	return cpl_tensor_hadamard(v, w);
}

scalar cpl_vector_dot(cpl_tensor *v, cpl_tensor *w) {
	scalar dot = 0;
	for (int i = 1; i <= cpl_vector_dim(v); ++i)
		dot += cpl_vector_get(v, i) * cpl_vector_get(w, i);

	return dot;
}

void cpl_vector_print(cpl_tensor *v) {
	cpl_matrix_print(v);
}

cpl_tensor *cpl_matrix_alloc(int rows, int cols) {
	return cpl_tensor_alloc(2, rows, cols);
}

cpl_tensor *cpl_matrix_zalloc(int rows, int cols) {
	cpl_tensor *O = cpl_tensor_alloc(2, rows, cols);
	for (int i = 1; i <= rows; ++i) {
		for (int j = 1; j <= cols; ++j) {
			cpl_matrix_set(O, i, j, 0.);
		}
	}

	return O;
}

cpl_tensor *cpl_matrix_ialloc(int dim) {
	cpl_tensor *Id = cpl_tensor_alloc(2, dim, dim);
	for (int i = 1; i <= dim; ++i) {
		for (int j = 1; j <= dim; ++j) {
			i == j ? cpl_matrix_set(Id, i, j, 1.) : cpl_matrix_set(Id, i, j, 0.);
		}
	}
	return Id;
}

int cpl_matrix_rows(cpl_tensor *M) {
	return cpl_tuple_get(cpl_tensor_shape(M), 1);
}

int cpl_matrix_cols(cpl_tensor *M) {
	return cpl_tuple_get(cpl_tensor_shape(M), 2);
}

scalar cpl_matrix_get(cpl_tensor *M, int i, int j) {
	return cpl_tensor_get(M, i, j);
}

scalar cpl_matrix_set(cpl_tensor *M, int i, int j, scalar Mij) {
	return cpl_tensor_set(M, i, j, Mij);
}

cpl_tensor *cpl_matrix_add(cpl_tensor *A, cpl_tensor *B) {
	return cpl_tensor_add(A, B);
}

scalar cpl_matrix_hadamard(cpl_tensor *A, cpl_tensor *B) {
	return cpl_tensor_hadamard(A, B);
}

cpl_tensor *cpl_matrix_mult(cpl_tensor *A, cpl_tensor *B) {
	cpl_tensor *C;

	if (cpl_tensor_rank(A) == 2 && cpl_tensor_rank(B) == 2) {
		int Crows = cpl_matrix_rows(A);
		int Ccols = cpl_matrix_cols(B);
		C = cpl_matrix_alloc(Crows, Ccols);
		for (int i = 1; i <= Crows; ++i) {
			for (int j = 1; j <= Ccols; ++j) {
				scalar Cij = 0;
				for (int k = 1; k <= cpl_matrix_cols(A); ++k) {
					Cij += cpl_matrix_get(A, i, k) * cpl_matrix_get(B, k, j);
				}
				cpl_matrix_set(C, i, j, Cij);
			}
		}

	} else if (cpl_tensor_rank(A) == 2 && cpl_tensor_rank(B) == 1) {
		int Cdim = cpl_matrix_rows(A);
		C = cpl_vector_alloc(Cdim);
		for (int i = 1; i <= Cdim; ++i) {
			scalar Ci = 0;
			for (int k = 1; k <= cpl_matrix_cols(A); ++k) {
				Ci += cpl_matrix_get(A, i, k) * cpl_vector_get(B, k);
			}
			cpl_vector_set(C, i, Ci);
		}
	}

	return C;
}

void cpl_matrix_print(cpl_tensor *M) {
	printf("\n");
	if (cpl_tensor_rank(M) == 1) {
		printf("%d-dimensional Vector:\n", cpl_vector_dim(M));
		for (int i = 1; i <= cpl_vector_dim(M); ++i)
			printf(" % .3e\n", cpl_vector_get(M, i));
	} else if (cpl_tensor_rank(M) == 2) {
		printf("%dx%d Matrix:\n", cpl_matrix_rows(M), cpl_matrix_cols(M));
		for (int i = 1; i <= cpl_matrix_rows(M); ++i) {
			printf(" ");
			for (int j = 1; j <= cpl_matrix_cols(M); ++j) {
				printf("% .3e  ", cpl_matrix_get(M, i, j));
			}
			printf("\n");
		}
	}
}
