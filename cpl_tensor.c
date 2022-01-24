#include "cpl_defines.h"
#include "cpl_includes.h"

#include "cpl.h"

cpl_tuple *cpl_tensor_shape(cpl_tensor *T) {
	return T->shape;
}

int cpl_tensor_rank(cpl_tensor *T) {
	return cpl_tuple_length(cpl_tensor_shape(T));
}

int cpl_tensor_length(cpl_tensor *T) {
	return cpl_tuple_pdt(cpl_tensor_shape(T), 1, cpl_tensor_rank(T));
}

int cpl_tensor_isequal(cpl_tensor *S, cpl_tensor *T) {
	if (!cpl_tuple_isequal(cpl_tensor_shape(S), cpl_tensor_shape(T)))
		return 0;

	for (int i = 0; i < cpl_tensor_length(S); ++i) {
		if (S->array[i] != T->array[i])
			return 0;
	}

	return 1;
}

void cpl_tensor_free(cpl_tensor* T) {
	cpl_tuple_free(T->shape);
	free(T->array);
	free(T);
}

cpl_tensor *cpl_tensor_alloc_shape(cpl_tuple *shape) {
	if (!cpl_tuple_ispos(shape)) {
		fprintf(stderr, "ERROR: Invalid shape for tensor. All elements of the shape tuple must be postive integers.\nAborting...\n");
		exit(EXIT_FAILURE);
	}

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
	for (int i = 1; i <= rank; ++i) {
		int ith_dim = va_arg(ap, int);

		if (ith_dim <= 0) {
			cpl_tuple_free(shape);
			fprintf(stderr, "ERROR: Invalid shape for tensor. All elements of the shape tuple must be postive integers.\nAborting...\n");
			exit(EXIT_FAILURE);
		}

		cpl_tuple_set(shape, i, ith_dim);
	}
	va_end(ap);

	return cpl_tensor_alloc_shape(shape);
}

cpl_tensor *cpl_tensor_copy(cpl_tensor *T) {
	cpl_tensor *S = cpl_tensor_alloc_shape(cpl_tuple_copy(cpl_tensor_shape(T)));
	for (int i = 0; i < cpl_tensor_length(T); ++i)
		S->array[i] = T->array[i];

	return S;
}

void cpl_tensor_reshape(cpl_tensor *T, cpl_tuple *new_shape) {
	if (cpl_tensor_length(T) 
			!= cpl_tuple_pdt(new_shape, 1, cpl_tuple_length(new_shape))) {
		fprintf(stderr, "ERROR: Cannot reshape, size mismatch.");
		exit(EXIT_FAILURE);
	}
	cpl_tuple_free(T->shape);
	T->shape = new_shape;
}

scalar cpl_tensor_get(cpl_tensor *T, ...) {
	va_list ap;
	va_start(ap, T);
	int pos = 0;
	for (int i = 1; i <= cpl_tensor_rank(T); ++i) {
		int steps = va_arg(ap, int) - 1;
		int stride = cpl_tuple_pdt(cpl_tensor_shape(T), 1, i - 1);
		pos += steps * stride;

		if (steps < 0 || steps > cpl_tuple_get(cpl_tensor_shape(T), i) - 1) {
			fprintf(stderr, "ERROR: Tensor index out of bounds.\nAborting...\n");
			exit(EXIT_FAILURE);
		}
	}

	va_end(ap);

	return T->array[pos];
}

scalar cpl_tensor_set(cpl_tensor *T, ...) {
	va_list ap;
	va_start(ap, T);
	int pos = 0;
	for (int i = 1; i <= cpl_tensor_rank(T); ++i) {
		int steps = va_arg(ap, int) - 1;
		int stride = cpl_tuple_pdt(cpl_tensor_shape(T), 1, i - 1);
		pos += steps * stride;

		if (steps < 0 || steps > cpl_tuple_get(cpl_tensor_shape(T), i) - 1) {
			fprintf(stderr, "ERROR: Tensor index out of bounds.\nAborting...\n");
			exit(EXIT_FAILURE);
		}
	}

	T->array[pos] = va_arg(ap, scalar);
	va_end(ap);

	return T->array[pos];
}

void cpl_tensor_set_all(cpl_tensor *T, scalar x) {
	for (int i = 0; i < cpl_tensor_length(T); ++i)
		T->array[i] = x;
}

void cpl_tensor_scale(cpl_tensor *T, scalar c) {
	for (int i = 0; i < cpl_tensor_length(T); ++i)
		T->array[i] = c * T->array[i];
}

cpl_tensor *cpl_tensor_add(cpl_tensor *R, cpl_tensor *S) {
	if (!cpl_tuple_isequal(cpl_tensor_shape(R), cpl_tensor_shape(S))) {
		fprintf(stderr, "ERROR: Cannot add tensors of different shapes.\nAborting...\n");
		exit(EXIT_FAILURE);
	}

	cpl_tensor *T;
	T = cpl_tensor_alloc_shape(cpl_tuple_copy(cpl_tensor_shape(R)));
	for (int i = 0; i < cpl_tensor_length(R); ++i)
		T->array[i] = R->array[i] + S->array[i];

	return T;
}

cpl_tensor* cpl_tensor_hadamard(cpl_tensor *S, cpl_tensor *T) {
	if (!cpl_tuple_isequal(cpl_tensor_shape(S), cpl_tensor_shape(T))) {
		fprintf(stderr, "ERROR: Cannot Hadamard tensors of different shapes.\nAborting...\n");
		exit(EXIT_FAILURE);
	}

	cpl_tensor* R = cpl_tensor_alloc_shape(cpl_tensor_shape(S));
	for (int i = 0; i < cpl_tensor_length(S); ++i)
		R->array[i] = S->array[i] * T->array[i];

	return R;
}

cpl_tensor *cpl_vector_alloc(int dim) {
	return cpl_tensor_alloc(1, dim);
}

int cpl_vector_dim(cpl_tensor *v) {
	if (cpl_tensor_rank(v) > 1) {
		fprintf(stderr, "ERROR: Argument passed to `cpl_vector_dim` is not a vector.\nAborting...\n");
		exit(EXIT_FAILURE);
	}

	return cpl_tensor_length(v);
}

int cpl_vector_swap(cpl_tensor *v, int i, int j) {
	if (i == j) 
		return 0;

	scalar vi = cpl_tensor_get(v, i);
	cpl_tensor_set(v, i, cpl_tensor_get(v, j));
	cpl_tensor_set(v, j, vi);
	return 1;
}

scalar cpl_vector_dot(cpl_tensor *v, cpl_tensor *w) {
	if (cpl_vector_dim(v) != cpl_vector_dim(w)) {
		fprintf(stderr, "ERROR: Cannot dot vectors of different dimensions.\nAborting...\n");
		exit(EXIT_FAILURE);
	}

	scalar dot = 0;
	for (int i = 1; i <= cpl_vector_dim(v); ++i)
		dot += cpl_tensor_get(v, i) * cpl_tensor_get(w, i);

	return dot;
}

cpl_tensor *cpl_matrix_alloc(int rows, int cols) {
	return cpl_tensor_alloc(2, rows, cols);
}

cpl_tensor *cpl_matrix_id(int dim) {
	cpl_tensor *Id = cpl_matrix_alloc(dim, dim);
	for (int i = 1; i <= dim; ++i) {
		for (int j = 1; j <= dim; ++j) {
			cpl_tensor_set(Id, i, j, (i == j ? 1. : 0.));
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

cpl_tensor *cpl_matrix_get_col(cpl_tensor *A, int i) {
	cpl_tensor *col = cpl_vector_alloc(cpl_matrix_rows(A));
	for (int j = 1; j < cpl_vector_dim(col); ++j) 
		cpl_tensor_set(col, j, cpl_tensor_get(A, j, i));

	return col;
}

cpl_tensor *cpl_matrix_set_col(cpl_tensor *A, int i, cpl_tensor *v) {
	for (int j = 1; j <= cpl_matrix_rows(A); ++j)
		cpl_tensor_set(A, j, i, cpl_tensor_get(v, j));

	return v;
}

int cpl_matrix_swap_cols(cpl_tensor *A, int i, int j) {
	if (i == j) 
		return 0;

	for (int k = 1; k <= cpl_matrix_rows(A); ++k) {
		scalar Aki = cpl_tensor_get(A, k, i);
		cpl_tensor_set(A, k, i, cpl_tensor_get(A, k, j));
		cpl_tensor_set(A, k, j, Aki);
	}

	return 1;
}

void cpl_matrix_scale_col(cpl_tensor *A, int i, scalar c) {
	for (int k = 1; k <= cpl_matrix_rows(A); ++k)
		cpl_tensor_set(A, k, i, c * cpl_tensor_get(A, k, i));
}

void cpl_matrix_add_to_col(cpl_tensor *A, int j, cpl_tensor *v) {
	for (int i = 1; i <= cpl_matrix_rows(A); ++i)
		cpl_tensor_set(A, i, j, cpl_tensor_get(A, i, j) + cpl_tensor_get(v, i));
}

cpl_tensor *cpl_matrix_get_row(cpl_tensor *A, int i) {
	cpl_tensor *row = cpl_vector_alloc(cpl_matrix_cols(A));
	for (int j = 1; j <= cpl_vector_dim(row); ++j)
		cpl_tensor_set(row, j, cpl_tensor_get(A, i, j));

	return row;
}

cpl_tensor *cpl_matrix_set_row(cpl_tensor *A, int i, cpl_tensor *v) {
	for (int j = 1; j <= cpl_matrix_cols(A); ++j)
		cpl_tensor_set(A, i, j, cpl_tensor_get(v, j));

	return v;
}

int cpl_matrix_swap_rows(cpl_tensor *A, int i, int j) {
	if (i == j) 
		return 0;

	for (int k = 1; k <= cpl_matrix_cols(A); ++k) {
		scalar Aik = cpl_tensor_get(A, i, k);
		cpl_tensor_set(A, i, k, cpl_tensor_get(A, j, k));
		cpl_tensor_set(A, j, k, Aik);
	}

	return 1;
}

void cpl_matrix_scale_row(cpl_tensor *A, int i, scalar c) {
	for (int k = 1; k <= cpl_matrix_cols(A); ++k)
		cpl_tensor_set(A, i, k, c * cpl_tensor_get(A, i, k));
}

void cpl_matrix_add_to_row(cpl_tensor *A, int i, cpl_tensor *v) {
	for (int j = 1; j <= cpl_matrix_cols(A); ++j)
		cpl_tensor_set(A, i, j, cpl_tensor_get(A, i, j) + cpl_tensor_get(v, j));
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
					Cij += cpl_tensor_get(A, i, k) * cpl_tensor_get(B, k, j);
				}
				cpl_tensor_set(C, i, j, Cij);
			}
		}

	} else if (cpl_tensor_rank(A) == 2 && cpl_tensor_rank(B) == 1) {
		int Cdim = cpl_matrix_rows(A);
		C = cpl_vector_alloc(Cdim);
		for (int i = 1; i <= Cdim; ++i) {
			scalar Ci = 0;
			for (int k = 1; k <= cpl_matrix_cols(A); ++k) {
				Ci += cpl_tensor_get(A, i, k) * cpl_tensor_get(B, k);
			}
			cpl_tensor_set(C, i, Ci);
		}
	}

	return C;
}

void cpl_matrix_print(cpl_tensor *M) {
	printf("\n");
	if (cpl_tensor_rank(M) == 1) {
		printf("%d-dimensional Vector:\n", cpl_vector_dim(M));
		for (int i = 1; i <= cpl_vector_dim(M); ++i)
			printf(" % .3e\n", cpl_tensor_get(M, i));
	} else if (cpl_tensor_rank(M) == 2) {
		printf("%dx%d Matrix:\n", cpl_matrix_rows(M), cpl_matrix_cols(M));
		for (int i = 1; i <= cpl_matrix_rows(M); ++i) {
			printf(" ");
			for (int j = 1; j <= cpl_matrix_cols(M); ++j) {
				printf("% .3e  ", cpl_tensor_get(M, i, j));
			}
			printf("\n");
		}
	}
}
