#include "cpl_includes.h"
#include "cpl_defines.h"

#include "cpl_commons.h"
#include "cpl_arrays.h"

/* Blocks */

cpl_block *cpl_block_alloc(int size) {
	cpl_block *block = malloc(sizeof(*block));

	block->array = malloc(sizeof(scalar) * size);
	block->size = size;

	return block;
}

cpl_block *cpl_block_calloc(int size) {
	cpl_block *block = malloc(sizeof(*block));

	block->array = calloc(size, sizeof(scalar));
	block->size = size;

	return block;
}

cpl_block *cpl_block_realloc(cpl_block *block, int size) {
	block->array = realloc(block->array, sizeof(scalar) * size);
	block->size = size;

	return block;
}

void cpl_block_free(cpl_block *block) {
	free(block->array);
	free(block);
}

cpl_block *cpl_block_copy(cpl_block *block) {
	cpl_block *block2 = cpl_block_alloc(block->size);
	for (int i = 0; i < block->size; ++i)
		block2->array[i] = block->array[i];

	return block2;
}

/* Vectors */

cpl_vector *cpl_vector_alloc(int dim) {
	cpl_vector *v = malloc(sizeof(*v));
	v->block = cpl_block_alloc(dim);
	v->dim = dim;
	return v;
}

cpl_vector *cpl_vector_calloc(int dim) {
	cpl_vector *v = malloc(sizeof(*v));
	v->block = cpl_block_calloc(dim);
	v->dim = dim;
	return v;
}

void cpl_vector_free(cpl_vector *v) {
	cpl_block_free(v->block);
	free(v);
}

cpl_vector *cpl_vector_copy(cpl_vector *v) {
	cpl_vector *u = malloc(sizeof(*u));
	u->dim = v->dim;
	u->block = cpl_block_copy(v->block);
	return u;
}

cpl_vector *cpl_vector_like(cpl_vector *v) {
	cpl_vector *u = cpl_vector_alloc(cpl_vector_dim(v));
	return u;
}

cpl_vector *cpl_vector_hot(int dim) {
	cpl_vector *v = cpl_vector_alloc(dim);
	for (int i = 1; i <= dim; ++i)
		cpl_set(v, i, 1.0);

	return v;
}

cpl_vector *cpl_vector_1hot(int dim, int i) {
	cpl_vector *v = cpl_vector_calloc(dim);
	cpl_vector_set(v, i, 1.0);

	return v;
}

int cpl_vector_dim(cpl_vector *v) {
	return v->dim;
}

int cpl_vector_size(cpl_vector *v) {
	return v->block->size;
}

int cpl_vector_isequal(cpl_vector *u, cpl_vector *v) {
	if (cpl_vector_dim(u) != cpl_vector_dim(v)) return 0;

	int size = cpl_vector_dim(u);
	for (int i = 0; i < size; ++i) {
		if (sabs(u->block->array[i] - v->block->array[i]) > TOL)
			return 0;
	}

	return 1;
}

scalar cpl_vector_get(cpl_vector *v, int i) {
	cpl_check(1 <= i && i <= cpl_vector_dim(v), 
			  "Vector index out of bounds");
	return v->block->array[i - 1];
}

scalar cpl_vector_set(cpl_vector *v, int i, scalar vi) {
	cpl_check(1 <= i && i <= cpl_vector_dim(v),
			  "Vector index out of bounds");

	v->block->array[i - 1] = vi;
	return vi;
}

void cpl_vector_set_all(cpl_vector *v, scalar c) {
	for (int i = 0; i < v->block->size; ++i)
		v->block->array[i] = c;
}

void cpl_vector_build(cpl_vector *v, ...) {
	va_list ap;
	va_start(ap, v);
	for (int i = 0; i < cpl_vector_dim(v); ++i)
		v->block->array[i] = va_arg(ap, scalar);

	va_end(ap);
}

scalar cpl_vector_inner(cpl_vector *u, cpl_vector *v) {
	scalar Σ = 0;
	for (int i = 1; i <= cpl_vector_dim(u); ++i)
		Σ += cpl_get(u, i) * cpl_get(v, i);

	return Σ;
}

scalar cpl_vector_l2norm(cpl_vector *v) {
	scalar Σ = 0;
	for (int i = 1; i <= cpl_vector_dim(v); ++i)
		Σ += pow(cpl_get(v, i), 2);

	return sqrt(Σ);
}

void cpl_vector_normalize(cpl_vector *v) {
	scalar norm = cpl_vector_l2norm(v);
	for (int i = 1; i <= cpl_vector_dim(v); ++i)
		cpl_set(v, i, cpl_get(v, i) / norm);
}

scalar cpl_vector_l2dist(cpl_vector *u, cpl_vector *v) {
	cpl_check(cpl_vector_dim(u) == cpl_vector_dim(v),
			  "Size mismatch for computing distance");

	scalar Σ = 0;
	for (int i = 1; i <= cpl_vector_dim(u); ++i)
		Σ += pow(cpl_get(u, i) - cpl_get(v, i), 2);

	return sqrt(Σ);
}

cpl_vector *cpl_vector_scale(cpl_vector *v, scalar c) {
	for (int i = 1; i <= cpl_vector_dim(v); ++i)
		cpl_vector_set(v, i, c * cpl_vector_get(v, i));

	return v;
}

cpl_vector *cpl_vector_push(cpl_vector *v, scalar vn) {
	cpl_check(cpl_vector_dim(v) <= cpl_vector_size(v),
			  "Vector dimension is greater than space allocated");

	if (cpl_vector_dim(v) == cpl_vector_size(v)) {
		v->block = cpl_block_realloc(v->block, cpl_vector_size(v) + 2);
	} 

	v->dim += 1;
	cpl_vector_set(v, cpl_vector_dim(v), vn);

	return v;
}

scalar cpl_vector_pop(cpl_vector *v) {
	cpl_check(cpl_vector_dim(v) <= cpl_vector_size(v),
			  "Vector dimension is greater than space allocated");

	if (2 * cpl_vector_dim(v) < cpl_vector_size(v))
		v->block = cpl_block_realloc(v->block, cpl_vector_dim(v));

	scalar vn = cpl_get(v, v->dim);
	v->dim -= 1;

	return vn;
}

void cpl_vector_overwrite(cpl_vector *u, cpl_vector *v) {
	cpl_check(cpl_vector_dim(u) == cpl_vector_dim(v),
			  "Size mismatch when overwriting vector");

	for (int i = 1; i <= cpl_vector_dim(u); ++i)
		cpl_set(u, i, cpl_get(v, i));
}

void cpl_vector_print(cpl_vector *v) {
	printf("%d-dimensional vector (" xstr(scalar) ")\n", cpl_vector_dim(v));
	for (int i = 1; i <= cpl_vector_dim(v); ++i)
		printf(" " xstr(sfmt) "\n", cpl_vector_get(v, i));
	printf("\n");
}

/* Matrices */

cpl_matrix *cpl_matrix_alloc(int rows, int cols) {
	cpl_matrix *M = malloc(sizeof(*M));
	M->rows = rows;
	M->cols = cols;
	M->block = cpl_block_alloc(rows * cols);
	return M;
}

cpl_matrix *cpl_matrix_calloc(int rows, int cols) {
	cpl_matrix *M = malloc(sizeof(*M));
	M->rows = rows;
	M->cols = cols;
	M->block = cpl_block_calloc(rows * cols);
	return M;
}

void cpl_matrix_free(cpl_matrix *M) {
	cpl_block_free(M->block);
	free(M);
}

cpl_matrix *cpl_matrix_copy(cpl_matrix *M) {
	cpl_matrix *N = malloc(sizeof(*N));
	N->rows = M->rows;
	N->cols = M->cols;
	N->block = cpl_block_copy(M->block);
	return N;
}

cpl_matrix *cpl_matrix_id(int dim) {
	cpl_matrix *id = cpl_matrix_calloc(dim, dim);
	for (int i = 1; i <= dim; ++i)
		cpl_matrix_set(id, i, i, 1.0);

	return id;
}

int cpl_matrix_rows(cpl_matrix *M) {
	return M->rows;
}

int cpl_matrix_cols(cpl_matrix *M) {
	return M->cols;
}

int cpl_matrix_issquare(cpl_matrix *M) {
	return cpl_matrix_rows(M) == cpl_matrix_cols(M);
}

int cpl_matrix_isequal(cpl_matrix *M, cpl_matrix *N) {
	if (cpl_matrix_rows(M) != cpl_matrix_rows(N)
			|| cpl_matrix_cols(M) != cpl_matrix_cols(N)) return 0;

	int size = cpl_matrix_rows(M) * cpl_matrix_cols(M);
	for (int i = 0; i < size; ++i) {
		if (sabs(M->block->array[i] - N->block->array[i]) > TOL)
			return 0;
	}

	return 1;
}

scalar cpl_matrix_get(cpl_matrix *M, int i, int j) {
	cpl_check(1 <= i && i <= cpl_matrix_rows(M) && 1 <= j && j <= cpl_matrix_cols(M),
			  "Matrix index out of bounds");

	return M->block->array[(i - 1) * cpl_matrix_cols(M) + (j - 1)];
}

scalar cpl_matrix_set(cpl_matrix *M, int i, int j, scalar Mij) {
	cpl_check(1 <= i && i <= cpl_matrix_rows(M) && 1 <= j && j <= cpl_matrix_cols(M),
			  "Matrix index out of bounds");

	M->block->array[(i - 1) * cpl_matrix_cols(M) + (j - 1)] = Mij;
	return Mij;
}

void cpl_matrix_set_all(cpl_matrix *M, scalar c) {
	for (int i = 0; i < M->block->size; ++i)
		M->block->array[i] = c;
}

void cpl_matrix_build(cpl_matrix *M, ...) {
	va_list ap;
	va_start(ap, M);
	for (int i = 0; i < cpl_matrix_rows(M) * cpl_matrix_cols(M); ++i)
		M->block->array[i] = va_arg(ap, scalar);

	va_end(ap);
}

scalar cpl_matrix_l2dist(cpl_matrix *M, cpl_matrix *N) {
	cpl_check(cpl_matrix_rows(M) == cpl_matrix_rows(N) 
				&& cpl_matrix_cols(M) == cpl_matrix_cols(N), 
			  "Size mismatch for computing distance");
	scalar Σ = 0;
	for (int i = 0; i < M->block->size; ++i)
		Σ += pow(M->block->array[i] - M->block->array[i], 2);

	return Σ;
}

void cpl_matrix_overwrite(cpl_matrix *M, cpl_matrix *N) {
	cpl_check(cpl_matrix_rows(M) == cpl_matrix_rows(N)
				&& cpl_matrix_cols(M) == cpl_matrix_cols(N),
			  "Size mismatch when overwriting matrix");

	for (int i = 0; i < M->block->size; ++i)
		M->block->array[i] = N->block->array[i];
}

cpl_matrix *cpl_matrix_scale(cpl_matrix *M, scalar c) {
	for (int i = 0; i < cpl_matrix_cols(M) * cpl_matrix_rows(M); ++i)
		M->block->array[i] = c * M->block->array[i];

	return M;
}

void cpl_matrix_print(cpl_matrix *M) {
	printf("%dx%d matrix (" xstr(scalar) ")\n", cpl_matrix_rows(M), 
												cpl_matrix_cols(M));
	for (int i = 1; i <= cpl_matrix_rows(M); ++i) {
		printf(" ");
		for (int j = 1; j <= cpl_matrix_cols(M); ++j) {
			printf(xstr(sfmt) "\t", cpl_matrix_get(M, i, j));
		}
		printf("\n");
	}

	printf("\n");
}

/* Matrix algebra */

scalar cpl_matrix_trace(cpl_matrix *M) {
	cpl_check(cpl_matrix_issquare(M), "Cannot take trace of non-square matrix");

	scalar TrM = 0;
	for (int i = 1; i <= cpl_matrix_rows(M); ++i)
		TrM += cpl_matrix_get(M, i, i);

	return TrM;
}

scalar cpl_matrix_frobnorm(cpl_matrix *M) {
	scalar Σ = 0;
	for (int i = 1; i <= cpl_matrix_rows(M); ++i) {
		for (int j = 1; j <= cpl_matrix_cols(M); ++j)
			Σ += pow(sabs(cpl_get(M, i, j)), 2);
	}

	return sqrt(Σ);
}

cpl_matrix *cpl_matrix_adjoint(cpl_matrix *M) {
	cpl_matrix *M_adj = cpl_matrix_alloc(cpl_matrix_cols(M), cpl_matrix_rows(M));

	for (int i = 1; i <= cpl_matrix_rows(M); ++i) {
		for (int j = 1; j <= cpl_matrix_cols(M); ++j) {
			cpl_set(M_adj, j, i, cpl_get(M, i, j));
		}
	}

	return M_adj;
}

cpl_vector *cpl_vector_add(cpl_vector *u, cpl_vector *v) {
	cpl_check(cpl_vector_dim(u) == cpl_vector_dim(v),
			  "Size mismatch while trying to add vectors");

	for (int i = 1; i <= cpl_vector_dim(u); ++i)
		cpl_set(v, i, cpl_get(u, i) + cpl_get(v, i));

	return v;
}

cpl_matrix *cpl_matrix_add(cpl_matrix *M, cpl_matrix *N) {
	cpl_check(cpl_matrix_rows(M) == cpl_matrix_rows(N) 
				&& cpl_matrix_cols(M) == cpl_matrix_cols(N),
			  "Size mismatch while trying to add matrices");
	for (int i = 1; i <= cpl_matrix_rows(N); ++i) {
		for (int j = 1; j <= cpl_matrix_cols(N); ++j) {
			cpl_set(N, i, j, cpl_get(M, i, j) + cpl_get(N, i, j));
		}
	}

	return N;
}

void cpl_mvector_mult_overwrite(cpl_matrix *M, cpl_vector *v, cpl_vector *x) {
	for (int i = 1; i <= cpl_vector_dim(x); ++i) {
		scalar xi = 0;
		for (int j = 1; j <= cpl_vector_dim(v); ++j)
			xi += cpl_get(M, i, j) * cpl_get(v, j);
		cpl_set(x, i, xi);
	}
}

cpl_vector *cpl_mvector_mult_alloc(cpl_matrix *M, cpl_vector *v) {
	cpl_vector *Mv = cpl_vector_alloc(cpl_matrix_rows(M));
	for (int i = 1; i <= cpl_vector_dim(Mv); ++i) {
		scalar Mvi = 0;
		for (int j = 1; j <= cpl_vector_dim(v); ++j)
			Mvi += cpl_get(M, i, j) * cpl_get(v, j);
		cpl_set(Mv, i, Mvi);
	}

	return Mv;
}

cpl_matrix *cpl_mmatrix_mult_alloc(cpl_matrix *M, cpl_matrix *N) {
	cpl_matrix *MN = cpl_matrix_alloc(cpl_matrix_rows(M), cpl_matrix_cols(N));
	for (int i = 1; i <= cpl_matrix_rows(MN); ++i) {
		for (int j = 1; j <= cpl_matrix_cols(MN); ++j) {
			scalar MNij = 0;
			for (int k = 1; k <= cpl_matrix_cols(M); ++k)
				MNij += cpl_get(M, i, k) * cpl_get(N, k, j);

			cpl_set(MN, i, j, MNij);
		}
	}

	return MN;
}

cpl_vector *cpl_mvector_mult(cpl_matrix *M, cpl_vector *v) {
	cpl_vector *v_copy = cpl_vector_copy(v);
	v->dim = M->rows;
	v->block = cpl_block_realloc(v->block, v->dim);

	for (int i = 1; i <= cpl_vector_dim(v); ++i) {
		scalar vi = 0;
		for (int j = 1; j <= cpl_vector_dim(v_copy); ++j)
			vi += cpl_matrix_get(M, i, j) * cpl_vector_get(v_copy, j);

		cpl_vector_set(v, i, vi);
	}

	cpl_free(v_copy);

	return v;
}

cpl_matrix *cpl_mmatrix_mult(cpl_matrix *M, cpl_matrix *N) {
	cpl_check(cpl_matrix_cols(M) == cpl_matrix_rows(N),
			  "Dimension mismatch in matrix multiplication");

	cpl_matrix *N_copy = cpl_matrix_copy(N);
	N->rows = cpl_matrix_rows(M);
	N->block = cpl_block_realloc(N->block, N->rows * N->cols);

	for (int i = 1; i <= cpl_matrix_rows(N); ++i) {
		for (int j = 1; j <= cpl_matrix_cols(N); ++j) {
			scalar Nij = 0;
			for (int k = 1; k <= cpl_matrix_cols(M); ++k)
				Nij += cpl_matrix_get(M, i, k) * cpl_matrix_get(N_copy, k, j);
			cpl_set(N, i, j, Nij);
		}
	}

	cpl_free(N_copy);

	return M;
}

/* Row & Column operations */

void cpl_matrix_get_row(cpl_matrix *M, int i, cpl_vector *v) {
	cpl_check(cpl_vector_dim(v) == cpl_matrix_cols(M),
			  "Size mismatch when getting row");

	for (int j = 1; j <= cpl_matrix_cols(M); ++j)
		cpl_vector_set(v, j, cpl_matrix_get(M, i, j));
}

void cpl_matrix_scale_row(cpl_matrix *M, int i, scalar c) {
	for (int j = 1; j <= cpl_matrix_cols(M); ++j)
		cpl_matrix_set(M, i, j, c * cpl_matrix_get(M, i, j));
}

int cpl_matrix_swap_rows(cpl_matrix *M, int i, int j) {
	if (i == j)
		return 0;

	for (int k = 1; k <= cpl_matrix_cols(M); ++k) {
		scalar Mik = cpl_matrix_get(M, i, k);
		cpl_matrix_set(M, i, k, cpl_matrix_get(M, j, k));
		cpl_matrix_set(M, j, k, Mik);
	}

	return 1;
}

void cpl_matrix_add_to_row(cpl_matrix *M, int i, cpl_vector *v) {
	for (int j = 1; j <= cpl_matrix_cols(M); ++j)
		cpl_matrix_set(M, i, j, cpl_vector_get(v, j) + cpl_matrix_get(M, i, j));
}

void cpl_matrix_get_col(cpl_matrix *M, int i, cpl_vector *v) {
	cpl_check(cpl_vector_dim(v) == cpl_matrix_rows(M),
			  "Size mismatch when getting column");

	for (int j = 1; j <= cpl_matrix_rows(M); ++j)
		cpl_vector_set(v, j, cpl_matrix_get(M, j, i));
}

cpl_vector *cpl_matrix_flatten(cpl_matrix *M) {
	int size = cpl_matrix_rows(M) * cpl_matrix_cols(M);
	cpl_vector *v = cpl_vector_alloc(size);
	for (int i = 0; i < size; ++i)
		v->block->array[i] = M->block->array[i];

	return v;
}

cpl_matrix *cpl_matrix_loadtxt(char *fpath, int start, int end, int cols) {
	int rows = end - start + 1;
	cpl_matrix *X = cpl_matrix_alloc(rows, cols);

	float x; /* Replace with scalar ASAP */
	char *item, line[MAXLINE];
	FILE *fp = fopen(fpath, "r");
	for (int i = 1; i <= end; ++i) {
		fgets(line, MAXLINE, fp);
		if (start <= i) {
			item = strtok(line, " \t");
			for (int j = 1; j <= cols; ++j) {
				sscanf(item, "%e", &x);
				cpl_set(X, i - start + 1, j, x);
				item = strtok(NULL, " \t");
			}
		}
	}

	return X;
}

cpl_vector *cpl_vector_loadtxt(char *fpath, int start, int end, int col) {
	int dim = end - start + 1;
	cpl_vector *v = cpl_vector_alloc(dim);

	float x;
	char *item, line[MAXLINE];
	FILE *fp = fopen(fpath, "r");
	for (int i = 1; i <= end; ++i) {
		fgets(line, MAXLINE, fp);
		if (start <= i) {
			item = strtok(line, " \t");
			for (int j = 0; j < col - 1; ++j)
				item = strtok(NULL, " \t");

			sscanf(item, "%e", &x);
			cpl_set(v, i - start + 1, x);
		}
	}

	return v;
}

/* For matrices generated on the fly */
scalar cpl_func_get(scalar (*func)(int, int), int i, int j) {
	return func(i, j);
}
cpl_vector *cpl_mvfly_mult_alloc(scalar (*fn)(int, int), cpl_vector *v) {
	int dim = cpl_vector_dim(v);

	cpl_vector *w = cpl_vector_alloc(dim);
	cpl_mvfly_mult_overwrite(fn, v, w);

	return w;
}

cpl_vector *cpl_mvfly_mult_overwrite(scalar (*fn)(int, int), cpl_vector *v, cpl_vector *w) {
	int dim = cpl_vector_dim(v);

	scalar Σ;
	for (int i = 1; i <= dim; ++i) {
		Σ = 0;
		for (int j = 1; j <= dim; ++j)
			Σ += fn(i, j) * cpl_get(v, j);
		cpl_set(w, i, Σ);
	}

	return w;
}
