/* Scalar Tuple */
typedef struct {
	scalar first;
	scalar second;
	scalar third;
	scalar fourth;
	scalar fifth;
	scalar sixth;
} cpl_tuple;

/* Blocks */

typedef struct {
	int size;
	scalar *array;
} cpl_block;

cpl_block *cpl_block_alloc(int);
cpl_block *cpl_block_calloc(int);
cpl_block *cpl_block_realloc(cpl_block*, int);
cpl_block *cpl_block_copy(cpl_block*);
void cpl_block_free(cpl_block*);

/* Vectors */

typedef struct {
	cpl_block *block;
	int dim;
} cpl_vector;

cpl_vector *cpl_vector_alloc(int);
cpl_vector *cpl_vector_calloc(int);
cpl_vector *cpl_vector_copy(cpl_vector*);
cpl_vector *cpl_vector_like(cpl_vector*);
void cpl_vector_free(cpl_vector*);

cpl_vector *cpl_vector_hot(int);
cpl_vector *cpl_vector_1hot(int, int);

int cpl_vector_size(cpl_vector*);
int cpl_vector_dim(cpl_vector*);

int cpl_vector_isequal(cpl_vector*, cpl_vector*);

scalar cpl_vector_get(cpl_vector*, int);
scalar cpl_vector_set(cpl_vector*, int, scalar);
void cpl_vector_set_all(cpl_vector*, scalar);
void cpl_vector_build(cpl_vector*, ...);
void cpl_vector_overwrite(cpl_vector*, cpl_vector*);

scalar cpl_vector_inner(cpl_vector*, cpl_vector*);
scalar cpl_vector_l2norm(cpl_vector*);
void cpl_vector_normalize(cpl_vector*);
scalar cpl_vector_l2dist(cpl_vector*, cpl_vector*);

cpl_vector *cpl_vector_scale(cpl_vector*, scalar);

cpl_vector *cpl_vector_push(cpl_vector*, scalar);
scalar cpl_vector_pop(cpl_vector*);

void cpl_vector_print(cpl_vector*);

/* Matrices */

typedef struct {
	int rows, cols;
	cpl_block *block;
} cpl_matrix;

cpl_matrix *cpl_matrix_alloc(int, int);
cpl_matrix *cpl_matrix_calloc(int, int);
cpl_matrix *cpl_matrix_copy(cpl_matrix*);
void cpl_matrix_free(cpl_matrix*);

cpl_matrix *cpl_matrix_id(int);

scalar cpl_matrix_get(cpl_matrix*, int, int);
scalar cpl_matrix_set(cpl_matrix*, int, int, scalar);
void cpl_matrix_set_all(cpl_matrix*, scalar);
void cpl_matrix_build(cpl_matrix*, ...);
void cpl_matrix_overwrite(cpl_matrix*, cpl_matrix*);
scalar cpl_matrix_l2dist(cpl_matrix*, cpl_matrix*);

int cpl_matrix_rows(cpl_matrix*);
int cpl_matrix_cols(cpl_matrix*);

int cpl_matrix_issquare(cpl_matrix*);
int cpl_matrix_isequal(cpl_matrix*, cpl_matrix*);

cpl_matrix *cpl_matrix_scale(cpl_matrix*, scalar);

void cpl_matrix_get_row(cpl_matrix*, int, cpl_vector*);
void cpl_matrix_scale_row(cpl_matrix*, int, scalar);
int cpl_matrix_swap_rows(cpl_matrix*, int ,int);
void cpl_matrix_add_to_row(cpl_matrix*, int, cpl_vector*);

void cpl_matrix_get_col(cpl_matrix*, int, cpl_vector*);
cpl_vector *cpl_matrix_flatten(cpl_matrix*);

void cpl_matrix_print(cpl_matrix*);

/* Matrix algebra */

cpl_vector *cpl_vector_add(cpl_vector*, cpl_vector*);
cpl_matrix *cpl_matrix_add(cpl_matrix*, cpl_matrix*);

scalar cpl_matrix_trace(cpl_matrix*);
scalar cpl_matrix_frobnorm(cpl_matrix*);
cpl_matrix *cpl_matrix_adjoint(cpl_matrix*);
cpl_vector *cpl_mvector_mult(cpl_matrix*, cpl_vector*);
cpl_matrix *cpl_mmatrix_mult(cpl_matrix*, cpl_matrix*);
void cpl_mvector_mult_overwrite(cpl_matrix*, cpl_vector*, cpl_vector*);
cpl_vector *cpl_mvector_mult_alloc(cpl_matrix*, cpl_vector*);
cpl_matrix *cpl_mmatrix_mult_alloc(cpl_matrix*, cpl_matrix*);

/* For matrices generated on the fly */

scalar cpl_func_get(scalar (*)(int, int), int, int);
cpl_vector *cpl_mvfly_mult_alloc(scalar (*)(int, int), cpl_vector*);
cpl_vector *cpl_mvfly_mult_overwrite(scalar (*)(int, int), cpl_vector*, cpl_vector*);

cpl_matrix *cpl_matrix_loadtxt(char *, int, int, int);
cpl_vector *cpl_vector_loadtxt(char *, int, int, int);

/* Generics */

#define cpl_get(X, ...) _Generic((X), \
	cpl_vector*: cpl_vector_get, \
	cpl_matrix*: cpl_matrix_get, \
	scalar (*)(int, int): cpl_func_get \
)(X, __VA_ARGS__)

#define cpl_set(X, ...) _Generic((X), \
	cpl_vector*: cpl_vector_set, \
	cpl_matrix*: cpl_matrix_set  \
)(X, __VA_ARGS__)

#define cpl_print(X) _Generic((X), \
	cpl_vector*: cpl_vector_print, \
	cpl_matrix*: cpl_matrix_print, \
	default: printf \
)(X)

#define cpl_free(X) _Generic ((X), \
	cpl_vector*: cpl_vector_free, \
	cpl_matrix*: cpl_matrix_free, \
	cpl_block*: cpl_block_free \
)(X)

#define cpl_overwrite(X, Y) _Generic ((X), \
	cpl_vector*: cpl_vector_overwrite, \
	cpl_matrix*: cpl_matrix_overwrite \
)(X, Y)

#define cpl_scale(X, c) _Generic ((X), \
	cpl_vector*: cpl_vector_scale, \
	cpl_matrix*: cpl_matrix_scale \
)(X, c)

#define cpl_add(X, Y) _Generic ((X), \
	cpl_vector*: cpl_vector_add, \
	cpl_matrix*: cpl_matrix_add \
)(X, Y);

#define cpl_mult(X, Y) _Generic ((Y), \
	cpl_vector*: cpl_mvector_mult, \
	cpl_matrix*: cpl_mmatrix_mult \
)(X, Y)

#define cpl_mult_alloc(X, Y) _Generic((Y), \
	cpl_vector*: cpl_mvector_mult_alloc, \
	cpl_matrix*: cpl_mmatrix_mult_alloc \
)(X, Y)

#define cpl_mult_overwrite(X, Y, Z) _Generic((Y), \
	cpl_vector*: cpl_mvector_mult_overwrite \
)(X, Y, Z)

