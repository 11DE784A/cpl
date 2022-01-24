#include <stdlib.h>

#include "cpl_defines.h"

/* Tuples */
typedef struct {
	int length;
	int *array;
} cpl_tuple;

int cpl_tuple_get(cpl_tuple*, int);
int cpl_tuple_set(cpl_tuple*, int, int);
void cpl_tuple_set_all(cpl_tuple*, int);

cpl_tuple* cpl_tuple_alloc(int);
cpl_tuple* cpl_tuple_alloc_assign(int, ...);
cpl_tuple* cpl_tuple_copy(cpl_tuple*);
void cpl_tuple_free(cpl_tuple*);

int cpl_tuple_length(cpl_tuple*);
void cpl_tuple_print(cpl_tuple*);

int cpl_tuple_pdt(cpl_tuple*, int, int);
int cpl_tuple_isequal(cpl_tuple*, cpl_tuple*);
int cpl_tuple_ispos(cpl_tuple*);

/* Tensors */
typedef struct {
	cpl_tuple *shape;
	scalar *array;
} cpl_tensor;

cpl_tuple* cpl_tensor_shape(cpl_tensor*);
int cpl_tensor_rank(cpl_tensor*);
int cpl_tensor_length(cpl_tensor*);
void cpl_tensor_reshape(cpl_tensor*, cpl_tuple*);

scalar cpl_tensor_get(cpl_tensor*, ...);
scalar cpl_tensor_set(cpl_tensor*, ...);
void cpl_tensor_set_all(cpl_tensor*, scalar);

void cpl_tensor_scale(cpl_tensor*, scalar);
cpl_tensor* cpl_tensor_add(cpl_tensor*, cpl_tensor*);
cpl_tensor* cpl_tensor_hadamard(cpl_tensor*, cpl_tensor*);

cpl_tensor* cpl_tensor_alloc(int, ...);
cpl_tensor* cpl_tensor_alloc_shape(cpl_tuple*);
cpl_tensor* cpl_tensor_copy(cpl_tensor*);
void cpl_tensor_free(cpl_tensor*);

cpl_tensor* cpl_vector_alloc(int);

int cpl_vector_dim(cpl_tensor*);
int cpl_vector_swap(cpl_tensor*, int, int);

scalar cpl_vector_dot(cpl_tensor*, cpl_tensor*);

cpl_tensor* cpl_matrix_alloc(int, int);
cpl_tensor* cpl_matrix_id(int);

cpl_tensor* cpl_matrix_get_row(cpl_tensor*, int);
cpl_tensor* cpl_matrix_set_row(cpl_tensor*, int, cpl_tensor*);
int cpl_matrix_swap_rows(cpl_tensor*, int, int);
void cpl_matrix_scale_row(cpl_tensor*, int, scalar);
void cpl_matrix_add_to_row(cpl_tensor*, int, cpl_tensor*);

cpl_tensor* cpl_matrix_get_col(cpl_tensor*, int);
cpl_tensor* cpl_matrix_set_col(cpl_tensor*, int, cpl_tensor*);
int cpl_matrix_swap_cols(cpl_tensor*, int, int);
void cpl_matrix_scale_col(cpl_tensor*, int, scalar);
void cpl_matrix_add_to_col(cpl_tensor*, int, cpl_tensor*);


int cpl_matrix_rows(cpl_tensor*);
int cpl_matrix_cols(cpl_tensor*);

void cpl_matrix_print(cpl_tensor*);

cpl_tensor* cpl_matrix_mult(cpl_tensor*, cpl_tensor*);

/* Linear Algebra */
cpl_tensor* cpl_linalg_gaussjordan(cpl_tensor*, cpl_tensor*);
