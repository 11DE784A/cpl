#include <stdlib.h>

#define scalar double

/* Tuples */
typedef struct {
	int length;
	int *array;
} cpl_tuple;

int cpl_tuple_get(cpl_tuple*, int);
int cpl_tuple_set(cpl_tuple*, int, int);

cpl_tuple* cpl_tuple_alloc(int);
cpl_tuple* cpl_tuple_salloc(int, ...);
cpl_tuple* cpl_tuple_ralloc(int, int);
cpl_tuple* cpl_tuple_zalloc(int);
cpl_tuple* cpl_tuple_1alloc(int);
void cpl_tuple_free(cpl_tuple*);

int cpl_tuple_length(cpl_tuple*);
void cpl_tuple_print(cpl_tuple*);
int cpl_tuple_mult(cpl_tuple*, int);
int cpl_tuple_equal(cpl_tuple*, cpl_tuple*);
cpl_tuple* cpl_tuple_copy(cpl_tuple*);

/* Tensors */
typedef struct {
	cpl_tuple *shape;
	scalar *array;
} cpl_tensor;

cpl_tuple* cpl_tensor_shape(cpl_tensor*);
int cpl_tensor_rank(cpl_tensor*);
int cpl_tensor_length(cpl_tensor*);
scalar cpl_tensor_get(cpl_tensor*, ...);
scalar cpl_tensor_set(cpl_tensor*, ...);
cpl_tensor* cpl_tensor_add(cpl_tensor*, cpl_tensor*);
scalar cpl_tensor_hadamard(cpl_tensor*, cpl_tensor*);

cpl_tensor* cpl_tensor_talloc(cpl_tuple*);
cpl_tensor* cpl_tensor_alloc(int, ...);
void cpl_tensor_free(cpl_tensor*);

cpl_tensor* cpl_vector_alloc(int);
scalar cpl_vector_get(cpl_tensor*, int);
scalar cpl_vector_set(cpl_tensor*, int, scalar);
cpl_tensor* cpl_vector_add(cpl_tensor*, cpl_tensor*);
scalar cpl_vector_hadamard(cpl_tensor*, cpl_tensor*);

cpl_tensor* cpl_matrix_alloc(int, int);
cpl_tensor* cpl_matrix_zalloc(int, int);
cpl_tensor* cpl_matrix_ialloc(int);
int cpl_matrix_rows(cpl_tensor*);
int cpl_matrix_cols(cpl_tensor*);
scalar cpl_matrix_get(cpl_tensor*, int, int);
scalar cpl_matrix_set(cpl_tensor*, int, int, scalar);
void cpl_matrix_print(cpl_tensor*);
cpl_tensor* cpl_matrix_add(cpl_tensor*, cpl_tensor*);
scalar cpl_matrix_hadamard(cpl_tensor*, cpl_tensor*);
cpl_tensor* cpl_matrix_mult(cpl_tensor*, cpl_tensor*);

