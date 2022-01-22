#include <stdio.h>
#include <stdarg.h>
#include <stdlib.h>

double cpl_rand() {
	return (rand() % 1000) / 1000.0;
}

typedef struct {
	int size;
	double *array;
} cpl_vector;

cpl_vector *cpl_vector_alloc(int size) {
	cpl_vector *v = malloc(sizeof(cpl_vector));
	v->size = size;
	double *array = malloc(size * sizeof(double));
	v->array = array;

	return v;
}

void cpl_vector_free(cpl_vector *v) {
	free(v->array);
	free(v);
}

double cpl_vector_get(cpl_vector *v, int i) {
	return v->array[i-1];
}

double cpl_vector_set(cpl_vector *v, int i, double v_i) {
	v->array[i-1] = v_i;
	return v_i;
};

cpl_vector *cpl_vector_ralloc(int size) {
	cpl_vector *v = cpl_vector_alloc(size);
	for (int i = 1; i <= v-> size; ++i)
		cpl_vector_set(v, i, cpl_rand());

	return v;
}

cpl_vector *cpl_vector_zalloc(int size) {
	cpl_vector *v = cpl_vector_alloc(size);
	for (int i = 1; i <= v-> size; ++i)
		cpl_vector_set(v, i, 0);

	return v;
}

cpl_vector *cpl_vector_halloc(int size, int i) {
	cpl_vector *v = cpl_vector_zalloc(size);
	cpl_vector_set(v, i, 1);
	return v;
}
void cpl_vector_print(cpl_vector *v) {
	printf("%d-dimensional vector of doubles:\n", v->size);
	for (int i = 1; i <= v->size; ++i)
		printf(" %f\n", cpl_vector_get(v, i));
	printf("\n");
}

cpl_vector *cpl_vector_add(cpl_vector *u, cpl_vector *v) {
	cpl_vector *w = cpl_vector_alloc(u->size);
	for (int i = 1; i <= u->size; ++i)
		cpl_vector_set(w, i, cpl_vector_get(u, i) + cpl_vector_get(v, i));

	return w;
}

double cpl_vector_dot(cpl_vector *u, cpl_vector *v) {
	double dot_product = 0;
	for (int i = 1; i <= u.size; ++i)
		dot_product += cpl_vector_get(u, i) * cpl_vector_get(v, i);

	return dot_product;
}

typedef struct {
	int rows;
	int cols;
	double *array;
} cpl_matrix;

cpl_matrix *cpl_matrix_alloc(int rows, int cols) {
	cpl_matrix *M = malloc(sizeof(cpl_matrix));
	M->rows = rows;
	M->cols = cols;
	double *array = malloc(rows*cols * sizeof(double));
	M->array = array;

	return M;
}

int main() {
	cpl_vector *u, *v, *w;
	u = cpl_vector_zalloc(4);
	cpl_vector_print(u);
	v = cpl_vector_ralloc(4);
	cpl_vector_print(v);
	w = cpl_vector_add(u, v);
	cpl_vector_print(w);

	cpl_vector_free(u);
	cpl_vector_free(v);
	cpl_vector_free(w);
}
