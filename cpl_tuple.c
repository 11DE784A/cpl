#include "cpl_defines.h"
#include "cpl_includes.h"

#include "cpl.h"

int cpl_tuple_length(cpl_tuple *t) {
	return t->length;
}

int cpl_tuple_get(cpl_tuple *t, int j) {
	if (1 > j || j > cpl_tuple_length(t)) {
		fprintf(stderr, "ERROR: Tuple index out of bounds.\nAborting...\n");
		exit(EXIT_FAILURE);
	}

	return t->array[j-1];
}

int cpl_tuple_set(cpl_tuple *t, int j, int x) {
	if (1 > j || j > cpl_tuple_length(t)) {
		fprintf(stderr, "ERROR: Tuple index out of bounds.\nAborting...\n");
		exit(EXIT_FAILURE);
	}

	t->array[j-1] = x;
	return x;
}

void cpl_tuple_set_all(cpl_tuple *t, int c) {
	for (int i = 1; i <= cpl_tuple_length(t); ++i)
		cpl_tuple_set(t, i, c);
}

void cpl_tuple_free(cpl_tuple *t) {
	free(t->array);
	free(t);
}

cpl_tuple *cpl_tuple_alloc(int length) {
	if (length < 1) {
		fprintf(stderr, "ERROR: Tuple length has to be at least 1.\nAborting...\n");
		exit(EXIT_FAILURE);
	}

	cpl_tuple *t = malloc(sizeof(cpl_tuple));
	t->length = length;

	int *array = malloc(length * sizeof(int));
	t->array = array;

	return t;
}

cpl_tuple *cpl_tuple_alloc_assign(int length, ...) {
	if (length < 1) {
		fprintf(stderr, "ERROR: Tuple length has to be at least 1.\nAborting...\n");
		exit(EXIT_FAILURE);
	}

	cpl_tuple *t = cpl_tuple_alloc(length);

	va_list ap;
	va_start(ap, length);
	for (int i = 1; i <= length; ++i) {
		cpl_tuple_set(t, i, va_arg(ap, int));
	}
	va_end(ap);

	return t;
}


void cpl_tuple_print(cpl_tuple *t) {
	if (cpl_tuple_length(t) == 1) {
		printf("(%d,)\n", cpl_tuple_get(t, 1));
		return;
	}

	for (int i = 1; i <= cpl_tuple_length(t); ++i) {
		if (i == 1) {
			printf("(%d, ", cpl_tuple_get(t, i));
		} else if (i == cpl_tuple_length(t)) {
			printf("%d)\n", cpl_tuple_get(t, i));
		} else {
			printf("%d, ", cpl_tuple_get(t, i));
		}
	}
}

int cpl_tuple_pdt(cpl_tuple *t, int i, int j) {
	int pdt = 1;
	for (int k = i; k <= j; ++k) {
		pdt *= cpl_tuple_get(t, k);
	}

	return pdt;
}

int cpl_tuple_isequal(cpl_tuple *s, cpl_tuple *t) {
	if (cpl_tuple_length(s) != cpl_tuple_length(t))
		return 0;

	for (int i = 1; i <= cpl_tuple_length(s); ++i) {
		if (cpl_tuple_get(s, i) != cpl_tuple_get(t, i))
			return 0;
	}

	return 1;
}

cpl_tuple* cpl_tuple_copy(cpl_tuple *t) {
	cpl_tuple *s = cpl_tuple_alloc(cpl_tuple_length(t));
	for (int i = 1; i <= cpl_tuple_length(t); ++i)
		cpl_tuple_set(s, i, cpl_tuple_get(t, i));

	return s;
}

int cpl_tuple_ispos(cpl_tuple *t) {
	for (int i = 1; i <= cpl_tuple_length(t); ++i) {
		if (cpl_tuple_get(t, i) <= 0)
			return 0;
	}

	return 1;
}
