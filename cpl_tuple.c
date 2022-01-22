#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>

#include "cpl.h"

int cpl_tuple_get(cpl_tuple *t, int j) {
	return t->array[j-1];
}

int cpl_tuple_set(cpl_tuple *t, int j, int x) {
	t->array[j-1] = x;
	return x;
}

int cpl_tuple_length(cpl_tuple *t) {
	return t->length;
}

void cpl_tuple_free(cpl_tuple *t) {
	free(t->array);
	free(t);
}

cpl_tuple *cpl_tuple_alloc(int length) {
	cpl_tuple *t = malloc(sizeof(cpl_tuple));
	t->length = length;

	int *array = malloc(length * sizeof(int));
	t->array = array;

	return t;
}

cpl_tuple *cpl_tuple_salloc(int length, ...) {
	cpl_tuple *t = cpl_tuple_alloc(length);

	va_list ap;
	va_start(ap, length);
	for (int i = 1; i <= length; ++i) {
		cpl_tuple_set(t, i, va_arg(ap, int));
	}
	va_end(ap);

	return t;
}

cpl_tuple *cpl_tuple_ralloc(int length, int c) {
	cpl_tuple *t = cpl_tuple_alloc(length);

	for (int i = 1; i <= cpl_tuple_length(t); ++i) {
		cpl_tuple_set(t, i, c);
	}

	return t;
};

cpl_tuple *cpl_tuple_zalloc(int length) {
	return cpl_tuple_ralloc(length, 0);
}

cpl_tuple *cpl_tuple_1alloc(int length) {
	return cpl_tuple_ralloc(length, 1);
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

int cpl_tuple_mult(cpl_tuple *t, int i) {
	int pdt = 1;
	for (int j = i; j <= cpl_tuple_length(t); ++j) {
		pdt *= cpl_tuple_get(t, j);
	}

	return pdt;
}

int cpl_tuple_equal(cpl_tuple *s, cpl_tuple *t) {
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
