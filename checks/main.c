#include <stdlib.h>
#include <check.h>
#include "../code/cpl.h"

START_TEST(test_vectors) {
	int dim = 6;

	/* Zero vector */
	cpl_vector *b = cpl_vector_calloc(dim);
	for (int i = 1; i <= dim; ++i)
		ck_assert(cpl_get(b, i) == 0.0);

	cpl_free(b);

	/* Creating vector */
	cpl_vector *u = cpl_vector_alloc(dim);
	cpl_vector_build(u, -5.0/2, 2.0/3, 3.0, -4.0/3, -1.0/3, 5.0/3);

	ck_assert(cpl_vector_dim(u) == dim);
	ck_assert(cpl_get(u, 1) == -5.0/2);
	ck_assert(cpl_get(u, 3) == 3.0);
	ck_assert(cpl_get(u, 6) == 5.0/3);

	/* Copying vector */
	cpl_vector *v = cpl_vector_copy(u);

	ck_assert(u != v);
	ck_assert(u->block != v->block);
	ck_assert(u->block->array != v->block->array);
	ck_assert(cpl_vector_dim(v) == dim);

	for (int i = 1; i <= dim; ++i)
		ck_assert(cpl_get(u, i) == cpl_get(v, i));

	/* Like vector */
	cpl_vector *w = cpl_vector_like(u);
	ck_assert(u != w);
	ck_assert(u->block != w->block);
	ck_assert(u->block->array != w->block->array);
	ck_assert(cpl_vector_dim(w) == dim);

	/* Set all */
	scalar x = 3.0;
	cpl_vector_set_all(w, x);
	ck_assert(cpl_vector_dim(w) == dim);
	for (int i = 1; i <= dim; ++i)
		ck_assert(cpl_get(w, i) == x);

	/* Hot vectors */
	cpl_vector *hot = cpl_vector_hot(dim);
	ck_assert(cpl_vector_dim(hot) == dim);
	for (int i = 1; i <= dim; ++i)
		ck_assert(cpl_get(hot, i) == 1.0);

	int j = 3;
	cpl_vector *onehot = cpl_vector_1hot(dim, j);
	ck_assert(cpl_vector_dim(onehot) == dim);
	for (int i = 1; i <= dim; ++i)
		ck_assert(cpl_get(onehot, i) == (i == j ? 1 : 0));

	/* Overwriting vector */
	cpl_overwrite(hot, onehot);
	for (int i = 1; i <= dim; ++i)
		ck_assert(cpl_get(onehot, i) == cpl_get(hot, i));

	/* Vector frees */
	cpl_free(onehot);
	cpl_free(hot);
	cpl_free(w);
	cpl_free(v);
	cpl_free(u);


} END_TEST

START_TEST(test_norms) {
	int dim = 4;

	cpl_vector *v = cpl_vector_hot(dim);
	cpl_vector *e1 = cpl_vector_1hot(dim, 1);
	cpl_vector *e2 = cpl_vector_1hot(dim, 2);

	ck_assert(cpl_vector_inner(v,  e1) == 1.0);
	ck_assert(cpl_vector_inner(e1, e1) == 1.0);
	ck_assert(cpl_vector_inner(e1, e2) == 0.0);
	ck_assert(cpl_vector_inner(e2, e2) == 1.0);

	ck_assert(cpl_vector_l2norm(v) == sqrt(dim));
	ck_assert(cpl_vector_l2norm(v) == sqrt(cpl_vector_inner(v, v)));

	ck_assert(cpl_vector_l2dist(v,  e1) == sqrt(dim - 1));
	ck_assert(cpl_vector_l2dist(e1, e2) == sqrt(2));

	cpl_vector_normalize(v);
	for (int i = 1; i <= dim; ++i)
		ck_assert(cpl_get(v, i) == 1.0 / sqrt(dim));

	cpl_free(v);
	cpl_free(e2);
	cpl_free(e1);
} END_TEST

START_TEST(test_stack) {
	int dim = 5, ops = 5;
	scalar x, y;
	cpl_vector *u, *v = cpl_vector_calloc(dim);

	/* Pushes */
	for (int i = 1; i <= ops; ++i) {
		u = cpl_vector_push(v, i);
		ck_assert(u == v);
		ck_assert(v->dim <= v->block->size);
		ck_assert(cpl_vector_dim(v) == dim + i);
		ck_assert(cpl_get(v, dim + i) == i);
	}

	/* Pops */
	dim = cpl_vector_dim(v);
	for (int i = 1; i <= ops; ++i) {
		x = cpl_get(v, cpl_vector_dim(v));
		y = cpl_vector_pop(v);
		ck_assert(x == y);
		ck_assert(v->dim <= v->block->size);
		ck_assert(cpl_vector_dim(v) == dim - i);
	}

	cpl_free(v);
} END_TEST

START_TEST(test_matrices) {
	int dim = 6;

	/* Zero matrix */
	int rows = dim, cols = dim - 2;
	cpl_matrix *C = cpl_matrix_calloc(rows, cols);

	ck_assert(cpl_matrix_rows(C) == rows);
	ck_assert(cpl_matrix_cols(C) == cols);

	for (int i = 1; i <= rows; ++i) {
		for (int j = 1; j <= cols; ++j) {
			ck_assert(cpl_get(C, i, j) == 0.0);
		}
	}

	/* Creating matrix */
	cpl_matrix *A = cpl_matrix_alloc(dim, dim);
	cpl_matrix_build(A,  1.0, -1.0,  4.0,  0.0,  2.0,  9.0,
						 0.0,  5.0, -2.0,  7.0,  8.0,  4.0,
						 1.0,  0.0,  5.0,  7.0,  3.0, -2.0,
						 6.0, -1.0,  2.0,  3.0,  0.0,  8.0,
						-4.0,  2.0,  0.0,  5.0, -5.0,  3.0,
						 0.0,  7.0, -1.0,  5.0,  4.0, -2.0);

	ck_assert(cpl_matrix_rows(A) == dim);
	ck_assert(cpl_matrix_cols(A) == dim);

	ck_assert(cpl_get(A, 1, 1) == 1.0);
	ck_assert(cpl_get(A, 3, 3) == 5.0);
	ck_assert(cpl_get(A, 3, 6) == -2.0);
	ck_assert(cpl_get(A, 6, 6) == -2.0);

	/* Copying matrix */
	cpl_matrix *B = cpl_matrix_copy(A);

	ck_assert(A != B);
	ck_assert(A->block != B->block);
	ck_assert(A->block->array != B->block->array);
	ck_assert(cpl_matrix_rows(A) == cpl_matrix_rows(B));
	ck_assert(cpl_matrix_cols(A) == cpl_matrix_cols(B));

	for (int i = 1; i <= dim; ++i) {
		for (int j = 1; j <= dim; ++j) {
			ck_assert(cpl_get(A, i, j) == cpl_get(B, i, j));
		}
	}

	/* Identity matrix */
	cpl_matrix *id = cpl_matrix_id(dim);

	ck_assert(cpl_matrix_rows(id) == dim);
	ck_assert(cpl_matrix_cols(id) == dim);

	for (int i = 1; i <= dim; ++i) {
		for (int j = 1; j <= dim; ++j) {
			ck_assert(cpl_get(id, i, j) == (i == j ? 1 : 0));
		}
	}

	/* Squareness */
	ck_assert(cpl_matrix_issquare(A));
	ck_assert(!cpl_matrix_issquare(C));

	/* Set all */
	scalar x = 4;
	cpl_matrix_set_all(A, x);
	for (int i = 1; i <= dim; ++i) {
		for (int j = 1; j <= dim; ++j) {
			ck_assert(cpl_get(A, i, j) == x);
		}
	}

	/* Overwrite */
	cpl_overwrite(A, B);
	for (int i = 1; i <= dim; ++i) {
		for (int j = 1; j <= dim; ++j) {
			ck_assert(cpl_get(A, i, j) == cpl_get(B, i, j));
		}
	}


	cpl_free(id);
	cpl_free(C);
	cpl_free(B);
	cpl_free(A);

} END_TEST

START_TEST(test_loadtxt) {
	char fpath[30] = "data/linear_system.txt";

	/* Vectors */
	cpl_vector *u = cpl_vector_loadtxt(fpath, 2, 7, 3);
	cpl_vector *v = cpl_vector_alloc(6);
	cpl_vector_build(v, 0.0, 0.5, 1.5, 0.0, 0.0, 0.0);
	ck_assert(cpl_vector_isequal(u, v));

	cpl_free(v);
	cpl_free(u);

	u = cpl_vector_loadtxt(fpath, 11, 16, 1);
	v = cpl_vector_alloc(6);
	cpl_vector_build(v, -1.0, 0.0, 2.75, 2.5, -3.0, 2.0); 
	ck_assert(cpl_vector_isequal(u, v));

	cpl_free(v);
	cpl_free(u);

	/* Matrices */
	cpl_matrix *A = cpl_matrix_loadtxt(fpath, 2, 7, 6);
	cpl_matrix *B = cpl_matrix_alloc(6, 6);
	cpl_matrix_build(B, -2.0,  0.0,  0.0, -1.0,  0.0,  0.5,
						 0.0,  4.0,  0.5,  0.0,  1.0,  0.0,
						 0.0,  0.5,  1.5,  0.0,  0.0,  0.0,
						-1.0,  0.0,  0.0, -2.0,  0.0,  1.0,
						 0.0,  1.0,  0.0,  0.0, -2.5,  0.0,
						 0.5,  0.0,  0.0,  1.0,  0.0, -3.75);
	ck_assert(cpl_matrix_isequal(A, B));

	cpl_matrix *b = cpl_matrix_loadtxt(fpath, 11, 16, 1);
	cpl_matrix *c = cpl_matrix_alloc(6, 1);
	cpl_matrix_build(c, -1.0, 0.0, 2.75, 2.5, -3.0, 2.0); 
	ck_assert(cpl_matrix_isequal(b, c));

	cpl_free(c);
	cpl_free(b);
	cpl_free(B);
	cpl_free(A);
} END_TEST

START_TEST(test_casting) {
	/* Flattening matrices */
	cpl_matrix *X = cpl_matrix_alloc(6, 1);
	cpl_matrix_build(X, -1.0, 0.0, 2.75, 2.5, -3.0, 2.0); 

	cpl_vector *x = cpl_matrix_flatten(X);
	ck_assert(cpl_vector_dim(x) == cpl_matrix_rows(X));
	for (int i = 1; i <= cpl_vector_dim(x); ++i)
		ck_assert(cpl_get(x, i) == cpl_get(X, i, 1));

	cpl_free(x);
	cpl_free(X);

} END_TEST

START_TEST(test_algebra) {
	int dim = 6;

	cpl_matrix *A = cpl_matrix_alloc(dim, dim);
	cpl_matrix_build(A,  1.0, -1.0,  4.0,  0.0,  2.0,  9.0,
						 0.0,  5.0, -2.0,  7.0,  8.0,  4.0,
						 1.0,  0.0,  5.0,  7.0,  3.0, -2.0,
						 6.0, -1.0,  2.0,  3.0,  0.0,  8.0,
						-4.0,  2.0,  0.0,  5.0, -5.0,  3.0,
						 0.0,  7.0, -1.0,  5.0,  4.0, -2.0);

	cpl_vector *b = cpl_vector_alloc(dim);
	cpl_vector_build(b, -5.0/2, 2.0/3, 3.0, -4.0/3, -1.0/3, 5.0/3);

	/* Vector scaling */
	cpl_vector *b_copy = cpl_vector_copy(b);
	scalar c = 3.0;
	cpl_vector_scale(b_copy, c);
	for (int i = 1; i <= dim; ++i)
		ck_assert(cpl_get(b_copy, i) == c * cpl_get(b, i));

	/* Matrix scaling */
	cpl_matrix *A_copy = cpl_matrix_copy(A);
	cpl_matrix_scale(A_copy, c);
	for (int i = 1; i <= dim; ++i) {
		for (int j = 1; j <= dim; ++j) {
			ck_assert(cpl_get(A_copy, i, j) == c * cpl_get(A, i, j));
		}
	}

	/* Row operations */
	cpl_vector *v = cpl_vector_alloc(dim);
	int row = 3;

	/* Get */
	cpl_matrix_get_row(A, row, v);
	for (int j = 1; j <= dim; ++j)
		ck_assert(cpl_get(v, j) == cpl_get(A, row, j));

	/* Scale */
	c = row;
	cpl_overwrite(A_copy, A);
	cpl_matrix_scale_row(A, row, c);
	for (int i = 1; i <= dim; ++i) {
		for (int j = 1; j <= dim; ++j) {
			ck_assert(cpl_get(A, i, j) == ((i == row ? c : 1) * cpl_get(A_copy, i, j)));
		}
	}

	/* Swap */
	cpl_overwrite(A, A_copy);
	int i1 = 2, i2 = 4;
	cpl_matrix_swap_rows(A, i1, i2);
	for (int i = 1; i <= dim; ++i) {
		if (i == i1 || i == i2) continue;
		for (int j = 1; j <= dim; ++j) {
			ck_assert(cpl_get(A, i, j) == cpl_get(A_copy, i, j));
		}
	}

	for (int j = 1; j <= dim; ++j) {
		ck_assert(cpl_get(A, i1, j) == cpl_get(A_copy, i2, j));
		ck_assert(cpl_get(A, i2, j) == cpl_get(A_copy, i1, j));
	}

	/* Add to row */
	cpl_overwrite(A, A_copy);
	cpl_matrix_add_to_row(A, row, b);
	for (int j = 1; j <= dim; ++j)
		ck_assert(cpl_get(A, row, j) == cpl_get(b, j) + cpl_get(A_copy, row, j));

	for (int i = 1; i <= dim; ++i) {
		if (i == row) continue;
		for (int j = 1; j <= dim; ++j) {
			ck_assert(cpl_get(A, i, j) == cpl_get(A_copy, i, j));
		}
	}

	/* Get cols */
	for (int col = 1; col <= dim; ++col) {
		cpl_matrix_get_col(A, col, v);
		for (int i = 1; i <= dim; ++i)
			ck_assert(cpl_get(v, i) == cpl_get(A, i, col));
	}

	cpl_free(v);
	cpl_free(A_copy);
	cpl_free(b_copy);
	cpl_free(b);
	cpl_free(A);

} END_TEST

START_TEST(test_addition) {
	char fpath[30] = "data/addition.txt";

	/* Vectors */
	cpl_vector *u, *v, *w;
	u = cpl_vector_loadtxt(fpath, 4, 9, 1);
	v = cpl_vector_loadtxt(fpath, 12, 17, 1);
	cpl_add(u, v);

	w = cpl_vector_loadtxt(fpath, 20, 25, 1);
	ck_assert(cpl_vector_isequal(w, v));

	cpl_free(w);
	cpl_free(v);
	cpl_free(u);

	/* Square Matrices */
	cpl_matrix *A, *B, *C;
	A = cpl_matrix_loadtxt(fpath, 30, 35, 6);
	B = cpl_matrix_loadtxt(fpath, 38, 43, 6);
	cpl_add(A, B);

	C = cpl_matrix_loadtxt(fpath, 46, 51, 6);
	ck_assert(cpl_matrix_isequal(C, B));

	cpl_free(C);
	cpl_free(B);
	cpl_free(A);

	/* Tall Matrices */
	A = cpl_matrix_loadtxt(fpath, 56, 63, 6);
	B = cpl_matrix_loadtxt(fpath, 66, 73, 6);
	cpl_add(A, B);

	C = cpl_matrix_loadtxt(fpath, 76, 83, 6);
	ck_assert(cpl_matrix_isequal(C, B));

	cpl_free(C);
	cpl_free(B);
	cpl_free(A);

} END_TEST

START_TEST(test_mult) {
	char fpath[30] = "data/multiplication.txt";

	/* Matrix-Vector */
	cpl_matrix *A = cpl_matrix_loadtxt(fpath, 4, 9, 8);
	cpl_vector *v = cpl_vector_loadtxt(fpath, 12, 19, 1);
	cpl_vector *Av = cpl_vector_loadtxt(fpath, 22, 27, 1);

	cpl_vector *Av_calc = cpl_mult_alloc(A, v);
	ck_assert(cpl_vector_isequal(Av_calc, Av));

	cpl_mult_overwrite(A, v, Av_calc);
	ck_assert(cpl_vector_isequal(Av_calc, Av));

	cpl_mult(A, v);
	ck_assert(cpl_vector_isequal(v, Av));

	cpl_free(Av_calc);
	cpl_free(Av);
	cpl_free(v);
	cpl_free(A);

	/* Matrix-Matrix */
	A = cpl_matrix_loadtxt(fpath, 32, 37, 8);
	cpl_matrix *B = cpl_matrix_loadtxt(fpath, 40, 47, 7);
	cpl_matrix *AB = cpl_matrix_loadtxt(fpath, 50, 55, 7);

	cpl_matrix *AB_calc = cpl_mult_alloc(A, B);
	ck_assert(cpl_matrix_isequal(AB_calc, AB));

	cpl_mult(A, B);
	ck_assert(cpl_matrix_isequal(B, AB));

	cpl_free(AB_calc);
	cpl_free(AB);
	cpl_free(B);
	cpl_free(A);

	/* Vector-Matrix */
	cpl_matrix *u = cpl_matrix_loadtxt(fpath, 60, 67, 1);
	cpl_matrix *wT = cpl_matrix_loadtxt(fpath, 70, 70, 6);
	cpl_matrix *uwT = cpl_matrix_loadtxt(fpath, 73, 80, 6);

	cpl_matrix *uwT_calc = cpl_mult_alloc(u, wT);
	ck_assert(cpl_matrix_isequal(uwT_calc, uwT));

	cpl_mult(u, wT);
	ck_assert(cpl_matrix_isequal(wT, uwT));

	cpl_free(uwT_calc);
	cpl_free(uwT);
	cpl_free(wT);
	cpl_free(u);

} END_TEST

/*
 * mult_alloc, matrix*vector
 * mult_overwrite, matrix*vector
 */

START_TEST(test_fly) {
	int dim = 6;
	char fpath[30] = "data/fly.txt";
	cpl_matrix *A = cpl_matrix_loadtxt(fpath, 3, 8, 6);
	cpl_vector *b = cpl_vector_loadtxt(fpath, 12, 17, 1);
	cpl_vector *Ab = cpl_vector_loadtxt(fpath, 21, 26, 1);

	scalar m = 0.2;
	scalar A_fn (int i, int j) { 
		return (cpl_delta(i + 1, j) + cpl_delta(i - 1, j)) / 2.0 
				+ (m*m - 1) * cpl_delta(i, j); }

	for (int i = 1; i <= dim; ++i) {
		for (int j = 1; j <= dim; ++j)
			ck_assert(sabs(A_fn(i, j) - cpl_get(A, i, j)) < TOL);
	}

	cpl_vector *Ab_calc = cpl_mvfly_mult_alloc(A_fn, b);
	ck_assert(cpl_vector_isequal(Ab_calc, Ab));

	cpl_mvfly_mult_overwrite(A_fn, b, Ab_calc);
	ck_assert(cpl_vector_isequal(Ab_calc, Ab));

	cpl_free(Ab_calc);
	cpl_free(Ab);
	cpl_free(b);
	cpl_free(A);

} END_TEST

Suite *array_suite(void) {
	Suite *s;
	TCase *tc_vectors, *tc_matrices, *tc_loadtxt, *tc_algebra;

	s = suite_create("Arrays");

	tc_vectors = tcase_create("Vectors");
	tcase_add_test(tc_vectors, test_vectors);
	tcase_add_test(tc_vectors, test_norms);
	tcase_add_test(tc_vectors, test_stack);

	tc_matrices = tcase_create("Matrices");
	tcase_add_test(tc_matrices, test_matrices);

	tc_loadtxt = tcase_create("Loading Data");
	tcase_add_test(tc_loadtxt, test_loadtxt);

	tc_algebra = tcase_create("Algebra");
	tcase_add_test(tc_algebra, test_casting);
	tcase_add_test(tc_algebra, test_algebra);
	tcase_add_test(tc_algebra, test_addition);
	tcase_add_test(tc_algebra, test_mult);
	tcase_add_test(tc_algebra, test_fly);

	suite_add_tcase(s, tc_vectors);
	suite_add_tcase(s, tc_matrices);
	suite_add_tcase(s, tc_loadtxt);
	suite_add_tcase(s, tc_algebra);

	return s;
}

START_TEST(test_fit_linear) {
	char fpath[30] = "data/fitting.txt";
	cpl_vector *xdata = cpl_vector_loadtxt(fpath, 34, 54, 1);
	cpl_vector *ydata = cpl_vector_loadtxt(fpath, 34, 54, 2);
	cpl_vector *params = cpl_vector_loadtxt(fpath, 58, 59, 1);

	cpl_vector *params_calc = cpl_vector_alloc(2);
	cpl_vector_set_all(params_calc, 1.0);
	cpl_stats_linfit(xdata, ydata, NULL, params_calc, NULL);

	ck_assert(cpl_vector_l2dist(params, params_calc) < TOL);

	cpl_free(params_calc);
	cpl_free(params);
	cpl_free(ydata);
	cpl_free(xdata);
}

START_TEST(test_fit_cubic) {
	char fpath[30] = "data/fitting.txt";
	cpl_vector *xdata = cpl_vector_loadtxt(fpath, 3, 23, 1);
	cpl_vector *ydata = cpl_vector_loadtxt(fpath, 3, 23, 2);
	cpl_vector *params = cpl_vector_loadtxt(fpath, 27, 30, 1);

	cpl_vector *params_calc = cpl_vector_alloc(4);
	cpl_vector_set_all(params_calc, 1.0);

	cpl_stats_linreg(4, polynomial, xdata, ydata, NULL, params_calc, NULL);

	ck_assert(cpl_vector_l2dist(params, params_calc) < TOL);

	cpl_free(params_calc);
	cpl_free(params);
	cpl_free(xdata);
	cpl_free(ydata);

} END_TEST

Suite *fit_suite(void) {
	Suite *s;
	TCase *tc_lin, *tc_cubic;

	s = suite_create("Least squares fitting");

	tc_cubic = tcase_create("Cubic fitting");
	tcase_add_test(tc_cubic, test_fit_cubic);
	suite_add_tcase(s, tc_cubic);

	tc_lin = tcase_create("Linear fitting");
	tcase_add_test(tc_lin, test_fit_linear);
	suite_add_tcase(s, tc_lin);

	return s;
}

int main(void) {
	int number_failed = 0;

	#define NUM_SUITES 2
	Suite *suites[NUM_SUITES] = {array_suite(), fit_suite()};
	SRunner *runner;

	for (int i = 0; i < NUM_SUITES; ++i) {
		runner = srunner_create(suites[i]);
		srunner_run_all(runner, CK_NORMAL);
		number_failed += srunner_ntests_failed(runner);
		srunner_free(runner);
	}

	return (number_failed == 0) ? EXIT_SUCCESS : EXIT_FAILURE;
}
