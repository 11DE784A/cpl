#include "cpl_defines.h"
#include "cpl_includes.h"

#include "cpl_commons.h"

#include "cpl_arrays.h"
#include "cpl_linalg.h"
#include "cpl_diffeq.h"

cpl_matrix *cpl_diffeq_rk4(cpl_vector *(*F)(cpl_vector *, double), 
						   cpl_matrix *Y, cpl_vector *t, double h) {

	int d = cpl_matrix_cols(Y),
		N = cpl_matrix_rows(Y);

	cpl_vector *dY = cpl_vector_alloc(d);
	cpl_vector *yi = cpl_vector_alloc(d);

	cpl_vector *k1 = cpl_vector_alloc(d);
	cpl_vector *k2 = cpl_vector_alloc(d);
	cpl_vector *k3 = cpl_vector_alloc(d);
	cpl_vector *k4 = cpl_vector_alloc(d);

	double ti;
	for (int i = 1; i <= N - 1; ++i) {
		ti = cpl_get(t, i);
		cpl_matrix_get_row(Y, i, yi);

		cpl_vector_overwrite(k1, yi);
		cpl_vector_overwrite(k2, yi);
		cpl_vector_overwrite(k3, yi);
		cpl_vector_overwrite(k4, yi);

		F(k1, ti);
		F(cpl_vector_axby(1, k2, h/2, k1), ti + h/2);
		F(cpl_vector_axby(1, k3, h/2, k2), ti + h/2);
		F(cpl_vector_axby(1, k4, h, k3), ti + h);

		cpl_vector_set_all(dY, 0);
		cpl_vector_add(k1, dY);
		cpl_vector_axby(1, dY, 2, k2);
		cpl_vector_axby(1, dY, 2, k3);
		cpl_vector_add(k4, dY);

		cpl_vector_axby(1, yi, h/6, dY);
		cpl_matrix_set_row(Y, i + 1, yi);
	}

	cpl_free(k4);
	cpl_free(k3);
	cpl_free(k2);
	cpl_free(k1);

	cpl_free(yi);
	cpl_free(dY);

	return Y;
}

cpl_matrix *cpl_diffeq_heat_forward(cpl_matrix *U, double κ, double dx, double dt) {
	int Nx = cpl_matrix_rows(U),
		Nt = cpl_matrix_cols(U);

	double α = κ*dt / pow(dx, 2), Uxt;

	for (int it = 2; it <= Nt; ++it) {
		for (int ix = 2; ix <= Nx - 1; ++ix) {
			Uxt = α * (cpl_get(U, ix+1, it-1) + cpl_get(U, ix-1, it-1))
					+ (1 - 2*α) * cpl_get(U, ix, it-1);
			cpl_set(U, ix, it, Uxt);
		}
	}

	return U;
}

cpl_matrix *cpl_diffeq_heat_backward(cpl_matrix *U, double κ, double dx, double dt) {
	int Nx = cpl_matrix_rows(U),
		Nt = cpl_matrix_cols(U);

	double α = κ*dt / pow(dx, 2), Aij;

	cpl_matrix *A = cpl_matrix_id(Nx);
	for (int i = 2; i <= Nx - 1; ++i) {
		for (int j = 2; j <= Nx - 1; ++j) {
			Aij = - α*(cpl_delta(i+1, j) + cpl_delta(i-1, j)) 
					+ (1 + 2*α) * cpl_delta(i, j);
			cpl_set(A, i, j, Aij);
		}
	}

	cpl_matrix *id = cpl_matrix_id(Nx);
	cpl_matrix *Ainv = cpl_matrix_alloc(Nx, Nx);
	cpl_linalg_seidel(A, id, Ainv, NULL);

	cpl_vector *v = cpl_vector_alloc(Nx);
	cpl_matrix_get_col(U, 1, v);
	for (int j = 2; j <= Nt; ++j) {
		cpl_mult(Ainv, v);
		cpl_matrix_set_col(U, j, v);
	}

	cpl_free(v);
	cpl_free(id);
	cpl_free(Ainv);
	cpl_free(A);

	return U;
}

cpl_matrix *cpl_diffeq_laplace_dirichlet(cpl_matrix *U, double dx, double dy) {
	int Nx = cpl_matrix_rows(U),
		Ny = cpl_matrix_cols(U);

	int iters = 0;
	double ϵ = 1;
	double Uxy, dx2 = dx*dx, dy2 = dy*dy;
	while (ϵ > TOL) {
		if (++iters > MAX_ITERS) {
			fprintf(stderr, "Failed to converge after %d iterations", MAX_ITERS);
			break;
		}

		ϵ = 0;
		for (int ix = 2; ix <= Nx - 1; ++ix) {
			for (int iy = 2; iy <= Ny - 1; ++iy) {
				Uxy = (dy2 * (cpl_get(U, ix - 1, iy) + cpl_get(U, ix + 1, iy))
						+ dx2 * (cpl_get(U, ix, iy - 1) + cpl_get(U, ix, iy + 1)))
							/ (2 * (dx2 + dy2));

				ϵ += sabs(Uxy - cpl_get(U, ix, iy));
				cpl_set(U, ix, iy, Uxy);
			}
		}
	}

	return U;
}
