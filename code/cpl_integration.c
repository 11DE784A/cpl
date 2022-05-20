#include "cpl_includes.h"
#include "cpl_defines.h"
#include "cpl_commons.h"
#include "cpl_arrays.h"
#include "cpl_integration.h"

#include "cpl_quadgl.h"

double cpl_integrate_quadgl(cpl_function F, double a, double b, int N) {
	quadg_vals vals = gl[N - 4];
	double Σ = 0;
	for (int i = 0; i < N; ++i)
		Σ += vals.wt[i] * cpl_fn_eval(F, vals.pt[i] * (b - a) / 2 + (a + b) / 2) * (b - a) / 2;

	return Σ;
}

double cpl_integrate_trapezoid(cpl_function F, double a, double b, int N) {
	_INT_SCAF(Fcurr + Fnext, h / 2);
}

double cpl_integrate_simpson(cpl_function F, double a, double b, int N) {
	_INT_SCAF(Fcurr + 4 * cpl_fn_eval(F, x + h/2) + Fnext, h / 6);
}


double cpl_integrate_box_iter(double (*G)(cpl_vector *), int dim, double a, double b, int N) {
	if (dim == 1) {

		double g(double x, void *params) {
			double xarr[1] = { x };
			cpl_block xblock = { .size = 1, .array = xarr };
			cpl_vector xvec = { .block = &xblock, .dim = 1 };
			return G(&xvec);
		}

		cpl_function F = { .function = g, .params = NULL };
		return cpl_integrate_trapezoid(F, a, b, N);

	} else if (dim > 1) {

		double g(double x, void *params) {

			double H(cpl_vector *y) {
				cpl_vector *z = cpl_vector_alloc(dim);

				cpl_set(z, 1, x);
				for (int i = 1; i <= dim - 1; ++i)
					cpl_set(z, i + 1, cpl_get(y, i));

				double Gz = G(z);

				cpl_free(z);

				return Gz;
			}

			return cpl_integrate_box_iter(H, dim - 1, a, b, N);
		}

		cpl_function F = { .function = g, .params = NULL };
		return cpl_integrate_trapezoid(F, a, b, N);
	}

	return 0;
}

