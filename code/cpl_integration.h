
#define _INT_SCAF(rule, factor) ({ \
	double Fnext, Fcurr = cpl_fn_eval(F, a); \
	double x, h, Σ = 0; \
	h = (b - a) / N; \
	for (int i = 0; i < N; ++i) { \
		x = a + i*h; \
		Fnext = cpl_fn_eval(F, x + h); \
		Σ += rule; \
		Fcurr = Fnext; \
	} \
	return factor * Σ; \
})

typedef struct {
	double *pt;
	double *wt;
	int degree;
} quadg_vals;

double cpl_integrate_trapezoid(cpl_function, double, double, int);
double cpl_integrate_simpson(cpl_function, double, double, int);
double cpl_integrate_quadgl(cpl_function F, double, double, int N);
double cpl_integrate_box_iter(double (*)(cpl_vector *), int, double, double, int);

