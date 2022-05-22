
cpl_matrix *cpl_diffeq_rk4(cpl_vector *(*F)(cpl_vector *, double), 
						   cpl_matrix *Y, cpl_vector *t, double h);

cpl_matrix *cpl_diffeq_heat_forward(cpl_matrix *, double, double, double);
cpl_matrix *cpl_diffeq_heat_backward(cpl_matrix *, double, double, double);
cpl_matrix *cpl_diffeq_laplace_dirichlet(cpl_matrix *, double, double);

