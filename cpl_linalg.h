
cpl_matrix *cpl_linalg_gaussjordan(cpl_matrix*, cpl_matrix*);
cpl_matrix *cpl_linalg_ludecomp(cpl_matrix*);
cpl_matrix *cpl_linalg_cholesky(cpl_matrix*);

int cpl_linalg_jacobi(cpl_matrix*, cpl_matrix*, cpl_matrix*);
int cpl_linalg_seidel(cpl_matrix*, cpl_matrix*, cpl_matrix*);
void cpl_linalg_conjgrad(cpl_matrix*, cpl_vector*, cpl_vector*);
