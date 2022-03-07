
cpl_matrix *cpl_linalg_gaussjordan(cpl_matrix*, cpl_matrix*);
cpl_matrix *cpl_linalg_ludecomp(cpl_matrix*);
cpl_matrix *cpl_linalg_cholesky(cpl_matrix*);

void cpl_linalg_luseparate(cpl_matrix*, cpl_matrix*, cpl_matrix*);
cpl_vector *cpl_linalg_backsub(cpl_matrix*, cpl_vector*);
cpl_vector *cpl_linalg_forwardsub(cpl_matrix*, cpl_vector*);
cpl_vector *cpl_linalg_lusolve(cpl_matrix*, cpl_vector*);
scalar cpl_linalg_ludet(cpl_matrix*);

int cpl_linalg_jacobi(cpl_matrix*, cpl_matrix*, cpl_matrix*, cpl_vector*);
int cpl_linalg_seidel(cpl_matrix*, cpl_matrix*, cpl_matrix*);
void cpl_linalg_conjgrad_solve(cpl_matrix*, cpl_vector*, cpl_vector*);
void cpl_linalg_conjgrad(cpl_matrix*, cpl_matrix*, cpl_matrix*);


