
void cpl_lcg_seed(int);
void cpl_lcg_modulus(int);
void cpl_lcg_multiplier(int);

int cpl_rand();
double cpl_rand_uniform(double, double);

double cpl_mcint1d_naive(double (*)(double), double, double, int);
double cpl_mcint_box_naive(double (*)(cpl_vector *), int, double, double, int);

