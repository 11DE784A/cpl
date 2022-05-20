
typedef struct {
	/*
	 * Woods-Saxon density parameters
	 * A : Number of nucleons
	 * R : Nuclear radius
	 * a : Skin depth
	 * w : Deformation
	 */

	double A, R, a, w, ρ0;
} density_params;

static density_params  O = { .A = 16 , .R = 2.608, .a = 0.513, .w = -0.051, .ρ0 = 0.1654};
static density_params Si = { .A = 27 , .R = 3.07,  .a = 0.519, .w = 0, .ρ0 = 0.1739};
static density_params Ca = { .A = 40 , .R = 3.766, .a = 0.586, .w = -0.061, .ρ0 = 0.1699};
static density_params Cu = { .A = 63 , .R = 4.214, .a = 0.586, .w = 0, .ρ0 = 0.1701};
static density_params  J = { .A = 110, .R = 5.33,  .a = 0.535, .w = 0, .ρ0 = 0.1577};
static density_params Au = { .A = 197, .R = 6.38 , .a = 0.535, .w = 0, .ρ0 = 0.1693};
static density_params Pb = { .A = 208, .R = 6.624, .a = 0.549, .w = 0, .ρ0 = 0.1600};

double hard_sphere_density(double, void *);
double woods_saxon_density(double, void *);

double nuclear_thickness_radial(cpl_function, double, int);
double nuclear_thickness(cpl_function, cpl_vector *, cpl_vector *, int);
double nuclear_overlap(cpl_function, cpl_function, cpl_vector *b, int);
double nuclear_overlap_radial(cpl_function, cpl_function, double, int);

double inelastic_cross_section(cpl_function, cpl_function, double, int);
double number_of_collisions(cpl_function, cpl_function, double, double, int);
double number_of_participants(cpl_function, cpl_function, double, double, int);

double rand_impact_parameter(double);
double rand_woods_saxon(cpl_function, double);

cpl_matrix *generate_coords_3d(cpl_function, cpl_matrix *, double);
cpl_matrix *generate_coords(cpl_function, cpl_matrix *, double);

cpl_tuple gmc_simulate(cpl_function, cpl_function, double, double);
cpl_tuple gmc_compute(cpl_function, cpl_function, double, double, int);

double gmc_collision_probability(cpl_function, cpl_function, double, double, int);
double gmc_cross_section(cpl_function, cpl_function, double, int);
