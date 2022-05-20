#include "../code/cpl.h"

#include "glauber.h"

double hard_sphere_density(double r, void *params) {
	/* Normalized */
	density_params Nu = * (density_params *) params;
	double ρ0 = 3.0 / (4*PI*pow(Nu.R, 3));
	return (r <= Nu.R ? ρ0 : 0);
}

double woods_saxon_density(double r, void *params) {
	/* Normalized */
	density_params Nu = * (density_params *) params;
	return Nu.ρ0 * (1 + Nu.w*pow(r/Nu.R, 2)) / (1 + exp((r - Nu.R) / Nu.a));
}

double nuclear_thickness_radial(cpl_function density, double b, int N) {
	/* Not normalized */

	double density_longitudinal(double z, void *params) {
		return cpl_fn_eval(density, sqrt(b*b + z*z));
	}

	cpl_function F = { density_longitudinal, NULL };

	double zmax = 1;
	while (cpl_fn_eval(F, zmax) > TOL)
		++zmax;

	return 2*cpl_integrate_simpson(F, 0, zmax, N);
}

double nuclear_thickness(cpl_function density, cpl_vector *b, cpl_vector *s, int N) {
	double r = sqrt(cpl_vector_inner(b, b) 
					+ (s ? cpl_vector_inner(s, s) - 2*cpl_vector_inner(s, b) : 0));
	return nuclear_thickness_radial(density, r, N);
}

double nuclear_overlap(cpl_function ρA, cpl_function ρB, cpl_vector *b, int N) {
	double overlap_integrand(cpl_vector *s) {
		return nuclear_thickness(ρA, s, NULL, N) * nuclear_thickness(ρB, s, b, N);
	}

	/* Finding bounds of integration */
	double rmax = 1;
	cpl_vector *r = cpl_vector_alloc(2);
	cpl_vector_build(r, rmax, 0.0);
	while (rmax*overlap_integrand(r) > TOL)
		cpl_set(r, 1, ++rmax);
	cpl_free(r);

	return cpl_integrate_box_iter(overlap_integrand, 2, -rmax, rmax, N);
}

double nuclear_overlap_radial(cpl_function ρA, cpl_function ρB, double b, int N) {
	cpl_vector *r = cpl_vector_alloc(2);
	cpl_vector_build(r, b, 0.0);
	double T_AB = nuclear_overlap(ρA, ρB, r, N);
	cpl_free(r);

	return T_AB;
}

double inelastic_cross_section(cpl_function ρA, cpl_function ρB, double σ_NN, int N) {
	double A = ((density_params *) ρA.params)->A;
	double B = ((density_params *) ρB.params)->A;

	double cs_integrand(double b, void *params) {
		double n = σ_NN * nuclear_overlap_radial(ρA, ρB, b, N) / (A*B);
		// printf("b = %g, n = %g\n", b, n);
		return b * (1 - pow(1 - n, A*B));
	}

	cpl_function F = { cs_integrand, NULL };

	double bmax = 16; // From all simulations
	return 2*PI*cpl_integrate_simpson(F, 0, bmax, N);
}

double number_of_collisions(cpl_function ρA, cpl_function ρB, 
							double b, double σ_NN, 
							int N) {
	return σ_NN * nuclear_overlap_radial(ρA, ρB, b, N);
}

double number_of_participants(cpl_function ρA, cpl_function	ρB,
							  double b, double σ_NN,
							  int N) {

	density_params NuA = * (density_params *) ρA.params;
	density_params NuB = * (density_params *) ρB.params;

	cpl_vector *r = cpl_vector_alloc(2);
	cpl_vector_build(r, b, 0.0);

	double Npart_integrand(cpl_vector *s) {
		double TA, TB;
		TA = nuclear_thickness(ρA, s, NULL, N);
		TB = nuclear_thickness(ρB, s, r, N);

		return TA * (1 - pow(1 - σ_NN * TB / NuB.A, NuB.A)) 
					+ TB * (1 - pow(1 - σ_NN * TA / NuA.A, NuA.A));
	}

	double smax = 1;
	cpl_vector *smax_vec = cpl_vector_alloc(2);
	cpl_vector_build(smax_vec, smax, 0.0);
	while (smax*Npart_integrand(smax_vec) > TOL)
		cpl_set(smax_vec, 1, ++smax);
	cpl_free(smax_vec);

	double Npart = cpl_integrate_box_iter(Npart_integrand, 2, -smax, smax, N);

	cpl_free(r);

	return Npart;
}

double rand_impact_parameter(double bmax) {
	return sqrt(2*cpl_rand_uniform(0, pow(bmax, 2) / 2));
}

double rand_woods_saxon(cpl_function ρ, double M) {
	density_params Nu = * (density_params *) ρ.params;
	double r, u;

	while (1) {
		r = cpl_rand_uniform(0, 2*Nu.R);
		u = cpl_rand_uniform(0, M);
		if (u <= 4*PI*r*r*cpl_fn_eval(ρ, r))
			return r;
	}
}

cpl_matrix *generate_coords_3d(cpl_function ρ, cpl_matrix *coords, double b) {
	density_params Nu = * (density_params *) ρ.params;

	cpl_check(cpl_matrix_rows(coords) == Nu.A, "Number of rows in coords must equal A");

	double r, θ, ϕ;

	/* finding ρmax */
	double ρmax = 0;
	for (double x = 0; x <= Nu.R; x += 0.002) {
		if (ρmax < 4*PI*x*x*cpl_fn_eval(ρ, x))
			ρmax = 4*PI*x*x*cpl_fn_eval(ρ, x);
	}

	for (int i = 1; i <= Nu.A; ++i) {
		r = rand_woods_saxon(ρ, ρmax);
//		r = cbrt(3 * cpl_rand_uniform(0, pow(Nu.R, 3) / 3.0));
		θ = acos(cpl_rand_uniform(-1, 1));
		ϕ = cpl_rand_uniform(0, 2*PI);

		cpl_set(coords, i, 1, r*sin(θ)*cos(ϕ) + b/2);
		cpl_set(coords, i, 2, r*sin(θ)*sin(ϕ));
		cpl_set(coords, i, 3, r*cos(θ));
	}

	return coords;
}

cpl_matrix *generate_coords(cpl_function ρ, cpl_matrix *coords, double b) {
	density_params Nu = * (density_params *) ρ.params;

	cpl_check(cpl_matrix_rows(coords) == Nu.A, "Number of rows in coords must equal A");

	double r, θ, ϕ;

	/* finding ρmax */
	double ρmax = 0;
	for (double x = 0; x <= Nu.R; x += 0.002) {
		if (ρmax < 4*PI*x*x*cpl_fn_eval(ρ, x))
			ρmax = 4*PI*x*x*cpl_fn_eval(ρ, x);
	}

	for (int i = 1; i <= Nu.A; ++i) {
		r = rand_woods_saxon(ρ, ρmax);
//		r = cbrt(3 * cpl_rand_uniform(0, pow(Nu.R, 3) / 3.0));
		θ = acos(cpl_rand_uniform(-1, 1));
		ϕ = cpl_rand_uniform(0, 2*PI);

		cpl_set(coords, i, 1, r*sin(θ)*cos(ϕ) + b/2);
		cpl_set(coords, i, 2, r*sin(θ)*sin(ϕ));
	}

	return coords;
}

cpl_tuple gmc_simulate(cpl_function	ρA, cpl_function ρB, double b, double σ_NN) {
	cpl_tuple result = { .first = 0, .second = 0 };

	density_params NuA = * (density_params *) ρA.params;
	density_params NuB = * (density_params *) ρB.params;

	cpl_matrix *coords_A = cpl_matrix_alloc(NuA.A, 2);
	cpl_matrix *coords_B = cpl_matrix_alloc(NuB.A, 2);

	generate_coords(ρA, coords_A, b);
	generate_coords(ρB, coords_B, -b);

	cpl_vector *v = cpl_vector_alloc(2);
	cpl_vector *w = cpl_vector_alloc(2);

	double D = sqrt(σ_NN / PI);

	for (int i = 1; i <= NuA.A; ++i) {
		cpl_matrix_get_row(coords_A, i, v);
		for (int j = 1; j <= NuB.A; ++j) {
			cpl_matrix_get_row(coords_B, j, w);
			if (cpl_vector_l2dist(v, w) < D) {
				result.first += 1;
			}
		}
	}

	for (int i = 1; i <= NuA.A; ++i) {
		cpl_matrix_get_row(coords_A, i, v);
		for (int j = 1; j <= NuB.A; ++j) {
			cpl_matrix_get_row(coords_B, j, w);
			if (cpl_vector_l2dist(v, w) < D) {
				result.second += 1;
				break;
			}
		}
	}

	for (int i = 1; i <= NuB.A; ++i) {
		cpl_matrix_get_row(coords_B, i, v);
		for (int j = 1; j <= NuA.A; ++j) {
			cpl_matrix_get_row(coords_A, j, w);
			if (cpl_vector_l2dist(v, w) < D) {
				result.second += 1;
				break;
			}
		}
	}

	cpl_free(v);
	cpl_free(w);
	cpl_free(coords_A);
	cpl_free(coords_B);

	return result;
}

cpl_tuple gmc_compute(cpl_function ρA, cpl_function ρB, double b, double σ_NN, int N) {
	cpl_tuple res;
	double N_coll = 0, N_part = 0;
	for (int i = 0; i < N; ++i) {
		printf("\rit's not stuck, progress is being made: %2.2f", (i + 1.0) / N);
		res = gmc_simulate(ρA, ρB, b, σ_NN);
		N_coll += res.first;
		N_part += res.second;
	}

	res.first = N_coll / N;
	res.second = N_part / N;

	return res;
}

double gmc_collision_probability(cpl_function ρA, cpl_function ρB, 
								 double b, double σ_NN, 
								 int N) {
	double D = sqrt(σ_NN / PI);

	density_params NuA = * (density_params *) ρA.params;
	cpl_matrix *coords_A = cpl_matrix_alloc(NuA.A, 2);

	density_params NuB = * (density_params *) ρB.params;
	cpl_matrix *coords_B = cpl_matrix_alloc(NuB.A, 2);

	cpl_vector *v = cpl_vector_alloc(2);
	cpl_vector *w = cpl_vector_alloc(2);

	int collisions = 0, collided = 0;
	for (int e = 0; e < N; ++e) {
		generate_coords(ρA, coords_A, b);
		generate_coords(ρB, coords_B, -b);
		for (int i = 1; i <= NuA.A; ++i) {
			cpl_matrix_get_row(coords_A, i, v);
			for (int j = 1; j <= NuB.A; ++j) {
				cpl_matrix_get_row(coords_B, j, w);
				if (cpl_vector_l2dist(v, w) < D) {
					collided = 1;
					break;
				}
			}

			if (collided) {
				collided = 0;
				++collisions;
				break;
			}
		}
	}

	cpl_free(v);
	cpl_free(w);
	cpl_free(coords_A);
	cpl_free(coords_B);

	return (double) collisions / N;
}

double gmc_cross_section(cpl_function ρA, cpl_function ρB, double σ_NN, int N) {
	density_params NuA = * (density_params *) ρA.params;
	density_params NuB = * (density_params *) ρB.params;

	cpl_matrix *coords_A = cpl_matrix_alloc(NuA.A, 2);
	cpl_matrix *coords_B = cpl_matrix_alloc(NuB.A, 2);

	cpl_vector *v = cpl_vector_alloc(2);
	cpl_vector *w = cpl_vector_alloc(2);

	double D = sqrt(σ_NN / PI);

	double b, bmax = 20; // fm
	double σ_total = 0;
	for (int c = 0; c < N; ++c) {
		b = rand_impact_parameter(bmax);
		σ_total += gmc_collision_probability(ρA, ρB, b, σ_NN, 100);
	}

	return PI*bmax*bmax*σ_total / N;
}

