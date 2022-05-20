#include "../code/cpl.h"
#include "../code/cpl_bench.h"

#include "glauber.h"

#include "../gnuplot/gnuplot_i.h"

#define gnuplot_setup(_plt, _name) ({ \
	gnuplot_cmd(_plt, "set term png font \"monospace, 15\" linewidth 2"); \
	gnuplot_cmd(_plt, "set pointsize 2"); \
	gnuplot_cmd(_plt, "set grid"); \
	gnuplot_cmd(_plt, "set output \"" _name ".png\""); \
})

int main(void) {
	int N;

	// cpl_function ρ_Si = { woods_saxon_density, &Si };
	// cpl_function ρ_Ca = { woods_saxon_density, &Ca };
	cpl_function ρ_O  = { woods_saxon_density, &O };
	cpl_function ρ_Cu = { woods_saxon_density, &Cu };
	cpl_function ρ_I  = { woods_saxon_density, &J };
	cpl_function ρ_Au = { woods_saxon_density, &Au };
	cpl_function ρ_Pb = { woods_saxon_density, &Pb };

	/* Plotting densities */
	N = 20;

	cpl_function ρ_fn[5] = { ρ_O, ρ_Cu, ρ_I, ρ_Au, ρ_Pb };

	int numcases = 5;
	double ρ_norm[numcases];
	for (int i = 0; i < numcases; ++i)
		ρ_norm[i] = ((density_params *) ρ_fn[i].params)->ρ0;

	double ri = 0, rf = 10, rx = ri;
	double h = (rf - ri) / N;

//	cpl_matrix *rρ_mat = cpl_matrix_alloc(N, numcases + 1);
//	for (int i = 1; i <= N; ++i) {
//		cpl_set(rρ_mat, i, 1, rx);
//
//		for (int j = 0; j < numcases; ++j)
//			cpl_set(rρ_mat, i, j + 2, ρ_norm[j] * cpl_fn_eval(ρ_fn[j], rx));
//
//		rx += h;
//	}
//
//	cpl_matrix_savetxt("densities.txt", rρ_mat, "r O Cu I Au Pb");
//
	gnuplot_ctrl *plt;
//	plt = gnuplot_init();
//	gnuplot_setup(plt, "densities");
//	gnuplot_cmd(plt, "set xlabel \"radius (fm)\"");
//	gnuplot_cmd(plt, "set xtics 2");
//	gnuplot_cmd(plt, "set ylabel \"density (fm^{-3})\"");
//	gnuplot_cmd(plt, "plot for [col=2:6] \"densities.txt\" using 1:col with lines title columnheader");
//	gnuplot_close(plt);
//
//	cpl_free(rρ_mat);
//
//	/* Plotting overlap, N_coll */
//	N = 20;
//	ri = 0; rf = 15, rx = ri;
//	h = (rf - ri) / N;
//
//	double σ_NN = 4.2; // fermi squared
//
//	cpl_matrix *T_mat = cpl_matrix_alloc(N, numcases + 1);
//	cpl_matrix *Ncoll = cpl_matrix_alloc(N, numcases + 1);
//	cpl_matrix *Npart = cpl_matrix_alloc(N, numcases + 1);
//
//	for (int i = 1; i <= N; ++i) {
//		cpl_set(T_mat, i, 1, rx);
//		cpl_set(Ncoll, i, 1, rx);
//		cpl_set(Npart, i, 1, rx);
//
//		for (int j = 0; j < numcases; ++j) {
//			cpl_set(T_mat, i, j + 2, nuclear_overlap_radial(ρ_fn[j], ρ_fn[j], rx, N));
//			cpl_set(Ncoll, i, j + 2, number_of_collisions(ρ_fn[j], ρ_fn[j], rx, σ_NN, N));
//			cpl_set(Npart, i, j + 2, number_of_participants(ρ_fn[j], ρ_fn[j], rx, σ_NN, N));
//		}
//
//		rx += h;
//	}
//
//	/* Overlap */
//	cpl_matrix_savetxt("overlap.txt", T_mat, "r O Cu I Au Pb");
//
//	plt = gnuplot_init();
//	gnuplot_setup(plt, "overlap");
//	gnuplot_cmd(plt, "set xlabel \"impact parameter (fm)\"");
//	gnuplot_cmd(plt, "set xtics 2.5");
//	gnuplot_cmd(plt, "set ylabel \"nuclear overlap, T_{AB} (fm^{-2})\"");
//	gnuplot_cmd(plt, "set yrange [0:325]");
//	gnuplot_cmd(plt, "plot for [col=2:6] \"overlap.txt\" using 1:col with lines title columnheader");
//
//	/* Collisions */
//	cpl_matrix_savetxt("collisions.txt", Ncoll, "r O Cu I Au Pb");
//
//	gnuplot_setup(plt, "collisions");
//	gnuplot_cmd(plt, "set ylabel \"N_{coll}\"");
//	gnuplot_cmd(plt, "set yrange [0:1400]");
//	gnuplot_cmd(plt, "plot for [col=2:6] \"collisions.txt\" using 1:col with lines title columnheader");
//
//	cpl_matrix_savetxt("participants.txt", Npart, "r O Cu I Au Pb");
//
//	gnuplot_setup(plt, "participants");
//	gnuplot_cmd(plt, "set ylabel \"N_{part}\"");
//	gnuplot_cmd(plt, "set yrange [0:450]");
//	gnuplot_cmd(plt, "plot for [col=2:6] \"participants.txt\" using 1:col with lines title columnheader");
//	gnuplot_close(plt);

	/* Cross section */
	N = 25;
	ri = 2e-3; rf = 20;
	h = pow(rf/ri, 1.0 / N);
	cpl_matrix *σmat = cpl_matrix_alloc(20, 3);
	for (int i = 1; i <= 20; ++i) {
		printf("point %d of 20\n", i);
		cpl_set(σmat, i, 1, ri * pow(h, i));
		cpl_set(σmat, i, 2, inelastic_cross_section(ρ_Au, ρ_Au, cpl_get(σmat, i, 1), N));
		cpl_set(σmat, i, 3, gmc_cross_section(ρ_Au, ρ_Au, cpl_get(σmat, i, 1), 1e3));
	}

	for (int i = 1; i <= 20; ++i) {
		cpl_set(σmat, i, 1, 10.0*cpl_get(σmat, i, 1));  // mb
		cpl_set(σmat, i, 2, 0.01*cpl_get(σmat, i, 2));  // b
		cpl_set(σmat, i, 3, 0.01*cpl_get(σmat, i, 3));  // b
	}

	cpl_matrix_savetxt("cross_section.txt", σmat, "σ_{NN} Optical Monte-Carlo");

	plt = gnuplot_init();
	gnuplot_setup(plt, "cross_section");
	gnuplot_cmd(plt, "set xlabel \"σ_{NN} (mb)\"");
	gnuplot_cmd(plt, "set logscale x 10");
	gnuplot_cmd(plt, "set xrange [0.01:100]");
//	gnuplot_cmd(plt, "set xtics 2.5");
	gnuplot_cmd(plt, "set ylabel \"σ_{AB} (barns)\"");
	gnuplot_cmd(plt, "set yrange [0:8]");
	gnuplot_cmd(plt, "plot for [col=2:3] \"cross_section.txt\" using 1:col with points title columnheader");
	gnuplot_close(plt);

	return 0;
}
