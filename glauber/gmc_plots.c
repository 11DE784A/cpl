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
	int N, Nsim = 1e2;

	cpl_function ρ_O  = { woods_saxon_density, &O };
	cpl_function ρ_Cu = { woods_saxon_density, &Cu };
	cpl_function ρ_I  = { woods_saxon_density, &J };
	cpl_function ρ_Au = { woods_saxon_density, &Au };
	cpl_function ρ_Pb = { woods_saxon_density, &Pb };


	/* Plotting a 2d sample */
	cpl_matrix *target = cpl_matrix_alloc(Au.A, 2);
	generate_coords(ρ_Au, target, 6);
	cpl_matrix_savetxt("target.txt", target, "x y");

	cpl_matrix *projectile = cpl_matrix_alloc(Au.A, 2);
	generate_coords(ρ_Au, projectile, -6);
	cpl_matrix_savetxt("projectile.txt", projectile, "x y");

	gnuplot_ctrl *plt;
	plt = gnuplot_init();

	gnuplot_setup(plt, "sample");
	gnuplot_cmd(plt, "set size square");
	gnuplot_cmd(plt, "set xrange [-12:12]");
	gnuplot_cmd(plt, "set xlabel 'x (fm)'");
	gnuplot_cmd(plt, "set yrange [-12:12]");
	gnuplot_cmd(plt, "set ylabel 'y (fm)'");
	gnuplot_cmd(plt, "plot 'target.txt' using 1:2 title 'Target' with points pt 19, 'projectile.txt' using 1:2 title 'Projectile' with points pt 19");

	gnuplot_close(plt);

	cpl_free(projectile);
	cpl_free(target);

	/* Plotting a 3d sample */
	/*target = cpl_matrix_alloc(Au.A, 3);
	generate_coords(ρ_Au, target, 6);
	cpl_matrix_savetxt("target3d.txt", target, NULL);

	projectile = cpl_matrix_alloc(Au.A, 3);
	generate_coords(ρ_Au, projectile, -6);
	cpl_matrix_savetxt("projectile3d.txt", projectile, NULL);

	gnuplot_ctrl *p3d;
	p3d = gnuplot_init();
	gnuplot_setup(p3d, "sample3d");
	gnuplot_cmd(p3d, "set xrange [-12:12]");
	gnuplot_cmd(p3d, "set xlabel 'x (fm)'");
	gnuplot_cmd(p3d, "set yrange [-12:12]");
	gnuplot_cmd(p3d, "set ylabel 'y (fm)'");
	gnuplot_cmd(p3d, "set zrange [-1:1]");
	gnuplot_cmd(p3d, "set zlabel 'z (fm)'");
	gnuplot_cmd(p3d, "splot 'target3d.txt' with points pt 19, 'projectile3d.txt' with points pt 19");
	gnuplot_close(p3d); */

	N = 25;

	cpl_function ρ_fn[5] = { ρ_O, ρ_Cu, ρ_I, ρ_Au, ρ_Pb };

	int numcases = 5;

	double ri = 0, rf = 15, rx = ri;
	double h = (rf - ri) / N;

	double σ_NN = 4.2;

	cpl_matrix *Ncoll = cpl_matrix_alloc(N, numcases + 1);
	cpl_matrix *Npart = cpl_matrix_alloc(N, numcases + 1);

	cpl_tuple res;
	for (int i = 1; i <= N; ++i) {
		cpl_set(Ncoll, i, 1, rx);
		cpl_set(Npart, i, 1, rx);

		for (int j = 0; j < numcases; ++j) {
			res = gmc_compute(ρ_fn[j], ρ_fn[j], rx, σ_NN, Nsim);
			cpl_set(Ncoll, i, j + 2, res.first);
			cpl_set(Npart, i, j + 2, res.second);
		}

		rx += h;
	}

	cpl_matrix_savetxt("collisions_gmc.txt", Ncoll, "r O Cu I Au Pb");

	plt = gnuplot_init();
	gnuplot_setup(plt, "collisions_gmc");
	gnuplot_cmd(plt, "set xlabel \"impact parameter (fm)\"");
	gnuplot_cmd(plt, "set xrange [0:15]");
	gnuplot_cmd(plt, "set xtics 2.5");
	gnuplot_cmd(plt, "set ylabel \"N_{coll}\"");
	gnuplot_cmd(plt, "set yrange [0:1400]");
	gnuplot_cmd(plt, "plot for [col=2:6] \"collisions_gmc.txt\" using 1:col with lines title columnheader");

	cpl_matrix_savetxt("participants_gmc.txt", Npart, "r O Cu I Au Pb");

	gnuplot_setup(plt, "participants_gmc");
	gnuplot_cmd(plt, "set ylabel \"N_{part}\"");
	gnuplot_cmd(plt, "set yrange [0:450]");
	gnuplot_cmd(plt, "plot for [col=2:6] \"participants_gmc.txt\" using 1:col with lines title columnheader");

	return 0;
}
