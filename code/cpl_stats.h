
scalar cpl_stats_mean(cpl_vector*);
scalar cpl_stats_var(cpl_vector*);
scalar cpl_stats_std(cpl_vector*);
scalar cpl_stats_varm(cpl_vector*, scalar);
scalar cpl_stats_stdm(cpl_vector*, scalar);

cpl_tuple cpl_stats_jackknife(scalar (*)(scalar), cpl_vector*);

scalar cpl_stats_linfit(cpl_vector*, cpl_vector*, cpl_vector*, cpl_vector*, cpl_matrix*);
scalar cpl_stats_polyfit(int, cpl_vector*, cpl_vector*, cpl_vector*, cpl_vector*, cpl_matrix*);
