#ifndef CPL_DEFINES

#define CPL_DEFINES

#define xstr(x) str(x)
#define str(x) #x

#define scalar double

#define sfmt % -#9.4g

#define TOL 1e-6
#define MAXLINE 150
#define MAX_ITERS 20000
#define PI 3.1415926535897

#define sabs(x) _Generic((x), \
	float: fabsf, \
	double: fabs, \
	default: abs \
)(x)

#endif
