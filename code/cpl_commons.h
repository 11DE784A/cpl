typedef struct {
	double (*function) (double, void *);
	void *params;
} cpl_function;

#define cpl_fn_eval(F, x) F.function(x, F.params)

void cpl_check(int, char*);

int cpl_delta(int, int);

