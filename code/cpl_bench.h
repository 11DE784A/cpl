#ifndef CPL_BENCH

#define CPL_BENCH

#include "cpl_format.h"

#define cpl_time(process) \
	({clock_t _tstart = clock(); \
	  int _nstart = NALLOC; \
	  int _bstart = BALLOC; \
	  process; \
	  double _tdelta = (double) (clock() - _tstart) / CLOCKS_PER_SEC; \
	  int _ndelta = NALLOC - _nstart; \
	  int _bdelta = BALLOC - _bstart; \
	  char _ftdelta[FMTWIDTH]; strfsecs(_ftdelta, _tdelta); \
	  char _fbdelta[FMTWIDTH]; strfbytes(_fbdelta, _bdelta); \
	  printf("%s:%d " #process ": %s (%d allocations, %s)\n\n", \
			  __FILE__, __LINE__, _ftdelta, _ndelta, _fbdelta);})

static int NALLOC = 0;
static int BALLOC = 0;

void* malloc(size_t sz) {
    void *(*libc_malloc)(size_t) = dlsym(RTLD_NEXT, "malloc");
	NALLOC += 1;
	BALLOC += sz;
    return libc_malloc(sz);
}

void* calloc(size_t n, size_t sz) {
    void *(*libc_calloc)(size_t, size_t) = dlsym(RTLD_NEXT, "calloc");
	NALLOC += 1;
	BALLOC += sz;
    return libc_calloc(n, sz);
}

void* realloc(void *ptr, size_t sz) {
    void *(*libc_realloc)(void*, size_t) = dlsym(RTLD_NEXT, "realloc");
	NALLOC += 1;
	BALLOC += sz;
    return libc_realloc(ptr, sz);
}

#endif
