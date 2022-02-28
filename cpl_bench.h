#ifndef CPL_BENCH

#define CPL_BENCH

#define time(process) \
	({clock_t _tstart = clock(); \
	  int _nstart = NALLOC; \
	  int _bstart = BALLOC; \
	  process; \
	  double _tdelta = (double) (clock() - _tstart) / CLOCKS_PER_SEC; \
	  int _ndelta = NALLOC - _nstart; \
	  int _bdelta = BALLOC - _bstart; \
	  printf("%s:%d " #process ": %f seconds (%d allocations, %d bytes)\n", \
			  __FILE__, __LINE__, _tdelta, _ndelta, _bdelta);})

int NALLOC = 0;
int BALLOC = 0;

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

#endif
