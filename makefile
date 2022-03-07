CC=gcc -ldl -lm -Wall -g -O0
HEADERS=cpl_defines.h cpl_includes.h cpl_bench.h cpl_commons.h cpl_arrays.h cpl_linalg.h cpl.h cpl_format.h cpl_eigen.h

LINKS=main.o cpl_arrays.o cpl_commons.o cpl_linalg.o cpl_format.o cpl_eigen.o
main: $(LINKS)
	$(CC) -o main $(LINKS)

main.o: main.c $(HEADERS)
	$(CC) -c main.c

cpl_commons.o: cpl_commons.c $(HEADERS)
	$(CC) -c cpl_commons.c

cpl_format.o: cpl_format.c $(HEADERS)
	$(CC) -c cpl_format.c

cpl_arrays.o: cpl_arrays.c $(HEADERS)
	$(CC) -c cpl_arrays.c

cpl_linalg.o: cpl_linalg.c $(HEADERS)
	$(CC) -c cpl_linalg.c

cpl_eigen.o: cpl_eigen.c $(HEADERS)
	$(CC) -c cpl_eigen.c

