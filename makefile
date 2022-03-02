CC=gcc -ldl -lm -Wall -g -O0
HEADERS=cpl_defines.h cpl_includes.h cpl_bench.h cpl_commons.h cpl_arrays.h cpl_linalg.h cpl.h

main: main.o cpl_arrays.o cpl_commons.o cpl_format.c cpl_linalg.o
	$(CC) -o main main.o cpl_arrays.o cpl_commons.o cpl_linalg.o

main.o: main.c $(HEADERS)
	$(CC) -c main.c

cpl_commons.o: cpl_commons.c $(HEADERS)
	$(CC) -c cpl_commons.c

cpl_arrays.o: cpl_arrays.c $(HEADERS)
	$(CC) -c cpl_arrays.c

cpl_linalg.o: cpl_linalg.c $(HEADERS)
	$(CC) -c cpl_linalg.c
