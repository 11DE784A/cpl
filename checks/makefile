LINKS=$(addprefix -l, m dl check)
CFLAGS=$(addprefix -, Wall g O0)
CC=gcc $(LINKS) $(CFLAGS)

VPATH=../code

OBJECTS=main.o cpl_arrays.o cpl_commons.o cpl_linalg.o cpl_format.o \
		cpl_eigen.o cpl_stats.o cpl_random.o

main: $(OBJECTS)
	$(CC) -o main $^

HEADERS=cpl_defines.h cpl_includes.h cpl_bench.h cpl_commons.h cpl_arrays.h \
		cpl_linalg.h cpl.h cpl_format.h cpl_eigen.h cpl_stats.h

main.o: main.c $(HEADERS)
	$(CC) -c $<

cpl_commons.o: cpl_commons.c $(HEADERS)
	$(CC) -c $<

cpl_format.o: cpl_format.c $(HEADERS)
	$(CC) -c $<

cpl_arrays.o: cpl_arrays.c $(HEADERS)
	$(CC) -c $<

cpl_linalg.o: cpl_linalg.c $(HEADERS)
	$(CC) -c $<

cpl_eigen.o: cpl_eigen.c $(HEADERS)
	$(CC) -c $<

cpl_stats.o: cpl_stats.c $(HEADERS)
	$(CC) -c $<

cpl_random.o: cpl_random.c $(HEADERS)
	$(CC) -c $<
