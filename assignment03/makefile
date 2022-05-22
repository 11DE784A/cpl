LINKS=$(addprefix -l, m dl)
CFLAGS=$(addprefix -, Wall g O3)
CC=gcc $(LINKS) $(CFLAGS)

VPATH=../code/

OBJECTS=main.o cpl_arrays.o cpl_commons.o cpl_linalg.o cpl_format.o cpl_random.o \
		cpl_stats.o cpl_diffeq.o cpl_integration.o \
		../gnuplot/gnuplot_i.o

main: $(OBJECTS)
	$(CC) -o main $^

HEADERS=cpl_defines.h cpl_includes.h cpl_bench.h cpl_commons.h cpl_arrays.h \
		cpl_linalg.h cpl.h cpl_format.h cpl_random.h cpl_diffeq.h cpl_stats.h \
		cpl_integration.h \
		../gnuplot/gnuplot_i.h

main.o: main.c $(HEADERS)
	$(CC) -c $<

