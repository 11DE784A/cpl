CC=gcc -Wall -g -O0
HEADERS=cpl_defines.h cpl_includes.h cpl.h

main: main.o cpl_tuple.o cpl_tensor.o cpl_linalg.o cpl_commons.o
	$(CC) -o main main.o cpl_tuple.o cpl_tensor.o cpl_linalg.o cpl_commons.o

main.o: main.c $(HEADERS)
	$(CC) -c main.c

cpl_commons.o: cpl_commons.c
	$(CC) -c cpl_commons.c

cpl_tuple.o: cpl_tuple.c cpl_commons.c $(HEADERS)
	$(CC) -c cpl_tuple.c

cpl_tensor.o: cpl_tensor.c cpl_tuple.c cpl_commons.c $(HEADERS)
	$(CC) -c cpl_tensor.c

cpl_linalg.o: cpl_linalg.c cpl_tensor.c cpl_tuple.c cpl_commons.c $(HEADERS)
	$(CC) -c cpl_linalg.c
