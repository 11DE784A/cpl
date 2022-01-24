CC=gcc -Wall -g -O0
COMMONS=cpl_defines.h cpl_includes.h cpl.h

main: main.o cpl_tuple.o cpl_tensor.o cpl_linalg.o
	$(CC) -o main main.o cpl_tuple.o cpl_tensor.o cpl_linalg.o

main.o: main.c $(COMMONS)
	$(CC) -c main.c

cpl_tuple.o: cpl_tuple.c $(COMMONS)
	$(CC) -c cpl_tuple.c

cpl_tensor.o: cpl_tensor.c cpl_tuple.c $(COMMONS)
	$(CC) -c cpl_tensor.c

cpl_linalg.o: cpl_linalg.c cpl_tensor.c cpl_tuple.c $(COMMONS)
	$(CC) -c cpl_linalg.c
