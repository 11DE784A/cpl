CC=gcc

main: main.o cpl_tuple.o cpl_tensor.o
	$(CC) -o main main.o cpl_tuple.o cpl_tensor.o

main.o: main.c cpl.h
	$(CC) -c main.c

cpl_tuple.o: cpl_tuple.c cpl.h
	$(CC) -c cpl_tuple.c

cpl_tensor.o: cpl_tensor.c cpl_tuple.c cpl.h
	$(CC) -c cpl_tensor.c
