# Makefile SoSe 2018 NumIng #

CC = gcc

CFLAGS = -Wall -Wextra -pedantic -std=c99

all: iteration_tridiag

matrix.o: matrix.c matrix.h 
	$(CC) -c $(CFLAGS) $< -o $@ 

iteration_tridiag: iteration_tridiag.o matrix.o
	$(CC) $(CFLAGS) $^ -o $@ -lm

iteration_tridiag.o: iteration_tridiag.c matrix.h
	$(CC) $(CFLAGS) -c $< -o $@



