# Makefile SoSe 2018 NumIng #

CC = gcc

CFLAGS = -Wall -Wextra -pedantic -std=c99

all: qr_solve

matrix.o: matrix.c matrix.h 
	$(CC) -c $(CFLAGS) $< -o $@ 

qr_solve: qr_solve.o matrix.o
	$(CC) $(CFLAGS) $^ -o $@ -lm

qr_solve.o: qr_solve.c matrix.h
	$(CC) $(CFLAGS) -c $< -o $@



