# Makefile SoSe 2018 NumIng #

CC = gcc

CFLAGS = -Wall -Wextra -pedantic -std=c99

all: iteration

matrix.o: matrix.c matrix.h 
	$(CC) -c $(CFLAGS) $< -o $@ 

iteration: iteration.o matrix.o
	$(CC) $(CFLAGS) $^ -o $@ -lm

iteration.o: iteration.c matrix.h
	$(CC) $(CFLAGS) -c $< -o $@



