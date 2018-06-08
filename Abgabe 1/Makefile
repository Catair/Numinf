# Makefile SoSe 2018 NumIng #
 
CC = gcc

CFLAGS = -Wall -Wextra -pedantic -std=c99


all: gnomesort mergesort

gnomesort: gnomesort.o 
	$(CC) -c $(CFLAGS) $< -o $@ -lm

gnomesort.o: gnomesort.c
	$(CC) -c $(CFLAGS) $< -o $@

mergesort: mergesort.o
	$(CC) -c $(CFLAGS) $< -o $@ -lm

mergesort.o: mergesort.c
	$(CC) -c $(CFLAGS) $< -o $@



