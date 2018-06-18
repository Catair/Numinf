# Makefile SoSe 2018 NumIng #
 
CC = gcc

CFLAGS = -Wall -Wextra -pedantic -std=c99

all: interpolation

interpolation: interpolation.o
	$(CC) $(CFLAGS) $^ -o $@ -lm

interpolation.o: interpolation.c 
	$(CC) $(CFLAGS) -c $< -o $@



