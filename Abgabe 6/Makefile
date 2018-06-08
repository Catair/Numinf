# Makefile SoSe 2018 NumIng #

CC = gcc

CFLAGS = -Wall -Wextra -pedantic -std=c99

all: newton_heron

newton_heron: newton_heron.o
	$(CC) -g $< -o $@ -lm

newton_heron.o: newton_heron.c
	$(CC) -g $(CFLAGS) -c $< -o $@



