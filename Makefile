# Makefile surface smoothing

CC      = gcc
INCLUDE = -I/usr/local/include 
LDFLAGS = -L/usr/lib -lm  
CFLAGS  = -O2 $(INCLUDE)
OBJ     = spriorsmooth.o ssmooth.o 
PROG    = spriorsmooth



all  : spriorsmooth

spriorsmooth  : $(OBJ)
	$(CC) $(CFLAGS)  $(OBJ) -o $(PROG)  $(LDFLAGS)

ssmooth.o              : ssmooth.c ssmooth.h
spriorsmooth.o         : ssmooth.h spriorsmooth.c
 


