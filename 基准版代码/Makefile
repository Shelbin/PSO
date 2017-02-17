#!/bin/bash

CC=icc
CFLAGS=  -O2  -g   
LDFLAGS=  -O2  -g 

EXEC=pso
OBJ1=pso.o 

all: $(EXEC)

pso : $(OBJ1)
	$(CC) -o $@ $^ $(LDFLAGS)

 
##########################################
# Generic rules
##########################################

%.o: %.c %.h
	$(CC) -o $@ -c $< $(CFLAGS)
	
%.o: %.c
	$(CC) -o $@ -c $< $(CFLAGS)
	
clean:
	rm -f *.o *~ $(EXEC)
