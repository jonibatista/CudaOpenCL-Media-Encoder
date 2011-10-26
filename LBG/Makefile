# MAKEFILE for mmp.pt
CC 	= gcc
# ORIGINAL
# CFLAGS 	= -g -pg -I$(LIBDIR) -I$(LIBDIR1) -Wall
# CFLAGS 	= -g -ffloat-store -O3 #-pg #-O3
CFLAGS 	= -g -ffloat-store -O3 -Wno-unused-result -Wall #-pg #-O3
LIBRARY = -lm #-lefence #-lmpr 
PROGRAM = vecquant 
MAIN	= main.c 
MAIN_OBJ= main.o
SOURCES1 = 
OBJECTS1= 

std:	$(MAIN_OBJ) 
	$(CC) $(CFLAGS) $(MAIN_OBJ) $(OBJECTS1) -o $(PROGRAM) $(LIBRARY)

all:	$(MAIN_OBJ)
	$(CC) $(CFLAGS) $(MAIN_OBJ) $(OBJECTS1) -o $(PROGRAM) $(LIBRARY)

clean:
	rm -f $(MAIN_OBJ)

lint:
	lint  $(MAIN) $(SOURCES1)

