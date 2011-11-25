# MAKEFILE for mmp.pt
CC 	= gcc
# ORIGINAL
# CFLAGS 	= -g -pg -I$(LIBDIR) -I$(LIBDIR1) -Wall
CFLAGS 	= -g -ffloat-store -O3 -Wno-unused-result -Wall #-pg #-O3
LIBRARY = -lm #-lefence #-lmpr 
PROGRAM = codvector 
MAIN	= main.c 
MAIN_OBJ = main.o
SOURCES1 = params.c
OBJECTS1 = params.o


params:	$(OBJECTS1)
	$(CC) -c ${SOURCES1}

std:	$(MAIN_OBJ) 
	$(CC) $(CFLAGS) $(MAIN_OBJ) $(OBJECTS1) -o $(PROGRAM) $(LIBRARY)

all:	$(MAIN_OBJ)
	$(CC) $(CFLAGS) $(MAIN_OBJ) $(OBJECTS1) -o $(PROGRAM) $(LIBRARY)

clean-all:
	rm -f $(MAIN_OBJ) $(PROGRAM) *~ *.o

cleano:
	rm -f *~ *.o #$(MAIN_OBJ)

docs: Doxyfile
	doxygen Doxyfile
lint:
	lint  $(MAIN) $(SOURCES1)

