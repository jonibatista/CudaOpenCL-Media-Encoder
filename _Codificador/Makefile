# MAKEFILE for mmp.pt
CC 	= gcc
# ORIGINAL
# CFLAGS 	= -g -pg -I$(LIBDIR) -I$(LIBDIR1) -Wall
#gcc params.c -c 
#gcc -o saida params.o Encoder.o
CFLAGS 	= -g -ffloat-store -O3 -Wno-unused-result -Wall #-pg #-O3
LIBRARY = -lm #-lefence #-lmpr 
PROGRAM = codvector 
MAIN	= main.c 
MAIN_OBJ = main.o
SOURCES1 = params.c
OBJECTS1 = params.o


###################### ENCODER ###################### 
# nome do executavel 

# Nome do ficheiro de opcoes do gengetopt (caso exista)
#SHOWSTATS_OPT=showstats_opt
#params.o: params.c params.h

#como compilar .o a partir de .c
#.c.o:
#	${CC} -c ${SOURCES1} # $<

params:	$(OBJECTS1)
	$(CC) -c ${SOURCES1}


std:	$(MAIN_OBJ) 
	$(CC) $(CFLAGS) $(MAIN_OBJ) $(OBJECTS1) -o $(PROGRAM) $(LIBRARY)

all:	$(MAIN_OBJ)
	$(CC) $(CFLAGS) $(MAIN_OBJ) $(OBJECTS1) -o $(PROGRAM) $(LIBRARY)

clean:
	rm -f $(MAIN_OBJ) $(PROGRAM) *~ *.o

lint:
	lint  $(MAIN) $(SOURCES1)

