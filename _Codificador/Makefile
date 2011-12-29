CUDA_INSTALL_PATH ?= /usr/local/cuda
#gcc -W -Wall -I/usr/local/cuda/include $1.cl -o $1.exe -lOpenCL
CPP := g++ 
CC := gcc 

# Includes
INCLUDES = -I. -I$(CUDA_INSTALL_PATH)/include

# Common flags
COMMONFLAGS += $(INCLUDES)
NVCCFLAGS += $(COMMONFLAGS)
CXXFLAGS += $(COMMONFLAGS)
CFLAGS += $(COMMONFLAGS)

SOURCES1 = params.c
OBJECTS1 = params.o


OBJS = main.o
TARGET = codvector
LINKLINE = $(CPP) -o $(TARGET) $(OBJS) $(OBJECTS1) $(INCLUDES) -lOpenCL

.SUFFIXES: .c .o

params:	params.c
	gcc -c params.c

%.c.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

$(TARGET): $(OBJS) Makefile
	$(LINKLINE)

all: $(OBJS) Makefile
	$(LINKLINE)

clean:
	rm -f $(BIN) *.o ${TARGET}

docs: Doxyfile
	doxygen Doxyfile

gengetopt:
	gengetopt --file-name=params --conf-parser --unamed-opts < params_opt.ggo 

indent:
	indent -bad -bap -bbb -br -bli1 -brf -brf -bs -ce -cdw -cs -i8 -l160 -lc160 -nsob -npsl -sc *.c *.cu


