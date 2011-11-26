CUDA_INSTALL_PATH ?= /usr/local/cuda

#LD = `$LD_LIBRARY_PATH | awk '{print substr($1,2)}'`

CXX := g++ 
CC := gcc 
LINK := g++ -fPIC 
NVCC  := nvcc --compiler-options -fno-strict-aliasing -O3

# Includes
INCLUDES = -I. -I$(CUDA_INSTALL_PATH)/include

# Common flags
COMMONFLAGS += $(INCLUDES)
NVCCFLAGS += $(COMMONFLAGS)
CXXFLAGS += $(COMMONFLAGS)
CFLAGS += $(COMMONFLAGS)

SOURCES1 = params.c
OBJECTS1 = params.o




LIB_CUDA := -L$(CUDA_INSTALL_PATH)/lib64 -lcudart
#LIB_CUDA := -L$(LD) -lcudart


OBJS = main.cu.o
TARGET = codvector
LINKLINE = $(LINK) -o $(TARGET) $(OBJS) $(OBJECTS1) $(LIB_CUDA)

.SUFFIXES: .c .cu .o

params:	$(OBJECTS1)
	$(CC) -c ${SOURCES1}

%.c.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

%.cu.o: %.cu
	$(NVCC) $(NVCCFLAGS) -c $< -o $@

$(TARGET): $(OBJS) Makefile
	$(LINKLINE)

clean:
	rm -f $(BIN) *.o *.cu_o ${TARGET}

docs: Doxyfile
	doxygen Doxyfile

indent:
	indent -bad -bap -bbb -br -bli1 -brf -brf -bs -ce -cdw -cs -i8 -l160 -lc160 -nsob -npsl -sc *.c *.cu


