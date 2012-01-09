CUDA_INSTALL_PATH ?= /usr/local/cuda

CXX := g++ 
CC := gcc 
LINK := g++ -fPIC 
NVCC  := nvcc --compiler-options -fno-strict-aliasing -O3 -arch sm_20

# Includes
INCLUDES = -I. -I$(CUDA_INSTALL_PATH)/include

# Common flags
COMMONFLAGS += $(INCLUDES)
NVCCFLAGS += $(COMMONFLAGS)
CXXFLAGS += $(COMMONFLAGS)
CFLAGS += $(COMMONFLAGS)

SOURCES1 = params.c
OBJECTS1 = params.o

LIB_CUDA := -L$(CUDA_INSTALL_PATH)/lib -lcudart


OBJS = main.cu.o params.o
TARGET = codvector
LINKLINE = $(LINK) -o $(TARGET) $(OBJS) $(LIB_CUDA)


std: $(OBJS)
	$(LINKLINE)

.SUFFIXES: .c .cu .o

params.o: params.c
	gcc -lm -c params.c

%.cu.o: %.cu
	$(NVCC) $(NVCCFLAGS) -c $< -o $@

docs: Doxyfile
	doxygen Doxyfile
clean:
	rm -f $(BIN) *.o ${TARGET}

gengetopt:
	gengetopt --file-name=params --conf-parser --unamed-opts < params_opt.ggo 

indent:
	indent -bad -bap -bbb -br -bli1 -brf -brf -bs -ce -cdw -cs -i8 -l160 -lc160 -nsob -npsl -sc *.c *.cu


