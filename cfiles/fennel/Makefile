CC=gcc
CXX=g++
RM=rm -f

SMC_ROOT=../external/sparse_matrix_converter
BBU_ROOT=../external/bebop_util

CPPFLAGS= -g -O3 $(BBI)
LDFLAGS= -g -Wl,-rpath $(BBL)
BBI= -I$(SMC_ROOT)/include -I$(BBU_ROOT)/include
BBL= -L$(SMC_ROOT) -L$(BBU_ROOT) -lbebop_util -lsparse_matrix_converter

SRCS= fennel.cc
OBJS= fennel.o

fennel: $(OBJS)
	$(CXX) $(LDFLAGS) -o fennel fennel.o $(LDLIBS) 

fennel.o: $(SRCS)
	$(CXX) $(CPPFLAGS) -c fennel.cc

clean:
	$(RM) $(OBJS) fennel
