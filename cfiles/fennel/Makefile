CC=gcc
CXX=gcc
RM=rm -f
MEX=mkoctfile

SMC_ROOT=../external/sparse_matrix_converter
BBU_ROOT=../external/bebop_util

CPPFLAGS= -std=c99 -g -O3 $(BBI)
LDFLAGS= -g -Wl,-rpath $(BBL)
BBI= -I$(SMC_ROOT)/include -I$(BBU_ROOT)/include
BBL= -L$(SMC_ROOT) -L$(BBU_ROOT) -lbebop_util -lsparse_matrix_converter

SRCS= fennel.c csr.c timer.c matvec.c
OBJS= fennel.o csr.o timer.o matvec.o

fennel: $(OBJS)
	$(CXX) $(LDFLAGS) -o fennel $(OBJS) $(LDLIBS) 

fennel.o: $(SRCS) $(OBJS)
	$(CXX) $(CPPFLAGS) -c fennel.c

csr.o: $(SRCS)
	$(CXX) $(CPPFLAGS) -c csr.c

timer.o: $(SRCS)
	$(CXX) $(CPPFLAGS) -c timer.c

matvec.o: $(SRCS)
	$(CXX) $(CPPFLAGS) -c matvec.c

clean:
	$(RM) $(OBJS) fennel
