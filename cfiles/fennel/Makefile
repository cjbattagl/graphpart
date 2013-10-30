CC=gcc
CXX=g++
RM=rm -f

CPPFLAGS= -g -O2 $(BBI)
LDFLAGS= -g $(BBL) 
BBI= -I../external/sparse_matrix_converter/include -I../external/bebop_util/include
BBL= -L../external/sparse_matrix_converter -lbebop_util -L../external/bebop_util

SRCS= fennel.cc
OBJS= fennel.o

fennel: $(OBJS)
	$(CXX) $(LDFLAGS) -o fennel fennel.o $(LDLIBS) 

fennel.o: $(SRCS)
	$(CXX) $(CPPFLAGS) -c fennel.cc

clean:
	$(RM) $(OBJS) fennel
