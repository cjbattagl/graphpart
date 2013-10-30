CC=gcc
CXX=g++
RM=rm -f

CPPFLAGS= -g -O2 
LDFLAGS= -g 

SRCS= fennel.cc
OBJS= fennel.o

fennel: $(OBJS)
	$(CXX) $(LDFLAGS) -o fennel fennel.o $(LDLIBS) 

fennel.o: $(SRCS)
	$(CXX) $(CPPFLAGS) -c fennel.cc
