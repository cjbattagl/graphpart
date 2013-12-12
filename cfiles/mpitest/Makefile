EXECS=test
MPICC?=mpicc
OPTS= -std=c99 -O3

all: ${EXECS}

test: test.c
        ${MPICC} ${OPTS} -o test test.c
clean:
        rm ${EXECS}