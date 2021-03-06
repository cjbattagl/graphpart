EXT_ROOT=../external
SMC_ROOT=$(EXT_ROOT)/sparse_matrix_converter
BBU_ROOT=$(EXT_ROOT)/bebop_util
BBM_ROOT=$(EXT_ROOT)/bebop_make
#MET_ROOT=$(EXT_ROOT)/metis-5.1.0

MPICC=mpicc

include $(BBM_ROOT)/options
include $(BBM_ROOT)/Makefile.include
include $(BBM_ROOT)/ecl
include $(BBM_ROOT)/setup

## Set up library and include paths
CPPFLAGS += -I$(BBU_ROOT)/include  -I$(SMC_ROOT)/include -std=c99 #-I$(MET_ROOT)/include
LIBS += -L$(BBU_ROOT) -lbebop_util -L$(SMC_ROOT) -lsparse_matrix_converter #-L$(MET_ROOT)/build/Darwin-i386/libmetis -lmetis

## Filename of static library
LIBSMC_STATIC = $(SMC_ROOT)/libsparse_matrix_converter$(STATIC_LIB_EXTENSION)

## Filename of dynamic library
LIBSMC_DYNAMIC = $(SMC_ROOT)/libsparse_matrix_converter$(DYNAMIC_LIB_EXTENSION)

EXEOUT = fennel mpifennel
FEXEOBJ = fennel.o fennel_kernel.o util.o randperm.o 
MEXEOBJ = mpifennel.o fennel_kernel.o util.o randperm.o 

## Location of SMC header files
HEADERDIR = $(SMC_ROOT)/include/bebop/smc

#########################################################################
## GENERAL BUILD RULES 
#########################################################################

.o: 
	$(CC) $(CPPFLAGS) $(CFLAGS) -c $< -o $@

.o.cc: 
	$(CPP) $(CPPFLAGS) $(CFLAGS) -c $< -o $@

#########################################################################
## OUTPUT RULES
#########################################################################

.PHONY: all

all: $(EXEOUT)

mpifennel.o: mpifennel.c $(HEADERDIR)/sparse_matrix.h $(HEADERDIR)/sparse_matrix_ops.h
	$(MPICC) $(CPPFLAGS) $(CFLAGS) -c $< -o $@

fennel.o: fennel.c $(HEADERDIR)/sparse_matrix.h $(HEADERDIR)/sparse_matrix_ops.h

fennel_kernel.o: fennel_kernel.c $(HEADERDIR)/sparse_matrix.h $(HEADERDIR)/sparse_matrix_ops.h

randperm.o: randperm.c

util.o: util.c

fennel: $(LIBSMC_STATIC) fennel.o fennel_kernel.o randperm.o util.o
	$(LINKER) $(CPPFLAGS) $(CFLAGS) -o fennel $(FEXEOBJ) \
		$(LIBSMC_STATIC) $(LIBS)

mpifennel: $(LIBSMC_STATIC) mpifennel.o fennel_kernel.o randperm.o util.o
	mpicc $(CPPFLAGS) $(CFLAGS) -o mpifennel $(MEXEOBJ) \
	$(LIBSMC_STATIC) $(LIBS)

	
##########################################################################
## UTILITY STEPS
##########################################################################

ctags: 
	ctags *.[ch]

etags:
	etags *.[ch]
	
clean:
	$(RM) $(LIBSMCOBJ) 
	$(RM) $(EXEOBJ)
	$(RM) $(EXEOUT)
