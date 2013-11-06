EXT_ROOT=../external
SMC_ROOT=$(EXT_ROOT)/sparse_matrix_converter
BBU_ROOT=$(EXT_ROOT)/bebop_util
BBM_ROOT=$(EXT_ROOT)/bebop_make

include $(BBM_ROOT)/options
include $(BBM_ROOT)/Makefile.include
include $(BBM_ROOT)/ecl
include $(BBM_ROOT)/setup

## Set up library and include paths
CPPFLAGS += -I$(BBU_ROOT)/include  -I$(SMC_ROOT)/include -std=c99
LIBS += -L$(BBU_ROOT) -lbebop_util -L$(SMC_ROOT) -lsparse_matrix_converter

## Filename of static library
LIBSMC_STATIC = $(SMC_ROOT)/libsparse_matrix_converter$(STATIC_LIB_EXTENSION)

## Filename of dynamic library
LIBSMC_DYNAMIC = $(SMC_ROOT)/libsparse_matrix_converter$(DYNAMIC_LIB_EXTENSION)

EXEOUT = fennel
EXEOBJ = fennel.o  

## Location of SMC header files
HEADERDIR = $(SMC_ROOT)/include/bebop/smc

#########################################################################
## GENERAL BUILD RULES 
#########################################################################

.o: 
	$(CC) $(CPPFLAGS) $(CFLAGS) -c $< -o $@

#########################################################################
## OUTPUT RULES
#########################################################################

.PHONY: all

all: $(EXEOUT)

fennel.o: fennel.c $(HEADERDIR)/sparse_matrix.h $(HEADERDIR)/sparse_matrix_ops.h 

$(EXEOUT): $(LIBSMC_STATIC) fennel.o
	$(LINKER) $(CPPFLAGS) $(CFLAGS) -o $(EXEOUT) $(EXEOBJ) \
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
