# Makefile
# Author: Mark Hoemmen
# Since: 05 Feb 2008
# Time-stamp: <2008-07-16 15:42:30 mhoemmen>
#
# For a sparse matrix utility that converts between various file formats 
# for sparse matrices, and also between various internal storage formats
# for sparse matrices.  Harwell-Boeing, MatrixMarket and Matlab file
# formats are supported.
#
# Copyright (c) 2008, Regents of the University of California 
# All rights reserved.
# Redistribution and use in source and binary forms, with or
# without modification, are permitted provided that the
# following conditions are met:
# 
# * Redistributions of source code must retain the above copyright
#   notice, this list of conditions and the following disclaimer.
#
# * Redistributions in binary form must reproduce the above copyright 
#   notice, this list of conditions and the following disclaimer in 
#   the documentation and/or other materials provided with the 
#   distribution.
#
# * Neither the name of the University of California, Berkeley, nor
#   the names of its contributors may be used to endorse or promote
#   products derived from this software without specific prior
#   written permission.  
# 
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
# "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
# LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
# FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
# COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
# INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
# (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
# SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
# HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
# STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED
# OF THE POSSIBILITY OF SUCH DAMAGE.
##########################################################################

include ../../bebop_make/options
include ../../bebop_make/Makefile.include
include ../../bebop_make/ecl
include ../../bebop_make/setup


## Set up library and include paths
CPPFLAGS += -I../../bebop_util/include -I../include
LIBS += -L../../bebop_util -lbebop_util


## Location of SMC header files
SMC_HEADERDIR = ../include/bebop/smc


## Where to put the generated header files
SMC_GENHEADER_INSTALLDIR=$(SMC_HEADERDIR)/lowlevel


## Source files (intermediate) to generate
GENSOURCE = csr_to_dense.c dense_setzero.c dense_to_csr.c dense_to_csr_nnz.c

## Header files (not intermediate) to generate
GENHEADERS = $(patsubst %.c,$(SMC_GENHEADER_INSTALLDIR)/%.h,$(GENSOURCE))

## Object files to generate
OBJ = $(patsubst %.c,%.o,$(GENSOURCE))

PYTHON = python
PYGENSCRIPT = ../scripts/sparsegen.py
GENERATE = $(PYTHON) $(PYGENSCRIPT)


#########################################################################
## GENERAL BUILD RULES 
#########################################################################

%.o : %.c
	$(CC) $(CPPFLAGS) $(CFLAGS) -c $< -o $@

% :: %.tmpl
	$(GENERATE) < $< > $@


#########################################################################
## OUTPUT RULES
#########################################################################

.PHONY: all clean install

all: $(OBJ) 

install: $(GENHEADERS)
	for f in "$(GENHEADERS)"; do cp $$f $(SMC_GENHEADER_INSTALLDIR); done


#########################################################################
## SOURCE FILE DEPENDENCIES
#########################################################################

csr_to_dense.o: csr_to_dense.c dense_setzero.h

dense_setzero.o: dense_setzero.c dense_setzero.h

dense_to_csr.o: dense_to_csr.c dense_to_csr.h

dense_to_csr_nnz.o: dense_to_csr_nnz.c dense_to_csr_nnz.h


##########################################################################
## UTILITY STEPS
##########################################################################

ctags: 
	ctags *.[ch]

etags:
	etags *.[ch]

clean:
	$(RM) $(GENSOURCE) $(GENHEADERS) $(OBJ)

