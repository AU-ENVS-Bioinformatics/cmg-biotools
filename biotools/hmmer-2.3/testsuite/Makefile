################################################################
# Makefile for HMMER testsuite
# CVS $Id: Makefile.in,v 1.17 2003/05/26 16:21:50 eddy Exp $
##########
# HMMER - Biological sequence analysis with profile HMMs
# Copyright (C) 1992-2003 Washington University School of Medicine
# All Rights Reserved
# 
#     This source code is distributed under the terms of the
#     GNU General Public License. See the files COPYING and LICENSE
#     for details.
###########

CC        = gcc
CFLAGS    = -g -O2
LIBS      =  -lm
MYLIBS    = -lhmmer -lsquid
DEFS      = -DHAVE_CONFIG_H

# Configuration for optional pthreads multiprocessor support
#
PTHREAD_LIBS   = 
PTHREAD_CFLAGS = 

SHELL     = /bin/sh

SHIVA = alignalign_test\
	evd_test\
	masks_test\
	parsingviterbi_test\
	tophits_test\
	trace_test\
	viterbi_exercise\
	weeviterbi_test

#######
## Targets defining how to make Shiva executables.
#######

.c.o:
	$(CC) $(CFLAGS) $(PTHREAD_CFLAGS) $(DEFS) -I../squid -I../src -c $<		

all: 	$(SHIVA)

$(SHIVA): %: %.o
	$(CC) $(CFLAGS) $(PTHREAD_CFLAGS) $(DEFS) -o $@ -L../squid -L../src $@.o $(MYLIBS) $(PTHREAD_LIBS) $(LIBS)

#######
## `make check` actually runs the tests.
#######

check: 
	@echo 
	@echo Running test suite exercises.
	@echo Warning: some tests may take several minutes to complete. 
	@echo 
	./sqc 2 exercises.sqc . ../src


#######
## Miscellaneous
#######

clean:
	-rm -f *.o *~ Makefile.bak core $(SHIVA) TAGS gmon.out

distclean:
	make clean
	-rm -f Makefile

binclean:
	-rm -f *.o *~ Makefile.bak core TAGS gmon.out

TAGS:
	etags -t *.c *.h Makefile.in


