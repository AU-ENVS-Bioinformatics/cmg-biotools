###############################################################
# Makefile for SQUID's testsuite
# CVS $Id: Makefile.in,v 1.6 2003/05/23 16:16:42 eddy Exp $
#
# Note: The autoconf variables in this file must be coordinated
#       with HMMER, if you change them, because HMMER will 
#       create a Makefile from this Makefile.in using its own
#       configure script, not SQUID's.
#
###########
# HMMER - Biological sequence analysis with profile HMMs
# Copyright (C) 1992-2003 Washington University School of Medicine
# All Rights Reserved
# 
#     This source code is distributed under the terms of the
#     GNU General Public License. See the files COPYING and LICENSE
#     for details.
###########

SHELL  = /bin/sh
CC     = gcc
CFLAGS = -g -O2
DEFS   = -DHAVE_CONFIG_H
LIBS   = -lsquid  -lm

TESTPROGS = iospeed rndspeed

X-BASE = x-base-afetch\
	x-base-alistat\
	x-base-seqstat\
	x-base-sfetch\
	x-base-shuffle\
	x-base-sindex\
	x-base-sreformat

BUGLIST = bug-1-sfetch-paths

################################################################
# Targets for building the test programs.
#
all: sqdconfig $(TESTPROGS) 

$(TESTPROGS): %: %_main.o 
	$(CC) $(CFLAGS) -L../ -I../ $(DEFS) -o $@ $@_main.o $(LIBS)

sqdconfig: %: %_main.o 
	$(CC) $(CFLAGS) -L../ -I../ $(DEFS) -o $@ $@_main.o $(LIBS)

.c.o:
	$(CC) $(CFLAGS) $(DEFS) -L../ -I../ -c $<		

################################################################


################################################################
# 'make check' actually runs the tests.
#
check: $(TESTPROGS) sqdconfig
	@for testprog in $(TESTPROGS); do\
	   if ./$$testprog; then\
	      echo $$testprog: ok;\
	   else\
	      echo $$testprog: FAILED;\
	   fi;\
        done
	@for xprog in $(X-BASE); do\
	   ./$$xprog;\
	done        
	@for bugprog in $(BUGLIST); do\
	   ./$$bugprog;\
	done        
  
#######
## Miscellaneous
#######

clean:
	-rm -f *.o *~ Makefile.bak core $(TESTPROGS) TAGS gmon.out sqdconfig

distclean:
	make clean
	-rm -f Makefile

binclean:
	-rm -f *.o *~ Makefile.bak core TAGS gmon.out sqdconfig

TAGS:
	etags -t *.c *.h Makefile.in


