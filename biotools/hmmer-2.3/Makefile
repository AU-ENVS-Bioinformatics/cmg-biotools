################################################################# 
# Makefile for HMMER: Main project directory: installation, documentation
# CVS $Id: Makefile.in,v 1.50 2003/05/27 18:30:07 eddy Exp $
###########
# HMMER - Biological sequence analysis with profile HMMs
# Copyright (C) 1992-2003 Washington University School of Medicine
# All Rights Reserved
# 
#     This source code is distributed under the terms of the
#     GNU General Public License. See the files COPYING and LICENSE
#     for details.
###########

# On most Unices, you can build the package without modifying anything 
#   by just typing "./configure; make".
#
# You may want to modify the following make variables:
#   BINDIR  - where the executables will be installed by a 'make install'
#   MANDIR  - where the man pages will be installed by a 'make install'
#   CC      - which compiler to use
#   CFLAGS  - compiler flags to use

# where you want things installed
# Sort of uses GNU coding standards. ${prefix} might be /usr/local.
# ${exec_prefix} gives you some flexibility for installing architecture
# dependent files (e.g. the programs): an example ${exec_prefix} might be
# /nfs/share/irix64/
#
prefix      = /usr/local
exec_prefix = ${prefix}
BINDIR      = ${exec_prefix}/bin
MANDIR      = ${prefix}/man

# how to install the man pages; 
# cp is generally fine, unless you preformat your pages.
#
INSTMAN   = cp
MANSUFFIX = 1

# your compiler and compiler flags
#
CC     = gcc
CFLAGS = -g -O2

#######
## You should not need to modify below this line.
## Much of it is concerned with maintenance of the development version
## and building the release (indeed, several commands will only work in 
## St. Louis)
#######
SHELL       = /bin/sh
PACKAGE     = HMMER
BASENAME    = hmmer
RELEASE     = 2.3
RELEASEDATE = "May 2003"
COPYRIGHT   = "Copyright (C) 1992-2003 HHMI/Washington University School of Medicine"
LICENSELINE = "Freely distributed under the GNU General Public License (GPL)"
RELCODE     = hmmer2_3
LICENSETAG  = gnu
COMPRESS    = gzip

# The program lists below for HMMER are not necessarily
# a complete manifest. They are the list of stable programs that the
# package will install. There must be a man page for each one of them
# in the appropriate places (documentation/man for HMMER)
#
PROGS = hmmalign\
	hmmbuild\
	hmmcalibrate\
	hmmconvert\
	hmmemit\
	hmmfetch\
	hmmindex\
	hmmpfam\
	hmmsearch

PVMPROGS = 


# all: Compile everything.
#
all: 
	@ if test -d binaries; then\
	  echo "This appears to be a binary distribution; nothing needs to be compiled." ; \
          echo "Use make install to install the precompiled binaries." ; \
	  echo "If you really want to rebuild from source, move or delete the" ;\
	  echo "binaries/ subdirectory."; \
	else\
	  (cd squid;     ${MAKE}; ${MAKE} module);\
	  (cd src;       ${MAKE}; ${MAKE} module);\
	  (cd testsuite; ${MAKE});\
	fi

# check: compiles and runs test suite in testsuite/
#       These are public tests, distributed with the package.
#
check:	
	(cd testsuite; make check)

squid/libsquid.a:
	(cd squid; make; make module)

src/libhmmer.a: 
	(cd src;   make; make module)


# install: installs the binaries in BINDIR/
#          installs man pages in MANDIR/man1/  (e.g. if MANSUFFIX is 1)
install: 
	for file in $(PROGS) $(PVMPROGS); do\
	   cp src/$$file $(BINDIR)/;\
	done
	for file in hmmer $(PROGS); do\
	   $(INSTMAN) documentation/man/$$file.man $(MANDIR)/man$(MANSUFFIX)/$$file.$(MANSUFFIX);\
	done

# uninstall: Reverses the steps of "make install".
#
uninstall: 
	for file in $(PROGS) $(PVMPROGS); do\
	   rm ${BINDIR}/$$file;\
	done
	for file in hmmer $(PROGS); do\
	   rm $(MANDIR)/man$(MANSUFFIX)/$$file.$(MANSUFFIX);\
	done

# bindist: after building, prep the whole directory for a binary
#          distribution: symlink supported binaries in binaries/
#          subdir, remove everything but binaries and Makefiles.
#
bindist:
	mkdir binaries
	for prog in $(PROGS); do\
	   (cd binaries; ln -s ../src/$$prog .);\
	done
	make binclean

# "make clean" removes almost everything except configuration files.
#
clean:
	(cd src;       make clean)
	(cd squid;     make clean)
	(cd testsuite; make clean)
	-rm -f *.o *~ Makefile.bak core TAGS gmon.out

# "make distclean" leaves a pristine source distribution.
#
distclean:
	-rm -rf binaries
	-rm config.log config.status
	-rm -rf autom4te.cache
	-rm -f *.o *~ Makefile.bak core TAGS gmon.out
	-rm Makefile
	(cd src;       make distclean)
	(cd squid;     make distclean)
	(cd testsuite; make distclean)


# "make binclean" is special: it cleans up and leaves only a binary
#       distribution behind, including configured Makefiles.
#
binclean:
	(cd src;       make binclean)
	(cd squid;     make binclean)
	(cd testsuite; make binclean)
	-rm -f *.o *~ Makefile.bak core TAGS gmon.out
	-rm config.log config.status
	-rm -rf autom4te.cache

# doc:  build the Userguide and on-line manual
#
doc:
	(cd documentation/userguide; make)

# dist: build a new distribution directory in hmmer-$RELEASE
#       Exports from the CVS repository.
#       tags RCS files with $(RELCODE) for later reconstruction
#       squid RCS files are tagged with hmmer$(RELCODE).
#       Adds a license statement to each file that has a @ LICENSE @ line.
#       Virtually identical to squid's make dist -- keep them in sync!
dist:	
	# Delete old versions of the same release
	#
	@if test -d ${BASENAME}-$(RELEASE);        then rm -rf ${BASENAME}-$(RELEASE);        fi
	@if test -e ${BASENAME}-$(RELEASE).tar;    then rm -f  ${BASENAME}-$(RELEASE).tar;    fi
	@if test -e ${BASENAME}-$(RELEASE).tar.Z;  then rm -f  ${BASENAME}-$(RELEASE).tar.Z;  fi
	@if test -e ${BASENAME}-$(RELEASE).tar.gz; then rm -f  ${BASENAME}-$(RELEASE).tar.gz; fi
	# 
	# CVS tag and extract.  -F: allow more than one "make dist" per rel
	# prep: must have done "cvs commit", and CVSROOT must be set
	# We also need the squid library, so tag and export it too.
	#
	cvs tag -F $(RELCODE)
	cvs export -r $(RELCODE) -d ${BASENAME}-${RELEASE} ${BASENAME}
	cvs rtag -F ${RELCODE} squid
	(cd ${BASENAME}-${RELEASE}; cvs export -r ${RELCODE} -d squid squid)
	# 
	# We won't include the ssdk in the distro, but we need a 
	# working copy of it temporarily, for rmanprocess.pl 
	# (when making the user guide); we also need to put a copy of sqc 
	# in testsuite. 
	#
	cp ssdk/sqc ${BASENAME}-${RELEASE}/testsuite/
	mkdir ${BASENAME}-${RELEASE}/ssdk
	cp ssdk/rmanprocess.pl ${BASENAME}-${RELEASE}/ssdk/
	# 
	# Make the configure scripts from configure.ac
	#
	(cd ${BASENAME}-${RELEASE};       autoconf)
	(cd ${BASENAME}-${RELEASE}/squid; autoconf)
	#
	# We need a basic squid Makefile, so we can later do
	# a "make implode" in squid. We'll delete it later.
	#
	(cd ${BASENAME}-${RELEASE}/squid; ./configure)
	# 
	# Include the appropriate license files and release notes
	#
	cp  Licenses/LICENSE.$(LICENSETAG)   ${BASENAME}-$(RELEASE)/LICENSE
	cp  Licenses/COPYING.$(LICENSETAG)   ${BASENAME}-$(RELEASE)/COPYRIGHT
	-cp Release-Notes/RELEASE-$(RELEASE) ${BASENAME}-$(RELEASE)/NOTES
	#
	# Attach license stamps on files that need 'em (replace LICENSE keyword everywhere)
	#
	find ${BASENAME}-${RELEASE} -type f -exec ssdk/sedition-pp LICENSE Licenses/$(LICENSETAG) {} \;
	#
	# All other automated keyword replacements are in the documentation
	# or README type files.
	# 
	find ${BASENAME}-${RELEASE}/documentation -type f -exec ssdk/sedition RELEASE ${RELEASE} RELEASEDATE ${RELEASEDATE} PACKAGE ${PACKAGE} COPYRIGHT ${COPYRIGHT} LICENSELINE ${LICENSELINE} {} \;
	ssdk/sedition RELEASE ${RELEASE} RELEASEDATE ${RELEASEDATE} COPYRIGHT ${COPYRIGHT} ${BASENAME}-${RELEASE}/00README
	ssdk/sedition RELEASE ${RELEASE} ${BASENAME}-${RELEASE}/INSTALL
	ssdk/sedition COPYRIGHT ${COPYRIGHT} ${BASENAME}-${RELEASE}/COPYRIGHT
	#
	# Compilation of the documentation.
	# Documentation is not provided in source form.
	# We make it from the distro version, because we've done the proper 
	# keyword substitutions there.
	#
	(cd ${BASENAME}-${RELEASE}/documentation/userguide; make)
	cp -f ${BASENAME}-${RELEASE}/documentation/userguide/Userguide.pdf $(BASENAME)-$(RELEASE)/
	# 
	# Finish the man pages off, adding necessary boilerplate.
	# This must be done *after* the user guide is created, since the
	# user guide steals and uses the pre-boilerplate man pages.
	#
	for prog in $(PROGS); do\
	   cat ${BASENAME}-${RELEASE}/documentation/man/boilerplate >> ${BASENAME}-${RELEASE}/documentation/man/$$prog.man;\
	done
	-rm -f ${BASENAME}-${RELEASE}/documentation/man/boilerplate
	#
	# Remove CVS-controlled files/directories that don't belong in 
	# the distro
	#
	-rm -rf ${BASENAME}-${RELEASE}/Bugs
	-rm -rf ${BASENAME}-${RELEASE}/CHECKLIST
	-rm -rf ${BASENAME}-${RELEASE}/Internal-Notes
	-rm -rf ${BASENAME}-${RELEASE}/Licenses
	-rm -rf ${BASENAME}-${RELEASE}/MAILING_LIST
	-rm -rf ${BASENAME}-${RELEASE}/Release-Notes
	-rm -rf ${BASENAME}-${RELEASE}/configure.ac
	-rm -rf ${BASENAME}-${RELEASE}/autom4te.cache 
	-rm -rf ${BASENAME}-${RELEASE}/checkoptions.pl
	-rm -rf ${BASENAME}-${RELEASE}/licenseadd.pl
	-rm -rf ${BASENAME}-${RELEASE}/nodebuild
	-rm -rf ${BASENAME}-${RELEASE}/rootbuild
	-rm -rf ${BASENAME}-${RELEASE}/buildall
	-rm -rf ${BASENAME}-${RELEASE}/documentation/userguide
	-rm -rf ${BASENAME}-${RELEASE}/squid/Licenses
	-rm -rf ${BASENAME}-${RELEASE}/squid/LOG
	-rm -rf ${BASENAME}-${RELEASE}/squid/configure.ac
	-rm -rf ${BASENAME}-${RELEASE}/ssdk
	(cd ${BASENAME}-${RELEASE}/squid; make implode)
	(cd ${BASENAME}-${RELEASE}/squid; make distclean)
	#
	# Set file permissions
	#
	find ${BASENAME}-${RELEASE} -type f -exec chmod +r {} \;
	chmod +x ${BASENAME}-${RELEASE}/testsuite/sqc
	#
	# Packaging commands
	#
	tar cvf ${BASENAME}-${RELEASE}.tar ${BASENAME}-${RELEASE}
	$(COMPRESS) ${BASENAME}-$(RELEASE).tar


# make ftpdist: install FTP distribution. Assumes a "make dist"
#               has already succeeded. Doesn't do the symlink;
#               doesn't install the on-line manual for the Web site.
#
ftpdist:
	cp -f $(READMES)                     $(FTPDIR)
	cp -f build/hmmer-$(RELEASE)/COPYING $(FTPDIR)
	cp -f build/hmmer-$(RELEASE)/LICENSETAG $(FTPDIR)
	cp -f build/hmmer-$(RELEASE)/NOTES   $(FTPDIR)
	cp -f build/hmmer-$(RELEASE).tar.Z   $(FTPDIR)
	cp -f build/hmmer-$(RELEASE)-*.tar.Z $(FTPDIR)
	cp -f Userguide/Userguide.ps         $(FTPDIR)
	cp -f Userguide/Userguide.pdf        $(FTPDIR)
	cp -f Userguide/hmmer-html.tar.Z     $(FTPDIR)


# make stable: Set up the FTP site symlink to the current stable HMMER release.
#
stable:
	(cd ${FTPDIR}; ln -sf hmmer-$(RELEASE).tar.Z hmmer.tar.Z)


# make optcheck:
#    Check that all program options are documented in man pages,
#    and tested in sqc scripts, by running ssdk's checkoptions.pl 
#    script on each program. 
#    
optcheck:
	@echo Checking options for consistency and documentation...
	@for prog in $(PROGS); do\
	   ssdk/checkoptions.pl src/$$prog documentation/man/$$prog.man testsuite/exercises.sqc testsuite/exercises-threaded.sqc;\
	done

