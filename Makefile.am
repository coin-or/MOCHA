# Copyright (C) 2009 David Haws.
# All Rights Reserved.
# This file is distributed under the Eclipse Public License.

## $Id$

# Author: David Haws, University of Kentucky, 10/26/09

AUTOMAKE_OPTIONS = foreign

EXTRA_DIST = doxydoc/doxygen.conf

########################################################################
#                          Subdirectories                              #
########################################################################

# subdirs is set by configure as the list of all subdirectories to recurse
# into
SUBDIRS = $(subdirs)

########################################################################
#                           Extra Targets                              #
########################################################################

test: all
	cd MOCHA; $(MAKE) test

tests: all
	for dir in $(subdirs); do \
          if test -r $$dir/test/Makefile; then \
            (cd $$dir; $(MAKE) test) \
          fi; \
        done

doxydoc:
	cd $(srcdir); doxygen doxydoc/doxygen.conf

.PHONY: test unitTest tests unitTests doxydoc cuter

########################################################################
#                         Maintainer Stuff                             #
########################################################################

# Files that are generated and should be cleaned with make distclean
DISTCLEANFILES =

include BuildTools/Makemain.inc
