# Copyright 2009 by Jesus De Loera, David Haws, Jon Lee, Allison
# O'Hair, University of California at Davis, IBM, Inc.
# All Rights Reserved.
# This file is distributed under the Eclips Public License v 1.0

## $Id: configure.ac 48 2009-10-27 19:06:48Z andreasw $

# Authors:  David Haws      University of California at Davis   2009-10-27

#############################################################################
#                       Names and other basic things                        #
#############################################################################

AC_PREREQ(2.59)

AC_INIT([MOCHA],[0.9],[dchaws+MOCHA@gmail.com])

AC_COPYRIGHT(
[Copyright 2009 by Jesus De Loera, David Haws, Jon Lee, Allison O'Hair, University of California at Davis, IBM, Inc.
All Rights Reserved.
This file is part of the open source package MOCHA which is distributed
under the Eclipse Public License.])

# List one file in the package so that the configure script can test
# whether the package is actually there
AC_CONFIG_SRCDIR(src/graph.h)

# Where should everything be installed by default?  Here, we want it
# to be installed directly in 'bin', 'lib', 'include' subdirectories
# of the directory where configure is run.  The default would be
# /usr/local.
AC_PREFIX_DEFAULT([`pwd`])

#############################################################################
#                         Standard build tool stuff                         #
#############################################################################

# Get the system type
AC_CANONICAL_BUILD

# If this project depends on external projects, the Externals file in
# the source root directory contains definition of where to find those
# externals.  The following macro ensures that those externals are
# retrieved by svn if they are not there yet.
AC_COIN_PROJECTDIR_INIT

# Check if user wants to produce debugging code
AC_COIN_DEBUG_COMPILE(MOCHA)

# We need egrep and ln -s somewherw
AC_PROG_EGREP
AC_PROG_LN_S

# Get the name of the C++ compiler and appropriate compiler options
AC_COIN_PROG_CXX

# Get the name of the Fortran compiler and appropriate compiler options
AC_COIN_PROG_F77

# Find out how to call Fortran from C and determine Fortran runtime libraries
AC_COIN_F77_WRAPPERS

# Initialize automake and libtool
AC_COIN_INIT_AUTO_TOOLS


#############################################################################
#                              COIN components                              #
#############################################################################

AC_COIN_HAS_PROJECT(MOCHA)


##########################################################################

# Additional libraries and packages

### We need to restore the value of LIBS before those test, otherwise libtool
### will hick-up
##LIBS_without_ADDLIBS="$LIBS"

#######################
#        BLAS         #
#######################

AC_COIN_HAS_BLAS
if test $coin_has_blas = no; then
  AC_MSG_ERROR([cannot find BLAS])
fi

if test "$use_blas" = BUILD; then
  THIRDBUILDLIBS="../../ThirdParty/Blas/libcoinblas.la $THIRDBUILDLIBS"
fi

#######################
#       LAPACK        #
#######################

AC_COIN_HAS_LAPACK
if test $coin_has_lapack = no; then
  AC_MSG_WARN([cannot find LAPACK])
fi

if test "$use_lapack" = BUILD; then
  THIRDBUILDLIBS="../../ThirdParty/Lapack/libcoinlapack.la $THIRDBUILDLIBS"
fi

AC_SUBST(THIRDBUILDLIBS)

#######################
#         GMP         #
#######################

AC_ARG_WITH([gmp],
  [AS_HELP_STRING([--with-gmp],
    [GMP Support. Optionally set to the location of gmp, or no to disable.])],
  [if test "$withval" != ""; then
    CPPFLAGS="${CPPFLAGS} -I$withval/include"
    LDFLAGS="${LDFLAGS} -L$withval/lib"
  fi
  ],
  [with_gmp=check])

#LIBGMP=
AS_IF([test "x$with_gmp" != xno],
  [AC_CHECK_LIB([gmp], [__gmpz_init], ,
    [if test "x$with_gmp" != xcheck; then
       AC_MSG_FAILURE(
         [--with-gmp was given, but test for gmp failed])
     fi
    ])
   AC_CHECK_LIB([gmpxx], [main], ,
    [if test "x$with_gmp" != xcheck; then
       AC_MSG_FAILURE([--with-gmp was given, but test for gmpxx failed])
     fi
    ])
  ])


AC_CONFIG_HEADERS([config.h])
AC_COIN_VPATH_LINK(Instances/Examples/k6.mo)
AC_COIN_VPATH_LINK(test/unitTestInput.txt)
AC_CONFIG_FILES([
 Makefile
 src/Makefile
 test/Makefile
 test/run_unitTests
])

AC_COIN_FINALIZE
