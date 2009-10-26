## $Rev: 42 $ $Date: 2009-10-26 15:13:29 -0400 (Mon, 26 Oct 2009) $
AC_INIT([MOCHA], [0.9], [dchaws+MOCHA@gmail.com])
AM_INIT_AUTOMAKE([-Wall -Werror foreign])
AC_LANG(C++)
AC_PROG_CXX
AC_CONFIG_HEADERS([config.h])
AC_CONFIG_FILES([
 Makefile
 src/Makefile
])

#m4_pattern_allow([AC_COIN_THIRDPARTY_SUBDIRS])

AC_COIN_HAS_BLAS
AC_COIN_HAS_LAPACK
#AC_COIN_THIRDPARTY_SUBDIRS([Blas],[ThirdParty/Blas],[daxpy.f])
#AC_COIN_THIRDPARTY_SUBDIRS([Lapack],[ThirdParty/Lapack],[LAPACK/SRC/dlarf.f])

AC_ARG_WITH(lapack-prefix,
[  --with-lapack-prefix=DIR   path to LAPACK],[
  if test "$withval" != "no"; then
    CPPFLAGS="${CPPFLAGS} -I$withval/include"
    LDFLAGS="${LDFLAGS} -L$withval/lib"
  fi
])

AC_ARG_WITH(blas-prefix,
[  --with-blas-prefix=DIR   path to BLAS],[
  if test "$withval" != "no"; then
    CPPFLAGS="${CPPFLAGS} -I$withval/include"
    LDFLAGS="${LDFLAGS} -L$withval/lib"
  fi
])

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


#AC_ARG_WITH(gmp,
#[  --with-gmp=DIR   path to GMP],[
#  if test "$withval" != ""; then
#    CPPFLAGS="${CPPFLAGS} -I$withval/include"
#    LDFLAGS="${LDFLAGS} -L$withval/lib"
#  fi
#  with_gmp=yes
#],[with_gmp=no])

#AC_CHECK_LIB(gmp, __gmpz_init, ,
#       [AC_MSG_ERROR([GNU MP not found, see http://gmplib.org/])])
#AC_CHECK_LIB(gmpxx, __gmpz_init, ,
#       [AC_MSG_ERROR([GNU MP with cxx compiled not found, see http://gmplib.org/])])
#AC_CHECK_HEADER([gmpxx.h], ,
#       [AC_MSG_ERROR([Header gmpxx.h not found])],[AC_INCLUDES_DEFAULT])

AC_CHECK_LIB(lapack, dgesvd_, ,
       [AC_MSG_ERROR([Lapack not found, see ])])
AC_CHECK_LIB(blas, main, ,
       [AC_MSG_ERROR([Blas not found, see ])])

AC_OUTPUT
