## @file c/configure.ac
## @author Tim van Werkhoven (werkhoven@strw.leidenuniv.nl)
## Copyright 2012 Tim van Werkhoven (werkhoven@strw.leidenuniv.nl)

AC_PREREQ(2.59)

AC_INIT([libunwrap], [0.1], [Tim van Werkhoven <werkhoven@strw.leidenuniv.nl>, Visa Korkiakoski <korkiakoski@strw.leidenuniv.nl>])
AC_CONFIG_AUX_DIR([build-aux])

AM_INIT_AUTOMAKE([AC_PACKAGE_NAME], [AC_PACKAGE_VERSION])

AC_COPYRIGHT([Copyright (C) 2012 Tim van Werkhoven <werkhoven@strw.leidenuniv.nl>, Visa Korkiakoski <korkiakoski@strw.leidenuniv.nl>])
AC_REVISION([0.1])

AC_CONFIG_SRCDIR(unwrap/libunwrap.cc)
AC_CONFIG_HEADERS(unwrap/autoconfig.h)

AC_CONFIG_FILES([
Makefile
unwrap/Makefile
test/Makefile
])

AC_GNU_SOURCE

AC_PROG_CPP
AC_PROG_CC
AC_PROG_CXXCPP
AC_PROG_CXX
AC_PROG_MAKE_SET
AC_PROG_INSTALL
AC_PROG_RANLIB

AC_LANG(C++)

PKG_PROG_PKG_CONFIG

# FOAM data & configuration directory

### These libraries are crucial for all targets
AC_SEARCH_LIBS([sin],
		[m],
		[],
		[AC_MSG_ERROR([Cannot build libunwrap! Missing libm!])])

### Check command-line options

AC_ARG_ENABLE([debug],
		AC_HELP_STRING([--enable-debug], [enable debug output]),
		[have_debug=$enableval],
		[have_debug=no])
AC_ARG_WITH([tests],
		AC_HELP_STRING([--without-tests], [do not build test programs]),
		[with_tests=$withval],
		[with_tests=yes])
AM_CONDITIONAL([WITH_TESTING], [test "x$with_tests" != "xno"])

AM_CONDITIONAL([HAVE_DEBUG], [test "x$have_debug" != "xno"])
AS_IF([test "x$have_debug" != "xno"],
		[AC_DEFINE([HAVE_DEBUG], [1], [Toggle debugging.])],
		[AC_DEFINE([HAVE_DEBUG], [0], [Toggle debugging.])])

### Finish

AC_OUTPUT
