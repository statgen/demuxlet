#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.63])
AC_INIT([demuxlet], [0.1.0], [hmkang@umich.edu])
AC_CONFIG_AUX_DIR(aux-dist)
AM_INIT_AUTOMAKE

AC_PROG_MAKE_SET
AC_PROG_CXX
AC_LANG(C++)

AC_CONFIG_HEADERS([config.h])
AC_CONFIG_MACRO_DIR([m4])

# Checks for programs.
AC_PROG_CXX

AC_CHECK_HEADERS([limits.h stddef.h stdint.h stdlib.h string.h unistd.h])

AC_HEADER_STDBOOL
AC_C_INLINE
AC_TYPE_INT32_T
AC_TYPE_INT64_T
AC_C_RESTRICT
AC_TYPE_SIZE_T
AC_TYPE_UINT16_T
AC_TYPE_UINT32_T
AC_TYPE_UINT64_T
AC_TYPE_UINT8_T
AC_CHECK_TYPES([ptrdiff_t])

# Checks for library functions.
AC_FUNC_ERROR_AT_LINE
AC_FUNC_MALLOC
AC_HEADER_MAJOR
AC_FUNC_REALLOC
AC_FUNC_STRTOD
AC_CHECK_FUNCS([floor gethostname memset pow select sqrt strchr strstr strtol])
AC_CHECK_LIB(z, gzopen, [], [AC_MSG_ERROR([libz.{so,a} was not found. Please install zlib at http://www.zlib.net/ first])])
AC_CHECK_LIB([m],[erf],[AC_DEFINE([HAVE_ERF],[1],[libm includes erf])])
AC_CHECK_LIB([m],[erfc],[AC_DEFINE([HAVE_ERFC],[1],[libm includes erfc])])

LT_PREREQ([2.2])
LT_INIT([dlopen shared])

AC_CONFIG_FILES([Makefile])

AC_OUTPUT
