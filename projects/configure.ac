# Process this file with autoconf to produce a configure script.
#AC_INIT
AC_CONFIG_SRCDIR([thor/src/field.cc])
AC_INIT([thor],[1.0])
AC_CONFIG_AUX_DIR(config)
AC_CONFIG_SRCDIR(src/thor.cc)
AM_INIT_AUTOMAKE()
AM_CONFIG_HEADER(config.h)
# Checks for programs.
AC_PROG_AWK
# Use GNU vs. Clang/LLVM.
AC_PROG_CC(gcc-15)
AC_PROG_INSTALL
AC_PROG_LN_S
AC_PROG_CPP
AC_PROG_RANLIB
# Use GNU gcc; vs. LLVM compilers.
AC_PROG_CXX(g++-15)
AC_PROG_F77(gfortran-15)

# Checks for libraries.

# Checks for header files.
AC_CHECK_HEADERS([stdlib.h string.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_C_CONST
AC_C_INLINE
AC_TYPE_SIZE_T

# Checks for library functions.
AC_CONFIG_FILES([Makefile src/Makefile src/Powell/Makefile src/Powell/src/Makefile])
AC_OUTPUT
