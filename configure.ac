# Process this file with autoconf to produce a configure script.
#AC_INIT(thor/src/field.cc)
AC_INIT(thor, 1.0)
AC_CONFIG_AUX_DIR(config)
AC_CONFIG_SRCDIR(thor/src/field.cc)
#AM_INIT_AUTOMAKE(thor, 1.0)
AM_CONFIG_HEADER(config.h)
#AM_INIT_AUTOMAKE
AM_INIT_AUTOMAKE([subdir-objects])
# Checks for programs.
AC_PROG_AWK
AC_PROG_CC
AC_PROG_INSTALL
AC_PROG_LN_S
AC_PROG_CPP
AC_PROG_RANLIB
# Use GNU gcc; vs. LLVM.
AC_PROG_CXX(gcc-11)
AC_PROG_F77(gfortran)

# Checks for libraries.

# Checks for header files.
AC_CHECK_HEADERS([stdlib.h string.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_C_CONST
AC_C_INLINE
AC_TYPE_SIZE_T

# Checks for library functions.
AC_OUTPUT(Makefile
          TPSA/Makefile
          thor/Makefile thor/src/Makefile)
