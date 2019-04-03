/* main.f -- translated by f2c (version 20160102).
   You must link the resulting object file with libf2c:
	on Microsoft Windows system, link with libf2c.lib;
	on Linux or Unix systems, link with .../path/to/libf2c.a -lm
	or, if you install libf2c.a in a standard place, with -lf2c -lm
	-- in that order, at the end of the command line, as in
		cc *.o -lf2c -lm
	Source for libf2c is in /netlib/f2c/libf2c.zip, e.g.,

		http://www.netlib.org/f2c/libf2c.zip
*/

#include "f2c.h"

/* Table of constant values */

static integer c__1 = 1;


/*     The Chebyquad test problem (Fletcher, 1965) for N = 2,4,6 and 8, */
/*     with NPT = 2N+1. */

/* Main program */ int MAIN__(void)
{
    /* Format strings */
    static char fmt_20[] = "(//4x,\002Results with N =\002,i2,\002 and NPT "
	    "=\002,i3)";

    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    integer s_wsfe(cilist *), do_fio(integer *, char *, ftnlen), e_wsfe(void);
    /* Subroutine */ int s_stop(char *, ftnlen);

    /* Local variables */
    static integer i__, n;
    static doublereal w[10000], x[10];
    static integer npt;
    static doublereal rhobeg, rhoend;
    static integer maxfun;
    extern /* Subroutine */ int newuoa_(integer *, integer *, doublereal *, 
	    doublereal *, doublereal *, integer *, integer *, doublereal *);
    static integer iprint;

    /* Fortran I/O blocks */
    static cilist io___9 = { 0, 6, 0, fmt_20, 0 };


    iprint = 2;
    maxfun = 5000;
    rhoend = 1e-6;
    for (n = 2; n <= 8; n += 2) {
	npt = (n << 1) + 1;
	i__1 = n;
	for (i__ = 1; i__ <= i__1; ++i__) {
/* L10: */
	    x[i__ - 1] = (doublereal) i__ / (doublereal) (n + 1);
	}
	rhobeg = x[0] * .2;
	s_wsfe(&io___9);
	do_fio(&c__1, (char *)&n, (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&npt, (ftnlen)sizeof(integer));
	e_wsfe();
	newuoa_(&n, &npt, x, &rhobeg, &rhoend, &iprint, &maxfun, w);
/* L30: */
    }
    s_stop("", (ftnlen)0);
    return 0;
} /* MAIN__ */

