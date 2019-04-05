// The Chebyquad test problem (Fletcher, 1965) for N = 2,4,6 and 8,
// with NPT = 2N+1.

#include <stdio.h>

#include "newuoa.h"


int main(int argc, char *argv[])
{
  long int n, npt, maxfun, iprint, n1;
  int      i;
  double   w[10000], x[10], rhobeg, rhoend;

  iprint = 2;
  maxfun = 5000;
  rhoend = 1e-6;
  for (n = 2; n <= 8; n += 2) {
    npt = 2*n + 1;
    for (i = 1; i <= n; ++i) {
      x[i-1] = (double)i/(double)(n+1);
    }
    rhobeg = 0.2e0*x[0];
    printf("Results with N = %1d and NPT = %2d\n", n, npt);
    newuoa_(n, npt, x, rhobeg, rhoend, iprint, maxfun, w);
  }
}

