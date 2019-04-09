// The Chebyquad test problem (Fletcher, 1965) for N = 2,4,6 and 8,
// with NPT = 2N+1.

#include <stdio.h>
#include <cstdlib>

#include "num_rec.h"

#include "newuoa.h"

extern "C" {
  void calfun_(int *n, double *x, double *f);
}



void calfun_(int *n, double *x, double *f)
{
  int    i, j, np, iw;
  double y[10][10], sum;

  for (j = 0; j < *n; ++j) {
    y[0][j] = 1e0;
    y[1][j] = 2e0*x[j] - 1e0;
  }
  for (i = 1; i < *n; ++i)
    for (j = 0; j < *n; ++j)
      y[i+1][j] = 2e0*y[1][j]*y[i][j] - y[i-1][j];
  *f = 0e0;
  np = *n + 1;
  iw = 1;
  for (i = 1; i <= np; ++i) {
    sum = 0e0;
    for (j = 0; j < *n; ++j)
      sum += y[i-1][j];
    sum /= (double)(*n);
    if (iw > 0) sum += 1e0/(double)(i*i-2*i);
    iw = -iw;
    *f += sum*sum;
  }
}


int main(int argc, char *argv[])
{
  long int n, npt, maxfun, iprint, n1;
  int    i;
  double w[10000], x[10], rhobeg, rhoend;

  iprint = 2;
  maxfun = 5000;
  rhoend = 1e-6;
  for (n = 2; n <= 8; n += 2) {
    npt = 2*n + 1;
    for (i = 1; i <= n; ++i)
      x[i-1] = (double)i/(double)(n+1);
    rhobeg = 0.2e0*x[0];
    printf("Results with N = %1d and NPT = %2d\n", n, npt);
    newuoa_(n, npt, x, rhobeg, rhoend, iprint, maxfun, w);
  }
}

