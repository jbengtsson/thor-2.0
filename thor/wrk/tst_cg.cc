#define NO 1

#include "thor_lib.h"
//#include "num_rec.h"

int  no_tps   = NO,
     ndpt_tps = 5;


int  n_iter = 0;


double f(double x[])
{
  const int n_prt = 1;

  n_iter++;

  if (n_iter % n_prt == 0)
    cout << scientific << setprecision(16)
	 << setw(4) << n_iter
	 << setw(24) << x[1] << setw(24) << x[2] << endl;

  return sin(x[1]+x[2]);
}


void df(double x[], double df[])
{
  df[1] = cos(x[1]+x[2]);
  df[2] = cos(x[1]+x[2]);
}


void find_min()
{
  int    iter;
  double *x, fret;

  const int    n_prm = 2;
  const double eps = 1e-20;

  x = dvector(1, n_prm);

  x[1] = 0e0; x[2] = 0e0;

  dfrprmn(x, n_prm, eps, &iter, &fret, f, df);

  cout << scientific << setprecision(16)
       << setw(3) << iter << setw(24) << x[1] << setw(24) << x[2] << endl;

  free_dvector(x, 1, n_prm);
}


int main(int argc, char *argv[])
{

  find_min();
}
