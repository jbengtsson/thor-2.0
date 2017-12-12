#include <stdio.h>
#include <math.h>

#include <num_rec.h>

#define sqr(x) ((x)*(x))

static double xsav;
static double (*nrfunc)(double, double);

double quad2d(double (*func)(double, double), double x1, double x2)
{
  double dqgaus(double (*func)(double), double a, double b);
  double f1(double x);

  nrfunc = func;
  return dqgaus(f1, x1, x2);
}

double f1(double x)
{
  double dqgaus(double (*func)(double), double a, double b);
  double f2(double y);
  double yy1(double), yy2(double);

  xsav = x;
  return dqgaus(f2, yy1(x), yy2(x));
}

double f2(double y)
{
  return (*nrfunc)(xsav, y);
}

double yy1(double x)
{
  return -sqrt(sqr(x)-1e0);
}

double yy2(double x)
{
  return sqrt(sqr(x)-1e0);
}

double f(double x, double y)
{
  return exp(-x);
}

int main(int argc, char *argv[])
{
  double intgrl;
  intgrl = quad2d(f, 1e0, 20e0);
  printf("integral = %11.3e\n", intgrl);
}
