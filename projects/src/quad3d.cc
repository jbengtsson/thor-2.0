#include <stdio.h>
#include <math.h>

#include <num_rec.h>

#define sqr(x) ((x)*(x))

static double xsav, ysav;
static double (*func_save)(double,double,double);


const double R = sqrt(1e0);


double dqgaus(double (*func)(double), double a, double b)
{
  int           j;
  double        xr, xm, dx, s;
  static double x[] =
    {0.0, 0.1488743389, 0.4333953941, 0.6794095682, 0.8650633666, 0.9739065285};
  static double w[] =
    {0.0, 0.2955242247, 0.2692667193, 0.2190863625, 0.1494513491, 0.0666713443};

  xm = 0.5*(b+a);
  xr = 0.5*(b-a);
  s = 0;
  for (j = 1; j <= 5; j++) {
    dx = xr*x[j];
    s += w[j]*((*func)(xm+dx)+(*func)(xm-dx));
  }
  return s *= xr;
}

double y0(double x) { return -sqrt(sqr(R)-sqr(x)); }

double y1(double x) { return sqrt(sqr(R)-sqr(x)); }

double z0(double x, double y) { return -sqrt(sqr(R)-sqr(x)-sqr(y)); }

double z1(double x, double y) { return sqrt(sqr(R)-sqr(x)-sqr(y)); }

double fz(double z) { return (*func_save)(xsav, ysav, z); }

double fy(double y)
{
  ysav = y;
  return dqgaus(fz, z0(xsav, y), z1(xsav, y));
}

double fx(double x)
{
  xsav = x;
  return dqgaus(fy, y0(x),y1(x));
}

double quad3d(double (*func)(double, double, double), double x1, double x2)
{
  func_save = func;
  return dqgaus(fx, x1, x2);
}

double f(double x, double y, double z)
{
  return 1e0;
}

int main(int argc, char *argv[])
{
  double intgrl;
  
  intgrl = quad3d(f, -R, R);
  printf("integral = %13.5e\n", intgrl*3e0/4e0);
}
