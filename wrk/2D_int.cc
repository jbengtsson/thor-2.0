static double  xsav;
static double  (*nrfunc)(const double, const double);


double f2(double y) { return (*nrfunc)(xsav, y); }


double f1(double x)
{

  xsav = x;

  return dqromb(f2, 0e0, sqrt(2e0*Jy));
}


double dqromb2D(double (*func)(const double, const double),
		const double x1, const double x2)
{

  nrfunc = func;

  return dqromb(f1, x1, x2);
}


static double f_Int(const double Jx, const double Jy)
{
  ss_vect<double>  ps;

  ps.zero();
  ps[x_] = sqrt(2e0*Jx); ps[px_] = sqrt(2e0*Jx);
  ps[y_] = sqrt(2e0*Jy); ps[py_] = sqrt(2e0*Jy);

  return fabs((PPUSH(map1, ps))[0]);
} 


double f_prm2(const int jj[])
{
  double  f;

  f = (jj[ss_dim-1] == 1)? 1e0 : 0e0;

 return f;
}


static double f_Int_prm(const double Jx, const double Jy)
{
  double           f, f_prm;
  ss_vect<double>  ps;

  ps.zero();
  ps[x_] = sqrt(2e0*Jx); ps[px_] = sqrt(2e0*Jx);
  ps[y_] = sqrt(2e0*Jy); ps[py_] = sqrt(2e0*Jy);

  f = (PPUSH(map1, ps))[0]; f_prm = (PPUSH(map2, ps))[0];

  return (f >= 0e0)? f_prm : -f_prm;
} 


void ini_nu2(const tps &K)
{
  tps           K1, dnudJ[2][3], jacob;
  ss_vect<tps>  nus, dnus[2];

  K1 = dacfu1(K, f_prm1);

  nus = dHdJ(K1); dnus[X_] = dHdJ(nus[0]); dnus[Y_] = dHdJ(nus[1]);

  dnudJ[X_][X_] = dnus[X_][3];
  dnudJ[X_][Y_] = dnus[X_][4];
  dnudJ[X_][2]  = Der(nus[0], 5);

  dnudJ[Y_][X_] = dnus[Y_][3];
  dnudJ[Y_][Y_] = dnus[Y_][4];
  dnudJ[Y_][2]  = Der(nus[1], 5);


  jacob = dnudJ[X_][X_]*dnudJ[Y_][Y_] - dnudJ[X_][Y_]*dnudJ[Y_][X_];

  map1.identity(); map1[0] = jacob; map1 = MTREE(map1);

  map2.identity(); map2[0] = Der(dacfu1(jacob, f_prm2), 7); map2 = MTREE(map2);
}


double get_nu2(void)
{
  double  nu2;

  cout << "calling dqromb2D" << endl;
  nu2 = dqromb2D(f_Int, 0e0, sqrt(2e0*Jy));
  cout << "dqromb2D completed" << endl;

  return nu2;
}


double get_nu2_prm(void)
{
  double  nu2_prm;

  cout << "calling dqromb2D" << endl;
  nu2_prm = dqromb2D(f_Int_prm, 0e0, sqrt(2e0*Jy));
  cout << "dqromb2D completed" << endl;

  return nu2_prm;
}
