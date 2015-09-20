#define NO 6

#include <limits>
#include "thor_lib.h"

int  no_tps   = NO,
     ndpt_tps = 5;


extern double        b2_max;
extern tps           K, g;
extern ss_vect<tps>  Map, A0, A1, Map_res;

int    dnn;      // Number of constraints.
double *dfvec;
void   (*dnrfuncv)(int n, double v[], double f[]);

const int max_ind = NO;

bool          fit_chrm;
int           n_prm, prm[n_prm_max], prm_n[n_prm_max];
int           adj_tune, adj_chrom, n_iter, sext_scheme;
long int      beta_loc1, beta_loc2, beta_loc3, beta_loc4, rseed;
double        Jx, Jy, delta, beta1[2], beta2[2], beta3[2], beta4[2];
double        nu0[2], eps_nu, ksi1[2], eps_ksi;
double        nu_x_min, nu_x_max, nu_y_min, nu_y_max;
double        bnL_max[mpole_max+1], scl_res, scl_ksi1, scl_ksi2, scl_dnu2[2];
double        scl_nl, scl_fold1, scl_fold2, scl_fold3, scl_flip;
ss_vect<tps>  Id_scl, map1, map2;
ofstream      sext_out;


// const double A_max[] = {20e-3, 7.5e-3};
const double A_max[] = {20e-3, 10e-3};
const double delta_max = 3.0e-2;

const int no_g = 5, no_dnu = 7;


tps dacfu1(const tps &a, double (*func)(const int []))
{
  char    name[11];
  int     j, n, jj[ss_dim], ibuf1[bufsize], ibuf2[bufsize];
  double  rbuf[bufsize];
  tps     b;

  a.exprt(rbuf, ibuf1, ibuf2, name); n = (int)rbuf[0];

  for (j = 0; j < n; j++) {
    dehash_(no_tps, ss_dim, ibuf1[j], ibuf2[j], jj);

    rbuf[j+1] *= (*func)(jj);
  }

  b.imprt(n, rbuf, ibuf1, ibuf2);

  // Remove zeroes.
  return 1e0*b;
}


tps Int(const tps &a, const int k)
{
  char    name[11];
  int     j, n, jj[ss_dim], ibuf1[bufsize], ibuf2[bufsize];
  double  rbuf[bufsize];
  tps     b;

  a.exprt(rbuf, ibuf1, ibuf2, name); n = (int)rbuf[0];

  for (j = 0; j < n; j++) {
    dehash_(no_tps, ss_dim, ibuf1[j], ibuf2[j], jj);

    rbuf[j+1] /= jj[k-1] + 1e0;
  }

  b.imprt(n, rbuf, ibuf1, ibuf2);

  return b*tps(0e0, k);
}


double f_fltr(const int jj[])
{

  return ((jj[x_] == jj[px_]) && (jj[y_] == jj[py_]))? 1e0 : 0e0;
}


double f_trunc_par(const int jj[])
{

 return (jj[ss_dim-1] > 1)? 0e0 : 1e0;
}


double f_g(const int jj[])
{
  int     k, n;
  double  f;

  n = 0;
  for (k = 0; k < ss_dim; k++)
    n += jj[k];

  f = ((n == no_g) && (jj[ss_dim-1] == 0))? 0e0 : 1e0;

 return f;
}


double f_dnu_(const int jj[])
{
  int     k, n;
  double  f;

  n = 0;
  for (k = 0; k < ss_dim; k++)
    n += jj[k];

  f = ((n > no_dnu-2) || (jj[ss_dim-1] > 1) ||
       ((n == no_dnu-2) && (jj[ss_dim-1] == 0)) ||
       (jj[delta_] > 3) ||
       ((jj[ss_dim-1] == 0) && (jj[delta_] != 0) && (n > 4)) ||
       ((jj[ss_dim-1] == 1) && (jj[delta_] != 0) && (n > 5)))?
    0e0 : 1e0;

 return f;
}


double f_dnu(const int jj[])
{
  int     k, n;
  double  f;

  n = 0;
  for (k = 0; k < ss_dim; k++)
    n += jj[k];

  // [h^+, h^-] basis.
  f = ((jj[ss_dim-1] > 1) || (n > no_dnu/2+2) ||
       ((n == no_dnu/2+3) && (jj[ss_dim-1] == 0)) ||
       (jj[delta_] > 4) ||
       ((jj[ss_dim-1] == 0) && (jj[delta_] != 0) && (n > 4)) ||
       ((jj[ss_dim-1] == 1) && (jj[delta_] != 0) && (n > 5)))?
    0e0 : 1e0;

 return f;
}


tps dfdJ(const tps &f, const int i)
{
  // df/dJ in resonance basis:
  // 
  // df/dJ = df/dh^+ * dh^+/dJ + df/dh^- * dh^-/dJ
  //       = f/dh^+ * 1/h^- + df/dh^- * 1/dh^+

  return pseudo_der(Der(f, 2*i-1), 2*i) + pseudo_der(Der(f, 2*i), 2*i-1);
}


void get_dnu2_(const tps &K_re, tps dnu2[])
{
  // Dnu^2 = Int{|dnu_x/dJ x dnu_y/dJ| dJ_x dJ_y}.
  //
  // [J, J^2, J^3]: [K_22000, K_33000, K_44000, ...]
  //                [nu_11000, nu_22000, nu_33000, ...].
  //
  // [d, d^2, d^3]: [K_11001, K_11002, K_11003, ...]
  //                [nu_00001, nu_00002, nu_00003, ...].

  int           i, j;
  tps           nu, dnudJ[2][3], dS[3];
  ss_vect<tps>  Id_tr, ps1, ps0;

  const bool   prt = false;

  for (i = 0; i < 2; i++) {
    nu = -dfdJ(K_re, i+1)/(2e0*M_PI);
    nu = dacfu1(nu, f_dnu);

    if (prt) cout << nu << endl;

    for (j = 0; j < 2; j++)
      dnudJ[i][j] = dfdJ(nu, j+1);

    dnudJ[i][Z_] = Der(nu, delta_+1);

    for (j = 0; j < 3; j++)
      dnudJ[i][j] = dacfu1(dnudJ[i][j], f_fltr);

    if (prt) cout << dnudJ[i][X_] << dnudJ[i][Y_] << dnudJ[i][Z_] << endl;
  }

  // [h^+, h^-] -> [J, phi], h^+ * h^- = 2J.
  Id_tr.identity(); Id_tr[px_] = 2e0; Id_tr[py_] = 2e0;

  // Scale to avoid numerical issues.
  ps1.identity(); ps1[x_] *= 2e0*Jx; ps1[y_] *= 2e0*Jy; ps1[delta_] *= delta;

  for (i = 0; i < 2; i++)
    for (j = 0; j < 3; j++) {
      dnudJ[i][j] = dnudJ[i][j]*Id_tr;
      dnudJ[i][j] = dnudJ[i][j]*ps1;
    }

  dS[X_] = dnudJ[X_][Y_]*dnudJ[Y_][Z_] - dnudJ[X_][Z_]*dnudJ[Y_][Y_];
  dS[Y_] = -dnudJ[X_][X_]*dnudJ[Y_][Z_] + dnudJ[X_][Z_]*dnudJ[Y_][X_];
  dS[Z_] = dnudJ[X_][X_]*dnudJ[Y_][Y_] - dnudJ[X_][Y_]*dnudJ[Y_][X_];

  if (prt)
    cout << dacfu1(dS[X_], f_trunc_par) << dacfu1(dS[Y_], f_trunc_par)
	 << dacfu1(dS[Z_], f_trunc_par) << endl;

  for (i = 0; i < 3; i++) {
    dnu2[i] = sqr(dS[i]);
    if (prt) cout << dacfu1(dnu2[i], f_trunc_par) << endl;
  }

  // Undo scaling.
  for (i = 0; i < 3; i++) {
    ps1.identity(); ps1[x_] /= 2e0*Jx; ps1[y_] /= 2e0*Jy; ps1[delta_] /= delta;
    dnu2[i] = dnu2[i]*ps1;
    if (prt) cout << dacfu1(dnu2[i], f_trunc_par) << endl;

    dnu2[i] = Int(dnu2[i], x_+1);
    if (prt) cout << dacfu1(dnu2[i], f_trunc_par) << endl;
    ps0.identity(); ps0[x_] = 0e0; ps1.identity(); ps1[x_] = 2e0*Jx;
    dnu2[i] = dnu2[i]*ps1 - dnu2[i]*ps0;
    if (prt) cout << dacfu1(dnu2[i], f_trunc_par) << endl;

    dnu2[i] = Int(dnu2[i], y_+1);
    if (prt) cout << dacfu1(dnu2[i], f_trunc_par) << endl;
    ps0.identity(); ps0[y_] = 0e0; ps1.identity(); ps1[y_] = 2e0*Jy;
    dnu2[i] = dnu2[i]*ps1 - dnu2[i]*ps0;
    if (prt) cout << dacfu1(dnu2[i], f_trunc_par) << endl;

    dnu2[i] = Int(dnu2[i], delta_+1);
    if (prt) cout << dacfu1(dnu2[i], f_trunc_par) << endl;
    ps0.identity(); ps0[delta_] = 0e0; ps1.identity(); ps1[delta_] = delta;
    dnu2[i] = dnu2[i]*ps1 - dnu2[i]*ps0;
  }

  if (prt) {
    for (i = 0; i < 3; i++)
      cout << dacfu1(dnu2[i], f_trunc_par) << endl;
    exit(0);
  }
}


void get_dnu2(const tps &K_re, tps dnu2[])
{
  // [J, J^2, J^3]: [K_22000, K_33000, K_44000, ...]
  //                [nu_11000, nu_22000, nu_33000, ...].
  //
  // [d, d^2, d^3]: [K_11001, K_11002, K_11003, ...]
  //                [nu_00001, nu_00002, nu_00003, ...].

  int           i;
  tps           nu[2], dnu[2];
  ss_vect<tps>  Id, Id_tr, ps1, ps0;

  const bool prt = false;

  for (i = 0; i < 2; i++) {
    nu[i] = -dfdJ(K_re, i+1)/(2e0*M_PI);
    if (prt) cout << nu[i] << endl;
    dnu[i] = nu[i] - nu[i].cst();
    dnu[i] = dacfu1(dnu[i], f_dnu);
    if (prt) cout << dnu[i] << endl;

    // [h^+, h^-] -> [J, phi], h^+ * h^- = 2J.
    Id_tr.identity(); Id_tr[px_] = 2e0; Id_tr[py_] = 2e0;
    dnu[i] = dnu[i]*Id_tr;

    dnu2[i] = sqr(dnu[i]);

    if (prt) cout << dacfu1(dnu2[i], f_trunc_par) << endl;

    if (false) {
      // Remove non-chromatic terms.
      Id.identity(); Id[delta_] = 0e0;
      dnu2[i] -= dnu2[i]*Id;
      if (prt) cout << dacfu1(dnu2[i], f_trunc_par) << endl;
    }

    dnu2[i] = Int(dnu2[i], x_+1);
    ps0.identity(); ps0[x_] = 0e0; ps1.identity(); ps1[x_] = Jx;
    dnu2[i] = dnu2[i]*ps1 - dnu2[i]*ps0;

    dnu2[i] = Int(dnu2[i], y_+1);
    ps0.identity(); ps0[y_] = 0e0; ps1.identity(); ps1[y_] = Jy;
    dnu2[i] = dnu2[i]*ps1 - dnu2[i]*ps0;

    dnu2[i] = Int(dnu2[i], delta_+1);
    if (prt) cout << dacfu1(dnu2[i], f_trunc_par) << endl;
    ps0.identity(); ps0[delta_] = -delta; ps1.identity(); ps1[delta_] = delta;
    dnu2[i] = dnu2[i]*ps1 - dnu2[i]*ps0;

    dnu2[i] /= Jx*Jy*2e0*delta;

    if (prt) cout << dacfu1(dnu2[i], f_trunc_par) << endl;
  }

  if (prt) exit(0);
}


tps get_abs2(const tps &f)
{
  // Compute abs2 and dabs2.
  char         name[11];
  int          j, n, jj[ss_dim], ibuf1[bufsize], ibuf2[bufsize];
  double       rbuf[bufsize];
  tps          abs2;

  f.exprt(rbuf, ibuf1, ibuf2, name); n = (int)rbuf[0];

  abs2 = 0e0;
  for (j = 0; j < n; j++) {
    dehash_(no_tps, ss_dim, ibuf1[j], ibuf2[j], jj);

    if (jj[ss_dim-1] == 0) {
      jj[ss_dim-1] = 1;
      abs2 += sqr(rbuf[j+1]) + 2e0*rbuf[j+1]*f[jj]*tps(0e0, ss_dim);
    }
  }

  return abs2;
}


void get_dyn_(tps &K_re, tps &g2, tps dnu2[])
{
  tps K_im;

  danot_(no_dnu-1);
  get_Map();
  danot_(no_dnu);

  K = MapNorm(Map, g, A1, A0, Map_res, 1);

  CtoR(K, K_re, K_im);
  // Filter numeric noise.
  K_re = dacfu1(K_re, f_fltr);

  danot_(no_g);
  g = dacfu1(g, f_g);
  danot_(no_tps);

  g2 = get_abs2(g*Id_scl); get_dnu2_(K_re, dnu2);
}


void get_dyn(tps &K_re, tps &g2, tps dnu2[])
{
  tps K_im;

  danot_(no_dnu-1);
  get_Map();
  danot_(no_dnu);

  K = MapNorm(Map, g, A1, A0, Map_res, 1); CtoR(K, K_re, K_im);

  danot_(no_g);
  g = dacfu1(g, f_g);
  danot_(no_tps);

  g2 = get_abs2(g*Id_scl); get_dnu2(K_re, dnu2);
}


template<typename T>
T f_step(const T &x)
{
  const double alpha = 1e0;

  return 1e0/(1e0+exp(-2e0*alpha*x));
}


void f_nl(int n, double bn[], double f[])
{
  // n is number of constraints, see dfmin_ Numerical Recipes.
  int    i, m;
  tps    K_re, K_re_scl, g2, dnu2[2], abs2;

  cout << endl;
  cout << "f_nl" << endl;
  for (i = 1; i <= n_prm; i++)
    set_bn(prm[i-1], prm_n[i-1], bn[i]);

  get_dyn(K_re, g2, dnu2);

  K_re_scl = K_re*Id_scl;

  m = 0;

  f[++m] = scl_ksi1*(ksi1[X_]*M_PI+h_ijklm(K_re, 1, 1, 0, 0, 1));
  f[++m] = scl_ksi1*(ksi1[Y_]*M_PI+h_ijklm(K_re, 0, 0, 1, 1, 1));

  f[++m] = scl_nl*h_ijklm(K_re_scl, 2, 2, 0, 0, 0);
  f[++m] = scl_nl*h_ijklm(K_re_scl, 0, 0, 2, 2, 0);
  f[++m] = scl_nl*h_ijklm(K_re_scl, 1, 1, 1, 1, 0);

  if (scl_flip != 0e0) {
    // f[++m] = scl_flip*sqr(h_ijklm(K_re_scl, 1, 1, 1, 1, 0)
    // 			-h_ijklm(K_re_scl, 2, 2, 0, 0, 0));
    f[++m] = scl_flip*f_step(-h_ijklm(K_re_scl, 1, 1, 1, 1, 0)
			     /h_ijklm(K_re_scl, 2, 2, 0, 0, 0));
    // f[++m] = scl_flip*sqr(h_ijklm(K_re_scl, 1, 1, 1, 1, 0)
    // 			-h_ijklm(K_re_scl, 0, 0, 2, 2, 0));
    f[++m] = scl_flip*f_step(-h_ijklm(K_re_scl, 1, 1, 1, 1, 0)
			     /h_ijklm(K_re_scl, 0, 0, 2, 2, 0));
  }

  if (scl_fold1 == 0e0) {
    f[++m] =
      scl_nl*h_ijklm(K_re_scl, 3, 3, 0, 0, 0);
    f[++m] =
      scl_nl*h_ijklm(K_re_scl, 0, 0, 3, 3, 0);
  } else {
    f[++m] =
      scl_fold1*(h_ijklm(K_re_scl, 2, 2, 0, 0, 0)/3e0
		 +h_ijklm(K_re_scl, 3, 3, 0, 0, 0));
    f[++m] =
      scl_fold1*(h_ijklm(K_re_scl, 0, 0, 2, 2, 0)/3e0
		 +h_ijklm(K_re_scl, 0, 0, 3, 3, 0));
  }

  if (scl_fold2 == 0e0) {
    f[++m] =
      scl_nl*h_ijklm(K_re_scl, 1, 1, 2, 2, 0);
    f[++m] =
      scl_nl*h_ijklm(K_re_scl, 2, 2, 1, 1, 0);
  } else {
    f[++m] =
      scl_fold2*(h_ijklm(K_re_scl, 1, 1, 1, 1, 0)/2e0
		 +h_ijklm(K_re_scl, 1, 1, 2, 2, 0));
    f[++m] =
      scl_fold2*(h_ijklm(K_re_scl, 1, 1, 1, 1, 0)/2e0
		 +h_ijklm(K_re_scl, 2, 2, 1, 1, 0));
  }

  f[++m] = scl_ksi2*h_ijklm(K_re_scl, 1, 1, 0, 0, 2);
  f[++m] = scl_ksi2*h_ijklm(K_re_scl, 0, 0, 1, 1, 2);

  if (scl_fold3 == 0e0) {
    f[++m] = scl_nl*h_ijklm(K_re_scl, 1, 1, 0, 0, 3);
    f[++m] = scl_nl*h_ijklm(K_re_scl, 0, 0, 1, 1, 3);
  } else {
    f[++m] =
      scl_fold3*(h_ijklm(K_re_scl, 1, 1, 0, 0, 1)/3e0
		 +h_ijklm(K_re_scl, 1, 1, 0, 0, 3));
    f[++m] =
      scl_fold3*(h_ijklm(K_re_scl, 0, 0, 1, 1, 1)/3e0
		 +h_ijklm(K_re_scl, 0, 0, 1, 1, 3));
  }

  // f[++m] = scl_dnu2[X_]*dnu2[X_].cst();
  // f[++m] = scl_dnu2[Y_]*dnu2[Y_].cst();

  f[++m] = scl_res*g2.cst();
}


void f_jacob(int n, double x[], double fvec[], double **df,
	     void (*vecfunc)(int, double [], double []))
{
  int    i, j;
  double h, temp, *f;

  const double eps = 1e-4;

  f = dvector(1, n);

  for (j = 1; j <= n; j++) {
    temp = x[j];
    h = eps*fabs(temp);
    if (h == 0e0) h = eps;
    x[j] = temp + h;
    h = x[j] - temp;
    (*vecfunc)(n, x, f);
    x[j] = temp;
    for (i = 1; i <= dnn; i++)
      df[i][j] = (f[i]-fvec[i])/h;
  }

  free_dvector(f, 1, n);
}


void H_zero(const double prm_tol, const int n_max, const bool prt_iter)
{
  const int m_max = 50;  // max no of constraints

  char   hs[m_max][max_str];
  int    i, j, n, i1, m1 = 0, check;
  double L, dprm_max, chi20, chi21, step, slope;
  double **A, *b, *dbn0, *dbn, *bn0, *bn, *bn_max, **Jacob;
  double *grad;
  tps    K_re, K_re_scl, dnu2[2], g2, abs2;

  const bool   prt      = true;
  const int    n_prt    = 9;
  const double stp_max  = 1e30, svd_eps = 1e-13;
  const double scl_bn[] = {0e0, 0e0, 1e0, 1e1, 1e2};

  // Interface to Numerical Recipes.
  // dnn      = 6; // Number of constraints.
  dnn      = 16; // Number of constraints.
  dnrfuncv = f_nl;
  // dfmin_: 1/2 f*f. Also computes dfvec.

  b = dvector(1, m_max); dbn0 = dvector(1, n_prm); dbn = dvector(1, n_prm);
  bn0 = dvector(1, n_prm); bn = dvector(1, n_prm); bn_max = dvector(1, n_prm);
  grad = dvector(1, n_prm); dfvec = dvector(1, dnn);
  A = dmatrix(1, m_max, 1, n_prm); Jacob = dmatrix(1, dnn, 1, n_prm);

  for (i = 1; i <= n_prm; i++) {
    L = get_L(prm[i-1], 1);
    if (L == 0e0) L = 1e0;
    bn_max[i] = bnL_max[prm_n[i-1]]/L; 

    bn[i] = get_bn(prm[i-1], 1, prm_n[i-1]);
  }

  // store initial sextupole strengths
  cout << endl;
  cout << "initial b3Ls:" << endl;
  for (i = 1; i <= n_prm; i++) {
    bn[i] = get_bn(prm[i-1], 1, prm_n[i-1]);
    cout << scientific << setprecision(3)
	 << setw(11) << get_bnL(prm[i-1], 1, prm_n[i-1])
	 << setw(2) << prm_n[i-1] << endl;
  }

  lieinit_(no_dnu, ss_dim, nd_tps, ndpt_tps, iref_tps, 0);

  n = 0;
  do {
    n++; n_iter++;
    cout << endl;
    for (i1 = 1; i1 <= n_prm; i1++) {
      set_bn_par(prm[i1-1], prm_n[i1-1], 7);

      get_dyn(K_re, g2, dnu2);

      K_re_scl = K_re*Id_scl;

      m1 = 0;

      A[++m1][i1] =
	scl_ksi1*h_ijklm_p(K_re, 1, 1, 0, 0, 1, 7);
      if (i1 == n_prm) {
	sprintf(hs[m1-1], "K_11001");
	b[m1] = -scl_ksi1*(ksi1[X_]*M_PI+h_ijklm(K_re, 1, 1, 0, 0, 1));
      }

      A[++m1][i1] =
	scl_ksi1*h_ijklm_p(K_re, 0, 0, 1, 1, 1, 7);
      if (i1 == n_prm) {
	sprintf(hs[m1-1], "K_00111");
	b[m1] = -scl_ksi1*(ksi1[Y_]*M_PI+h_ijklm(K_re, 0, 0, 1, 1, 1));
      }

      A[++m1][i1] =
	scl_nl*h_ijklm_p(K_re_scl, 2, 2, 0, 0, 0, 7);
      if (i1 == n_prm) {
	sprintf(hs[m1-1], "K_22000");
	b[m1] = -scl_nl*h_ijklm(K_re_scl, 2, 2, 0, 0, 0);
      }

      A[++m1][i1] =
	scl_nl*h_ijklm_p(K_re_scl, 0, 0, 2, 2, 0, 7);
      if (i1 == n_prm) {
	sprintf(hs[m1-1], "K_00220");
	b[m1] = -scl_nl*h_ijklm(K_re_scl, 0, 0, 2, 2, 0);
      }

      A[++m1][i1] =
	scl_nl*h_ijklm_p(K_re_scl, 1, 1, 1, 1, 0, 7);
      if (i1 == n_prm) {
	sprintf(hs[m1-1], "K_11110");
	b[m1] = -scl_nl*h_ijklm(K_re_scl, 1, 1, 1, 1, 0);
      }

      if (scl_flip != 0e0) {
	abs2 =
	  // sqr(h_ijklm(K_re_scl, 1, 1, 1, 1, 0)
	  //     +h_ijklm_p(K_re_scl, 1, 1, 1, 1, 0, 7)*tps(0e0, 7)
	  //     -h_ijklm(K_re_scl, 2, 2, 0, 0, 0)
	  //     -h_ijklm_p(K_re_scl, 2, 2, 0, 0, 0, 7)*tps(0e0, 7));
	  f_step(-(h_ijklm(K_re_scl, 1, 1, 1, 1, 0)
		   +h_ijklm_p(K_re_scl, 1, 1, 1, 1, 0, 7)*tps(0e0, 7))
		 /(h_ijklm(K_re_scl, 2, 2, 0, 0, 0)
		   +h_ijklm_p(K_re_scl, 2, 2, 0, 0, 0, 7)*tps(0e0, 7)));

	A[++m1][i1] = scl_flip*h_ijklm_p(abs2, 0, 0, 0, 0, 0, 7);
	if (i1 == n_prm) {
	  sprintf(hs[m1-1], "flip1  ");
	  b[m1] = -scl_flip*h_ijklm(abs2, 0, 0, 0, 0, 0);
	}

	abs2 =
	  // sqr(h_ijklm(K_re_scl, 1, 1, 1, 1, 0)
	  //     +h_ijklm_p(K_re_scl, 0, 0, 2, 2, 0, 7)*tps(0e0, 7)
	  //     -h_ijklm(K_re_scl, 1, 1, 1, 1, 0)
	  //     -h_ijklm_p(K_re_scl, 0, 0, 2, 2, 0, 7)*tps(0e0, 7));
	  f_step(-(h_ijklm(K_re_scl, 1, 1, 1, 1, 0)
		   +h_ijklm_p(K_re_scl, 1, 1, 1, 1, 0, 7)*tps(0e0, 7))
		 /(h_ijklm(K_re_scl, 0, 0, 2, 2, 0)
		   +h_ijklm_p(K_re_scl, 0, 0, 2, 2, 0, 7)*tps(0e0, 7)));

	A[++m1][i1] = scl_flip*h_ijklm_p(abs2, 0, 0, 0, 0, 0, 7);
	if (i1 == n_prm) {
	  sprintf(hs[m1-1], "flip2  ");
	  b[m1] = -scl_flip*h_ijklm(abs2, 0, 0, 0, 0, 0);
	}
      }

      // Use scaled K to avoid floating point issues.
      if (scl_fold1 == 0e0) {
	A[++m1][i1] = scl_nl*h_ijklm_p(K_re_scl, 3, 3, 0, 0, 0, 7);
	if (i1 == n_prm) {
	  sprintf(hs[m1-1], "K_33000");
	  b[m1] = -scl_nl*h_ijklm(K_re_scl, 3, 3, 0, 0, 0);
	}

	A[++m1][i1] = scl_nl*h_ijklm_p(K_re_scl, 0, 0, 3, 3, 0, 7);
	if (i1 == n_prm) {
	  sprintf(hs[m1-1], "K_00330");
	  b[m1] = -scl_nl*h_ijklm(K_re_scl, 0, 0, 3, 3, 0);
	}
      } else {
	A[++m1][i1] = scl_fold1*(h_ijklm_p(K_re_scl, 2, 2, 0, 0, 0, 7)/3e0
				 +h_ijklm_p(K_re_scl, 3, 3, 0, 0, 0, 7));
	if (i1 == n_prm) {
	  sprintf(hs[m1-1], "K_33000");
	  b[m1] =
	    -scl_fold1*(h_ijklm(K_re_scl, 2, 2, 0, 0, 0)/3e0
			+h_ijklm(K_re_scl, 3, 3, 0, 0, 0));
	}

	A[++m1][i1] = scl_fold1*(h_ijklm_p(K_re_scl, 0, 0, 2, 2, 0, 7)/3e0
				 +h_ijklm_p(K_re_scl, 0, 0, 3, 3, 0, 7));
	if (i1 == n_prm) {
	  sprintf(hs[m1-1], "K_00330");
	  b[m1] =
	    -scl_fold1*(h_ijklm(K_re_scl, 0, 0, 2, 2, 0)/3e0
			+h_ijklm(K_re_scl, 0, 0, 3, 3, 0));
	}
      }

      if (scl_fold2 == 0e0) {
	A[++m1][i1] = scl_nl*h_ijklm_p(K_re_scl, 1, 1, 2, 2, 0, 7);
	if (i1 == n_prm) {
	  sprintf(hs[m1-1], "K_11220");
	  b[m1] = -scl_nl*h_ijklm(K_re_scl, 1, 1, 2, 2, 0);
	}

	A[++m1][i1] = scl_nl*h_ijklm_p(K_re_scl, 2, 2, 1, 1, 0, 7);
	if (i1 == n_prm) {
	  sprintf(hs[m1-1], "K_22110");
	  b[m1] = -scl_nl*h_ijklm(K_re_scl, 2, 2, 1, 1, 0);
	}
      } else {
	A[++m1][i1] = scl_fold2*(h_ijklm_p(K_re_scl, 1, 1, 1, 1, 0, 7)/2e0
				 +h_ijklm_p(K_re_scl, 1, 1, 2, 2, 0, 7));
	if (i1 == n_prm) {
	  sprintf(hs[m1-1], "K_11220");
	  b[m1] =
	    -scl_fold2*(h_ijklm(K_re_scl, 1, 1, 1, 1, 0)/2e0
			+h_ijklm(K_re_scl, 1, 1, 2, 2, 0));
	}

	A[++m1][i1] = scl_fold2*(h_ijklm_p(K_re_scl, 1, 1, 1, 1, 0, 7)/2e0
				 +h_ijklm_p(K_re_scl, 2, 2, 1, 1, 0, 7));
	if (i1 == n_prm) {
	  sprintf(hs[m1-1], "K_22110");
	  b[m1] =
	    -scl_fold2*(h_ijklm(K_re_scl, 1, 1, 1, 1, 0)/2e0
			+h_ijklm(K_re_scl, 2, 2, 1, 1, 0));
	}
      }

      A[++m1][i1] = scl_ksi2*h_ijklm_p(K_re_scl, 1, 1, 0, 0, 2, 7);
      if (i1 == n_prm) {
	sprintf(hs[m1-1], "K_11002");
	b[m1] = -scl_ksi2*h_ijklm(K_re_scl, 1, 1, 0, 0, 2);
      }

      A[++m1][i1] = scl_ksi2*h_ijklm_p(K_re_scl, 0, 0, 1, 1, 2, 7);
      if (i1 == n_prm) {
	sprintf(hs[m1-1], "K_00112");
	b[m1] = -scl_ksi2*h_ijklm(K_re_scl, 0, 0, 1, 1, 2);
      }

      if (scl_fold3 == 0e0) {
	A[++m1][i1] = scl_nl*h_ijklm_p(K_re_scl, 1, 1, 0, 0, 3, 7);
	if (i1 == n_prm) {
	  sprintf(hs[m1-1], "K_11003");
	  b[m1] = -scl_nl*h_ijklm(K_re_scl, 1, 1, 0, 0, 3);
	}

	A[++m1][i1] = scl_nl*h_ijklm_p(K_re_scl, 0, 0, 1, 1, 3, 7);
	if (i1 == n_prm) {
	  sprintf(hs[m1-1], "K_00113");
	  b[m1] = -scl_nl*h_ijklm(K_re_scl, 0, 0, 1, 1, 3);
	}
      } else {
	A[++m1][i1] = scl_fold3*(h_ijklm_p(K_re_scl, 1, 1, 0, 0, 1, 7)/3e0
				 +h_ijklm_p(K_re_scl, 1, 1, 0, 0, 3, 7));
	if (i1 == n_prm) {
	  sprintf(hs[m1-1], "K_11003");
	  b[m1] =
	    -scl_fold3*(h_ijklm(K_re_scl, 1, 1, 0, 0, 1)/3e0
			+h_ijklm(K_re_scl, 1, 1, 0, 0, 3));
	}

	A[++m1][i1] = scl_fold3*(h_ijklm_p(K_re_scl, 0, 0, 1, 1, 1, 7)/3e0
				 +h_ijklm_p(K_re_scl, 0, 0, 1, 1, 3, 7));
	if (i1 == n_prm) {
	  sprintf(hs[m1-1], "K_00113");
	  b[m1] =
	    -scl_fold3*(h_ijklm(K_re_scl, 0, 0, 1, 1, 1)/3e0
			+h_ijklm(K_re_scl, 0, 0, 1, 1, 3));
	}
      }

      // A[++m1][i1] = scl_dnu2[X_]*h_ijklm_p(dnu2[X_], 0, 0, 0, 0, 0, 7);
      // if (i1 == n_prm) {
      // 	sprintf(hs[m1-1], "dnu2_x ");
      // 	b[m1] = -scl_dnu2[X_]*dnu2[X_].cst();
      // }

      // A[++m1][i1] = scl_dnu2[Y_]*h_ijklm_p(dnu2[Y_], 0, 0, 0, 0, 0, 7);
      // if (i1 == n_prm) {
      // 	sprintf(hs[m1-1], "dnu2_y ");
      // 	b[m1] = -scl_dnu2[Y_]*dnu2[Y_].cst();
      // }

      A[++m1][i1] = scl_res*h_ijklm_p(g2, 0, 0, 0, 0, 0, 7);
      if (i1 == n_prm) {
	sprintf(hs[m1-1], "g2     ");
	b[m1] = -scl_res*g2.cst();
      }

      clr_bn_par(prm[i1-1], prm_n[i1-1]);
    }

    cout << endl;
    cout << scientific << setprecision(3)
	 << "K_11110:" << setw(11) << h_ijklm(K_re_scl, 1, 1, 1, 1, 0) << endl;
    cout << scientific << setprecision(3)
	 << "K_22000:" << setw(11) << h_ijklm(K_re_scl, 2, 2, 0, 0, 0) << endl;
    cout << scientific << setprecision(3)
	 << "K_00220:" << setw(11) << h_ijklm(K_re_scl, 0, 0, 2, 2, 0) << endl;

    // Scale decapoles.
    for (j = 1; j <= n_prm; j++) {
      bn_max[j] /= scl_bn[prm_n[j-1]-1]; bn[j] /= scl_bn[prm_n[j-1]-1];
      for (i = 1; i <= m1; i++)
	A[i][j] *= scl_bn[prm_n[j-1]-1];
    }

    chi20 = 0e0;
    for (j = 1; j <= m1; j++)
      chi20 += sqr(b[j]);

    if (prt) {
      cout << endl;
      cout << n << " Ax = b:" << endl;
      cout << endl;
      for (i = 1; i <= m1; i++) {
	cout  << setw(3) << i << " " << hs[i-1];
	for (j = 1; j <= n_prm; j++)
	  cout << scientific << setprecision(2) << setw(10) << A[i][j];
	cout << scientific << setprecision(2) << setw(10) << b[i] << endl;
      }
    }

    if (false) {
      f_nl(m1, bn, b);

      f_jacob(n_prm, bn, b, Jacob, f_nl);

      cout << endl;
      cout << "Jacobian:" << endl;
      cout << endl;
      for (i = 1; i <= m1; i++) {
	cout  << setw(3) << i << " " << hs[i-1];
	for (j = 1; j <= n_prm; j++)
	  cout << scientific << setprecision(2) << setw(10) << Jacob[i][j];
	cout << endl;
      }

      exit(0);
    }

    // Check for floating point problem.
    for (i = 1; i <= n_prm; i++)
      if (fabs(bn[i]) > bn_max[i]) {
	cout << endl;
	cout << scientific << setprecision(3)
	     << "bn_max exceeded: " << fabs(bn[i])-bn_max[i]
	     << " (" << bn_max[i] << ")" << endl;
	bn[i] = sgn(bn[i])*bn_max[i];
      }

    SVD_lim(m1, n_prm, A, b, bn_max, svd_eps, bn, dbn);

    // Scale decapoles.
    for (j = 1; j <= n_prm; j++) {
      dbn[j] *= scl_bn[prm_n[j-1]-1];
      bn_max[j] *= scl_bn[prm_n[j-1]-1]; bn[j] *= scl_bn[prm_n[j-1]-1];
      for (i = 1; i <= m1; i++)
	A[i][j] /= scl_bn[prm_n[j-1]-1];
    }

    // Objective function for line search:
    //   g = f^T . f / 2 => nabla g = (df/dx)^T . f
    for (i = 1; i <= n_prm; i++) {
      grad[i] = 0e0;
      for (j = 1; j <= m1; j++)
	grad[i] += A[j][i]*(-b[j]); // b = -fvec.
    }

    // Check for floating point problem.
    slope = 0e0;
    for (i = 1; i <= n_prm; i++)
      slope += grad[i]*dbn[i];
    if (slope >= 0e0) {
      cout << endl;
      cout << scientific << setprecision(3)
	   << "Floating point problem: slope = " << slope << endl;
      exit(1);
    }

    for (i = 1; i <= n_prm; i++) {
       bn0[i] = bn[i]; dbn0[i] = dbn[i];
    }

    dlnsrch(n_prm, bn0, chi20/2e0, grad, dbn, bn, &chi21, stp_max, &check,
	    dfmin_);

    chi21 *= 2e0;

    // Find non-zero parameter change.
    j = 1;
    while ((j < n_prm) && (dbn[j] == 0e0))
      j++;

    step = (bn[j]-bn0[j])/dbn0[j];
 
    cout << endl;
    cout << "dcorr.:" << endl;
    for (i = 1; i <= n_prm; i++) {
      cout << scientific << setprecision(3)
	   << setw(11) << bn[i]-bn0[i];
      if (i % n_prt == 0) cout << endl;
    }
    if (n_prm % n_prt != 0) cout << endl;

    for (i = 1; i <= n_prm; i++)
      set_bn(prm[i-1], prm_n[i-1], bn[i]);

    cout << endl;
    cout << "dcorr. (int):" << endl;
    dprm_max = 0.0;
    for (i = 1; i <= n_prm; i++) {
      L = get_L(prm[i-1], 1);
      if (L == 0.0) L = 1e0;
      dprm_max = max(fabs((bn[i]-bn0[i])*L), dprm_max);
      cout << scientific << setprecision(3) << setw(11) << (bn[i]-bn0[i])*L;
      if (i % n_prt == 0) cout << endl;
    }
    if (n_prm % n_prt != 0) cout << endl;

    cout << endl;
    cout << "corr.:" << endl;
    for (i = 1; i <= n_prm; i++) {
      cout << scientific << setprecision(3)
	   << setw(11) << get_bnL(prm[i-1], 1, prm_n[i-1]);
      if (i % n_prt == 0) cout << endl;
    }
    if (n_prm % n_prt != 0) cout << endl;

    if (prt_iter) {
      sext_out << endl;
      sext_out << "n = " << n_iter << ":" << endl;
      for (i = 1; i <= n_prm; i++)
	for (j = 1; j <= get_n_Kids(abs(prm[i-1])); j++) {
	  sext_out << fixed << setprecision(7) 
		   << setw(9) << get_Name(abs(prm[i-1]))
		   << "(" << j << ") = "
		   << setw(11) << get_bnL(prm[i-1], 1, prm_n[i-1])
		   << setw(2) << prm_n[i-1] << endl;
	}
      sext_out.flush();
    }

    cout << endl;
    cout << scientific << setprecision(1)
	 << setw(2) << n_iter << ", step = " << step << ", slope = " << slope
	 << ", chi2: " << chi20 << " -> " << chi21 << endl;
  } while ((dprm_max > prm_tol) && (n < n_max));

  if (!prt_iter) {
    for (i = 1; i <= n_prm; i++)
      for (j = 1; j <= get_n_Kids(prm[i-1]); j++) {
	sext_out << fixed << setprecision(7) 
		 << setw(6) << get_Name(abs(prm[i-1])) << "(" << j << ") = "
		 << setw(11) << get_bnL(prm[i-1], 1, prm_n[i-1])
		 << setw(2) << prm_n[i-1] << endl;
      }
    sext_out.flush();
  }

  free_dvector(b, 1, m_max); free_dvector(dbn0, 1, n_prm);
  free_dvector(dbn, 1, n_prm); free_dvector(bn0, 1, n_prm);
  free_dvector(bn, 1, n_prm); free_dvector(bn_max, 1, n_prm);
  free_dvector(grad, 1, n_prm); free_dvector(dfvec, 1, dnn);
  free_dmatrix(A, 1, m_max, 1, n_prm); free_dmatrix(Jacob, 1, dnn, 1, n_prm);
}


void no_mpoles(void)
{
  int j, k;

  cout << endl;
  cout << "zeroing multipoles" << endl;
  cout << endl;
  for (j = 0; j < n_elem; j++)
    if (elem[j].kind == Mpole)
      for (k = Sext; k < mpole_max; k++) {
//	cout << "zeroing " << elem[j].Name << endl;
	set_bn(elem[j].Fnum, elem[j].Knum, k, 0.0);
      }
}


void get_prm(const char *file_name)
{
  char      line[max_str];      
  ifstream  prm_in;

  file_rd(prm_in, file_name);

  do
    prm_in.getline(line, max_str);
  while (strstr(line, "#") != NULL);

  sscanf(line, "%*s %d %lf %lf %lf", &adj_tune, &nu0[X_], &nu0[Y_], &eps_nu);

  prm_in.getline(line, max_str);
  sscanf(line, "%*s %d %lf %lf %lf",
	 &adj_chrom, &ksi1[X_], &ksi1[Y_], &eps_ksi);
  fit_chrm = eps_ksi < 0.0; eps_ksi = fabs(eps_ksi);

  prm_in.getline(line, max_str);
  sscanf(line, "%*s %lf", &ds_max);

  prm_in.getline(line, max_str);
  sscanf(line, "%*s %lf", &b2_max);

  prm_in.getline(line, max_str);
  sscanf(line, "%*s %lf", &bnL_max[Sext]);

  prm_in.getline(line, max_str);
  sscanf(line, "%*s %lf", &bnL_max[Oct]);

  prm_in.getline(line, max_str);
  sscanf(line, "%*s %lf", &bnL_max[Dec]);

  prm_in.getline(line, max_str);
  sscanf(line, "%*s %lf", &scl_res);

  prm_in.getline(line, max_str);
  sscanf(line, "%*s %lf", &scl_ksi1);

  prm_in.getline(line, max_str);
  sscanf(line, "%*s %lf", &scl_ksi2);

  prm_in.getline(line, max_str);
  sscanf(line, "%*s %lf", &scl_nl);
  prm_in.getline(line, max_str);
  sscanf(line, "%*s %lf", &scl_fold1);
  prm_in.getline(line, max_str);
  sscanf(line, "%*s %lf", &scl_fold2);
  prm_in.getline(line, max_str);
  sscanf(line, "%*s %lf", &scl_fold3);
  prm_in.getline(line, max_str);
  sscanf(line, "%*s %lf", &scl_flip);

  cout << endl;
  cout << fixed << setprecision(6)
       << "fit_tune      = " << adj_tune
       << ", nu0_x  = " << nu0[X_] << ", nu0_y  = " << nu0[Y_]
       << scientific << setprecision(1) << ", eps_nu = " << eps_nu << endl;
  cout << fixed << setprecision(6)
       << "fit_chrom     = " << adj_chrom
       << ", ksi0_x = " << ksi1[X_] << ", ksi0_y = " << ksi1[Y_]
       << scientific << setprecision(1) << ", eps_ksi = " << eps_ksi
       << ", fit_chrm = " << fit_chrm << endl;
  cout << endl;
  cout << fixed << setprecision(2)
       << "ds_max        = " << ds_max << endl;
  cout << fixed << setprecision(1)
       << "b2_max        = " << b2_max << endl;
  cout << fixed << setprecision(1)
       << "b3L_max       = " << bnL_max[Sext] << endl;
  cout << fixed << setprecision(1)
       << "b4L_max       = " << bnL_max[Oct] << endl;
  cout << fixed << setprecision(1)
       << "b5L_max       = " << bnL_max[Dec] << endl;
  cout << fixed << setprecision(1)
       << "scl_res       = " << scl_res << endl;
  cout << scientific << setprecision(1)
       << "scl_ksi1      = " << scl_ksi1 << endl;
  cout << scientific << setprecision(1)
       << "scl_ksi2      = " << scl_ksi2 << endl;
  cout << scientific << setprecision(1)
       << "scl_nl        = " << scl_nl << endl;
  cout << scientific << setprecision(1)
       << "scl_fold1     = " << scl_fold1 << endl;
  cout << scientific << setprecision(1)
       << "scl_fold2     = " << scl_fold2 << endl;
  cout << scientific << setprecision(1)
       << "scl_fold3     = " << scl_fold3 << endl;
  cout << scientific << setprecision(1)
       << "scl_flip      = " << scl_flip << endl;
}


void get_prm()
{

  n_prm = 0;

  switch (sext_scheme) {
  case 0:
  case 1:
  case 2:
  case 3:
  case 4:
  case 5:
  case 6:
  case 7:
  case 8:
  case 10:
    prm[n_prm] = get_Fnum("sl1"); prm_n[n_prm++] = Sext;
    prm[n_prm] = get_Fnum("sl2"); prm_n[n_prm++] = Sext;
    prm[n_prm] = get_Fnum("sl3"); prm_n[n_prm++] = Sext;

    prm[n_prm] = get_Fnum("sh1"); prm_n[n_prm++] = Sext;
    prm[n_prm] = get_Fnum("sh3"); prm_n[n_prm++] = Sext;
    prm[n_prm] = get_Fnum("sh4"); prm_n[n_prm++] = Sext;
    break;
  }

  switch (sext_scheme) {
  case 0:
    prm[n_prm] = get_Fnum("sm1"); prm_n[n_prm++] = Sext;
    prm[n_prm] = get_Fnum("sm2"); prm_n[n_prm++] = Sext;
    // prm[n_prm] = get_Fnum("sm3"); prm_n[n_prm++] = Sext;
    break;
  case 1:
    prm[n_prm] = get_Fnum("sm1a"); prm_n[n_prm++] = Sext;
    prm[n_prm] = get_Fnum("sm1b"); prm_n[n_prm++] = Sext;
    prm[n_prm] = get_Fnum("sm2");  prm_n[n_prm++] = Sext;
    break;
  case 2:
    prm[n_prm] = get_Fnum("sm1a"); prm_n[n_prm++] = Sext;
    prm[n_prm] = get_Fnum("sm1b"); prm_n[n_prm++] = Sext;
    prm[n_prm] = get_Fnum("sm2");  prm_n[n_prm++] = Sext;

    prm[n_prm] = get_Fnum("sm3");  prm_n[n_prm++] = Sext;
    break;
  case 3:
    prm[n_prm] = get_Fnum("sm1a"); prm_n[n_prm++] = Sext;
    prm[n_prm] = get_Fnum("sm1b"); prm_n[n_prm++] = Sext;
    prm[n_prm] = get_Fnum("sm2");  prm_n[n_prm++] = Sext;

    prm[n_prm] = get_Fnum("sm3");  prm_n[n_prm++] = Oct;
    break;
  case 4:
    prm[n_prm] = get_Fnum("sm1a"); prm_n[n_prm++] = Sext;
    prm[n_prm] = get_Fnum("sm1b"); prm_n[n_prm++] = Sext;
    prm[n_prm] = get_Fnum("sm2");  prm_n[n_prm++] = Sext;

    prm[n_prm] = get_Fnum("sm3");  prm_n[n_prm++] = Dec;
    break;
  case 5:
    prm[n_prm] = get_Fnum("sm1a"); prm_n[n_prm++] = Sext;
    prm[n_prm] = get_Fnum("sm1b"); prm_n[n_prm++] = Sext;
    prm[n_prm] = get_Fnum("sm2");  prm_n[n_prm++] = Sext;

    prm[n_prm] = get_Fnum("sm3a"); prm_n[n_prm++] = Oct;
    prm[n_prm] = get_Fnum("sm3b"); prm_n[n_prm++] = Oct;
    break;
  case 6:
    prm[n_prm] = get_Fnum("sm1a"); prm_n[n_prm++] = Sext;
    prm[n_prm] = get_Fnum("sm1b"); prm_n[n_prm++] = Sext;
    prm[n_prm] = get_Fnum("sm2");  prm_n[n_prm++] = Sext;

    prm[n_prm] = get_Fnum("sm3a"); prm_n[n_prm++] = Dec;
    prm[n_prm] = get_Fnum("sm3b"); prm_n[n_prm++] = Dec;
    break;
  case 7:
    prm[n_prm] = get_Fnum("sm1a"); prm_n[n_prm++] = Sext;
    prm[n_prm] = get_Fnum("sm1b"); prm_n[n_prm++] = Sext;
    prm[n_prm] = get_Fnum("sm2");  prm_n[n_prm++] = Sext;

    prm[n_prm] = get_Fnum("sm3a"); prm_n[n_prm++] = Oct;
    prm[n_prm] = get_Fnum("sm3b"); prm_n[n_prm++] = Oct;

    prm[n_prm] = get_Fnum("sm3a"); prm_n[n_prm++] = Dec;
    prm[n_prm] = get_Fnum("sm3b"); prm_n[n_prm++] = Dec;
    break;
  case 8:
    prm[n_prm] = get_Fnum("sm1a"); prm_n[n_prm++] = Sext;
    prm[n_prm] = get_Fnum("sm1b"); prm_n[n_prm++] = Sext;
    prm[n_prm] = get_Fnum("sm2");  prm_n[n_prm++] = Sext;

    prm[n_prm] = get_Fnum("sh1_dw"); prm_n[n_prm++] = Sext;
    prm[n_prm] = get_Fnum("sh3_dw"); prm_n[n_prm++] = Sext;
    prm[n_prm] = get_Fnum("sh4_dw"); prm_n[n_prm++] = Sext;
    break;
  case 9:
    // Carbon.
    prm[n_prm] = get_Fnum("s1a"); prm_n[n_prm++] = Sext;
    prm[n_prm] = get_Fnum("s2a"); prm_n[n_prm++] = Sext;
    prm[n_prm] = get_Fnum("s1b"); prm_n[n_prm++] = Sext;
    prm[n_prm] = get_Fnum("s2b"); prm_n[n_prm++] = Sext;
    prm[n_prm] = get_Fnum("s1c"); prm_n[n_prm++] = Sext;
    prm[n_prm] = get_Fnum("s2c"); prm_n[n_prm++] = Sext;
    prm[n_prm] = get_Fnum("s1d"); prm_n[n_prm++] = Sext;
    prm[n_prm] = get_Fnum("s2d"); prm_n[n_prm++] = Sext;
    break;
  case 10:
    prm[n_prm] = get_Fnum("sm1"); prm_n[n_prm++] = Sext;
    prm[n_prm] = get_Fnum("sm2");  prm_n[n_prm++] = Sext;

    prm[n_prm] = get_Fnum("sm3"); prm_n[n_prm++] = Oct;

    prm[n_prm] = get_Fnum("sm3"); prm_n[n_prm++] = Dec;
    break;
  }

  if (n_prm > n_prm_max) {
    cout << "get_prm: n_prm_max exceeded " << n_prm << "(" << n_prm_max
	 << ")" << endl;
    exit(0);
  }

  cout << endl;
  cout << "get_prm: no of multipole families " << n_prm << endl;
}


void chk_lat(double nu[], double ksi[])
{
  double        alpha1[2];
  ss_vect<tps>  nus;

//  get_Map();
  danot_(2);
  get_COD(10, 1e-10, 0.0, true);
  lieinit_(3, ss_dim, nd_tps, ndpt_tps, iref_tps, 0);
  K = MapNorm(Map, g, A1, A0, Map_res, 1);
  nus = dHdJ(K); get_nu_ksi(nus, nu, ksi); get_ab(alpha1, beta1, 0);
  danot_(no_tps);
  cout << endl;
  cout << fixed << setprecision(3)
       << "alpha_x  = " << alpha1[X_] << ", alpha_y = " << alpha1[Y_]
       << ", beta_x = " << beta1[X_] << ", beta_y  = " << beta1[Y_] << endl;
  prt_nu(nus);
}


void fit_chrom(void)
{
  int k, n_b3, b3s[n_prm_max];

  n_b3 = 0;
  for (k = 0; k < n_prm; k++) {
    if (prm_n[k] == Sext) {
      n_b3++; b3s[n_b3-1] = prm[k];
    }
  }

  // fit chromaticity
  cavity_on = false;
  fit_chrom(ksi1[X_], ksi1[Y_], n_b3, b3s, true);
//  fit_chrom1(0.0, 0.0, n_prm, prm, eps_ksi, true);
}


int main(int argc, char *argv[])
{
  double           nu[2], ksi[2];
  double           alpha[2], beta[2];
  tps              h_re, h_im, H, H_re, H_im, K_re, K_im, g_re, g_im;
  ss_vect<tps>     nus;
  ofstream         outf, K_out, nus_out;
  ifstream         inf;

  danot_(no_tps-1);

  rad_on    = false; H_exact        = false; totpath_on   = false;
  cavity_on = false; quad_fringe_on = false; emittance_on = false;
  IBS_on    = false;

  rd_mfile(argv[1], elem); rd_mfile(argv[1], elem_tps);

  // initialize the symplectic integrator after energy has been defined
  ini_si();

  // disable from TPSALib- and LieLib log messages
  idprset(-1);

  // Required for dnu2.
  // daeps_(1e-38);
  daeps_(1e-50);

  cout << endl;
  cout << scientific << setprecision(3)
       << "numeric_limits<double>::max_exponent = "
       << std::numeric_limits<double>::max_exponent<< endl;
  cout << scientific << setprecision(3)
       << "numeric_limits<double>::max_exponent = "
       << std::numeric_limits<double>::min_exponent << endl;

  sscanf(argv[2], "%d", &sext_scheme);

  danot_(3);

  if (true) chk_lat(nu, ksi);

  get_ab(alpha, beta, 0);
  cout << endl;
  cout << fixed << setprecision(3)
       << "alpha_x  = " << setw(6) << alpha[X_]
       << ", alpha_y = " << setw(6) << alpha[Y_]
       << ", beta_x = " << setw(6) << beta[X_]
       << ", beta_y  = " << setw(6) << beta[Y_] << endl;

  if (true) prt_alphac();

  if (false) {
    // no scaling
    Jx = 0.5, Jy = 0.5, delta = 1.0;
  } else {
    Jx = sqr(A_max[X_])/(2e0*beta1[X_]); Jy = sqr(A_max[Y_])/(2e0*beta1[Y_]);
    delta = delta_max;
  }

  Id_scl.identity();
  Id_scl[x_] *= sqrt(2e0*Jx); Id_scl[px_] *= sqrt(2e0*Jx);
  Id_scl[y_] *= sqrt(2e0*Jy); Id_scl[py_] *= sqrt(2e0*Jy);
  Id_scl[delta_] *= delta;

  get_prm("nl_ns.dat");

  if (adj_chrom) {
    n_prm = 0;
    switch (sext_scheme) {
    case 0:
    case 10:
      prm[n_prm] = get_Fnum("sm1"); prm_n[n_prm++] = Sext;
      prm[n_prm] = get_Fnum("sm2"); prm_n[n_prm++] = Sext;
//       prm[n_prm] = get_Fnum("sm3"); prm_n[n_prm++] = Sext;
      break;
    case 1:
    case 2:
    case 3:
    case 4:
    case 5:
    case 6:
    case 7:
    case 8:
      prm[n_prm] = get_Fnum("sm1a"); prm_n[n_prm++] = Sext;
      prm[n_prm] = get_Fnum("sm1b"); prm_n[n_prm++] = Sext;
      prm[n_prm] = get_Fnum("sm2");  prm_n[n_prm++] = Sext;
      break;
    case 9:
      // Carbon.
      prm[n_prm] = get_Fnum("s1a"); prm_n[n_prm++] = Sext;
      prm[n_prm] = get_Fnum("s2a"); prm_n[n_prm++] = Sext;
      prm[n_prm] = get_Fnum("s1b"); prm_n[n_prm++] = Sext;
      prm[n_prm] = get_Fnum("s2b"); prm_n[n_prm++] = Sext;
      prm[n_prm] = get_Fnum("s1c"); prm_n[n_prm++] = Sext;
      prm[n_prm] = get_Fnum("s2c"); prm_n[n_prm++] = Sext;
      prm[n_prm] = get_Fnum("s1d"); prm_n[n_prm++] = Sext;
      prm[n_prm] = get_Fnum("s2d"); prm_n[n_prm++] = Sext;
      break;
    }

    if (fit_chrm) {
      danot_(3);
      no_mpoles(); fit_chrom();
    }

    get_prm();

    file_wr(sext_out, "sext.dat");

    H_zero(eps_ksi, 25, true);

    sext_out.close();
  }


  danot_(no_tps-1);

  get_Map();

  file_wr(outf, "map.dat"); outf << Map; outf.close();

  danot_(no_tps);

  CtoR(get_h()*Id_scl, h_re, h_im);
  outf.open("h.dat"); outf << h_re; outf.close();

  K = MapNorm(Map, g, A1, A0, Map_res, no_tps); CtoR(K, K_re, K_im);

  outf.open("K.dat");
  outf << K_re << K_re*Id_scl;
  outf.close();

  nus = dHdJ(K);
  outf.open("nus.dat");
  outf << nus[3] << nus[3]*Id_scl << nus[4] << nus[4]*Id_scl;
  outf.close();

  CtoR(get_H()*Id_scl, H_re, H_im);
  outf.open("H.dat");
  outf << H_re << H_re*Id_scl;
  outf.close();

  CtoR(g*Id_scl, g_re, g_im);
  outf.open("g.dat");
  outf << g_im << g_im*Id_scl;
  outf.close();
}
