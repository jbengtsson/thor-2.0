#define NO 10

#include "thor_lib.h"

int  no_tps   = NO,
     ndpt_tps = 5;


extern double        b2_max;
extern tps           K, g;
extern ss_vect<tps>  Map, A0, A1, Map_res;

const int max_ind = NO;

bool          fit_chrm;
int           n_prm, prm[n_prm_max], prm_n[n_prm_max];
int           adj_tune, adj_chrom, n_iter, sext_scheme;
long int      beta_loc1, beta_loc2, beta_loc3, beta_loc4, rseed;
double        Jx, Jy, delta, beta1[2], beta2[2], beta3[2], beta4[2];
double        nu0[2], eps_nu, ksi1[2], eps_ksi;
double        nu_x_min, nu_x_max, nu_y_min, nu_y_max;
double        bnL_max[mpole_max], scl_res, scl_ksi1, scl_dnu2, chi2_min;
ss_vect<tps>  Id_scl, map1, map2;
ofstream      sext_out;


// const double A_max[] = {20e-3, 7.5e-3};
const double A_max[] = {25e-3, 7.5e-3};
const double delta_max = 3.0e-2, delta_A_max[] = {10e-3, 0e-3};

const int no_g = 5, no_dnu = 9;


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

  f = ((n == 5) && (jj[ss_dim-1] == 0))? 0e0 : 1e0;

 return f;
}


double f_dnu(const int jj[])
{
  int     k, n;
  double  f;

  n = 0;
  for (k = 0; k < ss_dim; k++)
    n += jj[k];

  f = ((n > 7) || (jj[ss_dim-1] > 1) || ((n == 7) && (jj[ss_dim-1] == 0)) ||
       (jj[delta_] > 3) ||
       ((jj[ss_dim-1] == 0) && (jj[delta_] != 0) && (n > 4)) ||
       ((jj[ss_dim-1] == 1) && (jj[delta_] != 0) && (n > 5)))?
    0e0 : 1e0;

 return f;
}


tps dfdJ(const tps &f, const int i)
{
  // df/dJ in resonance basis:
  // 
  // df/dJ = df/dh^+ dh^+/dJ + df/dh^- dh^-/dJ
  //       = f/dh^+ 1/h^- + df/dh^- 1/dh^+

  return pseudo_der(Der(f, 2*i-1), 2*i) + pseudo_der(Der(f, 2*i), 2*i-1);
}


void get_dnu2(const tps &K_re, tps &dnu2)
{
  // Dnu^2 = Int{|dnu_x/dJ x dnu_y/dJ| dJ_x dJ_y}.
  //
  // [J, J^2, J^3]: [K_22000, K_33000, K_44000, ...]
  //                [nu_11000, nu_22000, nu_33000, ...].
  //
  // [d, d^2, d^3]: [K_11001, K_11002, K_11003, ...]
  //                [nu_00001, nu_00002, nu_00003, ...].

  int           i, j;
  tps           nu[2], dnudJ[2][3], dS[3];
  ss_vect<tps>  Id_tr, ps1, ps0;

  const bool prt = false, threeD = true;

  for (i = 0; i < 2; i++) {
    nu[i] = -dfdJ(K_re, i+1)/(2e0*M_PI);
    nu[i] = dacfu1(nu[i], f_dnu);

    if (prt) cout << nu[i] << endl;

    for (j = 0; j < 2; j++)
      dnudJ[i][j] = dfdJ(nu[i], j+1);

    dnudJ[i][Z_] = Der(nu[i], delta_+1);

    for (j = 0; j < 3; j++)
      dnudJ[i][j] = dacfu1(dnudJ[i][j], f_fltr);

    if (prt) cout << dnudJ[i][X_] << dnudJ[i][Y_] << dnudJ[i][Z_] << endl;
  }

  // [h^+, h^-] -> [J, phi], h^+ h^- = 2J.
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

  if (threeD) {
    dnu2 = 0e0;
    for (i = 0; i < 3; i++)
      dnu2 += sqr(dS[i]);
  } else
    dnu2 = sqr(dS[Z_]);

  // Undo scaling.
  ps1.identity(); ps1[x_] /= 2e0*Jx; ps1[y_] /= 2e0*Jy; ps1[delta_] /= delta;
  dnu2 = dnu2*ps1;

  if (prt) {
    cout << dacfu1(dnu2, f_trunc_par) << endl;
    exit(0);
  }

  if (false) dnu2 = sqrt(dnu2);

  dnu2 = Int(dnu2, x_+1);
  ps0.identity(); ps0[x_] = 0e0; ps1.identity(); ps1[x_] = 2e0*Jx;
  dnu2 = dnu2*ps1 - dnu2*ps0;

  dnu2 = Int(dnu2, y_+1);
  ps0.identity(); ps0[y_] = 0e0; ps1.identity(); ps1[y_] = 2e0*Jy;
  dnu2 = dnu2*ps1 - dnu2*ps0;

  if (threeD) {
    dnu2 = Int(dnu2, delta_+1);
    ps0.identity(); ps0[delta_] = 0e0; ps1.identity(); ps1[delta_] = delta;
    dnu2 = dnu2*ps1 - dnu2*ps0;
  } else {
    ps0.identity(); ps0[delta_] = 0e0;
    dnu2 = dnu2*ps0;
  }
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


void get_dyn(tps &K_re, tps &g2, tps &dnu2)
{
  tps K_im;

  danot_(no_dnu-1);
  get_Map();
  danot_(no_dnu);

  lieinit_(no_dnu, ss_dim, nd_tps, ndpt_tps, iref_tps, 0);

  K = MapNorm(Map, g, A1, A0, Map_res, 1); CtoR(K, K_re, K_im);

  danot_(no_g);
  g = dacfu1(g, f_g);
  danot_(no_tps);

  g2 = get_abs2(g*Id_scl); get_dnu2(K_re, dnu2);
}


double f_nl(double bn[])
{
  int    i, j;
  double chi2;
  tps    K_re, g2, dnu2;

  const bool prt = true;

  n_iter++;

  for (i = 1; i <= n_prm; i++)
    set_bn(prm[i-1], prm_n[i-1], bn[i]);

  get_dyn(K_re, g2, dnu2);

  chi2 = 0e0;

  chi2 += scl_ksi1*sqr(ksi1[X_]*M_PI+h_ijklm(K_re, 1, 1, 0, 0, 1));
  chi2 += scl_ksi1*sqr(ksi1[Y_]*M_PI+h_ijklm(K_re, 0, 0, 1, 1, 1));
  chi2 += scl_dnu2*dnu2.cst();

  chi2 += scl_res*g2.cst();

  if (chi2 < chi2_min) {
    if (prt) {
      cout << scientific << setprecision(3)
	   << "ksi = "
	   << scl_ksi1*sqr(ksi1[X_]*M_PI+h_ijklm(K_re, 1, 1, 0, 0, 1))
	   << " " << scl_ksi1*sqr(ksi1[Y_]*M_PI+h_ijklm(K_re, 0, 0, 1, 1, 1))
	   << endl;
      cout << scientific << setprecision(3)
	   << "g2 = " << scl_res*g2.cst()
	   << ", dnu2 = " << scl_dnu2*dnu2.cst() << endl;
    }

    chi2_min = min(chi2, chi2_min);

    cout << endl;
    cout << "bnL:";
    for (i = 1; i <= n_prm; i++)
      cout << scientific << setprecision(3)
	   << setw(11) << get_bnL(prm[i-1], 1, prm_n[i-1]);
    cout << endl;

    cout << endl;
    cout << scientific << setprecision(1)
	 << setw(2) << n_iter << ", chi2_min: " << chi2_min << endl;

    sext_out << endl;
    sext_out << "n = " << n_iter << ":" << endl;
    for (i = 1; i <= n_prm; i++)
      for (j = 1; j <= get_n_Kids(prm[i-1]); j++) {
	sext_out << fixed << setprecision(7) 
		 << setw(9) << get_Name(prm[i-1])
		 << "(" << j << ") = "
		 << setw(11) << get_bnL(prm[i-1], 1, prm_n[i-1])
		 << setw(2) << prm_n[i-1] << endl;
      }

    sext_out.flush();
  }

  return chi2;
}


void df_nl(double bn[], double df[])
{
  int    k;
  tps    K_re, g2, dnu2;

  cout << "df_nl" << endl;

  for (k = 1; k <= n_prm; k++)
    set_bn(prm[k-1], prm_n[k-1], bn[k]);

  for (k = 1; k <= n_prm; k++) {
    set_bn_par(prm[k-1], prm_n[k-1], 7);

    get_dyn(K_re, g2, dnu2);

    df[k] = 0e0;

    df[k] +=
      scl_ksi1*2e0*h_ijklm_p(K_re, 1, 1, 0, 0, 1, 7)
      *(ksi1[X_]*M_PI+h_ijklm(K_re, 1, 1, 0, 0, 1));

    df[k] +=
      scl_ksi1*2e0*h_ijklm_p(K_re, 0, 0, 1, 1, 1, 7)
      *(ksi1[Y_]*M_PI+h_ijklm(K_re, 0, 0, 1, 1, 1));

    df[k] += scl_dnu2*h_ijklm_p(dnu2, 0, 0, 0, 0, 0, 7);

    df[k] += scl_res*h_ijklm_p(g2, 0, 0, 0, 0, 0, 7);

    clr_bn_par(prm[k-1], prm_n[k-1]);
  }
}


void H_zero(const double eps)
{
  int    i, iter;
  double *bn, fret;

  bn = dvector(1, n_prm);

  // store initial sextupole strengths
  cout << endl;
  cout << "initial bns:" << endl;
  for (i = 1; i <= n_prm; i++) {
    bn[i] = get_bn(prm[i-1], 1, prm_n[i-1]);
    cout << scientific << setprecision(3)
	 << setw(11) << get_bnL(prm[i-1], 1, prm_n[i-1])
	 << setw(2) << prm_n[i-1] << endl;
  }

  n_iter = 0; chi2_min = 1e30;

  dfrprmn(bn, n_prm, eps, &iter, &fret, f_nl, df_nl);

  cout << endl;
  cout << scientific << setprecision(1)
       << setw(2) << iter << ", chi2: " << fret << endl;

  free_dvector(bn, 1, n_prm);
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
  sscanf(line, "%*s %lf", &scl_dnu2);

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
       << "scl_dnu2      = " << scl_dnu2 << endl;
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
    prm[n_prm] = get_Fnum("sm3"); prm_n[n_prm++] = Sext;
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

    prm[n_prm] = get_Fnum("sm3");  prm_n[n_prm++] = Oct;

    prm[n_prm] = get_Fnum("sh1_dw"); prm_n[n_prm++] = Sext;
    prm[n_prm] = get_Fnum("sh3_dw"); prm_n[n_prm++] = Sext;
    prm[n_prm] = get_Fnum("sh4_dw"); prm_n[n_prm++] = Sext;
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

    prm[n_prm] = get_Fnum("sm3a"); prm_n[n_prm++] = Oct;
    prm[n_prm] = get_Fnum("sm3a"); prm_n[n_prm++] = Dec;
    prm[n_prm] = get_Fnum("sm3b"); prm_n[n_prm++] = Oct;
    prm[n_prm] = get_Fnum("sm3b"); prm_n[n_prm++] = Dec;
    break;
  case 7:
    prm[n_prm] = get_Fnum("sm1a"); prm_n[n_prm++] = Sext;
    prm[n_prm] = get_Fnum("sm1b"); prm_n[n_prm++] = Sext;
    prm[n_prm] = get_Fnum("sm2");  prm_n[n_prm++] = Sext;

    prm[n_prm] = get_Fnum("sm3a"); prm_n[n_prm++] = Sext;
    prm[n_prm] = get_Fnum("sm3a"); prm_n[n_prm++] = Oct;
    prm[n_prm] = get_Fnum("sm3a"); prm_n[n_prm++] = Dec;
    prm[n_prm] = get_Fnum("sm3b"); prm_n[n_prm++] = Sext;
    prm[n_prm] = get_Fnum("sm3b"); prm_n[n_prm++] = Oct;
    prm[n_prm] = get_Fnum("sm3b"); prm_n[n_prm++] = Dec;
    break;
  case 8:
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
  daeps_(1e-38);

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

  get_prm("nl_cg.dat");

  if (adj_chrom) {
    n_prm = 0;
    switch (sext_scheme) {
    case 0:
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
      prm[n_prm] = get_Fnum("sm1a"); prm_n[n_prm++] = Sext;
      prm[n_prm] = get_Fnum("sm1b"); prm_n[n_prm++] = Sext;
      prm[n_prm] = get_Fnum("sm2");  prm_n[n_prm++] = Sext;
      break;
    case 8:
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
      no_mpoles(); //fit_chrom();
    }

    get_prm();

    file_wr(sext_out, "sext.dat");

    H_zero(1e-20);

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
