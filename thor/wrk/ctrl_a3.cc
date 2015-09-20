#define NO 6

#include <complex>
#include "thor_lib.h"

int  no_tps   = NO,
     ndpt_tps = 5;

extern tps           K, g;
extern ss_vect<tps>  Map, A0, A1, Map_res;

ss_vect<tps> Id_scl;

const bool   cod_prt = false;
const double cod_eps = 1e-14;


double get_an(const int Fnum, const int Knum, const int n)
{
  return elem[get_loc(Fnum, Knum)-1].mpole->an[n-1];
}


void set_dan(const int Fnum, const int Knum, const int n, const double dan)
{
  int  k;

  k = get_loc(Fnum, Knum) - 1;
  elem[k].mpole->an[n-1] += dan; elem_tps[k].mpole->an[n-1] += dan;
  if (n > elem[k].mpole->order) {
    elem[k].mpole->order = n; elem_tps[k].mpole->order = n;
  }
}


void set_dan(const int Fnum, const int n, const double dan)
{
  int  j;

  for (j = 1; j <= get_n_Kids(Fnum); j++)
    set_dan(Fnum, j, n, dan);
}


void set_danL(const int Fnum, const int Knum, const int n, const double danL)
{
  int  k;

  k = get_loc(Fnum, Knum) - 1;
  if (elem[k].L == 0.0) {
    elem[k].mpole->an[n-1] += danL; elem_tps[k].mpole->an[n-1] += danL;
  } else {
    elem[k].mpole->an[n-1] += danL/elem[k].L;
    elem_tps[k].mpole->an[n-1] += danL/elem[k].L;
  }
  if (n > elem[k].mpole->order) {
    elem[k].mpole->order = n; elem_tps[k].mpole->order = n;
  }
}


void set_danL(const int Fnum, const int n, const double danL)
{
  int  j;

  for (j = 1; j <= get_n_Kids(Fnum); j++)
    set_danL(Fnum, j, n, danL);
}


void set_an_par(const int Fnum, const int Knum, const int n, const int j)
{
  // set parameter dependence
  int     k;
  double  an;

  k = get_loc(Fnum, Knum) - 1;
  an = elem_tps[k].mpole->an[n-1].cst();
  elem_tps[k].mpole->an[n-1] = tps(an, j);
  if (n > elem_tps[k].mpole->order) elem_tps[k].mpole->order = n;
}


void set_an_par(const int Fnum, const int n, const int j)
{
  // set parameter dependence
  int  k;

  for (k = 1; k <= get_n_Kids(Fnum); k++)
    set_an_par(Fnum, k, n, j);
}


void clr_an_par(const int Fnum, const int Knum, const int n)
{
  // clear parameter dependence
  int     k;
  double  an;

  k = get_loc(Fnum, Knum) - 1;
  an = elem_tps[k].mpole->an[n-1].cst(); elem_tps[k].mpole->an[n-1] = an;
  // clear order
}


void clr_an_par(const int Fnum, const int n)
{
  // set parameter dependence
  int  k;

  for (k = 1; k <= get_n_Kids(Fnum); k++)
    clr_an_par(Fnum, k, n);
}


void get_dH_dbn(string fam_name)
{
  int             Fnum;
  complex<double> Hp, Hm, dH_da3;
  tps             H, H_re, H_im;

  const double da3 = 1e-2;

  // sm3a - sm3d has a 90 degree design roll.

  Fnum = get_Fnum(fam_name.c_str());

  set_danL(Fnum, Sext, -da3);

  danot_(no_tps-1);
  if (false)
    get_Map();
  else
    get_COD(10, 1e-14, 0.0, cod_prt);
  danot_(no_tps);
  K = MapNorm(Map, g, A1, A0, Map_res, 1); H = get_H(); CtoR(H, H_re, H_im);
  Hm =
    complex<double>(h_ijklm(H_re, 2, 0, 0, 1, 0), h_ijklm(H_im, 2, 0, 0, 1, 0));

  set_danL(Fnum, Sext, 2e0*da3);

  danot_(no_tps-1);
  if (false)
    get_Map();
  else
    get_COD(10, 1e-14, 0.0, cod_prt);
  danot_(no_tps);
  K = MapNorm(Map, g, A1, A0, Map_res, 1); H = get_H(); CtoR(H, H_re, H_im);
  Hp =
    complex<double>(h_ijklm(H_re, 2, 0, 0, 1, 0), h_ijklm(H_im, 2, 0, 0, 1, 0));

  set_danL(Fnum, Sext, -da3);

  dH_da3 = (Hp-Hm)/(2e0*da3);

  cout << scientific << setprecision(3)
       << setw(11) << real(dH_da3) << setw(11) << imag(dH_da3) << endl;
}


void ctrl_bn(void)
{
  const int m = 2, n = 4;

  int      Fnum, Knum, prm[n], i, j, k, m1, n1;
  double   *h, *dbn, *bn, *bn_max, **A;
  tps      H, H_re, H_im;
  ofstream outf;

  const bool   h_20010 = false, h_00500 = true;
  const int    n_iter  = 2, COD_n_iter = 10;
  const double SVD_eps = 1e-11, COD_eps = 1e-14;

  // const bool   symm[]     = {false, false, false, false};
  // const int    prm_type[] = {-Sext, -Sext, -Sext, -Sext};
  const bool   symm[]     = {true, true, true, true};
  const int    prm_type[] = {-Dec, -Dec, -Dec, -Dec};
  const double bnL_max[]  = {0e0, 0e0, 0e0, 1e0, 0e0, 200e0};

  h = dvector(1, m); dbn = dvector(1, n); bn = dvector(1, n);
  bn_max = dvector(1, n); A = dmatrix(1, m, 1, n);

  n1 = 0;
  prm[n1++] = get_loc(get_Fnum("sm3a"), 1) - 1;
  prm[n1++] = get_loc(get_Fnum("sm3b"), 1) - 1;
  prm[n1++] = get_loc(get_Fnum("sm3c"), 1) - 1;
  prm[n1++] = get_loc(get_Fnum("sm3d"), 1) - 1;
  // prm[n1++] = get_loc(get_Fnum("sm3a"), 2) - 1;
  // prm[n1++] = get_loc(get_Fnum("sm3b"), 2) - 1;
  // prm[n1++] = get_loc(get_Fnum("sm3c"), 2) - 1;
  // prm[n1++] = get_loc(get_Fnum("sm3d"), 2) - 1;

  for (j = 1; j <= n; j++) {
    Fnum = elem[prm[j-1]].Fnum; Knum = elem[prm[j-1]].Knum;

    if (prm_type[j-1] > 0)
      bn[j] = get_bn(Fnum, Knum, prm_type[j-1]);
    else
      bn[j] = get_an(Fnum, Knum, abs(prm_type[j-1]));
    bn_max[j] = bnL_max[abs(prm_type[j-1])];
  }

  for (i = 1; i <= n_iter; i++) {
    for (j = 1; j <= n; j++) {
      Fnum = elem[prm[j-1]].Fnum; Knum = elem[prm[j-1]].Knum;

      if (prm_type[j-1] > 0)
	if (!symm[j-1])
	  set_bn_par(Fnum, Knum, prm_type[j-1], 7);
	else {
	  set_bn_par(Fnum, 1, prm_type[j-1], 7);
	  set_bn_par(Fnum, 6, prm_type[j-1], 7);
	  set_bn_par(Fnum, 11, prm_type[j-1], 7);
	}
      else
	if (!symm[j-1])
	  set_an_par(Fnum, Knum, abs(prm_type[j-1]), 7);
	else {
	  set_an_par(Fnum, 1, abs(prm_type[j-1]), 7);
	  set_an_par(Fnum, 6, abs(prm_type[j-1]), 7);
	  set_an_par(Fnum, 11, abs(prm_type[j-1]), 7);
	}

      danot_(no_tps-1);

      if (false)
	get_Map();
      else
	get_COD(COD_n_iter, COD_eps, 0.0, cod_prt);

      danot_(no_tps);

      K = MapNorm(Map, g, A1, A0, Map_res, 1); H = get_H(); CtoR(H, H_re, H_im);
      H_re = H_re*Id_scl; H_im = H_im*Id_scl;

      if (j == 1) {
	cout << endl;
	cout << scientific << setprecision(3)
	     << "h_20010 = [" << h_ijklm(H_re, 2, 0, 0, 1, 0) << ", "
	     << setw(10) << h_ijklm(H_im, 2, 0, 0, 1, 0) << "]" << endl;
	cout << scientific << setprecision(3)
	     << "h_00500 = [" << setw(10) << h_ijklm(H_re, 0, 0, 5, 0, 0)
	     << ", " << h_ijklm(H_im, 0, 0, 5, 0, 0) << "]" << endl;
      }

      m1 = 0;
      if (h_20010) {
	A[++m1][j] = h_ijklm_p(H_re, 2, 0, 0, 1, 0, 7);
	A[++m1][j] = h_ijklm_p(H_im, 2, 0, 0, 1, 0, 7);
      }
      if (h_00500) {
	A[++m1][j] = h_ijklm_p(H_re, 0, 0, 5, 0, 0, 7);
	A[++m1][j] = h_ijklm_p(H_im, 0, 0, 5, 0, 0, 7);
      }

      if (prm_type[j-1] > 0)
	if (!symm[j-1])
	  clr_bn_par(Fnum, Knum, prm_type[j-1]);
	else {
	  clr_bn_par(Fnum, 1, prm_type[j-1]);
	  clr_bn_par(Fnum, 6, prm_type[j-1]);
	  clr_bn_par(Fnum, 11, prm_type[j-1]);
	}
      else
	if (!symm[j-1])
	  clr_an_par(Fnum, Knum, abs(prm_type[j-1]));
	else {
	  clr_an_par(Fnum, 1, abs(prm_type[j-1]));
	  clr_an_par(Fnum, 6, abs(prm_type[j-1]));
	  clr_an_par(Fnum, 11, abs(prm_type[j-1]));
	}
    }

    m1 = 0;
    if (h_20010) {
      h[++m1] = -h_ijklm(H_re, 2, 0, 0, 1, 0);
      h[++m1] = -h_ijklm(H_im, 2, 0, 0, 1, 0);
    }
    if (h_00500) {
      h[++m1] = -h_ijklm(H_re, 0, 0, 5, 0, 0);
      h[++m1] = -h_ijklm(H_im, 0, 0, 5, 0, 0);
    }

    cout << endl;
    for (j = 1; j <= m; j++) {
      for (k = 1; k <= n; k++)
	cout << scientific << setprecision(3) << setw(11) << A[j][k];
      cout << scientific << setprecision(3) << setw(11) << h[j] << endl;
    }

    SVD_lim(m, n, A, h, bn_max, SVD_eps, bn, dbn);

    for (j = 1; j <= n; j++) {
      Fnum = elem[prm[j-1]].Fnum; Knum = elem[prm[j-1]].Knum;

      if (prm_type[j-1] > 0) {
	if (!symm[j-1])
	  set_dbn(Fnum, Knum, prm_type[j-1], dbn[j]);
	else {
	  set_dbn(Fnum, 1, prm_type[j-1], dbn[j]);
	  set_dbn(Fnum, 6, prm_type[j-1], dbn[j]);
	  set_dbn(Fnum, 11, prm_type[j-1], dbn[j]);
	}
	bn[j] = get_bn(Fnum, Knum, prm_type[j-1]);
      } else {
	if (!symm[j-1])
	  set_dan(elem[prm[j-1]].Fnum, elem[prm[j-1]].Knum, abs(prm_type[j-1]),
		  dbn[j]);
	else {
	  set_dan(Fnum, 1, abs(prm_type[j-1]), dbn[j]);
	  set_dan(Fnum, 6, abs(prm_type[j-1]), dbn[j]);
	  set_dan(Fnum, 11, abs(prm_type[j-1]), dbn[j]);
	}
	bn[j] = get_an(Fnum, Knum, abs(prm_type[j-1]));
      }
    }

    outf.open("bn.dat", ios::out);
    cout << endl;
    cout << i << " bn:" << endl;
    for (j = 1; j <= n; j++) {
      Fnum = elem[prm[j-1]].Fnum; Knum = elem[prm[j-1]].Knum;

      cout << scientific << setprecision(3) << setw(11) << bn[j];
      outf << scientific << setprecision(16)
	   << setw(10) << elem[prm[j-1]].Name
	   << setw(3) << Fnum << setw(3) << Knum << setw(3) << prm_type[j-1]
	   << setw(24) << bn[j] << endl;
    }
    cout << endl;
    outf.close();

    danot_(no_tps-1);
  }

  if (false)
    get_Map();
  else
    get_COD(COD_n_iter, COD_eps, 0.0, cod_prt);

  danot_(no_tps);

  K = MapNorm(Map, g, A1, A0, Map_res, 1); H = get_H(); CtoR(H, H_re, H_im);
  H_re = H_re*Id_scl; H_im = H_im*Id_scl;

  cout << endl;
  cout << scientific << setprecision(3)
       << "h_20010 = [" << setw(10) << h_ijklm(H_re, 2, 0, 0, 1, 0) << ", "
       << h_ijklm(H_im, 2, 0, 0, 1, 0) << "]" << endl;
  cout << scientific << setprecision(3)
       << "h_00500 = [" << setw(10) << h_ijklm(H_re, 0, 0, 5, 0, 0) << ", "
       << h_ijklm(H_im, 0, 0, 5, 0, 0) << "]" << endl;

  prt_mfile("flat_file_corr.dat");

  free_dvector(h, 1, m); free_dvector(dbn, 1, n); free_dvector(bn, 1, n);
  free_dvector(bn_max, 1, n); free_dmatrix(A, 1, m, 1, n);
}


int main(int argc, char *argv[])
{
  int          k;
  double       twoJ[2], alpha[2], beta[2];
  tps          H, H_re, H_im, g_re, g_im, K_re, K_im;
  ss_vect<tps> nus;
  ofstream     outf;

  const bool   scl = true;
  const double A_max[] = {15e-3, 5e-3}, delta_max = 3e-2;

  danot_(no_tps-1);

  rad_on    = false; H_exact        = false; totpath_on   = false;
  cavity_on = false; quad_fringe_on = false; emittance_on = false;
  IBS_on    = false;


  rd_mfile(argv[1], elem); rd_mfile(argv[1], elem_tps);

  // initialize the symplectic integrator after energy has been defined
  ini_si();

  // disable from TPSALib- and LieLib log messages
  idprset(-1);

  daeps_(1e-30);

  danot_(no_tps-1);

  if (false)
    get_Map();
  else
    get_COD(10, 1e-14, 0.0, cod_prt);

  danot_(no_tps);

  K = MapNorm(Map, g, A1, A0, Map_res, 1); nus = dHdJ(K);
  CtoR(K, K_re, K_im); CtoR(g, g_re, g_im); H = get_H(); CtoR(H, H_re, H_im);

  get_ab(alpha, beta, 0);

  cout << endl;
  cout << fixed << setprecision(3)
       << "beta = [" << beta[X_] << ", " << beta[Y_] << "]" << endl;
  cout << fixed << setprecision(3)
       << "nu   = [" << nus[3].cst() << ", " << nus[4].cst() << "]" << endl;

  if (scl) {
    Id_scl.identity();
    for (k = 0; k < 2; k++) {
      twoJ[k] = sqr(A_max[k])/beta[k];
      Id_scl[2*k] *= sqrt(twoJ[k]); Id_scl[2*k+1] *= sqrt(twoJ[k]);
    }
    Id_scl[delta_] *= delta_max;

    nus[3] = nus[3]*Id_scl; nus[4] = nus[4]*Id_scl;
    K_re = K_re*Id_scl;
    g_re = g_re*Id_scl; g_im = g_im*Id_scl;
    H_re = H_re*Id_scl; H_im = H_im*Id_scl;
  }

  danot_(no_tps-1);
  outf.open("H.dat", ios::out);
  outf << H_re << H_im;
  outf.close();

  outf.open("K.dat", ios::out);
  outf << K_re << K_im;
  outf.close();

  outf.open("g.dat", ios::out);
  outf << g_re << g_im;
  outf.close();

  outf.open("nus.dat", ios::out);
  outf << nus[3] << nus[4];
  outf.close();

  if (false) {
    cout << endl;
    get_dH_dbn("sm3a");
    get_dH_dbn("sm3b");
    get_dH_dbn("sm3c");
    get_dH_dbn("sm3d");
  }

  if (true) ctrl_bn();
}
