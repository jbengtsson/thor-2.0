#define NO 11

#include "thor_lib.h"

int  no_tps   = NO,
     ndpt_tps = 5;


extern double        b2_max;
extern tps           K, g;
extern ss_vect<tps>  Map, A0, A1, Map_res;

const int max_ind = NO;

bool          h[max_ind][max_ind][max_ind][max_ind][max_ind], fit_chrm;
int           n_prm, prms[n_prm_max], bns[n_prm_max];
int           adj_tune, adj_chrom, n_iter, sext_scheme;
long int      beta_loc1, beta_loc2, beta_loc3, beta_loc4, rseed;
double        Jx, Jy, delta, beta1[2], beta2[2], beta3[2], beta4[2];
double        b3s[n_prm_max], chi2 = 0e0;
double        nu0[2], eps_nu, ksi1[2], eps_ksi;
double        nu_x_min, nu_x_max, nu_y_min, nu_y_max;
double        bnL_max[mpole_max];
double        scl_res, scl_dnu, scl_ksi1, scl_ksi_nl, scl_dnu2, scl_dksi2;
double        step;
ss_vect<tps>  Id_scl, map1, map2;
ofstream      sext_out;


// const double A_max[] = {20e-3, 7.5e-3};
const double A_max[] = {25e-3, 7.5e-3};
const double delta_max = 3.0e-2, delta_A_max[] = {10e-3, 0e-3};


void select_h(void)
{
  int  i, j, k, l, m;

  for (i = 0; i <= no_tps-1; i++)
    for (j = 0; j <= no_tps-1; j++)
      for (k = 0; k <= no_tps-1; k++)
	for (l = 0; l <= no_tps-1; l++)
	  for (m = 0; m <= no_tps-1; m++)
	    if (false)
	      // enable all driving terms
	      h[i][j][k][l][m] = true;
	    else
	      if (i+j+k+l <= 4)
		// enable all non-chromatic driving terms
		h[i][j][k][l][m] = (m == 0)? true : false;

  // linear chromaticity
  h[1][1][0][0][1] = true; h[0][0][1][1][1] = true;

  // 2nd order amplitude dependent tune shift
  h[2][2][0][0][0] = true; h[1][1][1][1][0] = true; h[0][0][2][2][0] = true;

  // 4th order amplitude dependent tune shift
  h[3][3][0][0][0] = true; h[2][2][1][1][0] = true;
  h[0][0][3][3][0] = true; h[1][1][2][2][0] = true;

  // 6th order amplitude dependent tune shift
  h[4][4][0][0][0] = true; h[0][0][4][4][0] = true;
  h[3][3][1][1][0] = true; h[1][1][3][3][0] = true;
  h[2][2][2][2][0] = true;

  if (true) {
    // nonlinear chromaticity
    for (m = 2; m <= no_tps-3; m++) {
      h[1][1][0][0][m] = true; h[0][0][1][1][m] = true;
    }
  }

  if (true) {
    // amplitude dependent chromaticity

    for (m = 1; m <= no_tps-5; m++) {
      h[2][2][0][0][m] = true; h[0][0][2][2][m] = true;
      h[1][1][1][1][m] = true;
    }

    for (m = 1; m <= no_tps-7; m++) {
      h[3][3][0][0][m] = true; h[0][0][3][3][m] = true;
      h[2][2][1][1][m] = true; h[1][1][2][2][m] = true;
    }
  }

  // exclude tune
  h[1][1][0][0][0] = false; h[0][0][1][1][0] = false;

  // exclude delta dependence of pathlength
  for (m = 2; m <= no_tps-1; m++)
    h[0][0][0][0][m] = false;

  if (true) {
    // higher order dispersion
    for (m = 2; m <= no_tps-2; m++) {
      h[1][0][0][0][m] = true; h[0][1][0][0][m] = true;
    }

    // delta dependance of the beta functions
    for (m = 1; m <= no_tps-3; m++) {
      h[2][0][0][0][m] = true; h[0][2][0][0][m] = true;
      h[0][0][2][0][m] = true; h[0][0][0][2][m] = true;
    }
  }
}


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


double f_num_fltr(const int jj[])
{

  return ((jj[x_] == jj[px_]) && (jj[y_] == jj[py_]))? 1e0 : 0e0;
}


double f_prm1(const int jj[])
{
  int     k, n;
  double  f;

  n = 0;
  for (k = 0; k < ss_dim; k++)
    n += jj[k];

  f =
    ((jj[ss_dim-1] > 1) || (jj[delta_] > 3) ||
     ((n == no_tps-2) && (jj[ss_dim-1] == 0)))?
    0e0 : 1e0;

 return f;
}


double f_prm2(const int jj[])
{

 return (jj[ss_dim-1] > 1)? 0e0 : 1e0;
}


double f_prm3(const int jj[])
{
  int     k, n;
  double  f;

  n = 0;
  for (k = 0; k < ss_dim; k++)
    n += jj[k];

  f = ((n == 5) && (jj[ss_dim-1] == 0))? 0e0 : 1e0;
  // f = ((n == 3) && (jj[ss_dim-1] == 0))? 0e0 : 1e0;

 return f;
}


void get_dnu2(const tps &K, tps &dnu2, tps &dksi2)
{
  // no must be 8+1 for J^3.
  // no must be 5+1 for delta^3.
  // [J, J^2, J^3]: [K_22000, K_33000, K_44000, ...].
  //                [nu_11000, nu_22000, nu_33000, ...].
  // [d, d^2, d^3]: [K_11001, K_11002, K_11003, ...].
  //                [nu_00001, nu_00002, nu_00003, ...].

  int          i;
  tps          dnu[2], dksi[2];
  ss_vect<tps> Id_tr, Id_no_delta, dnus, ps0, ps1;

  const double twoJ[]       = {2e0*Jx, 2e0*Jy};
  const double twoJ_x_delta = sqr(delta_A_max[X_])/(2.0*beta1[X_]);

  dnus = dHdJ(K);

  Id_tr.identity(); Id_tr[px_] = 2e0; Id_tr[py_] = 2e0;
  Id_no_delta.identity(); Id_no_delta[delta_] = 0e0;
  for (i = 0; i < 2; i++) {
    // Remove numeric noise.
    dnus[3+i] = dacfu1(dnus[3+i], f_num_fltr);
    // [h^+, h^-] -> [J, phi].
    dnus[3+i] = dacfu1(dnus[3+i]-dnus[3+i].cst(), f_prm1)*Id_tr;

    dnu[i] = dnus[3+i]*Id_no_delta; dksi[i] = dnus[3+i] - dnu2[i];
  }

  dnu2 = dacfu1(sqr(dnu[X_])+sqr(dnu[Y_]), f_prm2);
  dksi2 = dacfu1(sqr(dksi[X_])+sqr(dksi[Y_]), f_prm2);

  dnu2 = Int(dnu2, x_+1);
  ps1.identity(); ps1[x_] = twoJ[X_]; ps0.identity(); ps0[x_] = 0e0;
  dnu2 = dnu2*ps1 - dnu2*ps0;

  dnu2 = Int(dnu2, y_+1);
  ps1.identity(); ps1[y_] = twoJ[Y_]; ps0.identity(); ps0[y_] = 0e0;
  dnu2 = dnu2*ps1 - dnu2*ps0;

  dnu2 /= twoJ[X_]*twoJ[Y_];

  dksi2 = Int(dksi2, x_+1);
  ps1.identity(); ps1[x_] = twoJ_x_delta; ps0.identity(); ps0[x_] = 0e0;
  dksi2 = dksi2*ps1 - dksi2*ps0;

  ps1.identity(); ps1[y_] = 0e0; dksi2 = dksi2*ps1;

  dksi2 = Int(dksi2, delta_+1);
  ps1.identity(); ps1[delta_] = delta_max;
  ps0.identity(); ps0[delta_] = -delta_max;
  dksi2 = dksi2*ps1 - dksi2*ps0;

  dksi2 /= twoJ[X_]*2e0*delta_max;
}


void get_dnu2(const tps &K, tps dnu2[], tps dksi2[])
{
  // no must be 8+1 for J^3.
  // no must be 5+1 for delta^3.
  // [J, J^2, J^3]: [K_22000, K_33000, K_44000, ...].
  //                [nu_11000, nu_22000, nu_33000, ...].
  // [d, d^2, d^3]: [K_11001, K_11002, K_11003, ...].
  //                [nu_00001, nu_00002, nu_00003, ...].

  int          i;
  tps          dnu[2], dksi[2];
  ss_vect<tps> Id_tr, Id_no_delta, dnus, ps0, ps1;

  const double twoJ[]       = {2e0*Jx, 2e0*Jy};
  const double twoJ_x_delta = sqr(delta_A_max[X_])/(2.0*beta1[X_]);

  dnus = dHdJ(K);

  Id_tr.identity(); Id_tr[px_] = 2e0; Id_tr[py_] = 2e0;
  Id_no_delta.identity(); Id_no_delta[delta_] = 0e0;
  for (i = 0; i < 2; i++) {
    // Remove numeric noise.
    dnus[3+i] = dacfu1(dnus[3+i], f_num_fltr);
    // [h^+, h^-] -> [J, phi].
    dnus[3+i] = dacfu1(dnus[3+i]-dnus[3+i].cst(), f_prm1)*Id_tr;

    dnu[i] = dnus[3+i]*Id_no_delta; dksi[i] = dnus[3+i] - dnu2[i];
  }

  for (i = 0; i < 2; i++) {
    dnu2[i] = dacfu1(sqr(dnu[i]), f_prm2);
    dksi2[i] = dacfu1(sqr(dksi[i]), f_prm2);

    dnu2[i] = Int(dnu2[i], 1);
    ps1.identity(); ps1[x_] = twoJ[X_]; ps0.identity(); ps0[x_] = 0e0;
    dnu2[i] = dnu2[i]*ps1 - dnu2[i]*ps0;

    dnu2[i] = Int(dnu2[i], 3);
    ps1.identity(); ps1[y_] = twoJ[Y_]; ps0.identity(); ps0[y_] = 0e0;
    dnu2[i] = dnu2[i]*ps1 - dnu2[i]*ps0;

    dnu2[i] /= twoJ[X_]*twoJ[Y_];

    dksi2[i] = Int(dksi2[i], 1);
    ps1.identity(); ps1[x_] = twoJ_x_delta; ps0.identity(); ps0[x_] = 0e0;
    dksi2[i] = dksi2[i]*ps1 - dksi2[i]*ps0;

    ps1.identity(); ps1[y_] = 0e0; dksi2[i] = dksi2[i]*ps1;

    dksi2[i] = Int(dksi2[i], 5);
    ps1.identity(); ps1[delta_] = delta_max;
    ps0.identity(); ps0[delta_] = -delta_max;
    dksi2[i] = dksi2[i]*ps1 - dksi2[i]*ps0;

    dksi2[i] /= twoJ[X_]*2e0*delta_max;
  }
}


void get_dnu2(const tps &K, tps &dnu2)
{
  // Dnu^2 = Int{||dnu/dJ|| dJ_x dJ_y} = Int{sqrt((dnu/dJ)^2) )dJ_x dJ_y}.
  //
  // [J, J^2, J^3]: [K_22000, K_33000, K_44000, ...]
  //                [nu_11000, nu_22000, nu_33000, ...].
  //
  // [d, d^2, d^3]: [K_11001, K_11002, K_11003, ...]
  //                [nu_00001, nu_00002, nu_00003, ...].

  int           i, j;
  tps           dnudJ[2][2], dS;
  ss_vect<tps>  Id_tr, nus, dnus[2], ps1, ps0;

  danot_(7);
  nus = dHdJ(K);
  danot_(no_tps);

  for (i = 0; i < 2; i++) {
    nus[i] = dacfu1(nus[i], f_prm3);
    dnus[i] = dHdJ(-2e0*M_PI*nus[i]);

    for (j = 0; j < 2; j++)
      dnudJ[i][j] = dnus[i][3+j];
  }

  dS = dnudJ[X_][X_]*dnudJ[Y_][Y_] - dnudJ[X_][Y_]*dnudJ[Y_][X_];

  // [h^+, h^-] -> [J, phi], h^+ h^- = 2J.
  Id_tr.identity(); Id_tr[px_] = 2e0; Id_tr[py_] = 2e0;

  dS = dS*Id_tr;

  dnu2 = sqr(dS);

  if (true) dnu2 = sqrt(dnu2);

  dnu2 = Int(dnu2, x_+1);

  ps0.identity(); ps0[x_] = 0e0;
  ps1.identity(); ps1[x_] = 2e0*Jx;

  dnu2 = dnu2*ps1 - dnu2*ps0;

  dnu2 = Int(dnu2, y_+1);

  ps0.identity(); ps0[y_] = 0e0;
  ps1.identity(); ps1[y_] = 2e0*Jy;

  dnu2 = dnu2*ps1 - dnu2*ps0;

  // dnu2 = Int(dnu2, delta_+1);

  // ps0.identity(); ps0[delta_] = 0e0;
  // ps1.identity(); ps1[delta_] = delta;

  // dnu2 = dnu2*ps1 - dnu2*ps0;
}


void get_dnu2_(const tps &K, tps &dnu2)
{
  // Dnu^2 = Int{|dnu_x/dJ x dnu_y/dJ| dJ_x dJ_y}.
  //
  // [J, J^2, J^3]: [K_22000, K_33000, K_44000, ...]
  //                [nu_11000, nu_22000, nu_33000, ...].
  //
  // [d, d^2, d^3]: [K_11001, K_11002, K_11003, ...]
  //                [nu_00001, nu_00002, nu_00003, ...].

  int           i, j;
  tps           dnudJ[2][3], dS[3];
  ss_vect<tps>  Id_tr, nus, dnus[2], ps1, ps0;

  const bool prt = false;

  if (true)
    danot_(7);
  else
    // Terms up to J^2 and delta^3.
    danot_(5);
  nus = dHdJ(K);
  danot_(NO);

  for (i = 0; i < 2; i++) {
    nus[i] = dacfu1(nus[i], f_prm3);
    dnus[i] = dHdJ(-2e0*M_PI*nus[i]);

    for (j = 0; j < 2; j++)
      dnudJ[i][j] = dnus[i][3+j];

    dnudJ[i][Z_] = Der(dacfu1(nus[i], f_prm3), delta_+1);
  }

  dS[X_] = dnudJ[X_][Y_]*dnudJ[Y_][Z_] - dnudJ[X_][Z_]*dnudJ[Y_][Y_];
  dS[Y_] = -dnudJ[X_][X_]*dnudJ[Y_][Z_] + dnudJ[X_][Z_]*dnudJ[Y_][X_];
  dS[Z_] = dnudJ[X_][X_]*dnudJ[Y_][Y_] - dnudJ[X_][Y_]*dnudJ[Y_][X_];

  // [h^+, h^-] -> [J, phi], h^+ h^- = 2J.
  Id_tr.identity(); Id_tr[px_] = 2e0; Id_tr[py_] = 2e0;

  if (prt)
    cout << dacfu1(dS[X_]*Id_tr, f_prm2) << dacfu1(dS[Y_]*Id_tr, f_prm2)
	 << dacfu1(dS[Z_]*Id_tr, f_prm2) << endl;

  dnu2 = 0e0;
  for (i = 0; i < 3; i++)
    dnu2 += sqr(dS[i]*Id_tr);

  if (false) dnu2 = sqrt(dnu2);

  if (prt) {
    cout << dacfu1(dnu2, f_prm2) << endl;
    exit(0);
  }

  dnu2 = Int(dnu2, x_+1);

  ps0.identity(); ps0[x_] = 0e0;
  ps1.identity(); ps1[x_] = 2e0*Jx;

  dnu2 = dnu2*ps1 - dnu2*ps0;

  dnu2 = Int(dnu2, y_+1);

  ps0.identity(); ps0[y_] = 0e0;
  ps1.identity(); ps1[y_] = 2e0*Jy;

  dnu2 = dnu2*ps1 - dnu2*ps0;

  dnu2 = Int(dnu2, delta_+1);

  ps0.identity(); ps0[delta_] = 0e0;
  ps1.identity(); ps1[delta_] = delta;

  dnu2 = dnu2*ps1 - dnu2*ps0;
}


void H_zero(const double prm_tol, const int n_max, const bool prt_iter)
{
  // Minimize tune footprint.
  const int     m_max = 500;  // max no of constraints

  bool          first = true;
  char          hs[m_max][max_str], str[max_str];
  int           i, j, k, l, m, n, i1, m1 = 0;
  double        L, dprm_max, scl_L = 0.0;
  double        **A, *b, *w, **U, **V, *dbn, *bn, *bn_max, **A_inv;
  tps           K_re, K_im, K_re_scl, g_re, g_im, g_re_scl, g_im_scl;
  // tps           dnu2[2], dksi2[2];
  tps           dnu2, dksi2;

  const bool    prt    = true, mirror_sym = false;
  const int     n_prt  = 9;

  b = dvector(1, m_max); w = dvector(1, n_prm); dbn = dvector(1, n_prm);
  bn = dvector(1, n_prm); bn_max = dvector(1, n_prm);
  A = dmatrix(1, m_max, 1, n_prm); U = dmatrix(1, m_max, 1, n_prm);
  V = dmatrix(1, n_prm, 1, n_prm);
  A_inv = dmatrix(1, n_prm, 1, m_max);

  for (i = 1; i <= n_prm; i++) {
    if (prms[i-1] > 0) {
      // Note, Jacobian is a function of multipole strengths
      L = get_L(prms[i-1], 1);
      if (L == 0.0) L = 1.0;
      if (false && (strcmp(get_Name(prms[i-1]), "om1") == 0)) {
	cout << endl;
	cout << get_Name(prms[i-1]) << ": bnL_max set to 1.0" << endl;
	bn_max[i] = 1.0/L; 
      } else
	bn_max[i] = bnL_max[bns[i-1]]/L; 
    } else
      bn_max[i] = ds_max;

    bn[i] = get_bn_s(prms[i-1], 1, bns[i-1]);
  }

  if (first) {
    // store initial sextupole strengths
    first = false;
    cout << endl;
    cout << "initial b3s:" << endl;
    for (i = 1; i <= n_prm; i++) {
      b3s[i-1] = get_bnL_s(prms[i-1], 1, bns[i-1]);
      cout << scientific << setprecision(3)
	   << setw(11) << b3s[i-1] << setw(2) << bns[i-1] << endl;
    }

    select_h();
  } else if (true) {
    // initialize sextupoles
    for (i = 1; i <= n_prm; i++)
      set_bnL_s(prms[i-1], bns[i-1], b3s[i-1]);
  }

  n = 0;
  do {
    n++; n_iter++;
    cout << endl;
    for (i1 = 1; i1 <= n_prm; i1++) {
      if (prms[i1-1] > 0)
	set_bn_par(prms[i1-1], bns[i1-1], 7);
      else
	set_s_par(abs(prms[i1-1]), 7);

      scl_L = (prms[i1-1] > 0)? 1e0 : scl_ds;

//      Jx = 0.5; Jy = 0.5; delta = 1.0;

      danot_(no_tps-1);

      get_Map();

      danot_(no_tps);

      K = MapNorm(Map, g, A1, A0, Map_res, 1);

      CtoR(K, K_re, K_im); CtoR(g, g_re, g_im);

      K_re_scl = K_re*Id_scl;

      g_re_scl = scl_res*(g_re*Id_scl);
      g_im_scl = scl_res*(g_im*Id_scl);

      // get_dnu2(K, dnu2, dksi2);

      // mirror symmetric cell => g_re = 0
      m1 = 0;

      A[++m1][i1] = scl_L*h_ijklm_p(scl_dnu2*dnu2, 0, 0, 0, 0, 0, 7);
      if (i1 == n_prm) {
	sprintf(hs[m1-1], "dnu2   ");
	b[m1] = -(h_ijklm(scl_dnu2*dnu2, 0, 0, 0, 0, 0));
      }

      A[++m1][i1] = scl_L*h_ijklm_p(scl_dksi2*dksi2, 0, 0, 0, 0, 0, 7);
      if (i1 == n_prm) {
      	sprintf(hs[m1-1], "dksi2  ");
      	b[m1] = -(h_ijklm(scl_dksi2*dksi2, 0, 0, 0, 0, 0));
      }

      for (i = 0; i <= no_tps-1; i++)
	for (j = 0; j <= no_tps-1; j++)
	  for (k = 0; k <= no_tps-1; k++)
	    for (l = 0; l <= no_tps-1; l++)
	      for (m = 0; m <= no_tps-1; m++) {
		if ((0 < i+j+k+l+m) && (i+j+k+l+m <= no_tps-1)) {
		  if (m1 == m_max-1) {
		    cout << "m_max reached " << m1 << "(" << m_max << ")"
			 << endl;
		    exit(1);
		  }

		  if (h[i][j][k][l][m] &&
		      (fabs(h_ijklm(scl_ksi1*K_re_scl, i, j, k, l, m))
		       > 0e0)) {
		    if (is_h_ijklm(1, 1, 0, 0, 1, i, j, k, l, m)) {
		      // horizontal linear chromaticity
		      A[++m1][i1] =
			scl_L*scl_ksi1*h_ijklm_p(K_re_scl, i, j, k, l, m, 7);
		      if (i1 == n_prm) {
			sprintf(hs[m1-1], "K_%d%d%d%d%d", i, j, k, l, m);
			b[m1] =
			  -scl_ksi1*(ksi1[X_]*M_PI*2.0*Jx*delta
				     +h_ijklm(K_re_scl, i, j, k, l, m));
		      }
		    } else if (is_h_ijklm(0, 0, 1, 1, 1, i, j, k, l, m)) {
		      // vertical linear chromaticity
		      A[++m1][i1] =
			scl_L*scl_ksi1
			*h_ijklm_p(K_re_scl, i, j, k, l, m, 7);
		      if (i1 == n_prm) {
			sprintf(hs[m1-1], "K_%d%d%d%d%d", i, j, k, l, m);
			b[m1] =
			  -scl_ksi1*(ksi1[Y_]*M_PI*2.0*Jy*delta
				     +h_ijklm(K_re_scl, i, j, k, l, m));
		      }
		    }
		  }

		  if (h[i][j][k][l][m] &&
		      (fabs(h_ijklm(scl_ksi_nl*K_re_scl, i, j, k, l, m))
			      > 0e0)) {
		    if (is_h_ijklm(1, 1, 0, 0, 2, i, j, k, l, m) ||
			is_h_ijklm(0, 0, 1, 1, 2, i, j, k, l, m) ||
			is_h_ijklm(1, 1, 0, 0, 3, i, j, k, l, m) ||
			is_h_ijklm(0, 0, 1, 1, 3, i, j, k, l, m)) {
		      // nonlinear chromaticity
		      A[++m1][i1] =
			scl_L*scl_ksi_nl
			*h_ijklm_p(K_re_scl, i, j, k, l, m, 7);
		      if (i1 == n_prm) {
			sprintf(hs[m1-1], "K_%d%d%d%d%d", i, j, k, l, m);
			b[m1] = -scl_ksi_nl*h_ijklm(K_re_scl, i, j, k, l, m);
		      }
		    }
		  }

		  if (h[i][j][k][l][m] &&
		      (fabs(h_ijklm(scl_dnu*K_re_scl, i, j, k, l, m))
		       > 0e0)) {
		    if (is_h_ijklm(2, 2, 0, 0, 0, i, j, k, l, m) ||
			is_h_ijklm(0, 0, 2, 2, 0, i, j, k, l, m) ||
			is_h_ijklm(1, 1, 1, 1, 0, i, j, k, l, m)) {
		      // amplitude dependent tune shift
		      A[++m1][i1] =
			scl_L*scl_dnu*h_ijklm_p(K_re_scl, i, j, k, l, m, 7);
		      if (i1 == n_prm) {
			sprintf(hs[m1-1], "K_%d%d%d%d%d", i, j, k, l, m);
			b[m1] = -scl_dnu*h_ijklm(K_re_scl, i, j, k, l, m);
		      }
		    }
		  }

		  if (h[i][j][k][l][m] &&
		      (fabs(h_ijklm(g_im_scl, i, j, k, l, m)) > 0e0)) {
		    A[++m1][i1] = scl_L*h_ijklm_p(g_im_scl, i, j, k, l, m, 7);
		    if (i1 == n_prm) {
		      sprintf(hs[m1-1], "i_%d%d%d%d%d", i, j, k, l, m);
		      b[m1] = -h_ijklm(g_im_scl, i, j, k, l, m);
		    }

		    if (!mirror_sym) {
		      A[++m1][i1] =
			scl_L*h_ijklm_p(g_re_scl, i, j, k, l, m, 7);
		      if (i1 == n_prm) {
			sprintf(hs[m1-1], "r_%d%d%d%d%d", i, j, k, l, m);
			b[m1] = -h_ijklm(g_re_scl, i, j, k, l, m);
		      }
		    }
		  }
		}
	      }

      if (prms[i1-1] > 0)
	clr_bn_par(prms[i1-1], bns[i1-1]);
      else
	clr_s_par(abs(prms[i1-1]));
    }

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
    
    SVD_lim(m1, n_prm, A, b, bn_max, 1e-13, bn, dbn);

    cout << endl;
    cout << "dcorr. (int):" << endl;
    dprm_max = 0.0;
    for (i = 1; i <= n_prm; i++) {
      if (prms[i-1] > 0) {
	L = get_L(prms[i-1], 1);
	if (L == 0.0) L = 1.0;
	dprm_max = max(fabs(step*dbn[i]*L), dprm_max);
	cout << scientific << setprecision(3) << setw(11) << step*dbn[i]*L;
      } else {
	dprm_max = max(fabs(step*scl_ds*dbn[i]), dprm_max);
	cout << scientific << setprecision(3)
	     << setw(11) << step*scl_ds*dbn[i];
      }
      if (i % n_prt == 0) cout << endl;
    }
    if (n_prm % n_prt != 0) cout << endl;

    cout << endl;
    cout << "corr.:" << endl;
    for (i = 1; i <= n_prm; i++) {
      set_dbn_s(prms[i-1], bns[i-1], step*dbn[i]);
      bn[i] = get_bn_s(prms[i-1], 1, bns[i-1]);
      cout << scientific << setprecision(3)
	   << setw(11) << get_bnL_s(prms[i-1], 1, bns[i-1]);
      if (i % n_prt == 0) cout << endl;
    }
    if (n_prm % n_prt != 0) cout << endl;

    if (prt_iter) {
      sext_out << endl;
      sext_out << "n = " << n_iter << ":" << endl;
      for (i = 1; i <= n_prm; i++)
	for (j = 1; j <= get_n_Kids(abs(prms[i-1])); j++) {
	  if (prms[i-1] > 0) 
	    sext_out << fixed << setprecision(7) 
		     << setw(9) << get_Name(abs(prms[i-1]))
		     << "(" << j << ") = "
		     << setw(11) << get_bnL_s(prms[i-1], 1, bns[i-1])
		     << setw(2) << bns[i-1] << endl;
	  else {
	    sprintf(str, "du_%s", get_Name(abs(prms[i-1])));
	    sext_out << fixed << setprecision(7) 
		     << setw(9) << str << "(" << j << ") = "
		     << setw(11) << get_bnL_s(prms[i-1], 1, bns[i-1])
		     << setw(2) << 0 << endl;
	    sprintf(str, "dd_%s", get_Name(abs(prms[i-1])));
	    sext_out << fixed << setprecision(7) 
		      << setw(9) << str << "(" << j << ") = "
		     << setw(11) << -get_bnL_s(prms[i-1], 1, bns[i-1])
		     << setw(2) << 0 << endl;
	  }
	}
      sext_out.flush();
    }

    cout << endl;
    cout << scientific << setprecision(1)
	 << setw(2) << n_iter << " chi2: " << chi2;

    chi2 = 0.0;
    for (i = 1; i <= m1; i++)
      chi2 += sqr(b[i]);

    cout << scientific << setprecision(1)
	 << " -> " << chi2 << endl;
  } while ((dprm_max > prm_tol) && (n < n_max));

  if (!prt_iter) {
    for (i = 1; i <= n_prm; i++)
      for (j = 1; j <= get_n_Kids(abs(prms[i-1])); j++) {
	sext_out << fixed << setprecision(7) 
		 << setw(6) << get_Name(abs(prms[i-1])) << "(" << j << ") = "
		 << setw(11) << get_bnL_s(prms[i-1], 1, bns[i-1])
		 << setw(2) << bns[i-1] << endl;
      }
    sext_out.flush();
  }

  free_dvector(b, 1, m_max); free_dvector(w, 1, n_prm);
  free_dvector(dbn, 1, n_prm); free_dvector(bn, 1, n_prm);
  free_dvector(bn_max, 1, n_prm);
  free_dmatrix(A, 1, m_max, 1, n_prm);
  free_dmatrix(U, 1, m_max, 1, n_prm); free_dmatrix(V, 1, n_prm, 1, n_prm);
  free_dmatrix(A_inv, 1, n_prm, 1, m_max);
}


void get_dyn(tps &K_re, tps &g2, tps dnu2[], tps dksi2[])
{
  tps          K_im, g_re, g_im;
  ss_vect<tps> ps;

  danot_(no_tps-1);

  get_Map();

  danot_(no_tps);

  K = MapNorm(Map, g, A1, A0, Map_res, 1);

  CtoR(K, K_re, K_im); CtoR(g, g_re, g_im);

  danot_(5);
  g_re = dacfu1(g_re, f_prm3); g_im = dacfu1(g_im, f_prm3);
  danot_(no_tps);

  ps.identity();
  ps[x_] = sqrt(2e0*Jx); ps[px_] = sqrt(2e0*Jx);
  ps[y_] = sqrt(2e0*Jy); ps[py_] = sqrt(2e0*Jy);
  ps[delta_] = delta_max;

  g2 = (sqr(g_re)+sqr(g_im))*ps;

  if (false)
    get_dnu2(K, dnu2, dksi2);
  else
    get_dnu2_(K, dnu2[X_]);
}


void H_zero2(const double prm_tol, const int n_max, const bool prt_iter)
{
  // Minimize tune footprint.
  const int     m_max = 50;  // max no of constraints

  bool          first = true;
  char          hs[m_max][max_str], str[max_str];
  int           i, j, n, i1, m1 = 0;
  double        L, dprm_max, scl_L = 0.0, chi20, chi21, stp;
  double        **A, *b, *w, **U, **V, *dbn, *bn, *bn_max, **A_inv;
  tps           K_re, dnu2[2], dksi2[2], g2;

  const bool    prt    = true;
  const int     n_prt  = 9;

  const double y_scl = 10e0;

  b = dvector(1, m_max); w = dvector(1, n_prm); dbn = dvector(1, n_prm);
  bn = dvector(1, n_prm); bn_max = dvector(1, n_prm);
  A = dmatrix(1, m_max, 1, n_prm); U = dmatrix(1, m_max, 1, n_prm);
  V = dmatrix(1, n_prm, 1, n_prm);
  A_inv = dmatrix(1, n_prm, 1, m_max);

  for (i = 1; i <= n_prm; i++) {
    if (prms[i-1] > 0) {
      L = get_L(prms[i-1], 1);
      if (L == 0.0) L = 1.0;
      if (false && (strcmp(get_Name(prms[i-1]), "om1") == 0)) {
	cout << endl;
	cout << get_Name(prms[i-1]) << ": bnL_max set to 1.0" << endl;
	bn_max[i] = 1.0/L; 
      } else
	bn_max[i] = bnL_max[bns[i-1]]/L; 
    } else
      bn_max[i] = ds_max;

    bn[i] = get_bn_s(prms[i-1], 1, bns[i-1]);
  }

  if (first) {
    // store initial sextupole strengths
    first = false;
    cout << endl;
    cout << "initial b3s:" << endl;
    for (i = 1; i <= n_prm; i++) {
      b3s[i-1] = get_bnL_s(prms[i-1], 1, bns[i-1]);
      cout << scientific << setprecision(3)
	   << setw(11) << b3s[i-1] << setw(2) << bns[i-1] << endl;
    }

    select_h();
  } else if (true) {
    // initialize sextupoles
    for (i = 1; i <= n_prm; i++)
      set_bnL_s(prms[i-1], bns[i-1], b3s[i-1]);
  }

  n = 0;
  do {
    n++; n_iter++;
    cout << endl;
    for (i1 = 1; i1 <= n_prm; i1++) {
      if (prms[i1-1] > 0)
	set_bn_par(prms[i1-1], bns[i1-1], 7);
      else
	set_s_par(abs(prms[i1-1]), 7);

      scl_L = (prms[i1-1] > 0)? 1e0 : scl_ds;

//      Jx = 0.5; Jy = 0.5; delta = 1.0;

      get_dyn(K_re, g2, dnu2, dksi2);

      // mirror symmetric cell => g_re = 0
      m1 = 0;

      A[++m1][i1] =
	scl_L*scl_ksi1*h_ijklm_p(K_re, 1, 1, 0, 0, 1, 7);
      if (i1 == n_prm) {
	sprintf(hs[m1-1], "K_11001");
	b[m1] = -scl_ksi1*(ksi1[X_]*M_PI+h_ijklm(K_re, 1, 1, 0, 0, 1));
      }

      A[++m1][i1] =
	scl_L*scl_ksi1
	*h_ijklm_p(K_re, 0, 0, 1, 1, 1, 7);
      if (i1 == n_prm) {
	sprintf(hs[m1-1], "K_00111");
	b[m1] = -scl_ksi1*(ksi1[Y_]*M_PI+h_ijklm(K_re, 0, 0, 1, 1, 1));
      }

      A[++m1][i1] = scl_L*scl_dnu2*h_ijklm_p(dnu2[X_], 0, 0, 0, 0, 0, 7);
      if (i1 == n_prm) {
	sprintf(hs[m1-1], "dnu2_x ");
	b[m1] = -scl_dnu2*h_ijklm(dnu2[X_], 0, 0, 0, 0, 0);
      }

      if (false) {
	A[++m1][i1] =
	  scl_L*y_scl*scl_dnu2*h_ijklm_p(dnu2[Y_], 0, 0, 0, 0, 0, 7);
	if (i1 == n_prm) {
	  sprintf(hs[m1-1], "dnu2_y ");
	  b[m1] = -y_scl*scl_dnu2*h_ijklm(dnu2[Y_], 0, 0, 0, 0, 0);
	}

	A[++m1][i1] = scl_L*scl_dksi2*h_ijklm_p(dksi2[X_], 0, 0, 0, 0, 0, 7);
	if (i1 == n_prm) {
	  sprintf(hs[m1-1], "dksi2_x");
	  b[m1] = -scl_dksi2*h_ijklm(dksi2[X_], 0, 0, 0, 0, 0);
	}

	A[++m1][i1] = scl_L*scl_dksi2*h_ijklm_p(dksi2[Y_], 0, 0, 0, 0, 0, 7);
	if (i1 == n_prm) {
	  sprintf(hs[m1-1], "dksi2_y");
	  b[m1] = -h_ijklm(scl_dksi2*dksi2[Y_], 0, 0, 0, 0, 0);
	}
      }

      A[++m1][i1] = scl_L*scl_res*h_ijklm_p(g2, 0, 0, 0, 0, 0, 7);
      if (i1 == n_prm) {
	sprintf(hs[m1-1], "g2     ");
	b[m1] = -scl_res*h_ijklm(g2, 0, 0, 0, 0, 0);
      }

      if (prms[i1-1] > 0)
	clr_bn_par(prms[i1-1], bns[i1-1]);
      else
	clr_s_par(abs(prms[i1-1]));
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
    
    SVD_lim(m1, n_prm, A, b, bn_max, 1e-13, bn, dbn);

    cout << endl;
    stp = step;
    do {
      for (i = 1; i <= n_prm; i++) {
	set_dbn_s(prms[i-1], bns[i-1], stp*dbn[i]);
	bn[i] = get_bn_s(prms[i-1], 1, bns[i-1]);
      }

      get_dyn(K_re, g2, dnu2, dksi2);

      m1 = 0;

      b[++m1] = -scl_ksi1*(ksi1[X_]*M_PI+h_ijklm(K_re, 1, 1, 0, 0, 1));
      b[++m1] = -scl_ksi1*(ksi1[Y_]*M_PI+h_ijklm(K_re, 0, 0, 1, 1, 1));
      b[++m1] = -scl_dnu2*(h_ijklm(dnu2[X_], 0, 0, 0, 0, 0));
      if (false) {
	b[++m1] = -y_scl*scl_dnu2*(h_ijklm(dnu2[Y_], 0, 0, 0, 0, 0));
	b[++m1] = -scl_dksi2*(h_ijklm(dksi2[X_], 0, 0, 0, 0, 0));
	b[++m1] = -scl_dksi2*(h_ijklm(dksi2[Y_], 0, 0, 0, 0, 0));
      }
      b[++m1] = -scl_res*h_ijklm(g2, 0, 0, 0, 0, 0);

      chi21 = 0e0;
      for (j = 1; j <= m1; j++)
	chi21 += sqr(b[j]);

      cout << fixed << setprecision(3)
	   << "step = " << stp << endl;

      if (chi21 >= chi20) {
	// Roll back.
	cout << "rolling back" << endl;
	for (i = 1; i <= n_prm; i++) {
	  set_dbn_s(prms[i-1], bns[i-1], -stp*dbn[i]);
	  bn[i] = get_bn_s(prms[i-1], 1, bns[i-1]);
	}
	stp *= 0.9;
	if (stp < 1e-3) {
	  cout << "H_zero2: failed to converge." << endl;
	  exit(1);
	}
      }
    } while (chi21 >= chi20);

    cout << endl;
    cout << "dcorr. (int):" << endl;
    dprm_max = 0.0;
    for (i = 1; i <= n_prm; i++) {
      if (prms[i-1] > 0) {
	L = get_L(prms[i-1], 1);
	if (L == 0.0) L = 1.0;
	dprm_max = max(fabs(stp*dbn[i]*L), dprm_max);
	cout << scientific << setprecision(3) << setw(11) << stp*dbn[i]*L;
      } else {
	dprm_max = max(fabs(stp*scl_ds*dbn[i]), dprm_max);
	cout << scientific << setprecision(3)
	     << setw(11) << stp*scl_ds*dbn[i];
      }
      if (i % n_prt == 0) cout << endl;
    }
    if (n_prm % n_prt != 0) cout << endl;

    cout << endl;
    cout << "corr.:" << endl;
    for (i = 1; i <= n_prm; i++) {
      cout << scientific << setprecision(3)
	   << setw(11) << get_bnL_s(prms[i-1], 1, bns[i-1]);
      if (i % n_prt == 0) cout << endl;
    }
    if (n_prm % n_prt != 0) cout << endl;

    if (prt_iter) {
      sext_out << endl;
      sext_out << "n = " << n_iter << ":" << endl;
      for (i = 1; i <= n_prm; i++)
	for (j = 1; j <= get_n_Kids(abs(prms[i-1])); j++) {
	  if (prms[i-1] > 0) 
	    sext_out << fixed << setprecision(7) 
		     << setw(9) << get_Name(abs(prms[i-1]))
		     << "(" << j << ") = "
		     << setw(11) << get_bnL_s(prms[i-1], 1, bns[i-1])
		     << setw(2) << bns[i-1] << endl;
	  else {
	    sprintf(str, "du_%s", get_Name(abs(prms[i-1])));
	    sext_out << fixed << setprecision(7) 
		     << setw(9) << str << "(" << j << ") = "
		     << setw(11) << get_bnL_s(prms[i-1], 1, bns[i-1])
		     << setw(2) << 0 << endl;
	    sprintf(str, "dd_%s", get_Name(abs(prms[i-1])));
	    sext_out << fixed << setprecision(7) 
		      << setw(9) << str << "(" << j << ") = "
		     << setw(11) << -get_bnL_s(prms[i-1], 1, bns[i-1])
		     << setw(2) << 0 << endl;
	  }
	}
      sext_out.flush();
    }

    cout << endl;
    cout << fixed << setprecision(3)
	 << setw(2) << n_iter << ", step = " << stp
	 << scientific << setprecision(1)
	 << ", chi2: " << chi20 << " -> " << chi21 << endl;
  } while ((dprm_max > prm_tol) && (n < n_max));

  if (!prt_iter) {
    for (i = 1; i <= n_prm; i++)
      for (j = 1; j <= get_n_Kids(abs(prms[i-1])); j++) {
	sext_out << fixed << setprecision(7) 
		 << setw(6) << get_Name(abs(prms[i-1])) << "(" << j << ") = "
		 << setw(11) << get_bnL_s(prms[i-1], 1, bns[i-1])
		 << setw(2) << bns[i-1] << endl;
      }
    sext_out.flush();
  }

  free_dvector(b, 1, m_max); free_dvector(w, 1, n_prm);
  free_dvector(dbn, 1, n_prm); free_dvector(bn, 1, n_prm);
  free_dvector(bn_max, 1, n_prm);
  free_dmatrix(A, 1, m_max, 1, n_prm);
  free_dmatrix(U, 1, m_max, 1, n_prm); free_dmatrix(V, 1, n_prm, 1, n_prm);
  free_dmatrix(A_inv, 1, n_prm, 1, m_max);
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
  sscanf(line, "%*s %lf", &scl_dnu);

  prm_in.getline(line, max_str);
  sscanf(line, "%*s %lf", &scl_ksi1);

  prm_in.getline(line, max_str);
  sscanf(line, "%*s %lf", &scl_ksi_nl);

  prm_in.getline(line, max_str);
  sscanf(line, "%*s %lf", &scl_dnu2);

  prm_in.getline(line, max_str);
  sscanf(line, "%*s %lf", &scl_dksi2);

  prm_in.getline(line, max_str);
  sscanf(line, "%*s %lf", &step);

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
  cout << fixed << setprecision(1)
       << "scl_dnu       = " << scl_dnu << endl;
  cout << scientific << setprecision(1)
       << "scl_ksi1      = " << scl_ksi1 << endl;
  cout << scientific << setprecision(1)
       << "scl_ksi_nl    = " << scl_ksi_nl << endl;
  cout << scientific << setprecision(1)
       << "scl_dnu2      = " << scl_dnu2 << endl;
  cout << scientific << setprecision(1)
       << "scl_dksi2     = " << scl_dksi2 << endl;
  cout << fixed << setprecision(2)
       << "step          = " << step << endl;
}


void get_prms()
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
    prms[n_prm] = get_Fnum("sl1"); bns[n_prm++] = Sext;
    prms[n_prm] = get_Fnum("sl2"); bns[n_prm++] = Sext;
    prms[n_prm] = get_Fnum("sl3"); bns[n_prm++] = Sext;

    prms[n_prm] = get_Fnum("sh1"); bns[n_prm++] = Sext;
    prms[n_prm] = get_Fnum("sh3"); bns[n_prm++] = Sext;
    prms[n_prm] = get_Fnum("sh4"); bns[n_prm++] = Sext;
    break;
  }

  switch (sext_scheme) {
  case 0:
    prms[n_prm] = get_Fnum("sm1"); bns[n_prm++] = Sext;
    prms[n_prm] = get_Fnum("sm2"); bns[n_prm++] = Sext;
    prms[n_prm] = get_Fnum("sm3"); bns[n_prm++] = Sext;
    break;
  case 1:
    prms[n_prm] = get_Fnum("sm1a"); bns[n_prm++] = Sext;
    prms[n_prm] = get_Fnum("sm1b"); bns[n_prm++] = Sext;
    prms[n_prm] = get_Fnum("sm2");  bns[n_prm++] = Sext;
    break;
  case 2:
    prms[n_prm] = get_Fnum("sm1a"); bns[n_prm++] = Sext;
    prms[n_prm] = get_Fnum("sm1b"); bns[n_prm++] = Sext;
    prms[n_prm] = get_Fnum("sm2");  bns[n_prm++] = Sext;

    prms[n_prm] = get_Fnum("sm3");  bns[n_prm++] = Sext;
    break;
  case 3:
    prms[n_prm] = get_Fnum("sm1a"); bns[n_prm++] = Sext;
    prms[n_prm] = get_Fnum("sm1b"); bns[n_prm++] = Sext;
    prms[n_prm] = get_Fnum("sm2");  bns[n_prm++] = Sext;

    prms[n_prm] = get_Fnum("sm3");  bns[n_prm++] = Oct;
    break;
  case 4:
    prms[n_prm] = get_Fnum("sm1a"); bns[n_prm++] = Sext;
    prms[n_prm] = get_Fnum("sm1b"); bns[n_prm++] = Sext;
    prms[n_prm] = get_Fnum("sm2");  bns[n_prm++] = Sext;

    prms[n_prm] = get_Fnum("sm3");  bns[n_prm++] = Oct;

    prms[n_prm] = get_Fnum("sh1_dw"); bns[n_prm++] = Sext;
    prms[n_prm] = get_Fnum("sh3_dw"); bns[n_prm++] = Sext;
    prms[n_prm] = get_Fnum("sh4_dw"); bns[n_prm++] = Sext;
    break;
  case 5:
    prms[n_prm] = get_Fnum("sm1a"); bns[n_prm++] = Sext;
    prms[n_prm] = get_Fnum("sm1b"); bns[n_prm++] = Sext;
    prms[n_prm] = get_Fnum("sm2");  bns[n_prm++] = Sext;

    prms[n_prm] = get_Fnum("sm3a"); bns[n_prm++] = Oct;
    prms[n_prm] = get_Fnum("sm3b"); bns[n_prm++] = Oct;
    break;
  case 6:
    prms[n_prm] = get_Fnum("sm1a"); bns[n_prm++] = Sext;
    prms[n_prm] = get_Fnum("sm1b"); bns[n_prm++] = Sext;
    prms[n_prm] = get_Fnum("sm2");  bns[n_prm++] = Sext;

    prms[n_prm] = get_Fnum("sm3a"); bns[n_prm++] = Oct;
    prms[n_prm] = get_Fnum("sm3a"); bns[n_prm++] = Dec;
    prms[n_prm] = get_Fnum("sm3b"); bns[n_prm++] = Oct;
    prms[n_prm] = get_Fnum("sm3b"); bns[n_prm++] = Dec;
    break;
  case 7:
    prms[n_prm] = get_Fnum("sm1a"); bns[n_prm++] = Sext;
    prms[n_prm] = get_Fnum("sm1b"); bns[n_prm++] = Sext;
    prms[n_prm] = get_Fnum("sm2");  bns[n_prm++] = Sext;

    prms[n_prm] = get_Fnum("sm3a"); bns[n_prm++] = Sext;
    prms[n_prm] = get_Fnum("sm3a"); bns[n_prm++] = Oct;
    prms[n_prm] = get_Fnum("sm3a"); bns[n_prm++] = Dec;
    prms[n_prm] = get_Fnum("sm3b"); bns[n_prm++] = Sext;
    prms[n_prm] = get_Fnum("sm3b"); bns[n_prm++] = Oct;
    prms[n_prm] = get_Fnum("sm3b"); bns[n_prm++] = Dec;
    break;
  case 8:
    // Carbon.
    prms[n_prm] = get_Fnum("s1a"); bns[n_prm++] = Sext;
    prms[n_prm] = get_Fnum("s2a"); bns[n_prm++] = Sext;
    prms[n_prm] = get_Fnum("s1b"); bns[n_prm++] = Sext;
    prms[n_prm] = get_Fnum("s2b"); bns[n_prm++] = Sext;
    prms[n_prm] = get_Fnum("s1c"); bns[n_prm++] = Sext;
    prms[n_prm] = get_Fnum("s2c"); bns[n_prm++] = Sext;
    prms[n_prm] = get_Fnum("s1d"); bns[n_prm++] = Sext;
    prms[n_prm] = get_Fnum("s2d"); bns[n_prm++] = Sext;
    break;
  }

  if (n_prm > n_prm_max) {
    cout << "get_prms: n_prm_max exceeded " << n_prm << "(" << n_prm_max
	 << ")" << endl;
    exit(0);
  }

  cout << endl;
  cout << "get_prms: no of multipole families " << n_prm << endl;
}


void chk_lat(double nu[], double ksi[])
{
  double        alpha1[2];
  ss_vect<tps>  nus;

//  get_Map();
  get_COD(10, 1e-10, 0.0, true);
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
    if (bns[k] == Sext) {
      n_b3++; b3s[n_b3-1] = prms[k];
    }
  }

  // fit chromaticity
  cavity_on = false;
  fit_chrom(ksi1[X_], ksi1[Y_], n_b3, b3s, true);
//  fit_chrom1(0.0, 0.0, n_prm, prms, eps_ksi, true);
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
    // RHIC
//    Jx = sqr(0.5e-3)/(2e0*beta1[X_]); Jy = sqr(0.5e-3)/(2e0*beta1[Y_]);
//    delta = 1e-2;
  }

  Id_scl.identity();
  Id_scl[x_] *= sqrt(2e0*Jx); Id_scl[px_] *= sqrt(2e0*Jx);
  Id_scl[y_] *= sqrt(2e0*Jy); Id_scl[py_] *= sqrt(2e0*Jy);
  Id_scl[delta_] *= delta;

  get_prm("nl_opt.dat");

  if (adj_chrom) {
    n_prm = 0;
    switch (sext_scheme) {
    case 0:
      prms[n_prm] = get_Fnum("sm1"); bns[n_prm++] = Sext;
      prms[n_prm] = get_Fnum("sm2"); bns[n_prm++] = Sext;
//       prms[n_prm] = get_Fnum("sm3"); bns[n_prm++] = Sext;
      break;
    case 1:
    case 2:
    case 3:
    case 4:
    case 5:
    case 6:
    case 7:
      prms[n_prm] = get_Fnum("sm1a"); bns[n_prm++] = Sext;
      prms[n_prm] = get_Fnum("sm1b"); bns[n_prm++] = Sext;
      prms[n_prm] = get_Fnum("sm2");  bns[n_prm++] = Sext;
      break;
    case 8:
      // Carbon.
      prms[n_prm] = get_Fnum("s1a"); bns[n_prm++] = Sext;
      prms[n_prm] = get_Fnum("s2a"); bns[n_prm++] = Sext;
      prms[n_prm] = get_Fnum("s1b"); bns[n_prm++] = Sext;
      prms[n_prm] = get_Fnum("s2b"); bns[n_prm++] = Sext;
      prms[n_prm] = get_Fnum("s1c"); bns[n_prm++] = Sext;
      prms[n_prm] = get_Fnum("s2c"); bns[n_prm++] = Sext;
      prms[n_prm] = get_Fnum("s1d"); bns[n_prm++] = Sext;
      prms[n_prm] = get_Fnum("s2d"); bns[n_prm++] = Sext;
      break;
    }

    if (fit_chrm) {
      danot_(3);
      no_mpoles(); fit_chrom();
    }

    get_prms();

    file_wr(sext_out, "sext.dat");

    H_zero2(eps_ksi, 50, true);

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
