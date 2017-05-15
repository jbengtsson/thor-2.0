#define NO 6

#include "thor_lib.h"

int  no_tps   = NO,
     ndpt_tps = 5;


extern double        b2_max;
extern tps           K, g;
extern ss_vect<tps>  Map, A0, A1, Map_res;

const bool  tune_scn = false;

const int   max_ind = 10;

bool          h[max_ind][max_ind][max_ind][max_ind][max_ind], fit_chrm;
int           n_bn, bns_fam[n_prm_max], bns_n[n_prm_max];
int           check_range, adj_tune, adj_chrom, n_steps, n_cell;
long int      beta_loc1, beta_loc2, beta_loc3;
double        Jx, Jy, delta, beta1[2], beta2[2], beta3[2], b3s[n_prm_max];
double        nu0[2], eps_nu, ksi1[2], eps_ksi;
double        nu_x_min, nu_x_max, nu_y_min, nu_y_max;
double        bnL_max[mpole_max];
double        scl_dnu, scl_ksi1, scl_ksi2, scl_ksi_nl;
double        scl_dnuddelta, scl_dnudJ, step;
ss_vect<tps>  Id_scl;
ofstream      sext_out;


const bool    sext_prm = true, oct_prm = false, dodec_prm = false;

const double  max_Ax = 20e-3, max_Ay = 10e-3, max_delta = 3e-2;


void select_h(void)
{
  int  i, j, k, l, m;

  for (i = 0; i <= no_tps-1; i++)
    for (j = 0; j <= no_tps-1; j++)
      for (k = 0; k <= no_tps-1; k++)
	for (l = 0; l <= no_tps-1; l++)
	  for (m = 0; m <= no_tps-1; m++)
	    if (false)
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
  // balance nonlinear terms
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

    // balance nonlinear terms
//    h[1][1][0][0][3] = false; h[0][0][1][1][3] = false;
//    h[1][1][0][0][5] = false; h[0][0][1][1][5] = false;
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


void H_zero(const double bn_tol, const int n_max, const bool prt_iter)
{
  const int     m_max = 500;  // max no of constraints

  bool          first = true;
  char          hs[m_max][max_str], str[max_str];
  int           i, j, k, l, m, i1, m1 = 0, n;
  double        L, dbn_max, chi2 = 0.0, scl_L = 0.0;
  double        **A, *b, *w, **U, **V, *dbn, *bn, *bn_max, **A_inv;
  tps           r, K_re, K_im, g_re, g_im;
  ss_vect<tps>  Id_Fl, A1_inv, A_nl_inv;

  const bool    prt             = true, mirror_sym = false;
  const int     n_prt           = 9;
  const double  dnudJ_twoJx     = sqr(20e-3)/beta1[X_];
  const double  dnudJ_twoJy     = sqr(7.5e-3)/beta1[Y_];
  const double  dnuddelta_delta = 2.5e-2;

  b = dvector(1, m_max); w = dvector(1, n_bn); dbn = dvector(1, n_bn);
  bn = dvector(1, n_bn); bn_max = dvector(1, n_bn);
  A = dmatrix(1, m_max, 1, n_bn); U = dmatrix(1, m_max, 1, n_bn);
  V = dmatrix(1, n_bn, 1, n_bn);
  A_inv = dmatrix(1, n_bn, 1, m_max);

  for (i = 1; i <= n_bn; i++) {
    if (bns_fam[i-1] > 0) {
      // Note, Jacobian is a function of multipole strengths
      L = get_L(bns_fam[i-1], 1);
      if (L == 0.0) L = 1.0;
      bn_max[i] = bnL_max[bns_n[i-1]]/L; 
    } else
      bn_max[i] = ds_max;

    bn[i] = get_bn_s(bns_fam[i-1], 1, bns_n[i-1]);
 }

  if (first) {
    // store initial sextupole strengths
    first = false;
    cout << endl;
    cout << "initial b3s:" << endl;
    for (i = 1; i <= n_bn; i++) {
      b3s[i-1] = get_bnL_s(bns_fam[i-1], 1, bns_n[i-1]);
      cout << scientific << setprecision(3)
	   << setw(11) << b3s[i-1] << setw(2) << bns_n[i-1] << endl;
    }

    select_h();
  } else if (true) {
    // initialize sextupoles
    for (i = 1; i <= n_bn; i++)
      set_bnL_s(bns_fam[i-1], bns_n[i-1], b3s[i-1]);
  }

  n = 0;
  do {
    n++;
    cout << endl;
    for (i1 = 1; i1 <= n_bn; i1++) {
      if (bns_fam[i1-1] > 0)
	set_bn_par(bns_fam[i1-1], bns_n[i1-1], 7);
      else
	set_s_par(abs(bns_fam[i1-1]), 7);

//      Jx = 0.5; Jy = 0.5; delta = 1.0;

      danot_(no_tps-1);
      get_Map();
      danot_(no_tps);
      K = MapNorm(Map, g, A1, A0, Map_res, 1);
      CtoR(K*Id_scl, K_re, K_im); CtoR(g*Id_scl, g_re, g_im);

      // mirror symmetric cell => g_re = 0
      m1 = 0;
      for (i = 0; i <= no_tps-1; i++)
	for (j = 0; j <= no_tps-1; j++)
	  for (k = 0; k <= no_tps-1; k++)
	    for (l = 0; l <= no_tps-1; l++)
	      for (m = 0; m <= no_tps-1; m++) {
		if ((0 < i+j+k+l+m) && (i+j+k+l+m <= no_tps-1) &&
		    h[i][j][k][l][m] &&
		    ((fabs(h_ijklm(g_im, i, j, k, l, m)) > 0.0) ||
		     (fabs(h_ijklm(K_re, i, j, k, l, m)) > 0.0))) {

		  scl_L = (bns_fam[i1-1] > 0)? 1.0 : scl_ds;

		  A[++m1][i1] = scl_L*h_ijklm_p(g_im, i, j, k, l, m, 7);

		  if (m1 >= m_max) {
		    cout << "m_max reached " << m1 << "(" << m_max << ")"
			 << endl;
		    exit(1);
		  }

		  if (is_h_ijklm(1, 1, 0, 0, 1, i, j, k, l, m))
		    // horizontal linear chromaticity
		    A[m1][i1] =
		      scl_L*scl_ksi1*h_ijklm_p(K_re, i, j, k, l, m, 7);
		  else if (is_h_ijklm(0, 0, 1, 1, 1, i, j, k, l, m))
		    // vertical linear chromaticity
		    A[m1][i1] =
		      scl_L*scl_ksi1*h_ijklm_p(K_re, i, j, k, l, m, 7);
		  else if (is_h_ijklm(2, 2, 0, 0, 0, i, j, k, l, m) ||
			   is_h_ijklm(0, 0, 2, 2, 0, i, j, k, l, m) ||
			   is_h_ijklm(1, 1, 1, 1, 0, i, j, k, l, m) ||
			   is_h_ijklm(3, 3, 0, 0, 0, i, j, k, l, m) ||
			   is_h_ijklm(2, 2, 1, 1, 0, i, j, k, l, m) ||
			   is_h_ijklm(1, 1, 2, 2, 0, i, j, k, l, m) ||
			   is_h_ijklm(0, 0, 3, 3, 0, i, j, k, l, m) ||
			   is_h_ijklm(4, 4, 0, 0, 0, i, j, k, l, m) ||
			   is_h_ijklm(0, 0, 4, 4, 0, i, j, k, l, m) ||
			   is_h_ijklm(3, 3, 1, 1, 0, i, j, k, l, m) ||
			   is_h_ijklm(1, 1, 3, 3, 0, i, j, k, l, m) ||
			   is_h_ijklm(2, 2, 2, 2, 0, i, j, k, l, m))
		    // amplitude dependent tune shift
		    A[m1][i1] =
		      scl_L*scl_dnu*h_ijklm_p(K_re, i, j, k, l, m, 7);
		  else if (is_h_ijklm(1, 1, 0, 0, 2, i, j, k, l, m) ||
			   is_h_ijklm(0, 0, 1, 1, 2, i, j, k, l, m))
		    // 2nd order chromaticity
		    A[m1][i1] =
		      scl_L*scl_ksi2*h_ijklm_p(K_re, i, j, k, l, m, 7);
		  else if (is_h_ijklm(1, 1, 0, 0, 2, i, j, k, l, m) ||
			   is_h_ijklm(0, 0, 1, 1, 2, i, j, k, l, m) ||
			   is_h_ijklm(1, 1, 0, 0, 3, i, j, k, l, m) ||
			   is_h_ijklm(0, 0, 1, 1, 3, i, j, k, l, m) ||
			   is_h_ijklm(1, 1, 0, 0, 4, i, j, k, l, m) ||
			   is_h_ijklm(0, 0, 1, 1, 4, i, j, k, l, m) ||
			   is_h_ijklm(1, 1, 0, 0, 5, i, j, k, l, m) ||
			   is_h_ijklm(0, 0, 1, 1, 5, i, j, k, l, m) ||
			   is_h_ijklm(1, 1, 0, 0, 6, i, j, k, l, m) ||
			   is_h_ijklm(0, 0, 1, 1, 6, i, j, k, l, m))
		    // nonlinear chromaticity
		    A[m1][i1] =
		      scl_L*scl_ksi_nl*h_ijklm_p(K_re, i, j, k, l, m, 7);
		  else if (is_h_ijklm(1, 1, 1, 1, 1, i, j, k, l, m) ||
			   is_h_ijklm(2, 2, 0, 0, 2, i, j, k, l, m) ||
			   is_h_ijklm(0, 0, 2, 2, 2, i, j, k, l, m) ||
			   is_h_ijklm(1, 1, 1, 1, 2, i, j, k, l, m) ||
			   is_h_ijklm(2, 2, 0, 0, 3, i, j, k, l, m) ||
			   is_h_ijklm(0, 0, 2, 2, 3, i, j, k, l, m) ||
			   is_h_ijklm(1, 1, 1, 1, 3, i, j, k, l, m) ||
			   is_h_ijklm(2, 2, 0, 0, 4, i, j, k, l, m) ||
			   is_h_ijklm(0, 0, 2, 2, 4, i, j, k, l, m) ||
			   is_h_ijklm(1, 1, 1, 1, 4, i, j, k, l, m) ||
			   is_h_ijklm(3, 3, 0, 0, 1, i, j, k, l, m) ||
			   is_h_ijklm(2, 2, 1, 1, 1, i, j, k, l, m) ||
			   is_h_ijklm(1, 1, 2, 2, 1, i, j, k, l, m) ||
			   is_h_ijklm(0, 0, 3, 3, 1, i, j, k, l, m) ||
			   is_h_ijklm(3, 3, 0, 0, 2, i, j, k, l, m) ||
			   is_h_ijklm(2, 2, 1, 1, 2, i, j, k, l, m) ||
			   is_h_ijklm(1, 1, 2, 2, 2, i, j, k, l, m) ||
			   is_h_ijklm(0, 0, 3, 3, 2, i, j, k, l, m))
		    // amplitude dependent chromaticity
		    A[m1][i1] =
		      scl_L*scl_ksi_nl*h_ijklm_p(K_re, i, j, k, l, m, 7);
		  else if (!mirror_sym)
		    A[++m1][i1] = scl_L*h_ijklm_p(g_re, i, j, k, l, m, 7);

		  if (i1 == n_bn) {
                    sprintf(hs[m1-1], "h_%d%d%d%d%d", i, j, k, l, m);
                    b[m1] = -h_ijklm(g_im, i, j, k, l, m);

		    if (is_h_ijklm(1, 1, 0, 0, 1, i, j, k, l, m)) {
		      // horizontal linear chromaticity
		      b[m1] =
			-scl_ksi1*(n_cell*ksi1[X_]*M_PI*2.0*Jx*delta
				   +h_ijklm(K_re, i, j, k, l, m));
		    } else if (is_h_ijklm(0, 0, 1, 1, 1, i, j, k, l, m)) {
		      // vertical linear chromaticity
		      b[m1] =
			-scl_ksi1*(n_cell*ksi1[Y_]*M_PI*2.0*Jy*delta
			  +h_ijklm(K_re, i, j, k, l, m));
		    } else if (is_h_ijklm(2, 2, 0, 0, 0, i, j, k, l, m) ||
			       is_h_ijklm(0, 0, 2, 2, 0, i, j, k, l, m) ||
			       is_h_ijklm(1, 1, 1, 1, 0, i, j, k, l, m) ||
			       is_h_ijklm(3, 3, 0, 0, 0, i, j, k, l, m) ||
			       is_h_ijklm(2, 2, 1, 1, 0, i, j, k, l, m) ||
			       is_h_ijklm(1, 1, 2, 2, 0, i, j, k, l, m) ||
			       is_h_ijklm(0, 0, 3, 3, 0, i, j, k, l, m) ||
			       is_h_ijklm(4, 4, 0, 0, 0, i, j, k, l, m) ||
			       is_h_ijklm(0, 0, 4, 4, 0, i, j, k, l, m) ||
			       is_h_ijklm(3, 3, 1, 1, 0, i, j, k, l, m) ||
			       is_h_ijklm(1, 1, 3, 3, 0, i, j, k, l, m) ||
			       is_h_ijklm(2, 2, 2, 2, 0, i, j, k, l, m)) {
		      // amplitude dependent tune shift
		      b[m1] = -scl_dnu*h_ijklm(K_re, i, j, k, l, m);
		    } else if (is_h_ijklm(1, 1, 0, 0, 2, i, j, k, l, m) ||
			       is_h_ijklm(0, 0, 1, 1, 2, i, j, k, l, m)) {
		      // 2nd order chromaticity
		      b[m1] = -scl_ksi2*h_ijklm(K_re, i, j, k, l, m);
		    } else if (is_h_ijklm(1, 1, 0, 0, 3, i, j, k, l, m) ||
			       is_h_ijklm(0, 0, 1, 1, 3, i, j, k, l, m) ||
			       is_h_ijklm(1, 1, 0, 0, 4, i, j, k, l, m) ||
			       is_h_ijklm(0, 0, 1, 1, 4, i, j, k, l, m) ||
			       is_h_ijklm(1, 1, 0, 0, 5, i, j, k, l, m) ||
			       is_h_ijklm(0, 0, 1, 1, 5, i, j, k, l, m) ||
			       is_h_ijklm(1, 1, 0, 0, 6, i, j, k, l, m) ||
			       is_h_ijklm(0, 0, 1, 1, 6, i, j, k, l, m)) {
		      // nonlinear chromaticity
		      b[m1] = -scl_ksi_nl*h_ijklm(K_re, i, j, k, l, m);
		    } else if (is_h_ijklm(2, 2, 0, 0, 1, i, j, k, l, m) ||
			       is_h_ijklm(0, 0, 2, 2, 1, i, j, k, l, m) ||
			       is_h_ijklm(1, 1, 1, 1, 1, i, j, k, l, m) ||
			       is_h_ijklm(2, 2, 0, 0, 2, i, j, k, l, m) ||
			       is_h_ijklm(0, 0, 2, 2, 2, i, j, k, l, m) ||
			       is_h_ijklm(1, 1, 1, 1, 2, i, j, k, l, m) ||
			       is_h_ijklm(2, 2, 0, 0, 3, i, j, k, l, m) ||
			       is_h_ijklm(0, 0, 2, 2, 3, i, j, k, l, m) ||
			       is_h_ijklm(1, 1, 1, 1, 3, i, j, k, l, m) ||
			       is_h_ijklm(2, 2, 0, 0, 4, i, j, k, l, m) ||
			       is_h_ijklm(0, 0, 2, 2, 4, i, j, k, l, m) ||
			       is_h_ijklm(1, 1, 1, 1, 4, i, j, k, l, m) ||
			       is_h_ijklm(3, 3, 0, 0, 1, i, j, k, l, m) ||
			       is_h_ijklm(2, 2, 1, 1, 1, i, j, k, l, m) ||
			       is_h_ijklm(1, 1, 2, 2, 1, i, j, k, l, m) ||
			       is_h_ijklm(0, 0, 3, 3, 1, i, j, k, l, m) ||
			       is_h_ijklm(3, 3, 0, 0, 2, i, j, k, l, m) ||
			       is_h_ijklm(2, 2, 1, 1, 2, i, j, k, l, m) ||
			       is_h_ijklm(1, 1, 2, 2, 2, i, j, k, l, m) ||
			       is_h_ijklm(0, 0, 3, 3, 2, i, j, k, l, m)) {
		      // amplitude dependent chromaticity
		      b[m1] = -scl_ksi_nl*h_ijklm(K_re, i, j, k, l, m);
		    } else if (!mirror_sym) {
		      sprintf(hs[m1-2], "i_%d%d%d%d%d", i, j, k, l, m);
		      b[m1-1] = -h_ijklm(g_im, i, j, k, l, m);
		      sprintf(hs[m1-1], "r_%d%d%d%d%d", i, j, k, l, m);
		      b[m1] = -h_ijklm(g_re, i, j, k, l, m);
		    }
		  }
		}
	      }

      // balance nonlinear terms
      A[++m1][i1] =
	scl_L*scl_dnudJ*3.0*h_ijklm_p(K_re, 3, 3, 0, 0, 0, 7)*dnudJ_twoJx;
      if (i1 == n_bn)
	b[m1] =
	  -scl_dnudJ*(h_ijklm(K_re, 2, 2, 0, 0, 0)
		      +3.0*h_ijklm(K_re, 3, 3, 0, 0, 0)*dnudJ_twoJx);
      sprintf(hs[m1-1], "dnuxdJx");

      A[++m1][i1] =
	scl_L*scl_dnudJ*2.0*h_ijklm_p(K_re, 2, 2, 1, 1, 0, 7)*dnudJ_twoJx;
      if (i1 == n_bn)
	b[m1] =
	  -scl_dnudJ*(h_ijklm(K_re, 1, 1, 1, 1, 0)
		      +2.0*h_ijklm(K_re, 2, 2, 1, 1, 0)*dnudJ_twoJx);
      sprintf(hs[m1-1], "dnuydJx");

     A[++m1][i1] =
	scl_L*scl_dnuddelta
	*(3.0*h_ijklm_p(K_re, 1, 1, 0, 0, 3, 7)*sqr(dnuddelta_delta)
	  +5.0*h_ijklm_p(K_re, 1, 1, 0, 0, 5, 7)*pow(dnuddelta_delta, 4));
      if (i1 == n_bn)
	b[m1] = -scl_dnuddelta
	*(h_ijklm(K_re, 1, 1, 0, 0, 1)
	  +3.0*h_ijklm(K_re, 1, 1, 0, 0, 3)*sqr(dnuddelta_delta)
	  +5.0*h_ijklm(K_re, 1, 1, 0, 0, 5)*pow(dnuddelta_delta, 4));
      sprintf(hs[m1-1], "dnuxdd ");

      A[++m1][i1] =
	scl_L*scl_dnuddelta
	*(3.0*h_ijklm_p(K_re, 0, 0, 1, 1, 3, 7)*sqr(dnuddelta_delta)
	  +5.0*h_ijklm_p(K_re, 0, 0, 1, 1, 5, 7)*pow(dnuddelta_delta, 4));
      if (i1 == n_bn)
	b[m1] = -scl_dnuddelta
	*(h_ijklm(K_re, 0, 0, 1, 1, 1)
	  +3.0*h_ijklm(K_re, 0, 0, 1, 1, 3)*sqr(dnuddelta_delta)
	  +5.0*h_ijklm(K_re, 0, 0, 1, 1, 5)*pow(dnuddelta_delta, 4));
      sprintf(hs[m1-1], "dnuydd ");

      A[++m1][i1] =
	scl_L*scl_dnudJ*2.0*h_ijklm_p(K_re, 1, 1, 2, 2, 0, 7)*dnudJ_twoJy;
      if (i1 == n_bn)
	b[m1] =
	  -scl_dnudJ*(h_ijklm(K_re, 1, 1, 1, 1, 0)
		      +2.0*h_ijklm(K_re, 1, 1, 2, 2, 0)*dnudJ_twoJy);
      sprintf(hs[m1-1], "dnuxdJy");

      A[++m1][i1] =
	scl_L*scl_dnudJ*3.0*h_ijklm_p(K_re, 0, 0, 3, 3, 0, 7)*dnudJ_twoJy;
      if (i1 == n_bn)
	b[m1] =
	  -scl_dnudJ*(h_ijklm(K_re, 0, 0, 2, 2, 0)
		      +3.0*h_ijklm(K_re, 0, 0, 3, 3, 0)*dnudJ_twoJy);
      sprintf(hs[m1-1], "dnuydJy");

      if (bns_fam[i1-1] > 0)
	clr_bn_par(bns_fam[i1-1], bns_n[i1-1]);
      else
	clr_s_par(abs(bns_fam[i1-1]));
    }

    // Exclude balancing of nonlinear terms
    m1 -= 6;
    // Exclude balancing of dnudJy
//    m1 -= 2;

    if (prt) {
      cout << endl;
      cout << n << " Ax = b:" << endl;
      cout << endl;
      for (i = 1; i <= m1; i++) {
	cout  << setw(3) << i << " " << hs[i-1];
	for (j = 1; j <= n_bn; j++)
	  cout << scientific << setprecision(2) << setw(10) << A[i][j];
	cout << scientific << setprecision(2) << setw(10) << b[i] << endl;
      }
    }
    
    cout << endl;
    cout << scientific << setprecision(3)
	 << "dnux/dJx:    "
	 << setw(10) << h_ijklm(K_re, 2, 2, 0, 0, 0)
	 << " + " << setw(10) << 3.0*h_ijklm(K_re, 3, 3, 0, 0, 0)*dnudJ_twoJx
	 << " + "
	 << setw(10) << 6.0*h_ijklm(K_re, 4, 4, 0, 0, 0)*sqr(dnudJ_twoJx)
	 << " = " << setw(10) << -b[m1-5]/scl_dnudJ << endl;
    cout << scientific << setprecision(3)
	 << "dnuy/dJx:    "
	 << setw(10) << h_ijklm(K_re, 1, 1, 1, 1, 0)
	 << " + " << setw(10) << 2.0*h_ijklm(K_re, 2, 2, 1, 1, 0)*dnudJ_twoJx
	 << " + "
	 << setw(10) << 3.0*h_ijklm(K_re, 3, 3, 1, 1, 0)*sqr(dnudJ_twoJx)
	 << " = " << setw(10) << -b[m1-4]/scl_dnudJ << endl;

    cout << scientific << setprecision(3)
	 << "dnux/ddelta: "
	 << setw(10) << h_ijklm(K_re, 1, 1, 0, 0, 1)
	 << " + "
	 << setw(10) << 3.0*h_ijklm(K_re, 1, 1, 0, 0, 3)*sqr(dnuddelta_delta)
	 << " + "
	 << setw(10)
	 << 5.0*h_ijklm(K_re, 1, 1, 0, 0, 5)*pow(dnuddelta_delta, 4)
	 << " = " << setw(10) << -b[m1-3]/scl_dnuddelta << endl;
    cout << scientific << setprecision(3)
	 << "dnuy/ddelta: "
	 << setw(10) << h_ijklm(K_re, 0, 0, 1, 1, 1)
	 << " + "
	 << setw(10) << 3.0*h_ijklm(K_re, 0, 0, 1, 1, 3)*sqr(dnuddelta_delta)
	 << " + "
	 << setw(10)
	 << 5.0*h_ijklm(K_re, 0, 0, 1, 1, 5)*pow(dnuddelta_delta, 4)
	 << " = " << setw(10) << -b[m1-2]/scl_dnuddelta << endl;

    cout << "excluded" << endl;

    cout << scientific << setprecision(3)
	 << "dnux/dJy:    "
	 << setw(10) << h_ijklm(K_re, 1, 1, 1, 1, 0)
	 << " + " << setw(10) << 2.0*h_ijklm(K_re, 1, 1, 2, 2, 0)*dnudJ_twoJy
	 << " + "
	 << setw(10) << 3.0*h_ijklm(K_re, 1, 1, 3, 3, 0)*sqr(dnudJ_twoJy)
	 << " = " << setw(10) << -b[m1-1]/scl_dnudJ << endl;
    cout << scientific << setprecision(3)
	 << "dnuy/dJy:    "
	 << setw(10) << h_ijklm(K_re, 0, 0, 2, 2, 0)
	 << " + " << setw(10) << 3.0*h_ijklm(K_re, 0, 0, 3, 3, 0)*dnudJ_twoJy
	 << " + "
	 << setw(10) << 6.0*h_ijklm(K_re, 0, 0, 4, 4, 0)*sqr(dnudJ_twoJy)
	 << " = " << setw(10) << -b[m1]/scl_dnudJ << endl;

    SVD_lim(m1, n_bn, A, b, bn_max, 1e-11, bn, dbn);

    cout << endl;
    cout << "dcorr. (int):" << endl;
    dbn_max = 0.0;
    for (i = 1; i <= n_bn; i++) {
      if (bns_fam[i-1] > 0) {
	L = get_L(bns_fam[i-1], 1);
	if (L == 0.0) L = 1.0;
	dbn_max = max(fabs(step*dbn[i]*L), dbn_max);
	cout << scientific << setprecision(3) << setw(11) << step*dbn[i]*L;
      } else {
	dbn_max = max(fabs(step*scl_ds*dbn[i]), dbn_max);
	cout << scientific << setprecision(3)
	     << setw(11) << step*scl_ds*dbn[i];
      }
      if (i % n_prt == 0) cout << endl;
    }
    if (n_bn % n_prt != 0) cout << endl;

    cout << endl;
    cout << "corr.:" << endl;
    for (i = 1; i <= n_bn; i++) {
      set_dbn_s(bns_fam[i-1], bns_n[i-1], step*dbn[i]);
      bn[i] = get_bn_s(bns_fam[i-1], 1, bns_n[i-1]);
      cout << scientific << setprecision(3)
	   << setw(11) << get_bnL_s(bns_fam[i-1], 1, bns_n[i-1]);
      if (i % n_prt == 0) cout << endl;
    }
    if (n_bn % n_prt != 0) cout << endl;

    if (prt_iter) {
      sext_out << endl;
      sext_out << "n = " << n << ":" << endl;
      for (i = 1; i <= n_bn; i++)
	for (j = 1; j <= get_n_Kids(abs(bns_fam[i-1])); j++) {
	  if (bns_fam[i-1] > 0) 
	    sext_out << fixed << setprecision(7) 
		     << setw(9) << get_Name(abs(bns_fam[i-1]))
		     << "(" << j << ") = "
		     << setw(11) << get_bnL_s(bns_fam[i-1], 1, bns_n[i-1])
		     << setw(2) << bns_n[i-1] << endl;
	  else {
	    sprintf(str, "du_%s", get_Name(abs(bns_fam[i-1])));
	    sext_out << fixed << setprecision(7) 
		     << setw(9) << str << "(" << j << ") = "
		     << setw(11) << get_bnL_s(bns_fam[i-1], 1, bns_n[i-1])
		     << setw(2) << 0 << endl;
	    sprintf(str, "dd_%s", get_Name(abs(bns_fam[i-1])));
	    sext_out << fixed << setprecision(7) 
		      << setw(9) << str << "(" << j << ") = "
		     << setw(11) << -get_bnL_s(bns_fam[i-1], 1, bns_n[i-1])
		     << setw(2) << 0 << endl;
	  }
	}
      sext_out.flush();
    }

    cout << endl;
    cout << scientific << setprecision(1)
	 << setw(2) << n << " chi2: " << chi2;

    chi2 = 0.0;
    for (i = 1; i <= m1; i++)
      chi2 += sqr(b[i]);

    cout << scientific << setprecision(1)
	 << " -> " << chi2 << endl;
  } while ((dbn_max > bn_tol) && (n < n_max));

  if (!prt_iter) {
    for (i = 1; i <= n_bn; i++)
      for (j = 1; j <= get_n_Kids(abs(bns_fam[i-1])); j++) {
	sext_out << fixed << setprecision(7) 
		 << setw(6) << get_Name(abs(bns_fam[i-1]))
		 << "(" << j << ") = "
		 << setw(11) << get_bnL_s(bns_fam[i-1], 1, bns_n[i-1])
		 << setw(2) << bns_n[i-1] << endl;
      }
    sext_out.flush();
  }

  free_dvector(b, 1, m_max); free_dvector(w, 1, n_bn);
  free_dvector(dbn, 1, n_bn); free_dvector(bn, 1, n_bn);
  free_dvector(bn_max, 1, n_bn);
  free_dmatrix(A, 1, m_max, 1, n_bn);
  free_dmatrix(U, 1, m_max, 1, n_bn); free_dmatrix(V, 1, n_bn, 1, n_bn);
  free_dmatrix(A_inv, 1, n_bn, 1, m_max);
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


void get_prm(char *file_name)
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
  sscanf(line, "%*s %d", &check_range);

  prm_in.getline(line, max_str);
  sscanf(line, "%*s %lf", &nu_x_min);

  prm_in.getline(line, max_str);
  sscanf(line, "%*s %lf", &nu_x_max);

  prm_in.getline(line, max_str);
  sscanf(line, "%*s %lf", &nu_y_min);

  prm_in.getline(line, max_str);
  sscanf(line, "%*s %lf", &nu_y_max);

  prm_in.getline(line, max_str);
  sscanf(line, "%*s %d", &n_steps);

  prm_in.getline(line, max_str);
  sscanf(line, "%*s %d", &n_cell);

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
  sscanf(line, "%*s %lf", &scl_dnu);

  prm_in.getline(line, max_str);
  sscanf(line, "%*s %lf", &scl_ksi1);

  prm_in.getline(line, max_str);
  sscanf(line, "%*s %lf", &scl_ksi2);

  prm_in.getline(line, max_str);
  sscanf(line, "%*s %lf", &scl_ksi_nl);

  prm_in.getline(line, max_str);
  sscanf(line, "%*s %lf", &scl_dnuddelta);

  prm_in.getline(line, max_str);
  sscanf(line, "%*s %lf", &scl_dnudJ);

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
       << scientific << setprecision(1) << ", eps_ksi = " << eps_ksi << endl;
  cout << "check_range   = " << check_range << endl;
  cout << endl;
  cout << fixed << setprecision(6)
       << "n_steps   = " << n_steps
       << ", nu_x_min = " << nu_x_min << ", nu_x_max = " << nu_x_max
       << ", nu_y_min = " << nu_y_min << ", nu_y_max = " << nu_y_max << endl;
  cout << endl;
  cout << "n_cell        = " << n_cell << endl;
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
       << "scl_dnu       = " << scl_dnu << endl;
  cout << fixed << setprecision(1)
       << "scl_ksi1      = " << scl_ksi1 << endl;
  cout << fixed << setprecision(1)
       << "scl_ksi2      = " << scl_ksi2 << endl;
  cout << fixed << setprecision(1)
       << "scl_ksi_nl    = " << scl_ksi_nl << endl;
  cout << scientific << setprecision(1)
       << "scl_dnuddelta = " << scl_dnuddelta << endl;
  cout << scientific << setprecision(1)
       << "scl_dnuddJ    = " << scl_dnudJ << endl;
  cout << fixed << setprecision(2)
       << "step          = " << step << endl;
}


void get_bns()
{

  n_bn = 0;
  bns_fam[n_bn] = get_Fnum("sl1");  bns_n[n_bn++] = Sext;
  bns_fam[n_bn] = get_Fnum("sl2");  bns_n[n_bn++] = Sext;
  bns_fam[n_bn] = get_Fnum("sl3");  bns_n[n_bn++] = Sext;

  bns_fam[n_bn] = get_Fnum("sm1a"); bns_n[n_bn++] = Sext;
  bns_fam[n_bn] = get_Fnum("sm1b"); bns_n[n_bn++] = Sext;
  bns_fam[n_bn] = get_Fnum("sm2");  bns_n[n_bn++] = Sext;

  bns_fam[n_bn] = get_Fnum("sh1");  bns_n[n_bn++] = Sext;
  bns_fam[n_bn] = get_Fnum("sh3");  bns_n[n_bn++] = Sext;
  bns_fam[n_bn] = get_Fnum("sh4");  bns_n[n_bn++] = Sext;

  if (n_bn > n_prm_max) {
    cout << "get_bns: n_prm_max exceeded " << n_bn << "(" << n_prm_max
	 << ")" << endl;
    exit(0);
  }

  cout << endl;
  cout << "get_bns: no of multipole families " << n_bn << endl;
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


void get_locs()
{
  double  alpha1[2], alpha2[2], alpha3[2];

  beta_loc1 = get_loc(get_Fnum("mp"), 1); get_ab(alpha1, beta1, beta_loc1);
  beta_loc2 = get_loc(get_Fnum("ss"), 1); get_ab(alpha2, beta2, beta_loc2);
//  beta2[X_] = 1.0; beta2[Y_] = 1.0;
  beta_loc3 = get_loc(get_Fnum("ls"), 1); get_ab(alpha3, beta3, beta_loc3);
//  beta3[X_] = 15.0; beta3[Y_] = 3.0;

  cout << endl;
  cout << fixed << setprecision(3)
       << "alpha1_x  = " << setw(6) << alpha1[X_]
       << ", alpha1_y = " << setw(6) << alpha1[Y_]
       << ", beta1_x = " << setw(6) << beta1[X_]
       << ", beta1_y  = " << setw(6) << beta1[Y_]
       << endl;
  cout << fixed << setprecision(3)
       << "alpha2_x  = " << setw(6) << alpha2[X_]
       << ", alpha2_y = " << setw(6) << alpha2[Y_]
       << ", beta2_x = " << setw(6) << beta2[X_]
	   << ", beta2_y  = " << setw(6) << beta2[Y_]
       << endl;
  cout << fixed << setprecision(3)
       << "alpha3_x  = " << setw(6) << alpha3[X_]
       << ", alpha3_y = " << setw(6) << alpha3[Y_]
       << ", beta3_x = " << setw(6) << beta3[X_]
       << ", beta3_y  = " << setw(6) << beta3[Y_]
       << endl;
}


int main()
{
  double           nu[2], ksi[2];
  double           alpha[2], beta[2];
  tps              Hnl, H2, gn, h, h_re, h_im, H_num, dH, H, H_re, H_im;
  tps              K_re, K_im, g_re, g_im;
  ss_vect<tps>     nus;
  ofstream         outf, K_out, nus_out, A1_out, J_out;
  ifstream         inf;

  danot_(no_tps-1);

  rad_on    = false; H_exact        = false; totpath_on   = false;
  cavity_on = false; quad_fringe_on = false; emittance_on = false;
  IBS_on    = false;


  rd_mfile("flat_file.dat", elem); rd_mfile("flat_file.dat", elem_tps);

  // initialize the symplectic integrator after energy has been defined
  ini_si();

  // disable from TPSALib- and LieLib log messages
  idprset(-1);

  danot_(3);

  if (true) chk_lat(nu, ksi);

  get_ab(alpha, beta, 0);
  cout << endl;
  cout << fixed << setprecision(3)
       << "alpha_x  = " << setw(6) << alpha[X_]
       << ", alpha_y = " << setw(6) << alpha[Y_]
       << ", beta_x = " << setw(6) << beta[X_]
       << ", beta_y  = " << setw(6) << beta[Y_] << endl;

  Jx = sqr(max_Ax)/(2.0*beta1[X_]); Jy = sqr(max_Ay)/(2.0*beta1[Y_]);
  delta = max_delta;

  Id_scl.identity();
  Id_scl[x_] *= sqrt(2.0*Jx); Id_scl[px_] *= sqrt(2.0*Jx);
  Id_scl[y_] *= sqrt(2.0*Jy); Id_scl[py_] *= sqrt(2.0*Jy);
  Id_scl[delta_] *= delta;

  get_prm("H_zero.prm");

  if (adj_chrom) {
    switch (no_tps) {
    case 4: file_wr(sext_out, "sext_1.dat");
      break;
    case 5:
    case 6: file_wr(sext_out, "sext_2.dat");
      break;
    case 7:
    case 8: file_wr(sext_out, "sext_3.dat"); 
      break;
    case 9: file_wr(sext_out, "sext_4.dat");
      break;
    default: ;
    }
  }
 
 if (adj_chrom) {
    n_bn = 0;
    bns_fam[n_bn++] = get_Fnum("sm1a"); bns_fam[n_bn++] = get_Fnum("sm1b");
    bns_fam[n_bn++] = get_Fnum("sm2");

    if (fit_chrm) {
      danot_(3);
      no_mpoles();
      fit_chrom(ksi1[X_], ksi1[Y_], n_bn, bns_fam, true);
    }

    get_bns(); H_zero(eps_ksi, 150, true);
  }

  danot_(no_tps-1);
  get_Map();
  danot_(no_tps);
  K = MapNorm(Map, g, A1, A0, Map_res, no_tps); nus = dHdJ(K);
  CtoR(K*Id_scl, K_re, K_im); CtoR(get_h()*Id_scl, h_re, h_im);

  file_wr(outf, "map.dat"); outf << Map; outf.close();

  if (false) Id_scl.identity();

  file_wr(outf, "h.dat"); outf << h_re; outf.close();

  file_wr(outf, "K.dat"); outf << K_re; outf.close();

  file_wr(outf, "nus.dat"); outf << nus[3] << nus[4]; outf.close();

  CtoR(get_H()*Id_scl, H_re, H_im);
  file_wr(outf, "H.dat"); outf << H_re; outf.close();

  CtoR(g*Id_scl, g_re, g_im);
  file_wr(outf, "g.dat"); outf << g_im; outf.close();
}
