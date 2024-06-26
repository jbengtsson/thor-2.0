#define NO 5

#include "thor_lib.h"

int no_tps   = NO,
    ndpt_tps = 5;


/* SLS:

   I [A] = | (b_3L)* 13.72 |,    pos = { sf, s*b }, net = { se, sd, s*a }

*/


extern double        b2_max;
extern tps           K, g;
extern ss_vect<tps>  Map, A0, A1, Map_res;

const int
  max_n_femto = 7,
  n_track     = 512,
  n_aper      = 25,
  max_ind     = 10,
  n_prm_max   = 20;

const bool  tune_scn = false;


bool          h[max_ind][max_ind][max_ind][max_ind][max_ind], fit_chrm;
int           n_prm, prms[n_prm_max], bns[n_prm_max], n_cell, n_iter;
int           check_range, adj_tune, adj_chrom, n_steps;
int           n_b2_femto, b2_femto[max_n_femto], k1_femto, k2_femto;
long int      beta_loc1, beta_loc2, beta_loc3, beta_loc4, rseed;
double        Jx, Jy, delta, beta1[2], beta2[2], beta3[2], beta4[2];
double        alpha_femto[2], beta_femto[2], nu_femto[2];
double        b3s[n_prm_max];
double        nu0[2], eps_nu, ksi1[2], eps_ksi;
double        nu_x_min, nu_x_max, nu_y_min, nu_y_max;
double        bnL_max[mpole_max];
double        scl_dnu, scl_ksi1, scl_ksi2, scl_ksi_nl;
double        scl_dnuddelta, scl_dnudJ, step;
ss_vect<tps>  Id_scl;
std::ofstream      quad_out, sext_out, oct_out, tune_out;


const char    lattice[]  = "NSLS-II_DBA";
//const char    lattice[]  = "RHIC";
//const char    lattice[]  = "ALBA_SCL";
//const char    lattice[]  = "NSLS-II_BOOSTER";
//const char    lattice[]  = "SLS";
//const char    lattice[]  = "DIAMOND";
//const char    lattice[]  = "ESRF";
//const char    lattice[]  = "SOLEIL";
//const char    lattice[]  = "ALBA";
//const char    lattice[]  = "MAX-IV";


const bool    sext_prm = true, oct_prm = false, dodec_prm = false;

const double  max_Ax = 20e-3, max_Ay = 10e-3, max_delta = 3e-2;
// SLS
//const double  max_Ax = 20e-3, max_Ay = 10e-3, max_delta = 2.5e-2;


void ini_ranf(const int i) { rseed = i; }


double ranf(void)
{

  const int       k = 19;
  const long int  c = 656329, m = 100000001;

  rseed = (k*rseed+c) % m;

  return rseed/1e8;
}


void get_stats(long n, double sum, double sum2, double &m, double &s)
{
  if (n > 0)
    m = sum/n;
  else
    printf("average not defined: n = %ld\n", n);
  if (n < 2) {
    printf("sigma not defined: n = %ld\n", n);
    s = 0.0;
  } else {
    s = (sum2-n*sqr(m))/(n-1);
    if (s >= 0.0)
      s = sqrt(s);
    else
      printf("sqrt of neg. value\n");
  }
}


ss_vect<tps> get_A_nl(const tps g)
{
  int           j;
  tps           gn;
  ss_vect<tps>  Id, A_nl;

  Id.identity();  A_nl = Id;
  for (j = 3; j <= no_tps; j++) {
    gn = Take(g, j); A_nl = A_nl*LieExp(gn, Id);
  }

  //  return FExpo(g, Id, 3, no_tps, -1);

  return A_nl;
}


ss_vect<tps> get_A_nl_inv(const tps g)
{
  int           j;
  tps           gn;
  ss_vect<tps>  Id, A_nl_inv;

  Id.identity();  A_nl_inv = Id;
  for (j = no_tps; j >= 3; j--) {
    gn = Take(g, j); A_nl_inv = A_nl_inv*LieExp(-gn, Id);
  }

  //  return FExpo(-g, Id, 3, no_tps, 1);

  return A_nl_inv;
}


void prt_h(const tps &H, const double J_max, const double delta)
{
  int           i, j;
  ss_vect<tps>  ps;
  std::ofstream      of;
  
  file_wr(of, "h.dat");
  of << "#  J_x    J_y    ||h||" << std::endl;
  ps.zero(); ps[delta_] = delta;
  for (i = -n_steps; i <= n_steps; i++) {
    ps[x_] = i*sqrt(2.0*J_max)/n_steps;
    for (j = -n_steps; j <= n_steps; j++) {
      ps[y_] = j*sqrt(2.0*J_max)/n_steps;
      of << std::scientific << std::setprecision(2)
	 << std::setw(10) << ps[x_].cst() << std::setw(10) << ps[y_].cst()
	//	 << std::setw(10) << (H*ps).cst() << std::endl;
	 << std::setw(10) << ((K-Take(K, 2))*ps).cst() << std::endl;
    }
    of << std::endl;
  }
  of.close();
}


void get_Map_N(const int n)
{
  ss_vect<tps>  Map2, Map4;

  get_Map();

  switch (n) {
  case 1:
    break;
  case 2:
    Map = Map*Map;
    break;
  case 3:
    Map2 = Map*Map; Map = Map2*Map;
    break;
  case 4:
    Map2 = Map*Map; Map = Map2*Map2;
    break;
  case 5:
    Map2 = Map*Map; Map = Map2*Map2*Map;
    break;
  case 6:
    Map2 = Map*Map; Map = Map2*Map2*Map2;
    break;
  case 7:
    Map2 = Map*Map; Map4 = Map2*Map2; Map = Map4*Map2*Map;
    break;
  case 8:
    Map2 = Map*Map; Map4 = Map2*Map2; Map = Map4*Map4;
    break;
  case 9:
    Map2 = Map*Map; Map4 = Map2*Map2; Map = Map4*Map4*Map;
    break;
  case 10:
    Map2 = Map*Map; Map4 = Map2*Map2; Map = Map4*Map4*Map2;
    break;
  case 11:
    Map2 = Map*Map; Map4 = Map2*Map2; Map = Map4*Map4*Map2*Map;
    break;
  case 12:
    Map2 = Map*Map; Map4 = Map2*Map2; Map = Map4*Map4*Map4;
    break;
  case 13:
    Map2 = Map*Map; Map4 = Map2*Map2; Map = Map4*Map4*Map4*Map;
    break;
  case 14:
    Map2 = Map*Map; Map4 = Map2*Map2; Map = Map4*Map4*Map4*Map2;
    break;
  case 15:
    Map2 = Map*Map; Map4 = Map2*Map2; Map = Map4*Map4*Map4*Map2*Map;
    break;
  default:
    std::cout << "get_Map_N: n not defined " << n << std::endl;
    exit(1);
    break;
  }
}


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


void H_zero(const double prm_tol, const int n_max, const bool prt_iter)
{
  const int     m_max = 500;  // max no of constraints

  bool          first = true;
  char          hs[m_max][max_str], str[max_str];
  int           i, j, k, l, m, i1, m1 = 0, n;
  double        L, dprm_max, chi2 = 0.0, scl_L = 0.0;
  double        **A, *b, *w, **U, **V, *dbn, *bn, *bn_max, **A_inv;
  tps           r, K_re, K_im, K_re_scl, g_re, g_im;
  ss_vect<tps>  Id_Fl, A1_inv, A_nl_inv;

  const bool    prt             = true, mirror_sym = false;
  const int     n_prt           = 9;
  const double  dnudJ_twoJx     = sqr(20e-3)/beta1[X_];
  const double  dnudJ_twoJy     = sqr(7.5e-3)/beta1[Y_];
  const double  dnuddelta_delta = 2.5e-2;

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
	std::cout << std::endl;
	std::cout << get_Name(prms[i-1]) << ": bnL_max set to 1.0" << std::endl;
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
    std::cout << std::endl;
    std::cout << "initial b3s:" << std::endl;
    for (i = 1; i <= n_prm; i++) {
      b3s[i-1] = get_bnL_s(prms[i-1], 1, bns[i-1]);
      std::cout << std::scientific << std::setprecision(3)
	   << std::setw(11) << b3s[i-1] << std::setw(2) << bns[i-1] << std::endl;
    }

    select_h();
  } else if (true) {
    // initialize sextupoles
    for (i = 1; i <= n_prm; i++)
      set_bnL_s(prms[i-1], bns[i-1], b3s[i-1]);
  }

  n = 0;
  do {
    n++;
    std::cout << std::endl;
    for (i1 = 1; i1 <= n_prm; i1++) {
      if (prms[i1-1] > 0)
	set_bn_par(prms[i1-1], bns[i1-1], 7);
      else
	set_s_par(abs(prms[i1-1]), 7);

      //      Jx = 0.5; Jy = 0.5; delta = 1.0;

      danot_(no_tps-1); get_Map_N(n_cell);
      danot_(no_tps); K = MapNorm(Map, g, A1, A0, Map_res, 1);
      CtoR(K, K_re, K_im); K_re_scl = K_re*Id_scl;

      g = g*Id_scl; CtoR(g, g_re, g_im);

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
		     (fabs(h_ijklm(K_re_scl, i, j, k, l, m)) > 0.0))) {

		  scl_L = (prms[i1-1] > 0)? 1.0 : scl_ds;

		  A[++m1][i1] = scl_L*h_ijklm_p(g_im, i, j, k, l, m, 7);

		  if (m1 >= m_max) {
		    std::cout << "m_max reached " << m1 << "(" << m_max << ")"
			 << std::endl;
		    exit(1);
		  }

		  if (is_h_ijklm(1, 1, 0, 0, 1, i, j, k, l, m))
		    // horizontal linear chromaticity
		    A[m1][i1] =
		      scl_L*scl_ksi1
		      *h_ijklm_p(K_re_scl, i, j, k, l, m, 7);
		  else if (is_h_ijklm(0, 0, 1, 1, 1, i, j, k, l, m))
		    // vertical linear chromaticity
		    A[m1][i1] =
		      scl_L*scl_ksi1
		      *h_ijklm_p(K_re_scl, i, j, k, l, m, 7);
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
		      scl_L*scl_dnu
		      *h_ijklm_p(K_re_scl, i, j, k, l, m, 7);
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
		      scl_L*scl_ksi_nl
		      *h_ijklm_p(K_re_scl, i, j, k, l, m, 7);
		  else if (is_h_ijklm(2, 2, 0, 0, 1, i, j, k, l, m) ||
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
			   is_h_ijklm(0, 0, 3, 3, 2, i, j, k, l, m))
		    // amplitude dependent chromaticity
		    A[m1][i1] =
		      scl_L*scl_ksi_nl
		      *h_ijklm_p(K_re_scl, i, j, k, l, m, 7);
		  else if (!mirror_sym)
		    A[++m1][i1] = scl_L
		      *h_ijklm_p(g_re, i, j, k, l, m, 7);

		  if (i1 == n_prm) {
                    sprintf(hs[m1-1], "h_%d%d%d%d%d", i, j, k, l, m);
                    b[m1] = -h_ijklm(g_im, i, j, k, l, m);

		    if (is_h_ijklm(1, 1, 0, 0, 1, i, j, k, l, m)) {
		      // horizontal linear chromaticity
		      b[m1] =
			-scl_ksi1
			*(n_cell*ksi1[X_]*M_PI*2.0*Jx*delta
			  +h_ijklm(K_re_scl, i, j, k, l, m));
		    } else if (is_h_ijklm(0, 0, 1, 1, 1, i, j, k, l, m)) {
		      // vertical linear chromaticity
		      b[m1] =
			-scl_ksi1
			*(n_cell*ksi1[Y_]*M_PI*2.0*Jy*delta
			  +h_ijklm(K_re_scl, i, j, k, l, m));
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
		      b[m1] = -scl_dnu*h_ijklm(K_re_scl, i, j, k, l, m);
		    } else if (is_h_ijklm(1, 1, 0, 0, 2, i, j, k, l, m) ||
			       is_h_ijklm(0, 0, 1, 1, 2, i, j, k, l, m) ||
			       is_h_ijklm(1, 1, 0, 0, 3, i, j, k, l, m) ||
			       is_h_ijklm(0, 0, 1, 1, 3, i, j, k, l, m) ||
			       is_h_ijklm(1, 1, 0, 0, 4, i, j, k, l, m) ||
			       is_h_ijklm(0, 0, 1, 1, 4, i, j, k, l, m) ||
			       is_h_ijklm(1, 1, 0, 0, 5, i, j, k, l, m) ||
			       is_h_ijklm(0, 0, 1, 1, 5, i, j, k, l, m) ||
			       is_h_ijklm(1, 1, 0, 0, 6, i, j, k, l, m) ||
			       is_h_ijklm(0, 0, 1, 1, 6, i, j, k, l, m)) {
		      // nonlinear chromaticity
		      b[m1] =
			-scl_ksi_nl*h_ijklm(K_re_scl, i, j, k, l, m);
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
		      b[m1] =
			-scl_ksi_nl*h_ijklm(K_re_scl, i, j, k, l, m);
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
      if (i1 == n_prm)
	b[m1] =
	  -scl_dnudJ
	  *(h_ijklm(K_re, 2, 2, 0, 0, 0)
	    +3.0*h_ijklm(K_re, 3, 3, 0, 0, 0)*dnudJ_twoJx);
      sprintf(hs[m1-1], "dnuxdJx");

      A[++m1][i1] =
	scl_L*scl_dnudJ*2.0*h_ijklm_p(K_re, 2, 2, 1, 1, 0, 7)*dnudJ_twoJx;
      if (i1 == n_prm)
	b[m1] =
	  -scl_dnudJ
	  *(h_ijklm(K_re, 1, 1, 1, 1, 0)
	    +2.0*h_ijklm(K_re, 2, 2, 1, 1, 0)*dnudJ_twoJx);
      sprintf(hs[m1-1], "dnuydJx");

      A[++m1][i1] =
	scl_L*scl_dnuddelta
	*(3.0*h_ijklm_p(K_re, 1, 1, 0, 0, 3, 7)*sqr(dnuddelta_delta)
	  +5.0*h_ijklm_p(K_re, 1, 1, 0, 0, 5, 7)*pow(dnuddelta_delta, 4));
      if (i1 == n_prm)
	b[m1] = -scl_dnuddelta
	  *(h_ijklm(K_re, 1, 1, 0, 0, 1)
	    +3.0*h_ijklm(K_re, 1, 1, 0, 0, 3)*sqr(dnuddelta_delta)
	    +5.0*h_ijklm(K_re, 1, 1, 0, 0, 5)*pow(dnuddelta_delta, 4));
      sprintf(hs[m1-1], "dnuxdd ");

      A[++m1][i1] =
	scl_L*scl_dnuddelta
	*(3.0*h_ijklm_p(K_re, 0, 0, 1, 1, 3, 7)*sqr(dnuddelta_delta)
	  +5.0*h_ijklm_p(K_re, 0, 0, 1, 1, 5, 7)*pow(dnuddelta_delta, 4));
      if (i1 == n_prm)
	b[m1] = -scl_dnuddelta
	  *(h_ijklm(K_re, 0, 0, 1, 1, 1)
	    +3.0*h_ijklm(K_re, 0, 0, 1, 1, 3)*sqr(dnuddelta_delta)
	    +5.0*h_ijklm(K_re, 0, 0, 1, 1, 5)*pow(dnuddelta_delta, 4));
      sprintf(hs[m1-1], "dnuydd ");

      A[++m1][i1] =
	scl_L*scl_dnudJ*2.0*h_ijklm_p(K_re, 1, 1, 2, 2, 0, 7)*dnudJ_twoJy;
      if (i1 == n_prm)
	b[m1] =
	  -scl_dnudJ
	  *(h_ijklm(K_re, 1, 1, 1, 1, 0)
	    +2.0*h_ijklm(K_re, 1, 1, 2, 2, 0)*dnudJ_twoJy);
      sprintf(hs[m1-1], "dnuxdJy");

      A[++m1][i1] =
	scl_L*scl_dnudJ*3.0*h_ijklm_p(K_re, 0, 0, 3, 3, 0, 7)*dnudJ_twoJy;
      if (i1 == n_prm)
	b[m1] =
	  -scl_dnudJ
	  *(h_ijklm(K_re, 0, 0, 2, 2, 0)
	    +3.0*h_ijklm(K_re, 0, 0, 3, 3, 0)*dnudJ_twoJy);
      sprintf(hs[m1-1], "dnuydJy");

      if (prms[i1-1] > 0)
	clr_bn_par(prms[i1-1], bns[i1-1]);
      else
	clr_s_par(abs(prms[i1-1]));
    }

    // Exclude balancing of nonlinear terms
    //    m1 -= 6;
    // Exclude balancing of dnudJy
    m1 -= 2;

    if (prt) {
      std::cout << std::endl;
      std::cout << n << " Ax = b:" << std::endl;
      std::cout << std::endl;
      for (i = 1; i <= m1; i++) {
	std::cout  << std::setw(3) << i << " " << hs[i-1];
	for (j = 1; j <= n_prm; j++)
	  std::cout << std::scientific << std::setprecision(2) << std::setw(10) << A[i][j];
	std::cout << std::scientific << std::setprecision(2) << std::setw(10) << b[i] << std::endl;
      }
    }
    
    std::cout << std::endl;
    std::cout << std::scientific << std::setprecision(3)
	 << "dnux/dJx:    "
	 << std::setw(10) << h_ijklm(K_re, 2, 2, 0, 0, 0)
	 << " + " << std::setw(10) << 3.0*h_ijklm(K_re, 3, 3, 0, 0, 0)*dnudJ_twoJx
	 << " + "
	 << std::setw(10) << 6.0*h_ijklm(K_re, 4, 4, 0, 0, 0)*sqr(dnudJ_twoJx)
	 << " = " << std::setw(10) << -b[m1-5]/scl_dnudJ << std::endl;
    std::cout << std::scientific << std::setprecision(3)
	 << "dnuy/dJx:    "
	 << std::setw(10) << h_ijklm(K_re, 1, 1, 1, 1, 0)
	 << " + " << std::setw(10) << 2.0*h_ijklm(K_re, 2, 2, 1, 1, 0)*dnudJ_twoJx
	 << " + "
	 << std::setw(10) << 3.0*h_ijklm(K_re, 3, 3, 1, 1, 0)*sqr(dnudJ_twoJx)
	 << " = " << std::setw(10) << -b[m1-4]/scl_dnudJ << std::endl;

    std::cout << std::scientific << std::setprecision(3)
	 << "dnux/ddelta: "
	 << std::setw(10) << h_ijklm(K_re, 1, 1, 0, 0, 1)
	 << " + "
	 << std::setw(10) << 3.0*h_ijklm(K_re, 1, 1, 0, 0, 3)*sqr(dnuddelta_delta)
	 << " + "
	 << std::setw(10)
	 << 5.0*h_ijklm(K_re, 1, 1, 0, 0, 5)*pow(dnuddelta_delta, 4)
	 << " = " << std::setw(10) << -b[m1-3]/scl_dnuddelta << std::endl;
    std::cout << std::scientific << std::setprecision(3)
	 << "dnuy/ddelta: "
	 << std::setw(10) << h_ijklm(K_re, 0, 0, 1, 1, 1)
	 << " + "
	 << std::setw(10) << 3.0*h_ijklm(K_re, 0, 0, 1, 1, 3)*sqr(dnuddelta_delta)
	 << " + "
	 << std::setw(10)
	 << 5.0*h_ijklm(K_re, 0, 0, 1, 1, 5)*pow(dnuddelta_delta, 4)
	 << " = " << std::setw(10) << -b[m1-2]/scl_dnuddelta << std::endl;

    std::cout << "excluded" << std::endl;

    std::cout << std::scientific << std::setprecision(3)
	 << "dnux/dJy:    "
	 << std::setw(10) << h_ijklm(K_re, 1, 1, 1, 1, 0)
	 << " + " << std::setw(10) << 2.0*h_ijklm(K_re, 1, 1, 2, 2, 0)*dnudJ_twoJy
	 << " + "
	 << std::setw(10) << 3.0*h_ijklm(K_re, 1, 1, 3, 3, 0)*sqr(dnudJ_twoJy)
	 << " = " << std::setw(10) << -b[m1-1]/scl_dnudJ << std::endl;
    std::cout << std::scientific << std::setprecision(3)
	 << "dnuy/dJy:    "
	 << std::setw(10) << h_ijklm(K_re, 0, 0, 2, 2, 0)
	 << " + " << std::setw(10) << 3.0*h_ijklm(K_re, 0, 0, 3, 3, 0)*dnudJ_twoJy
	 << " + "
	 << std::setw(10) << 6.0*h_ijklm(K_re, 0, 0, 4, 4, 0)*sqr(dnudJ_twoJy)
	 << " = " << std::setw(10) << -b[m1]/scl_dnudJ << std::endl;

    SVD_lim(m1, n_prm, A, b, bn_max, 1e-11, bn, dbn);

    std::cout << std::endl;
    std::cout << "dcorr. (int):" << std::endl;
    dprm_max = 0.0;
    for (i = 1; i <= n_prm; i++) {
      if (prms[i-1] > 0) {
	L = get_L(prms[i-1], 1);
	if (L == 0.0) L = 1.0;
	dprm_max = max(fabs(step*dbn[i]*L), dprm_max);
	std::cout << std::scientific << std::setprecision(3) << std::setw(11) << step*dbn[i]*L;
      } else {
	dprm_max = max(fabs(step*scl_ds*dbn[i]), dprm_max);
	std::cout << std::scientific << std::setprecision(3)
	     << std::setw(11) << step*scl_ds*dbn[i];
      }
      if (i % n_prt == 0) std::cout << std::endl;
    }
    if (n_prm % n_prt != 0) std::cout << std::endl;

    std::cout << std::endl;
    std::cout << "corr.:" << std::endl;
    for (i = 1; i <= n_prm; i++) {
      set_dbn_s(prms[i-1], bns[i-1], step*dbn[i]);
      bn[i] = get_bn_s(prms[i-1], 1, bns[i-1]);
      std::cout << std::scientific << std::setprecision(3)
	   << std::setw(11) << get_bnL_s(prms[i-1], 1, bns[i-1]);
      if (i % n_prt == 0) std::cout << std::endl;
    }
    if (n_prm % n_prt != 0) std::cout << std::endl;

    if (prt_iter) {
      sext_out << std::endl;
      sext_out << "n = " << n << ":" << std::endl;
      for (i = 1; i <= n_prm; i++)
	for (j = 1; j <= get_n_Kids(abs(prms[i-1])); j++) {
	  if (prms[i-1] > 0) 
	    sext_out << std::fixed << std::setprecision(7) 
		     << std::setw(9) << get_Name(abs(prms[i-1]))
		     << "(" << j << ") = "
		     << std::setw(11) << get_bnL_s(prms[i-1], 1, bns[i-1])
		     << std::setw(2) << bns[i-1] << std::endl;
	  else {
	    sprintf(str, "du_%s", get_Name(abs(prms[i-1])));
	    sext_out << std::fixed << std::setprecision(7) 
		     << std::setw(9) << str << "(" << j << ") = "
		     << std::setw(11) << get_bnL_s(prms[i-1], 1, bns[i-1])
		     << std::setw(2) << 0 << std::endl;
	    sprintf(str, "dd_%s", get_Name(abs(prms[i-1])));
	    sext_out << std::fixed << std::setprecision(7) 
		     << std::setw(9) << str << "(" << j << ") = "
		     << std::setw(11) << -get_bnL_s(prms[i-1], 1, bns[i-1])
		     << std::setw(2) << 0 << std::endl;
	  }
	}
      sext_out.flush();
    }

    std::cout << std::endl;
    std::cout << std::scientific << std::setprecision(1)
	 << std::setw(2) << n << " chi2: " << chi2;

    chi2 = 0.0;
    for (i = 1; i <= m1; i++)
      chi2 += sqr(b[i]);

    std::cout << std::scientific << std::setprecision(1)
	 << " -> " << chi2 << std::endl;
  } while ((dprm_max > prm_tol) && (n < n_max));

  if (!prt_iter) {
    for (i = 1; i <= n_prm; i++)
      for (j = 1; j <= get_n_Kids(abs(prms[i-1])); j++) {
	sext_out << std::fixed << std::setprecision(7) 
		 << std::setw(6) << get_Name(abs(prms[i-1])) << "(" << j << ") = "
		 << std::setw(11) << get_bnL_s(prms[i-1], 1, bns[i-1])
		 << std::setw(2) << bns[i-1] << std::endl;
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


double get_fixed_points(const double Ax, const double Ay, const double delta,
			const int n, const int N, const bool prt)
{
  int              i, j, k, l, n1;
  double           A[2], dJ2, dr2, phi0[2], phi1[2];
  double           dphi[2], phi[2], D[2], twoJ0, twoJ1, r_max, S;
  ss_vect<double>  ps0, ps1;
  ss_vect<tps>     ps0_Fl, ps1_Fl;
  std::ofstream         outf;

  const bool    Floq    = false;
  const int     n_turn  = 10;
  const double  D_max   = 10.0, dtwoJ_max = 100e-2;

  if (prt) file_wr(outf, "std::fixed_points.out");

  r_max = sqrt(sqr(Ax)+sqr(Ay));
  n1 = 0; S = 0.0;
  for (i = -n; i <= n; i++) {
    A[X_] = i*Ax/n; 
    for (j = -n; j <= n; j++) {
      n1++; A[Y_] = j*Ay/n;
      ps0.zero(); ps0[x_] = A[X_]; ps0[y_] = A[Y_]; ps0[delta_] = delta;
      ps0_Fl = get_FS(ps0);
      twoJ0 = 0.0;
      for (l = 0; l <= 3; l++)
	twoJ0 += sqr(ps0_Fl[l].cst());
      for (l = 0; l <= 1; l++) {
	phi0[l] = -atan2(ps0_Fl[2*l+1].cst(), ps0_Fl[2*l].cst()); phi[l] = 0.0;
      }
      ps1 = ps0; dJ2 = 0.0; dr2 = 0.0;
      for (k = 1; k <= n_turn; k++) {
	if (ps1.propagate(1, n_elem)) {
	  ps1_Fl = get_FS(ps1);
	  for (l = 0; l <= 1; l++) {
	    phi1[l] = -atan2(ps1_Fl[2*l+1].cst(), ps1_Fl[2*l].cst());
	    dphi[l] = phi1[l] - phi0[l];
	    if (dphi[l] < 0.0)
	      dphi[l] += 2.0*M_PI;
	    else if (dphi[l] > 2.0*M_PI)
	      dphi[l] -= 2.0*M_PI;
	    phi[l] += dphi[l]; phi0[l] = phi1[l];
	  }
	  twoJ1 = 0.0;
	  for (l = 0; l <= 3; l++) {
	    twoJ1 += sqr(ps1_Fl[l].cst());
	    dr2 += sqr(ps1_Fl[l].cst()-ps0_Fl[l].cst());
	  }
	  dJ2 += sqr(twoJ1-twoJ0);
	} else
	  break;
      }

      if (!lost) {
	if (twoJ0 != 0.0)
	  dJ2 /= n_turn*sqr(twoJ0);
	else
	  dJ2 = 0.0;
	dJ2 = min(dJ2, sqr(dtwoJ_max));
	dr2 /= n_turn; dr2 = min(dr2, sqr(r_max));
	for (l = 0; l <= 1; l++) {
	  phi[l] /= n_turn; D[l] = min(1.0/fabs(sin(N*phi[l])), D_max);
	}
      } else {
	dJ2 = sqr(dtwoJ_max); dr2 = sqr(r_max);
	for (l = 0; l <= 1; l++) {
	  phi[l] = 0.0; D[l] = D_max;
	}
      }

      S += D[X_] + D[Y_];

      if (prt) {
	if (!Floq)
	  outf << std::fixed << std::setprecision(3) << std::setw(3) << n1
	       << std::setw(8) << 1e3*ps0[x_] << std::setw(8) << 1e3*ps0[y_];
	else
	  outf << std::fixed << std::setprecision(3) << std::setw(3) << n1
	       << std::setw(8) << 1e3*ps0_Fl[x_].cst()
	       << std::setw(8) << 1e3*ps0_Fl[y_].cst();
	outf << std::scientific << std::setw(10)
	     << std::setw(10) << log(1.0+1e3*sqrt(dr2))
	     << std::setw(10) << log(1.0+sqrt(dJ2))
	     << std::fixed << std::setprecision(5)
	     << std::setw(8) << phi[0]/(2.0*M_PI) << std::setw(8) << phi[1]/(2.0*M_PI)
	     << std::scientific << std::setprecision(3)
	     << std::setw(10) << D[0] << std::setw(10) << D[1] << std::endl;
      }
    }
    if (prt) outf << std::endl;
  }
  if (prt) outf.close();

  return S;
}


void prt_scan(const int n, const double nu_x, const double nu_y,
	      const double DA, const double x_min[], const double x_max[])
{

  tune_out << std::fixed << std::setprecision(4)
	   << "n = " << std::setw(4) << n
	   << ", nu_x= " << nu_x << ", nu_y= " << nu_y
	   << std::setprecision(1)
	   << ", DA^= " << std::setw(6) << 1e6*DA
	   << ", x^= " << std::setw(5) << 1e3*x_min[X_] << " - "
	   << std::setw(4) << 1e3*x_max[X_] << " mm"
	   << ", y^= " << std::setw(5) << 1e3*x_min[Y_] << " - "
	   << std::setw(4) << 1e3*x_max[Y_] << " mm" << std::endl;
}


void get_DA(const int n, const double nu_x, const double nu_y,
	    const double delta, const bool prt)
{
  double  DA, DA_m, x_min[2], x_max[2];

  DA_m = 0.0;

  if (delta != 0.0) {
    std::cout << std::endl;
    std::cout << std::fixed << std::setprecision(4)
	 << "#n = " << std::setw(4) << n
	 << ", nu_x = " << nu_x << ", nu_y = " << nu_y << ", "; 

    DA = get_dynap(10e-3, -delta, n_track, 0.1e-3, n_aper,
		   x_min, x_max, false);
    DA_m += DA;
    if (prt) {
      tune_out << "#";
      prt_scan(n, nu_x, nu_y, DA, x_min, x_max);
    }

    std::cout << std::fixed << std::setprecision(3)
	 << "#n = " << std::setw(4) << n
	 << ", nu_x = " << nu_x << ", nu_y = " << nu_y << ", "; 
    DA = get_dynap(10e-3, delta, n_track, 0.1e-3, n_aper,
		   x_min, x_max, false); 
    DA_m += DA;
    if (prt) {
      tune_out << "#";
      prt_scan(n, nu_x, nu_y, DA, x_min, x_max);
    }
  }
	
  std::cout << std::fixed << std::setprecision(3)
       << " n = " << std::setw(4) << n
       << ", nu_x = " << nu_x << ", nu_y = " << nu_y << ", "; 
  DA = get_dynap(10e-3, 0.0, n_track, 0.1e-3, n_aper, x_min, x_max, false); 
  DA_m += DA;
  DA /= (delta != 0.0)? 3.0 : 1.0;
  if (prt) {
    tune_out << " ";
    prt_scan(n, nu_x, nu_y, DA, x_min, x_max);
  }
}


void get_FP(const int n, const double nu_x, const double nu_y)
{
  double  DA, DA_sum;

  const double  delta = 2.5e-2;

  A1_inv = Inv(A1);
  DA = get_fixed_points(max_Ax, max_Ay, -delta, 25, 15, false);

  DA_sum = DA;

  tune_out << std::fixed << std::setprecision(5)
	   << "#n = " << std::setw(4) << n
	   << ", nu_x= " << nu_x << ", nu_y= " << nu_y
	   << std::setprecision(3)
	   << ", S= " << std::setw(8) << DA << std::endl;

  A1_inv = Inv(A1);
  DA = get_fixed_points(max_Ax, max_Ay, delta, 25, 15, false);

  DA_sum += DA;

  tune_out << std::fixed << std::setprecision(5)
	   << "#n = " << std::setw(4) << n
	   << ", nu_x= " << nu_x << ", nu_y= " << nu_y
	   << std::setprecision(3)
	   << ", S= " << std::setw(8) << DA << std::endl;

  A1_inv = Inv(A1);
  DA = get_fixed_points(max_Ax, max_Ay, 0.0, 25, 15, false);

  DA_sum += DA; DA = DA_sum/3.0;

  tune_out << std::fixed << std::setprecision(5)
	   << " n = " << std::setw(4) << n
	   << ", nu_x= " << nu_x << ", nu_y= " << nu_y
	   << std::setprecision(3)
	   << ", S= " << std::setw(8) << DA << std::endl;
}


void tune_scan(const int n_b2, const int b2s[], const bool opt)
{
  int           i, j, k, l, n;
  long int      loc;
  double        nu_x, nu_y, nu[3], ksi[2];
  double        dnu_x, dnu_y;
  ss_vect<tps>  nus;

  file_wr(tune_out, "tune_scan.dat");
  file_wr(quad_out, "quad.dat");
  file_wr(sext_out, "sext.dat");

  std::cout << std::endl;
  danot_(2);
  get_Map();
  //  get_COD(10, 1e-10, 0.0, true);
  K = MapNorm(Map, g, A1, A0, Map_res, 1);
  nus = dHdJ(K); get_nu_ksi(nus, nu, ksi); prt_nu(nus);

  n = 0;
  dnu_x = (nu_x_max-nu_x_min)/(n_steps-1);
  dnu_y = (nu_y_max-nu_y_min)/(n_steps-1);
  for (i = 0; i < n_steps; i++) {
    nu_x = nu_x_min + i*dnu_x;
    for (j = 0; j < n_steps; j++) {
      n++; nu_y = nu_y_min + j*dnu_y;
      std::cout << std::endl;
      std::cout << std::fixed << std::setprecision(5)
	   << "n = " << std::setw(4) << n
	   << ", nu_x = " << nu_x << ", nu_y = " << nu_y << ", "; 
      if (strcmp(lattice, "NSLS-II_BOOSTER") == 0)
	fit_tune(nu_x, nu_y, n_b2, b2s, eps_nu, true);
      else if (strstr(lattice, "NSLS-II") != NULL)
	fit_tune(nu_x, nu_y,
		 beta1[X_], beta1[Y_], beta_loc1,
		 beta2[X_], beta2[Y_], beta_loc2,
		 beta3[X_], beta3[Y_], beta_loc3,
		 n_b2, b2s, eps_nu, true);

      quad_out << std::endl;
      quad_out << "n = " << n << ":" << std::endl;
      for (k = 1; k <= n_b2; k++)
	for (l = 1; l <= get_n_Kids(abs(b2s[k-1])); l++)
	  if (b2s[k-1] > 0)
	    quad_out << std::fixed << std::setprecision(7) 
		     << std::setw(6) << get_Name(b2s[k-1]) << "(" << l << ") = "
		     << std::setw(11) << get_bnL(b2s[k-1], l, Quad)
		     << std::setw(2) << Quad << std::endl;
	  else {
	    loc = get_loc(abs(b2s[k-1]), l) - 1;
	    quad_out << std::fixed << std::setprecision(7) 
		     << std::setw(6) << get_Name(elem[loc-1].Fnum)
		     << "(" << l << ") = "
		     << std::setw(11) << get_L(elem[loc-1].Fnum, l)
		     << std::setw(3) << -Quad << std::endl;
	    quad_out << std::fixed << std::setprecision(7) 
		     << std::setw(6) << get_Name(elem[loc+1].Fnum)
		     << "(" << l << ") = "
		     << std::setw(11) << get_L(elem[loc+1].Fnum, l)
		     << std::setw(3) << -Quad << std::endl;
	  }

 
      sext_out << std::endl;
      sext_out << "n = " << n << ":" << std::endl;

      if (opt)
	//	h_zero(no_tps, 5, nu_x, nu_y, false);
	H_zero(eps_ksi, 10, false);

      if (true)
	get_DA(n, nu_x, nu_y, 2.5e-2, true);
      else
	get_FP(n, nu_x, nu_y);
    }
    tune_out << std::endl;
    std::cout << std::endl;
  }
  tune_out.close(); quad_out.close(); sext_out.close();
}


void tune_scan(const bool opt)
{
  const int  b2s_max = 10;

  char          line[max_str];
  int           k, n;
  double        nu[2];
  double        b2s[b2s_max];
  ss_vect<tps>  nus;
  std::ifstream      inf1, inf2;
  std::ofstream      outf;

  char  file_name1[] = "/home/bengtsson/projects/src/quad-K1_1.txt";
  char  file_name2[] = "/home/bengtsson/projects/src/quad-K1_2.txt";

  const int  n_skip = 4, n_grid = 11;

  file_wr(tune_out, "tune_scan.dat"); file_wr(sext_out, "sext.dat");

  file_rd(inf1, file_name1); file_rd(inf2, file_name2);

  for (k = 1; k <= n_skip; k++) {
    inf1.getline(line, max_str); inf2.getline(line, max_str);
  }

  std::cout << std::endl;
  n = 0;
  while (inf1.getline(line, max_str)) {
    n++;

    sscanf(line, "%d %lf %lf %lf %lf %lf %lf",
	   &k, &nu[X_], &nu[Y_], &b2s[0], &b2s[1], &b2s[2], &b2s[3]);
    inf2.getline(line, max_str);
    sscanf(line, "%d %lf %lf %lf %lf %lf %lf",
	   &k, &nu[X_], &nu[Y_], &b2s[4], &b2s[5], &b2s[6], &b2s[7]);

    if (false)
      printf("%3d %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f\n",
	     k,
	     b2s[0], b2s[1], b2s[2], b2s[3], b2s[4], b2s[5], b2s[6], b2s[7]);

    set_bn(get_Fnum("qh1"), Quad, b2s[0]);
    set_bn(get_Fnum("qh2"), Quad, b2s[1]);
    set_bn(get_Fnum("qh3"), Quad, b2s[2]);
    set_bn(get_Fnum("ql1"), Quad, b2s[3]);
    set_bn(get_Fnum("ql2"), Quad, b2s[4]);
    set_bn(get_Fnum("ql3"), Quad, b2s[5]);
    set_bn(get_Fnum("qm1"), Quad, b2s[6]);
    set_bn(get_Fnum("qm2"), Quad, b2s[7]);

    //    danot_(2);
    //  get_Map();
    //    get_COD(10, 1e-10, 0.0, true);
    //    K = MapNorm(Map, g, A1, A0, Map_res, 1);
    //    nus = dHdJ(K); get_nu_ksi(nus, nu, ksi);

    std::cout << std::endl;
    std::cout << std::fixed << std::setprecision(4)
	 << "n = " << std::setw(4) << n
	 << ", nu_x = " << nu[X_] << ", nu_y = " << nu[Y_] << ", "; 

    sext_out << std::endl;
    sext_out << "n = " << n << ":" << std::endl;

    if (opt) H_zero(eps_ksi, 10, false);

    if (true)
      get_DA(n, nu[X_], nu[Y_], 2.5e-2, true);
    else
      get_FP(n, nu[X_], nu[Y_]);

    if (n % n_grid == 0) {
      tune_out << std::endl;
      std::cout << std::endl;
    }
  }

  inf1.close(); inf2.close();

  tune_out.close(); sext_out.close();
}


tps g_renorm(const double nu0_x, const double nu0_y,
	     const double nu1_x, const double nu1_y,
	     const tps &g)
{
  // Renormalize g: (1-R^-1)^-1 * h 

  long int      jj1[ss_dim], jj2[ss_dim];
  int           i, j, k, l, m;
  double        re, im, cotan0, cotan1, cotan0_sqr;
  tps           h_re, h_im, g_re, g_im, G_re, G_im, mn1, mn2;
  ss_vect<tps>  Id;

  CtoR(g, g_re, g_im);

  for (k = 0; k < ss_dim; k++) {
    jj1[k] = 0; jj2[k] = 0;
  }

  Id.identity(); G_re = 0.0; G_im = 0.0;
  for (i = 0; i <= no_tps; i++) {
    jj1[x_] = i; jj2[px_] = i;
    for (j = 0; j <= i; j++) {
      jj1[px_] = j; jj2[x_] = j;
      for (k = 0; k <= no_tps; k++) {
	jj1[y_] = k; jj2[py_] = k;
	for (l = 0; l <= no_tps; l++) {
	  jj1[py_] = l; jj2[y_] = l;

	  if (i+j+k+l <= no_tps) {
	    cotan0 = 1.0/tan(((i-j)*nu0_x+(k-l)*nu0_y)*M_PI);
	    cotan0_sqr = sqr(cotan0);
	    cotan1 = 1.0/tan(((i-j)*nu1_x+(k-l)*nu1_y)*M_PI);
	    mn1 =
	      pow(Id[x_], i)*pow(Id[px_], j)*pow(Id[y_], k)*pow(Id[py_], l);
	    mn2 =
	      pow(Id[x_], j)*pow(Id[px_], i)*pow(Id[y_], l)*pow(Id[py_], k);

	    for (m = 0; m <= no_tps; m++) {
	      if (i+j+k+l+m <= no_tps) {
		jj1[delta_] = m; jj2[delta_] = m;
		if ((i != j) || (k != l)) {
		  re = g_re[jj1]; im = g_im[jj1];

		  // compute h
		  h_re = (re+cotan0*im)*2.0/(1.0+cotan0_sqr);
		  h_im = (im-cotan0*re)*2.0/(1.0+cotan0_sqr);

		  // renormalize g
		  G_re += (h_re-cotan1*h_im)*(mn1+mn2)*pow(Id[delta_], m)/2.0;
		  G_im += (h_im+cotan1*h_re)*(mn1-mn2)*pow(Id[delta_], m)/2.0;
		  g_re.pook(jj2, 0.0); g_im.pook(jj2, 0.0);
		}
	      }
	    }
	  }
	}
      }
    }
  }

  return RtoC(G_re, G_im);
}


tps get_H(const tps &g)
{
  int           i;
  tps           H, gn;
  ss_vect<tps>  Id, Mn;

  // Construct generator.
  // K is in Dragt-Finn form but the generators commute.
  Id.identity(); H = K;
  for (i = no_tps; i >= 3; i--) {
    gn = Take(g, i); H = H*LieExp(-gn, Id);
  }

  return H;
}


void track_H(const char *file_name, const double Ax, const double Ay)
{
  int              k;
  double           h, h0, J0[2], J[2], phi[2], nu[2], scl, J2[2], phi2[2];
  tps              H, g_r, scl_tps;
  ss_vect<double>  ps, ps_Fl, ps_nl;
  ss_vect<tps>     nus, A1_inv, A_nl_inv, Id;
  std::ofstream         H_out;

  const bool  relative = false;

  file_wr(H_out, file_name);

  danot_(no_tps);
  get_Map(); K = MapNorm(Map, g, A1, A0, Map_res, 1); nus = dHdJ(K);
  A1_inv = Inv(A1); A_nl_inv = get_A_nl_inv(g);

  ps.zero();
  ps[x_] = Ax; ps[px_] = 0.0; ps[y_] = Ay; ps[py_] = 0.0;
  ps[delta_] = 0.0; ps[ct_] = 0.0;

  ps_Fl = (A1_inv*ps).cst(); ps_nl = (A_nl_inv*ps_Fl).cst();

  // scale action-angle variables
  for (k = 0; k < 2; k++) {
    nu[k] = (nus[k]*ps_nl).cst();
    scl = sqrt(nus[k].cst()/nu[k]);
    ps_Fl[2*k] *= scl; ps_Fl[2*k+1] *= scl;
  }

  if (false)
    g_r = g_renorm(nus[0].cst(), nus[1].cst(), nu[X_], nu[Y_], g);
  else
    g_r = g;

  H = get_H(g_r); A_nl_inv = get_A_nl_inv(g_r);

  J0[X_] = (sqr(ps_Fl[x_])+sqr(ps_Fl[px_]))/2.0;
  J0[Y_] = (sqr(ps_Fl[y_])+sqr(ps_Fl[py_]))/2.0;
  h0 = (H*ps_Fl).cst();

  for (k = 1; k <= 500; k++) {
    ps.propagate(1, n_elem);
    ps_Fl = (A1_inv*ps).cst();
    J[X_] = (sqr(ps_Fl[x_])+sqr(ps_Fl[px_]))/2.0;
    phi[X_] = atan2(ps_Fl[px_], ps_Fl[x_]);
    J[Y_] = (sqr(ps_Fl[y_])+sqr(ps_Fl[py_]))/2.0;
    phi[Y_] = atan2(ps_Fl[py_], ps_Fl[y_]);

    h = (H*ps_Fl).cst();

    ps_nl = (A_nl_inv*ps_Fl).cst();
    J2[X_] = (sqr(ps_nl[x_])+sqr(ps_nl[px_]))/2.0;
    phi2[X_] = atan2(ps_nl[px_], ps_nl[x_]);
    J2[Y_] = (sqr(ps_nl[y_])+sqr(ps_nl[py_]))/2.0;
    phi2[Y_] = atan2(ps_nl[py_], ps_nl[y_]);

    if (!relative)
      H_out << std::scientific << std::setprecision(5)
	    << std::setw(5) << k
	    << std::setw(13) << 1e3*ps_Fl[x_] << std::setw(13) << 1e3*ps_Fl[px_]
	    << std::setw(13) << 1e3*ps_Fl[y_] << std::setw(13) << 1e3*ps_Fl[py_]
	    << std::setw(13) << 1e6*J[X_] << std::setw(13) << phi[X_]
	    << std::setw(13) << 1e6*J[Y_] << std::setw(13) << phi[Y_]
	    << std::setw(13) << 1e6*J2[X_] << std::setw(13) << phi2[X_]
	    << std::setw(13) << 1e6*J2[Y_] << std::setw(13) << phi2[Y_]
	    << std::setw(13) << 1e6*fabs(h) << std::endl;
    else
      H_out << std::fixed << std::setprecision(1)
	    << std::setw(5) << k << std::setw(6) << 1e2*(J[X_]-J0[X_])/J0[X_]
	    << std::setw(6) << 1e2*(h-h0)/h0 << std::endl;
  }

  H_out.close();
}


void track_map(const double J_max)
{
  int              i, j, k;
  double           phi, A, x, px;
  ss_vect<double>  ps;
  ss_vect<tps>     A1_inv, PS;
  std::ofstream         ps_out;

  const int     n_r = 10, n_phi = 50, n = 10;

  file_wr(ps_out, "track_map.dat");

  danot_(no_tps); A1_inv = Inv(A1); Map = A1_inv*Map*A1; ps.zero();
  for (i = 1; i <= n_r; i++) {
    A = i*sqrt(2.0*J_max)/n_r;
    for (j = 0; j < n_phi; j++) {
      phi = j*2.0*pi/n_phi; x = A*cos(phi); px = A*sin(phi);

      ps[x_] = x; ps[px_] = px; ps[y_] = 0.0; ps[py_] = 0.0;
      ps[delta_] = 0.0; ps[ct_] = 0.0;
      ps = (A1*ps).cst();
      for (k = 1; k <= n; k++)
	ps.propagate(1, n_elem);
      ps = (A1_inv*ps).cst();
 
      ps_out << std::scientific << std::setprecision(5)
	     << std::setw(13) << ps[x_] << std::setw(13) << ps[px_]
	     << std::setw(13) << ps[y_] << std::setw(13) << ps[py_]
	     << std::setw(13) << ps[delta_] << std::setw(13) << ps[ct_];

      PS[x_] = x; PS[px_] = px; PS[y_] = 0.0; PS[py_] = 0.0;
      PS[delta_] = 0.0; PS[ct_] = 0.0;
      for (k = 1; k <= n; k++)
	PS = Map*PS;

      ps_out << std::scientific << std::setprecision(5)
	     << std::setw(13) << PS[x_].cst()-ps[x_]
	     << std::setw(13) << PS[px_].cst()-ps[px_]
	     << std::setw(13) << PS[y_].cst()-ps[y_]
	     << std::setw(13) << PS[py_].cst()-ps[py_]
	     << std::setw(13) << PS[delta_].cst()-ps[delta_]
	     << std::setw(13) << PS[ct_].cst()-ps[ct_]
	     << std::endl;
    }
  }

  ps_out.close();
}


void get_dnu(const int i0, const int i1, double dnu[])
{
  double  nu0[2], nu1[2];

  get_nu(nu0, i0); A1.propagate(i0, i1); get_nu(nu1, 0);
  dnu[X_] = nu1[X_] - nu0[X_];
  if (dnu[X_] < 0.0) dnu[X_] += 1.0;
  dnu[Y_] = nu1[Y_] - nu0[Y_];
  if (dnu[Y_] < 0.0) dnu[Y_] += 1.0;
}


void no_mpoles(void)
{
  int j, k;

  std::cout << std::endl;
  std::cout << "zeroing multipoles" << std::endl;
  std::cout << std::endl;
  for (j = 0; j < n_elem; j++)
    if (elem[j].kind == Mpole)
      for (k = Sext; k < mpole_max; k++) {
	//	std::cout << "zeroing " << elem[j].Name << std::endl;
	set_bn(elem[j].Fnum, elem[j].Knum, k, 0.0);
      }
}


void get_prm(const char *file_name)
{
  char      line[max_str];      
  std::ifstream  prm_in;

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

  std::cout << std::endl;
  std::cout << std::fixed << std::setprecision(6)
       << "fit_tune      = " << adj_tune
       << ", nu0_x  = " << nu0[X_] << ", nu0_y  = " << nu0[Y_]
       << std::scientific << std::setprecision(1) << ", eps_nu = " << eps_nu << std::endl;
  std::cout << std::fixed << std::setprecision(6)
       << "fit_chrom     = " << adj_chrom
       << ", ksi0_x = " << ksi1[X_] << ", ksi0_y = " << ksi1[Y_]
       << std::scientific << std::setprecision(1) << ", eps_ksi = " << eps_ksi << std::endl;
  std::cout << "check_range   = " << check_range << std::endl;
  std::cout << std::endl;
  std::cout << std::fixed << std::setprecision(6)
       << "n_steps   = " << n_steps
       << ", nu_x_min = " << nu_x_min << ", nu_x_max = " << nu_x_max
       << ", nu_y_min = " << nu_y_min << ", nu_y_max = " << nu_y_max << std::endl;
  std::cout << std::endl;
  std::cout << "n_cell        = " << n_cell << std::endl;
  std::cout << std::fixed << std::setprecision(2)
       << "ds_max        = " << ds_max << std::endl;
  std::cout << std::fixed << std::setprecision(1)
       << "b2_max        = " << b2_max << std::endl;
  std::cout << std::fixed << std::setprecision(1)
       << "b3L_max       = " << bnL_max[Sext] << std::endl;
  std::cout << std::fixed << std::setprecision(1)
       << "b4L_max       = " << bnL_max[Oct] << std::endl;
  std::cout << std::fixed << std::setprecision(1)
       << "b5L_max       = " << bnL_max[Dec] << std::endl;
  std::cout << std::fixed << std::setprecision(1)
       << "scl_dnu       = " << scl_dnu << std::endl;
  std::cout << std::fixed << std::setprecision(1)
       << "scl_ksi1      = " << scl_ksi1 << std::endl;
  std::cout << std::fixed << std::setprecision(1)
       << "scl_ksi2      = " << scl_ksi2 << std::endl;
  std::cout << std::fixed << std::setprecision(1)
       << "scl_ksi_nl    = " << scl_ksi_nl << std::endl;
  std::cout << std::scientific << std::setprecision(1)
       << "scl_dnuddelta = " << scl_dnuddelta << std::endl;
  std::cout << std::scientific << std::setprecision(1)
       << "scl_dnuddJ    = " << scl_dnudJ << std::endl;
  std::cout << std::fixed << std::setprecision(2)
       << "step          = " << step << std::endl;
}


void sxt_h1(void)
{
  long int      jj[ss_dim];
  int           i;
  tps           h, h_re, h_im;

  danot_(no_tps-1);

  get_Map();

  danot_(no_tps);

  h = get_h(); h = h*Id_scl; CtoR(h, h_re, h_im);

  for (i = 0; i < ss_dim; i++)
    jj[i] = 0;

  // linear chromaticity
  std::cout << std::endl;
  jj[x_] = 1; jj[px_] = 1; jj[y_] = 0; jj[py_] = 0; jj[delta_] = 1;
  std::cout << std::scientific << std::setprecision(16)
	    << "h_11001:" << std::setw(24) << 2.0*h_re[jj]
	    << std::setw(24) << 2.0*h_im[jj] << std::endl;
  jj[x_] = 0; jj[px_] = 0; jj[y_] = 1; jj[py_] = 1; jj[delta_] = 1;
  std::cout << std::scientific << std::setprecision(16)
	    << "h_00111:" << std::setw(24) << 2.0*h_re[jj]
	    << std::setw(24) << 2.0*h_im[jj] << std::endl;

  // first order chromatic terms
  std::cout << std::endl;
  jj[x_] = 2; jj[px_] = 0; jj[y_] = 0; jj[py_] = 0; jj[delta_] = 1;
  std::cout << std::scientific << std::setprecision(16)
	    << "h_20001:" << std::setw(24) << 2.0*h_re[jj]
	    << std::setw(24) << 2.0*h_im[jj] << std::endl;
  jj[x_] = 0; jj[px_] = 0; jj[y_] = 2; jj[py_] = 0;
  std::cout << std::scientific << std::setprecision(16)
	    << "h_00201:" << std::setw(24) << 2.0*h_re[jj]
	    << std::setw(24) << 2.0*h_im[jj] << std::endl;
  jj[x_] = 1; jj[px_] = 0; jj[y_] = 0; jj[py_] = 0; jj[delta_] = 2;
  std::cout << std::scientific << std::setprecision(16)
	    << "h_10002:" << std::setw(24) << 2.0*h_re[jj]
	    << std::setw(24) << 2.0*h_im[jj] << std::endl;

  // normal sextupoles
  std::cout << std::endl;
  jj[x_] = 2; jj[px_] = 1; jj[y_] = 0; jj[py_] = 0; jj[delta_] = 0;
  std::cout << std::scientific << std::setprecision(16)
	    << "h_21000:" << std::setw(24) << 2.0*h_re[jj]
	    << std::setw(24) << 2.0*h_im[jj] << std::endl;
  jj[x_] = 3; jj[px_] = 0; jj[y_] = 0; jj[py_] = 0;
  std::cout << std::scientific << std::setprecision(16)
	    << "h_30000:" << std::setw(24) << 2.0*h_re[jj]
	    << std::setw(24) << 2.0*h_im[jj] << std::endl;
  jj[x_] = 1; jj[px_] = 0; jj[y_] = 1; jj[py_] = 1;
  std::cout << std::scientific << std::setprecision(16)
	    << "h_10110:" << std::setw(24) << 2.0*h_re[jj]
	    << std::setw(24) << 2.0*h_im[jj] << std::endl;
  jj[x_] = 1; jj[px_] = 0; jj[y_] = 0; jj[py_] = 2;
  std::cout << std::scientific << std::setprecision(16)
	    << "h_10020:" << std::setw(24) << 2.0*h_re[jj]
	    << std::setw(24) << 2.0*h_im[jj] << std::endl;
  jj[x_] = 1; jj[px_] = 0; jj[y_] = 2; jj[py_] = 0;
  std::cout << std::scientific << std::setprecision(16)
	    << "h_10200:" << std::setw(24) << 2.0*h_re[jj]
	    << std::setw(24) << 2.0*h_im[jj] << std::endl;

  if (false) {
    // skew sextupoles
    std::cout << std::endl;
    jj[x_] = 0; jj[px_] = 0; jj[y_] = 2; jj[py_] = 1;
    std::cout << std::scientific << std::setprecision(2)
	      << "h_00210:" << std::setw(10) << 2.0*h_re[jj]
	      << std::setw(10) << 2.0*h_im[jj]<< std::endl;
    jj[x_] = 0; jj[px_] = 0; jj[y_] = 3; jj[py_] = 0;
    std::cout << std::scientific << std::setprecision(2)
	      << "h_00300:" << std::setw(10) << 2.0*h_re[jj]
	      << std::setw(10) << 2.0*h_im[jj] << std::endl;
    jj[x_] = 1; jj[px_] = 1; jj[y_] = 1; jj[py_] = 0;
    std::cout << std::scientific << std::setprecision(2)
	      << "h_11100:" << std::setw(10) << 2.0*h_re[jj]
	      << std::setw(10) << 2.0*h_im[jj] << std::endl;
    std::cout << std::endl;
    jj[x_] = 0; jj[px_] = 2; jj[y_] = 1; jj[py_] = 0;
    std::cout << std::scientific << std::setprecision(2)
	      << "h_02100:" << std::setw(10) << 2.0*h_re[jj]
	      << std::setw(10) << 2.0*h_im[jj] << std::endl;
    jj[x_] = 2; jj[px_] = 0; jj[y_] = 1; jj[py_] = 0;
    std::cout << std::scientific << std::setprecision(2)
	      << "h_20100:" << std::setw(10) << 2.0*h_re[jj]
	      << std::setw(10) << 2.0*h_im[jj] << std::endl;
  }
}


void sxt_h2(void)
{
  long int      jj[ss_dim];
  int           i;
  tps           K_re, K_im, h, h_re, h_im;

  danot_(no_tps-1);

  get_Map();
  K = MapNorm(Map, g, A1, A0, Map_res, 1); CtoR(K, K_re, K_im);

  danot_(no_tps);

  h = get_h(); h = h*Id_scl; CtoR(h, h_re, h_im);

  for (i = 0; i < ss_dim; i++)
    jj[i] = 0;

  std::cout << std::endl;
  jj[x_] = 4; jj[px_] = 0; jj[y_] = 0; jj[py_] = 0;
  std::cout << std::scientific << std::setprecision(16)
	    << "h_40000:" << std::setw(24) << 2.0*h_re[jj]
	    << std::setw(24) << 2.0*h_im[jj] << std::endl;
  jj[x_] = 3; jj[px_] = 1; jj[y_] = 0; jj[py_] = 0;
  std::cout << std::scientific << std::setprecision(16)
	    << "h_31000:" << std::setw(24) << 2.0*h_re[jj]
	    << std::setw(24) << 2.0*h_im[jj] << std::endl;
  jj[x_] = 2; jj[px_] = 0; jj[y_] = 1; jj[py_] = 1;
  std::cout << std::scientific << std::setprecision(16)
	    << "h_20110:" << std::setw(24) << 2.0*h_re[jj]
	    << std::setw(24) << 2.0*h_im[jj] << std::endl;
  jj[x_] = 2; jj[px_] = 0; jj[y_] = 0; jj[py_] = 2;
  std::cout << std::scientific << std::setprecision(16)
	    << "h_20020:" << std::setw(24) << 2.0*h_re[jj]
	    << std::setw(24) << 2.0*h_im[jj] << std::endl;
  jj[x_] = 2; jj[px_] = 0; jj[y_] = 2; jj[py_] = 0;
  std::cout << std::scientific << std::setprecision(16)
	    << "h_20200:" << std::setw(24) << 2.0*h_re[jj]
	    << std::setw(24) << 2.0*h_im[jj] << std::endl;
  jj[x_] = 1; jj[px_] = 1; jj[y_] = 2; jj[py_] = 0;
  std::cout << std::scientific << std::setprecision(16)
	    << "h_11200:" << std::setw(24) << 2.0*h_re[jj]
	    << std::setw(24) << 2.0*h_im[jj] << std::endl;
  jj[x_] = 0; jj[px_] = 0; jj[y_] = 3; jj[py_] = 1;
  std::cout << std::scientific << std::setprecision(16)
	    << "h_00310:" << std::setw(24) << 2.0*h_re[jj]
	    << std::setw(24) << 2.0*h_im[jj] << std::endl;
  jj[x_] = 0; jj[px_] = 0; jj[y_] = 4; jj[py_] = 0;
  std::cout << std::scientific << std::setprecision(16)
	    << "h_00400:" << std::setw(24) << 2.0*h_re[jj]
	    << std::setw(24) << 2.0*h_im[jj] << std::endl;

  std::cout << std::endl;
  jj[x_] = 2; jj[px_] = 2; jj[y_] = 0; jj[py_] = 0;
  std::cout << std::scientific << std::setprecision(16)
	    << "h_22000:" << std::setw(24) << 2.0*h_re[jj]
	    << std::setw(24) << 2.0*h_im[jj]
	    << std::endl;
  jj[x_] = 1; jj[px_] = 1; jj[y_] = 1; jj[py_] = 1;
  std::cout << std::scientific << std::setprecision(16)
	    << "h_11110:" << std::setw(24) << 2.0*h_re[jj]
	    << std::setw(24) << 2.0*h_im[jj] << std::endl;
  jj[x_] = 0; jj[px_] = 0; jj[y_] = 2; jj[py_] = 2;
  std::cout << std::scientific << std::setprecision(16)
	    << "h_00220:" << std::setw(24) << 2.0*h_re[jj]
	    << std::setw(24) << 2.0*h_im[jj] << std::endl;

  std::cout << std::endl;
  jj[x_] = 2; jj[px_] = 2; jj[y_] = 0; jj[py_] = 0;
  std::cout << std::scientific << std::setprecision(16) << "a_xx ="
	    << std::setw(24) << K_re[jj]/M_PI << std::endl;
  jj[x_] = 1; jj[px_] = 1; jj[y_] = 1; jj[py_] = 1;
  std::cout << std::scientific << std::setprecision(16) << "a_xy ="
	    << std::setw(24) << K_re[jj]/M_PI << std::endl;
  jj[x_] = 0; jj[px_] = 0; jj[y_] = 2; jj[py_] = 2;
  std::cout << std::scientific << std::setprecision(16) << "a_yy ="
	    << std::setw(24) << K_re[jj]/M_PI << std::endl;
}


void get_b2s(int &n_b2, int b2_Fams[])
{

  n_b2 = 0;
  if (strcmp(lattice, "NSLS-II_TBA") == 0) {
    b2_Fams[n_b2++] = get_Fnum("qi1");

    b2_Fams[n_b2++] = get_Fnum("qi2a");
    b2_Fams[n_b2++] = get_Fnum("qi2b");
    b2_Fams[n_b2++] = get_Fnum("qi3");

    if (false) {
      b2_Fams[n_b2++] = get_Fnum("qi1");
      b2_Fams[n_b2++] = get_Fnum("qi2a"); b2_Fams[n_b2++] = get_Fnum("qi2b");
      b2_Fams[n_b2++] = get_Fnum("qi3");

      b2_Fams[n_b2++] = -get_Fnum("twl1"); b2_Fams[n_b2++] = -get_Fnum("twl2");
      b2_Fams[n_b2++] = -get_Fnum("twl3"); b2_Fams[n_b2++] = -get_Fnum("twl4");
    }
  } else if (strcmp(lattice, "NSLS-II_DBA") == 0) {
    b2_Fams[n_b2++] = get_Fnum("ql1");
    b2_Fams[n_b2++] = get_Fnum("ql2");
    b2_Fams[n_b2++] = get_Fnum("ql3");

    b2_Fams[n_b2++] = get_Fnum("qh1");
    b2_Fams[n_b2++] = get_Fnum("qh2");
    b2_Fams[n_b2++] = get_Fnum("qh3");

    if (false) {
      b2_Fams[n_b2++] = -get_Fnum("q1");  b2_Fams[n_b2++] = -get_Fnum("q2");
      b2_Fams[n_b2++] = -get_Fnum("q3");  b2_Fams[n_b2++] = -get_Fnum("q4");

      //      b2_Fams[n_b2++] = -get_Fnum("q11");
      b2_Fams[n_b2++] = -get_Fnum("q22");
      b2_Fams[n_b2++] = -get_Fnum("q33"); b2_Fams[n_b2++] = -get_Fnum("q44");
    }

    if (false) {
      b2_Fams[n_b2++] = get_Fnum("qd3"); b2_Fams[n_b2++] = get_Fnum("qf2");
    }
  } else if (strcmp(lattice, "NSLS-II_BOOSTER") == 0) {
    b2_Fams[n_b2++] = get_Fnum("qf");
    b2_Fams[n_b2++] = get_Fnum("qg");
    //      b2_Fams[n_b2++] = get_Fnum("qd");
  } else if (strcmp(lattice, "SLS") == 0) {
    b2_Fams[n_b2++] = get_Fnum("qle"); b2_Fams[n_b2++] = get_Fnum("qlf");
    b2_Fams[n_b2++] = get_Fnum("qlg"); b2_Fams[n_b2++] = get_Fnum("qlh");
    b2_Fams[n_b2++] = get_Fnum("qse"); b2_Fams[n_b2++] = get_Fnum("qsf");
    b2_Fams[n_b2++] = get_Fnum("qsg");
    b2_Fams[n_b2++] = get_Fnum("qme"); b2_Fams[n_b2++] = get_Fnum("qmf");
    b2_Fams[n_b2++] = get_Fnum("qmg");
  } else if (strcmp(lattice, "ESRF") == 0) {
    b2_Fams[n_b2++] = get_Fnum("qd4");   b2_Fams[n_b2++] = get_Fnum("qf3");
    b2_Fams[n_b2++] = get_Fnum("qd5");
    b2_Fams[n_b2++] = -get_Fnum("twk4"); b2_Fams[n_b2++] = -get_Fnum("twk5");
    b2_Fams[n_b2++] = -get_Fnum("twk6");
  } else if (strcmp(lattice, "ALBA") == 0) {
    b2_Fams[n_b2++] = get_Fnum("qd1"); b2_Fams[n_b2++] = get_Fnum("qd2");
    b2_Fams[n_b2++] = get_Fnum("qd3"); b2_Fams[n_b2++] = get_Fnum("qf1");
    b2_Fams[n_b2++] = get_Fnum("qf2"); b2_Fams[n_b2++] = get_Fnum("qf3");
    b2_Fams[n_b2++] = get_Fnum("qf4"); b2_Fams[n_b2++] = get_Fnum("qf5");
    b2_Fams[n_b2++] = get_Fnum("qf6"); b2_Fams[n_b2++] = get_Fnum("qf7");
    b2_Fams[n_b2++] = get_Fnum("qf9");
  } else if (strcmp(lattice, "ESRF") == 0) {
    b2_Fams[n_b2++] = get_Fnum("qf2"); b2_Fams[n_b2++] = get_Fnum("qd3");
    b2_Fams[n_b2++] = get_Fnum("qd6"); b2_Fams[n_b2++] = get_Fnum("qf7");
  }
}


void get_prms()
{

  n_prm = 0;
  if (strcmp(lattice, "NSLS-II_DBA") == 0) {
    std::cout << std::endl;
    std::cout << "NSLS-II_DBA:" << std::endl;

    prms[n_prm] = get_Fnum("sl1"); bns[n_prm++] = Sext;
    prms[n_prm] = get_Fnum("sl2"); bns[n_prm++] = Sext;
    prms[n_prm] = get_Fnum("sl3"); bns[n_prm++] = Sext;

    prms[n_prm] = get_Fnum("sm1a"); bns[n_prm++] = Sext;
    prms[n_prm] = get_Fnum("sm1b"); bns[n_prm++] = Sext;
    prms[n_prm] = get_Fnum("sm2"); bns[n_prm++] = Sext;

    prms[n_prm] = get_Fnum("sh1"); bns[n_prm++] = Sext;
    prms[n_prm] = get_Fnum("sh3"); bns[n_prm++] = Sext;
    prms[n_prm] = get_Fnum("sh4"); bns[n_prm++] = Sext;

  } else if (strcmp(lattice, "DIAMOND") == 0) {
    prms[n_prm] = get_Fnum("ss1a"); bns[n_prm++] = Sext;
    prms[n_prm] = get_Fnum("ss2a"); bns[n_prm++] = Sext;
    prms[n_prm] = get_Fnum("ss1b"); bns[n_prm++] = Sext;
    prms[n_prm] = get_Fnum("ss2b"); bns[n_prm++] = Sext;
    prms[n_prm] = get_Fnum("ss1c"); bns[n_prm++] = Sext;
    prms[n_prm] = get_Fnum("ss2c"); bns[n_prm++] = Sext;
    prms[n_prm] = get_Fnum("ss1d"); bns[n_prm++] = Sext;
    prms[n_prm] = get_Fnum("ss2d"); bns[n_prm++] = Sext;
  }

  if (n_prm > n_prm_max) {
    std::cout << "get_prms: n_prm_max exceeded " << n_prm << "(" << n_prm_max
	 << ")" << std::endl;
    exit(0);
  }

  std::cout << std::endl;
  std::cout << "get_prms: no of multipole families " << n_prm << std::endl;
}


void get_lin_map(void)
{
  long int      jj[ss_dim];
  int           i;
  ss_vect<tps>  ps;

  rad_on = true; cavity_on = true;

  get_COD(10, 1e-14, 0.0, true);

  danot_(1);

  for (i = 0; i < ss_dim; i++)
    ps[i] = tps(fixed_point[i], i+1);

  ps.propagate(1, n_elem);
    
  std::cout << std::endl;
  std::cout << "std::fixed point:" << std::endl;
  std::cout << std::endl;
  std::cout << std::scientific << std::setprecision(5) << std::setw(24) << ps.cst() << std::endl;
  
  std::cout << std::endl;
  std::cout << "linear map:" << std::endl;
  std::cout << std::endl;

  for (i = 0; i < ss_dim; i++)
    jj[i] = 0;
  
  for (i = 0; i < 2*nd_tps; i++) {
    jj[x_] = 1;
    std::cout << std::scientific << std::setprecision(5) << std::setw(24) << ps[i][jj];
    jj[x_] = 0; jj[px_] = 1;
    std::cout << std::scientific << std::setprecision(5) << std::setw(24) << ps[i][jj];
    jj[px_] = 0; jj[y_] = 1;
    std::cout << std::scientific << std::setprecision(5) << std::setw(24) << ps[i][jj];
    jj[y_] = 0; jj[py_] = 1;
    std::cout << std::scientific << std::setprecision(5) << std::setw(24) << ps[i][jj];
    jj[py_] = 0; jj[delta_] = 1;
    std::cout << std::scientific << std::setprecision(5) << std::setw(24) << ps[i][jj];
    jj[delta_] = 0; jj[ct_] = 1;
    std::cout << std::scientific << std::setprecision(5) << std::setw(24) << ps[i][jj] << std::endl;
    jj[ct_] = 0;
  }

  exit(0);
}


void get_ampl(const double nux, const double nuy, const int n,
	      int &nx, int &ny, double &A)
{
  // n must be even
  int     i, j;
  double  r;

  // special case for nuy
  A = fabs(1.0/sin(M_PI*nuy)); nx = 0; ny = 1;
  for (i = 0; i <= n; i++) {
    for (j = -n; j <= n; j += 2) {
      // only even coeffs in the vertical plane for normal sextupoles
      if (((i != 0) || (j != 0)) && (i+j <= n) &&
	  ((i > 0) || ((i == 0) && (j > 0)))) {
	r = fabs(1.0/sin(M_PI*(i*nux+j*nuy)));
	if (r > A) {
	  nx = i; ny = j; A = r;
	}
      }
    }
  }
}


void chk_lat(double nu[], double ksi[])
{
  double        alpha1[2];
  ss_vect<tps>  nus;

  //  get_Map();
  get_COD(10, 1e-10, 0.0, true);
  K = MapNorm(Map, g, A1, A0, Map_res, 1);
  nus = dHdJ(K); get_nu_ksi(nus, nu, ksi); get_ab(alpha1, beta1, 0);
  std::cout << std::endl;
  std::cout << std::fixed << std::setprecision(3)
       << "alpha_x  = " << alpha1[X_] << ", alpha_y = " << alpha1[Y_]
       << ", beta_x = " << beta1[X_] << ", beta_y  = " << beta1[Y_] << std::endl;
  prt_nu(nus);
}


void get_locs()
{
  double  alpha1[2], alpha2[2], alpha3[2], alpha4[2];

  beta_loc1 = get_loc(get_Fnum("mp"), 1); get_ab(alpha1, beta1, beta_loc1);
  beta_loc2 = get_loc(get_Fnum("ss"), 1); get_ab(alpha2, beta2, beta_loc2);
  //  beta2[X_] = 1.0; beta2[Y_] = 1.0;
  beta_loc3 = get_loc(get_Fnum("ls"), 1); get_ab(alpha3, beta3, beta_loc3);
  //  beta3[X_] = 15.0; beta3[Y_] = 3.0;
  if (strcmp(lattice, "SLS") == 0)
    beta_loc4 = get_loc(get_Fnum("ms"), 1); get_ab(alpha4, beta4, beta_loc4);

  std::cout << std::endl;
  std::cout << std::fixed << std::setprecision(3)
       << "alpha1_x  = " << std::setw(6) << alpha1[X_]
       << ", alpha1_y = " << std::setw(6) << alpha1[Y_]
       << ", beta1_x = " << std::setw(6) << beta1[X_]
       << ", beta1_y  = " << std::setw(6) << beta1[Y_]
       << std::endl;
  std::cout << std::fixed << std::setprecision(3)
       << "alpha2_x  = " << std::setw(6) << alpha2[X_]
       << ", alpha2_y = " << std::setw(6) << alpha2[Y_]
       << ", beta2_x = " << std::setw(6) << beta2[X_]
       << ", beta2_y  = " << std::setw(6) << beta2[Y_]
       << std::endl;
  std::cout << std::fixed << std::setprecision(3)
       << "alpha3_x  = " << std::setw(6) << alpha3[X_]
       << ", alpha3_y = " << std::setw(6) << alpha3[Y_]
       << ", beta3_x = " << std::setw(6) << beta3[X_]
       << ", beta3_y  = " << std::setw(6) << beta3[Y_]
       << std::endl;
  if (strcmp(lattice, "SLS") == 0)
    std::cout << std::fixed << std::setprecision(3)
	 << "alpha4_x  = " << std::setw(6) << alpha4[X_]
	 << ", alpha4_y = " << std::setw(6) << alpha4[Y_]
	 << ", beta4_x = " << std::setw(6) << beta4[X_]
	 << ", beta4_y  = " << std::setw(6) << beta4[Y_]
	 << std::endl;
}


void chk_range(const int n_b2, const int b2_Fams[],
	       const double nu_x_min, const double nu_x_max,
	       const double nu_y_min, const double nu_y_max)
{

  if (strcmp(lattice, "SLS") == 0) {
    fit_tune_SLS(nu_x_min, nu_y_min,
		 beta1[X_], beta1[Y_], beta_loc1,
		 beta2[X_], beta2[Y_], beta_loc2,
		 beta3[X_], beta3[Y_], beta_loc3,
		 beta4[X_], beta4[Y_], beta_loc4,
		 n_b2, b2_Fams, eps_nu, true);
    fit_tune_SLS(nu_x_max, nu_y_max,
		 beta1[X_], beta1[Y_], beta_loc1,
		 beta2[X_], beta2[Y_], beta_loc2,
		 beta3[X_], beta3[Y_], beta_loc3,
		 beta4[X_], beta4[Y_], beta_loc4,
		 n_b2, b2_Fams, eps_nu, true);
  } else if (strcmp(lattice, "NSLS-II_BOOSTER") == 0) {
    fit_tune(nu_x_min, nu_y_min, n_b2, b2_Fams, eps_nu, true);
    fit_tune(nu_x_max, nu_y_max, n_b2, b2_Fams, eps_nu, true);
  } else {
    fit_tune(nu_x_min, nu_y_min,
	     beta1[X_], beta1[Y_], beta_loc1,
	     beta2[X_], beta2[Y_], beta_loc2,
	     beta3[X_], beta3[Y_], beta_loc3,
	     n_b2, b2_Fams, eps_nu, true);
    fit_tune(nu_x_max, nu_y_max,
	     beta1[X_], beta1[Y_], beta_loc1,
	     beta2[X_], beta2[Y_], beta_loc2,
	     beta3[X_], beta3[Y_], beta_loc3,
	     n_b2, b2_Fams, eps_nu, true);
  }
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


void get_dnu(const int n_nu, const double max_Ax, const double max_Ay)
{
  int              i, j, nx, ny;
  double           Ax, Ay, twoJx, twoJy, nux, nuy, A;
  ss_vect<double>  x_Floq;
  ss_vect<tps>     x,  A_inv, nus;
  std::ofstream         nus_out;

  file_wr(nus_out, "dnu_dAx.dat");
  x.zero();
  for (i = -n_nu; i <= n_nu; i++) {
    Ax = i*max_Ax/n_nu; x[x_] = Ax; x[px_] = 0.0;
    x_Floq = (A_inv*x).cst(); twoJx = sqr(x_Floq[0]) + sqr(x_Floq[1]);
    x[0] = sqrt(twoJx); x[1] = sqrt(twoJx);
    nux = (nus[3]*x).cst(); nuy = (nus[4]*x).cst();
    nus_out << std::fixed << std::setprecision(3)
	    << std::setw(7) << 1e3*Ax << " " << std::setw(7) << 1e6*twoJx/2.0
	    << std::setprecision(5)
	    << " " << std::setw(8) << nux << " " << std::setw(8) << nuy
	    << std::setprecision(3)
	    << " " << std::setw(8) << -(dHdJ(nus[0])[3]*x).cst()
	    << " " << std::setw(8) << -(dHdJ(nus[1])[4]*x).cst()
	    << std::endl;
  }
  nus_out.close();
  
  file_wr(nus_out, "dnu_dAy.dat");
  x.zero();
  for (i = -n_nu; i <= n_nu; i++) {
    Ay = i*max_Ay/n_nu; x[y_] = Ay; x[py_] = 0.0;
    x_Floq = (A_inv*x).cst(); twoJy = sqr(x_Floq[2]) + sqr(x_Floq[3]);
    x[2] = sqrt(twoJy); x[3] = sqrt(twoJy);
    nux = (nus[3]*x).cst(); nuy = (nus[4]*x).cst();
    nus_out << std::fixed << std::setprecision(3)
	    << std::setw(7) << 1e3*Ay << " " << std::setw(7) << 1e6*twoJy/2.0
	    << std::setprecision(5)
	    << " " << std::setw(8) << nux << " " << std::setw(8) << nuy
	    << std::setprecision(3)
	    << " " << std::setw(8) << -(dHdJ(nus[0])[3]*x).cst()
	    << " " << std::setw(8) << -(dHdJ(nus[1])[4]*x).cst()
	    << std::endl;
  }
  nus_out.close();
  
  file_wr(nus_out, "dnu_dA.dat");
  x.zero();
  for (i = -n_nu; i <= n_nu; i++) {
    Ax = i*max_Ax/n_nu;
    for (j = 0; j <= n_nu; j++) {
      Ay = j*max_Ay/n_nu;
      x[x_] = Ax; x[px_] = 0.0; x[y_] = Ay; x[py_] = 0.0;
      x_Floq = (A_inv*x).cst();
      twoJx = sqr(x_Floq[0]) + sqr(x_Floq[1]);
      twoJy = sqr(x_Floq[2]) + sqr(x_Floq[3]);
      x[0] = sqrt(twoJx); x[1] = sqrt(twoJx);
      x[2] = sqrt(twoJy); x[3] = sqrt(twoJy);
      nux = (nus[3]*x).cst(); nuy = (nus[4]*x).cst();
      get_ampl(nux, nuy, 6, nx, ny, A);
      nus_out << std::fixed << std::setprecision(3)
	      << std::setw(7) << 1e3*Ax << " " << std::setw(7) << 1e3*Ay
	      << " " << std::setw(7) << 1e6*twoJx/2.0
	      << " " << std::setw(7) << 1e6*twoJy/2.0
	      << std::setprecision(5)
	      << " " << std::setw(8) << nux << " " << std::setw(8) << nuy
	      << std::setprecision(3)
	      << " " << std::setw(8) << log(A) << std::setw(3) << nx << std::setw(3) << ny
	      << std::endl;
    }
    nus_out << std::endl;
  }
  nus_out.close();
}


void get_track(void)
{
  char             line[max_str];
  int              i, j;
  double           twoJx, twoJy;
  ss_vect<double>  x, x_Floq;
  ss_vect<tps>     A_inv;
  std::ifstream         inf;
  std::ofstream         J_out;

  file_rd(inf, "/home/bengtsson/tracy_JB/NSLS-II/track.out");
  file_wr(J_out, "J.dat");
  
  do
    inf.getline(line, max_str);
  while (strstr(line, "#") != NULL);
  
  // skip initial conditions
  inf.getline(line, max_str);
  do {
    sscanf(line, "%d %lf %lf %lf %lf %lf %lf",
	   &i, &x[x_], &x[px_], &x[y_], &x[py_], &x[delta_], &x[ct_]);
    
    for (j = 0; j <= 3; j++)
      x[j] *= 1e-3;
    x[delta_] = 0.0; x[ct_] = 0.0;
    
    x_Floq = (A_inv*x).cst();
    twoJx = sqr(x_Floq[0]) + sqr(x_Floq[1]);
    twoJy = sqr(x_Floq[2]) + sqr(x_Floq[3]);
    
    J_out << std::scientific << std::setprecision(3)
	  << std::setw(4) << i
	  << " " << std::setw(10) << 1e6*twoJx << " " << std::setw(10) << 1e6*twoJy
	  << std::endl;
  } while (inf.getline(line, max_str));
  inf.close();
  J_out.close();
}


float H_fun(float p[])
{
  char          hs[max_str];
  int           i, j, k, l, m, m1;
  double        H2, dH;
  tps           H, H_re, H_im, K_re, K_im;
  ss_vect<tps>  nus;

  const bool    use_K    = true;
  const int     n_prt    = 10, n_fold = 8;
  const double  scl_ksi1 = 1e4;

  n_iter++;

  for (i = 1; i <= n_prm; i++)
    if (fabs(p[i]) < bnL_max[Sext])
      set_bnL(prms[i-1], Sext, p[i]);
    else {
      H2 = 1e30;
      std::cout << std::scientific << std::setprecision(3)
	   << std::setw(4) << n_iter << std::setw(10) << H2 << std::endl;
      return H2;
    }

  danot_(no_tps-1); get_Map_N(n_cell);
  danot_(no_tps); K = MapNorm(Map, g, A1, A0, Map_res, 1);
  CtoR(K, K_re, K_im);
  K_re = K_re*Id_scl;

  // requires output from MapNorm
  H = get_H();
  //  H = get_h();
  CtoR(H, H_re, H_im);
  H_re = H_re*Id_scl;

  if (n_iter % n_prt == 0) std::cout << std::endl;
  H2 = 0.0; m1 = 0;
  for (i = 0; i <= no_tps-1; i++)
    for (j = 0; j <= no_tps-1; j++)
      for (k = 0; k <= no_tps-1; k++)
	for (l = 0; l <= no_tps-1; l++)
	  for (m = 0; m <= no_tps-1; m++) {
	    if ((0 < i+j+k+l+m) && (i+j+k+l+m <= no_tps-1) &&
		h[i][j][k][l][m] &&
		(fabs(h_ijklm(H_re, i, j, k, l, m)) > 0.0)) {
	      m1++;

	      dH = h_ijklm(H_re, i, j, k, l, m);

	      if (is_h_ijklm(1, 1, 0, 0, 1, i, j, k, l, m))
		// horizontal linear chromaticity
		dH = scl_ksi1*(n_cell*ksi1[X_]*M_PI*2.0*Jx*delta
			       +h_ijklm(H_re, i, j, k, l, m));
	      else if (is_h_ijklm(0, 0, 1, 1, 1, i, j, k, l, m))
		// vertical linear chromaticity
		dH = scl_ksi1*(n_cell*ksi1[Y_]*M_PI*2.0*Jy*delta
			       +h_ijklm(H_re, i, j, k, l, m));
	      else if (is_h_ijklm(2, 2, 0, 0, 0, i, j, k, l, m) ||
		       is_h_ijklm(0, 0, 2, 2, 0, i, j, k, l, m) ||
		       is_h_ijklm(1, 1, 1, 1, 0, i, j, k, l, m) ||
		       is_h_ijklm(3, 3, 0, 0, 0, i, j, k, l, m) ||
		       is_h_ijklm(2, 2, 1, 1, 0, i, j, k, l, m) ||
		       is_h_ijklm(1, 1, 2, 2, 0, i, j, k, l, m) ||
		       is_h_ijklm(0, 0, 3, 3, 0, i, j, k, l, m))
		// amplitude dependent tune shift
		if (!use_K)
		  dH *= scl_dnu;
		else
		  dH = scl_dnu*h_ijklm(K_re, i, j, k, l, m);
	      else if (is_h_ijklm(1, 1, 0, 0, 2, i, j, k, l, m) ||
		       is_h_ijklm(0, 0, 1, 1, 2, i, j, k, l, m) ||
		       is_h_ijklm(1, 1, 0, 0, 3, i, j, k, l, m) ||
		       is_h_ijklm(0, 0, 1, 1, 3, i, j, k, l, m) ||
		       is_h_ijklm(1, 1, 0, 0, 4, i, j, k, l, m) ||
		       is_h_ijklm(0, 0, 1, 1, 4, i, j, k, l, m) ||
		       is_h_ijklm(1, 1, 0, 0, 5, i, j, k, l, m) ||
		       is_h_ijklm(0, 0, 1, 1, 5, i, j, k, l, m))
		// nonlinear chromaticity
		if (!use_K)
		  dH *= scl_ksi_nl;
		else
		  dH = scl_ksi_nl*h_ijklm(K_re, i, j, k, l, m);
	      else if (is_h_ijklm(2, 2, 0, 0, 1, i, j, k, l, m) ||
		       is_h_ijklm(0, 0, 2, 2, 1, i, j, k, l, m) ||
		       is_h_ijklm(1, 1, 1, 1, 1, i, j, k, l, m) ||
		       is_h_ijklm(2, 2, 0, 0, 2, i, j, k, l, m) ||
		       is_h_ijklm(0, 0, 2, 2, 2, i, j, k, l, m) ||
		       is_h_ijklm(1, 1, 1, 1, 2, i, j, k, l, m) ||
		       is_h_ijklm(2, 2, 0, 0, 3, i, j, k, l, m) ||
		       is_h_ijklm(0, 0, 2, 2, 3, i, j, k, l, m) ||
		       is_h_ijklm(1, 1, 1, 1, 3, i, j, k, l, m) ||
		       is_h_ijklm(3, 3, 0, 0, 1, i, j, k, l, m) ||
		       is_h_ijklm(2, 2, 1, 1, 1, i, j, k, l, m) ||
		       is_h_ijklm(1, 1, 2, 2, 1, i, j, k, l, m) ||
		       is_h_ijklm(0, 0, 3, 3, 1, i, j, k, l, m))
		// amplitude dependent chromaticity
		if (!use_K) {
		  dH *= scl_ksi_nl;
		} else {
		  dH = scl_ksi_nl*h_ijklm(K_re, i, j, k, l, m);
		}

	      H2 += sqr(dH);

	      if (n_iter % n_prt == 0) {
		sprintf(hs, "h_%d%d%d%d%d", i, j, k, l, m);
		std::cout << std::setw(3) << m1 << " " << hs
		     << std::scientific << std::setprecision(2)
		     << std::setw(10) << dH << std::endl;
	      }
	    }
	  }

  if (n_iter % n_prt == 0) {
    std::cout << std::endl;
    std::cout << std::scientific << std::setprecision(3)
	 << std::setw(4) << n_iter << std::setw(10) << H2;

    for (i = 1; i <= n_prm; i++) {
      std::cout << std::fixed << std::setprecision(3) 
	   << std::setw(7) << get_bnL(prms[i-1], 1, bns[i-1]);
      if (i % n_fold == 0) std::cout << std::endl;
    }
    if (n_prm % n_fold != 0) std::cout << std::endl;

    sext_out << std::endl;
    sext_out << "n = " << n_iter << ":" << std::endl;
    for (i = 1; i <= n_prm; i++)
      for (j = 1; j <= get_n_Kids(prms[i-1]); j++)
	sext_out << std::fixed << std::setprecision(7)
		 << std::setw(6) << get_Name(prms[i-1]) << "(" << j << ") = "
		 << std::setw(11) << get_bnL(prms[i-1], 1, bns[i-1])
		 << std::setw(2) << bns[i-1] << std::endl;
  }

  return H2;
}


void H_min(void)
{
  int       i, j, iter;
  float     *p, **xi, fret;

  const float  ftol = 1e-8;
  p = vector(1, n_prm); xi = matrix(1, n_prm, 1, n_prm);

  for (i = 1; i <= n_prm; i++) {
    p[i] = get_bnL(prms[i-1], 1, Sext);
    for (j = 1; j <= n_prm; j++)
      xi[i][j] = (i == j)? 0.1 : 0.0;
  }

  std::cout << std::endl;
  select_h(); n_iter = 0;
  powell(p, xi, n_prm, ftol, &iter, &fret, H_fun);
  n_iter = -1; H_fun(p);

  free_vector(p, 1, n_prm); free_matrix(xi, 1, n_prm, 1, n_prm);
}


int main()
{
  char             line[max_str];
  long int         jj[ss_dim];
  int              i, j, h_rf;
  int              b2_Fams[n_prm_max], n_b2 = 0;
  long int         k1, k2;
  double           r, x_min[2], x_max[2], DA, nu[2], ksi[2], V_rf, f_rf;
  double           nu_x, nu_y, b3L;
  double           alpha[2], beta[2], A[2];
  Matrix           M;
  tps              Hnl, H2, gn, h, h_re, h_im, H_num, dH, H, H_re, H_im;
  tps              K_re, K_im, ab1[2], dnu1[2], ab2[2], dnu2[2], dnu[2];
  tps              g_re, g_im;
  ss_vect<double>  x;
  ss_vect<tps>     A_inv, ps, ps_re, ps_im, nus, R, R_inv, Map2;
  ss_vect<tps>     Map_Fl, Mn, Mk, J, eta, M1, M2, M3, M4, Mp;
  std::ofstream         outf, K_out, nus_out, A1_out, J_out;
  std::ifstream         inf;

  const int     seed = 1005, m = 6, n = 4;
  const int     n_nu = 25;

  danot_(no_tps-1);

  rad_on    = false; H_exact        = false; totpath_on   = false;
  cavity_on = false; quad_fringe_on = false; emittance_on = false;
  IBS_on    = false;


  rd_mfile("flat_file.dat", elem); rd_mfile("flat_file.dat", elem_tps);

  // initialize the symplectic integrator after energy has been defined
  ini_si();

  // disable from TPSALib- and LieLib log messages
  idprset(-1);

  if (false) {
    // change to 3 degrees of freedom, i.e., set ndpt_tps = 0 
    idprset(1);
    danot_(2);
    get_emittance();

    exit(0);
  }

  //  if (H_exact) bend_cal(get_Fnum("bd1"));

  danot_(3);

  if (true) chk_lat(nu, ksi);

  get_ab(alpha, beta, 0);
  std::cout << std::endl;
  std::cout << std::fixed << std::setprecision(3)
       << "alpha_x  = " << std::setw(6) << alpha[X_]
       << ", alpha_y = " << std::setw(6) << alpha[Y_]
       << ", beta_x = " << std::setw(6) << beta[X_]
       << ", beta_y  = " << std::setw(6) << beta[Y_] << std::endl;

  if (false &&
      (strcmp(lattice, "SLS") != 0) &&
      (strcmp(lattice, "NSLS-II_BOOSTER") != 0)) {
    get_ab(ab1, dnu1, get_loc(get_Fnum("mp"), 1));
    alpha[X_] = -h_ijklm(ab1[X_], 0, 1, 0, 0, 0);
    alpha[Y_] = -h_ijklm(ab1[Y_], 0, 0, 0, 1, 0);
    beta[X_]  =  h_ijklm(ab1[X_], 1, 0, 0, 0, 0);
    beta[Y_]  =  h_ijklm(ab1[Y_], 0, 0, 1, 0, 0);
    std::cout << std::fixed << std::setprecision(3)
	 << "alpha_x  = " << std::setw(6) << alpha[X_]
	 << ", alpha_y = " << std::setw(6) << alpha[Y_]
	 << ", beta_x = " << std::setw(6) << beta[X_]
	 << ", beta_y  = " << std::setw(6) << beta[Y_] << std::endl;
  }

  /*  Id.identity(); Map2.identity();
      Map2[x_]  += sqr(Id[x_]) + Id[x_]*Id[px_];
      Map2[px_] += Id[x_]*Id[px_] + sqr(Id[px_]);
      Map2[y_]  += sqr(Id[y_]) + Id[y_]*Id[py_];
      Map2[py_] += Id[y_]*Id[py_] + sqr(Id[py_]);

      Map2.propagate(1, n_elem);
      std::cout << Map2;*/

  if (false) {
    n_b2 = 0;
    b2_Fams[n_b2++] = get_Fnum("q1_jb");
    b2_Fams[n_b2++] = get_Fnum("q2_jb");

    k1 = get_loc(get_Fnum("q22"), 1); get_ab(alpha, beta, k1);
    k2 = get_loc(get_Fnum("ss"), 1);

    fit_beta(alpha[X_], beta[X_], alpha[Y_], beta[Y_],
	     1.0, 1.0, k1, k2, n_b2, b2_Fams, 0.3, true);

    b2_Fams[n_b2++] = get_Fnum("q3_jb");
    b2_Fams[n_b2++] = get_Fnum("q4_jb");
    n_b2 = 0;

    //    fit_alpha(alpha[X_], beta[X_], alpha[Y_], beta[Y_],
    //         k1, k2, n_b2, b2_Fams, 0.15, true);

    n_b2 = 0;
    b2_Fams[n_b2++] = get_Fnum("q1_jb");
    b2_Fams[n_b2++] = get_Fnum("q2_jb");
    b2_Fams[n_b2++] = get_Fnum("q3_jb");
    b2_Fams[n_b2++] = get_Fnum("q4_jb");
 
    fit_alpha_beta(alpha[X_], beta[X_], alpha[Y_], beta[Y_],
		   1.0, 1.0, k1, k2, n_b2, b2_Fams, 1e-4, true);

    exit(0);
  }

  if (true) prt_alphac();

  if (false) {
    // no scaling
    Jx = 0.5, Jy = 0.5, delta = 1.0;
  } else {
    Jx = sqr(max_Ax)/(2.0*beta1[X_]); Jy = sqr(max_Ay)/(2.0*beta1[Y_]);
    delta = max_delta;
    // RHIC
    //    Jx = sqr(0.5e-3)/(2.0*beta1[X_]); Jy = sqr(0.5e-3)/(2.0*beta1[Y_]);
    //    delta = 1e-2;
  }

  Id_scl.identity();
  // Id_scl[x_] *= sqrt(2.0*Jx); Id_scl[px_] *= sqrt(2.0*Jx);
  // Id_scl[y_] *= sqrt(2.0*Jy); Id_scl[py_] *= sqrt(2.0*Jy);
  // Id_scl[delta_] *= delta;

  if (false) { get_lin_map(); exit(0); }

  if (!false) {
    sxt_h1();
    sxt_h2();
    exit(0);
  }

  get_prm("tune_scan.prm");

  if (tune_scn) {
    if (adj_chrom) no_mpoles();

    if (strcmp(lattice, "NSLS-II_BOOSTER") != 0) get_locs();

    get_b2s(n_b2, b2_Fams);

    if (adj_tune) {
      if (strcmp(lattice, "SLS") == 0)
	fit_tune_SLS(nu0[X_], nu0[Y_],
		     beta1[X_], beta1[Y_], beta_loc1,
		     beta2[X_], beta2[Y_], beta_loc2,
		     beta3[X_], beta3[Y_], beta_loc3,
		     beta4[X_], beta4[Y_], beta_loc4,
		     n_b2, b2_Fams, eps_nu, true);
      else if (strcmp(lattice, "NSLS-II_BOOSTER") == 0)
	fit_tune(nu0[X_], nu0[Y_], n_b2, b2_Fams, eps_nu, true);
      else
	fit_tune(nu0[X_], nu0[Y_],
		 beta1[X_], beta1[Y_], beta_loc1,
		 beta2[X_], beta2[Y_], beta_loc2,
		 beta3[X_], beta3[Y_], beta_loc3,
		 n_b2, b2_Fams, eps_nu, true);

      //      exit(0);
    }

    // check tuning range (linear stability)
    if (check_range) {
      chk_range(n_b2, b2_Fams, nu_x_min, nu_x_max, nu_y_min, nu_y_max);
      exit(0);
    }

    if (adj_chrom) {
      n_prm = 0;
      prms[n_prm] = get_Fnum("sm1");  bns[n_prm++] = Sext;
      prms[n_prm] = get_Fnum("sm2"); bns[n_prm++] = Sext;

      danot_(3);
      fit_chrom();

      get_prms();
    }

    cavity_on = true;

    //    tune_scan(n_b2, b2_Fams, true);
    tune_scan(true);
  } else {
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
 
    //    get_fixed_points(50e-3, 50e-3, 25, 15, true);
 
    if (false) {
      track_H("track_H_1.dat",  1.0e-3,  1.0e-3);
      track_H("track_H_2.dat",  2.5e-3,  2.5e-3);
      track_H("track_H_3.dat",  5.0e-3,  3.0e-3);
      track_H("track_H_4.dat",  7.5e-3,  4.0e-3);
      track_H("track_H_5.dat", 10.0e-3,  4.0e-3);
      track_H("track_H_6.dat", 15.0e-3,  4.0e-3);
      //      track_H("track_H_7.dat", 12.5e-3,  8.0e-3);
      //      track_H("track_H_8.dat", 19.9e-3,  8.0e-3);

      exit(0);
    }

    //    H_min();
    
    //    sext_out.close();

    if (false)
      DA = get_dynap(10e-3, 0.0, n_track, 0.1e-3, n_aper, x_min, x_max, false);

    if (true) {
      if (adj_tune) {
	get_locs();

	get_b2s(n_b2, b2_Fams);

	fit_tune(nu0[X_], nu0[Y_],
		 beta1[X_], beta1[Y_], beta_loc1,
		 beta2[X_], beta2[Y_], beta_loc2,
		 beta3[X_], beta3[Y_], beta_loc3,
		 n_b2, b2_Fams, eps_nu, true);
      }

      if (adj_chrom) {
	n_prm = 0;
	if (strcmp(lattice, "NSLS-II_DBA") == 0) {
	  prms[n_prm] = get_Fnum("sm1a");  bns[n_prm++] = Sext;
	  prms[n_prm] = get_Fnum("sm1b");  bns[n_prm++] = Sext;
	  prms[n_prm] = get_Fnum("sm2"); bns[n_prm++] = Sext;
	} else if (strcmp(lattice, "DIAMOND") == 0) {
	  prms[n_prm] = get_Fnum("ss1a"); bns[n_prm++] = Sext;
	  prms[n_prm] = get_Fnum("ss2a"); bns[n_prm++] = Sext;
	  prms[n_prm] = get_Fnum("ss1b"); bns[n_prm++] = Sext;
	  prms[n_prm] = get_Fnum("ss2b"); bns[n_prm++] = Sext;
	  prms[n_prm] = get_Fnum("ss1c"); bns[n_prm++] = Sext;
	  prms[n_prm] = get_Fnum("ss2c"); bns[n_prm++] = Sext;
	  prms[n_prm] = get_Fnum("ss1d"); bns[n_prm++] = Sext;
	  prms[n_prm] = get_Fnum("ss2d"); bns[n_prm++] = Sext;
	}

	if (fit_chrm) {
	  danot_(3);
	  no_mpoles(); fit_chrom();
	}

	get_prms();

	if (false) {
	  // Random starting point for the non-chromatic sextupoles
	  std::cout << std::endl;
	  std::cout << "Random starting point for the non-chromatic sextupoles:"
	       << std::endl;
	  for (i = 0; i < 3; i++) {
	    b3L = 2.0*bnL_max[Sext]*(ranf()-0.5);
	    set_bnL(prms[i], bns[i], b3L);
	    std::cout << std::scientific << std::setprecision(3) << std::setw(11) << b3L;
	  }
	  for (i = 6; i < n_prm; i++) {
	    b3L = 2.0*bnL_max[Sext]*(ranf()-0.5);
	    set_bnL(prms[i], bns[i], b3L);
	    std::cout << std::scientific << std::setprecision(3) << std::setw(11) << b3L;
	  }
	  std::cout << std::endl;
	}

	H_zero(eps_ksi, 25, true);
	//   H_min();
      }

      //      danot_(no_tps-1); get_COD(10, 1e-10, 0.0, true);
      //      file_wr(outf, "map.dat"); outf << Map; outf.close();

      danot_(no_tps-1);

      //      if (true)
      get_Map();
      //      else
      //   get_Map_N(n_cell);

      //      get_Map_N(15);
      //      get_Map_N(3);

      file_wr(outf, "map.dat"); outf << Map; outf.close();

      danot_(no_tps);

      if (true) Id_scl.identity();

      CtoR(get_h()*Id_scl, h_re, h_im);
      file_wr(outf, "h.dat"); outf << h_re; outf.close();

      K = MapNorm(Map, g, A1, A0, Map_res, no_tps); CtoR(K*Id_scl, K_re, K_im);
      file_wr(outf, "K.dat"); outf << K_re; outf.close();

      nus = dHdJ(K);
      file_wr(outf, "nus.dat"); outf << nus[3] << nus[4]; outf.close();

      CtoR(get_H()*Id_scl, H_re, H_im);
      file_wr(outf, "H.dat"); outf << H_re; outf.close();

      CtoR(g*Id_scl, g_re, g_im);
      file_wr(outf, "g.dat"); outf << g_im; outf.close();
    }

    if (false) {
      danot_(no_tps-1); get_Map();
      danot_(no_tps); K = MapNorm(Map, g, A1, A0, Map_res, 1);
      H = get_H();
      CtoR(H, H_re, H_im);
      H_re = H_re*Id_scl;
      file_wr(outf, "H.dat"); outf << H_re; outf.close();

      /*      // single Lie exponent, i.e. not Dragt-Finn factored
	      h = get_h();
	      Id.identity();
	      Map_Fl = LieExp(h, Id)*R;
	      danot_(no_tps-1);
	      Map_Fl = Inv(A0*A1)*Map*A0*A1;
	      danot_(1); R = Map_Fl;
	      danot_(no_tps);
	      std::cout << h*R-h*Map_Fl;
	      std::cout << std::endl;
	      std::cout << std::scientific << std::setprecision(14) << std::setw(22) << abs2(h) << std::endl;
	      std::cout << std::scientific << std::setprecision(14)
	      << std::setw(22) << abs2(h*Map_Fl) << std::endl;
	      std::cout << std::scientific << std::setprecision(14)
	      << std::setw(22) << abs2(h*Map_Fl*Map_Fl) << std::endl;*/

      //      h = h*Id;
      //      CtoR(h, h_re, h_im);
      //      file_wr(outf, "h.dat"); outf << h_re << h_im; outf.close();

      exit(0);
    }

    if (false) {
      danot_(no_tps-1); get_Map();
      danot_(no_tps);
      K = MapNorm(Map, g, A1, A0, Map_res, 1); CtoR(K, K_re, K_im);
      nus = dHdJ(K); get_nu_ksi(nus, nu, ksi);
      //      K_re = K_re*Id; K_im = K_im*Id;

      file_wr(A1_out, "A1.dat"); A1_out << A1; A1_out.close();

      file_wr(K_out, "K.dat"); K_out << K_re << K_im; K_out.close();

      file_wr(nus_out, "nus.dat"); nus_out << nus[3] << nus[4];
      nus_out.close();

      A1_inv = Inv(A1); A_inv = get_A_inv();

      h = 0.5*tps(0.0, 1) + 2.0*tps(0.0, 2);
      CtoR(h, h_re, h_im);
      //      std::cout << h << h_re << h_im;
      exit(0);
    }

    if (false) get_track();

    if (false) { get_dnu(n_nu, max_Ax, max_Ay); }

    //    DA = get_dynap(10e-3, 0.0, n_track, 0.1e-3, n_aper, x_min, x_max, false); 
  }
  
  //  danot_(no_tps-1); get_Map(); danot_(no_tps);
  //  Map2 = Map*Map; Map = Map2*Map2*Map;
  // requires map normal form
  //  H = get_H();
  //  prt_h(H, pow(30e-3, 2)/(2.0*beta[X_]), 0e-2, 10);

  if (false) {
    prt_alphac();

    //    prt_H_long(10, 180*pi/180.0, 4e-2, -932.5e3);

    // MAX-IV
    //    prt_H_long(10, 180*pi/180.0, 6.2e-2, -359.4e3);

    // DIAMOND
    prt_H_long(10, 180*M_PI/180.0, 11e-2, -1005.2e-3, false);
  }

  exit(0);

  //  prt_H_long(10, 180*pi/180.0, 5e-2, 0.0);
  //  prt_H_long(10, 150*pi/180.0, 4e-2, -26.7*pi/180.0);
  //  prt_H_long(10, 200*pi/180.0, 10e-2, 0.0);
  //  prt_H_long(10, 150*pi/180.0, 4e-2, -28.3*pi/180.0);

  //  get_cav("cav", h_rf, V_rf, f_rf);
  //  track(0e-3, 0e-3, 0e-3, 0e-3, 2.4e-2, n_track, f_rf, true);

  //  get_dynap(5e-3, 0.0, n_track, 0.1e-3, n_aper, 0.0, 0.0,
  //       x_min, x_max, false);

  //  get_dynap(5e-3, 0.0, n_track, 0.1e-3, n_aper, 0.491, 0.604,
  //       x_min, x_max, true);

  if (false) {
    danot_(1); R = Map; danot_(no_tps); h = LieFact(Map*Inv(R));
    for (i = 0; i < ss_dim; i++)
      jj[i] = 0;
    jj[x_] = 2; jj[px_] = 2; h.pook(jj, 0.0);
    jj[x_] = 1; jj[px_] = 1; jj[y_] = 1; jj[py_] = 1; h.pook(jj, 0.0);
    jj[x_] = 0; jj[px_] = 0; jj[y_] = 2; jj[py_] = 2; h.pook(jj, 0.0);
    //    Id.identity(); Map = LieExp(h, Id)*R;
    danot_(no_tps-1); Map = Map; danot_(no_tps);
  }

  // requires map normal form
  //  H = get_H();
  //  std::cout << H*Id;

  //  get_dynap(5e-3, 0.0, n_track, 0.1e-3, n_aper, 0.491, 0.604,
  //       x_min, x_max, true);

  //  std::cout << get_h()*Id;

  //  danot_(3);
  // requires map normal form
  //  H = get_H();
  //  std::cout << H;

  exit(0);

  std::cout << "O.K.? "; std::cin.get();  

  //  std::cout << K; prt_dnu(nus);

  //  prt_h(H, 20e-3, 30e-3, 10);


  get_Map();
  // requires map normal form
  //  H = get_H();
  std::cout << H;

  /*  gn = Take(g, 4);
      for (i = 5; i <= no_tps; i++) {
      pq = Take(g, i); gn = BCH(gn, pq, 3);
      }
      std::cout << (H-K*LieExp(-gn, I()));*/

  /*  // Transform map to Floquet space
      Map = Inv(A1*A2)*Map*A1*A2;
      // Factor out the linear part.
      danot_(1); R = Map; danot_(no_tps-1); Map2 = Map*Inv(R);
      // Extract the Lie generator.
      danot_(no_tps); Hnl = LieFact(Map2, eps);*/

  /*  get_matrix(R, M); prt_matrix(M);
      nu_x = GetAngle(M[0][0], M[0][1])/(2.0*pi);
      nu_y = GetAngle(M[2][2], M[2][3])/(2.0*pi);
      std::cout << std::endl;
      std::cout << std::fixed << std::setprecision(5)
      << "nux = " << std::setw(7) << nu_x
      << " nuy = " << std::setw(7) << nu_y << std::endl;

      ps.identity();
      H2 = - 2.0*pi*nu_x*(sqr(ps[x_])+sqr(ps[px_]))/2.0
      - 2.0*pi*nu_y*(sqr(ps[y_])+sqr(ps[py_]))/2.0;*/

  /*  P = H2 + Hnl; Q = H2 - Hnl;
      pq = PB(P, Q);
      pqq = PB(pq, Q); pqp = PB(pq, P);
      pqpp = PB(pqp, P); pqqq = PB(pqq, Q); pqqp = PB(pqq, P);
      H = P - pq/4.0 + pqq/24.0 + (pqpp-pqqq)/192.0;
      // 5th order appears to be incorrect
      //  H = P - pq/4.0 + pqq/24.0 + (pqpp-pqqq)/192.0
      //      + (8.0*PB(pqqp, P)-15.0*PB(pqpp, Q)+3.0*PB(pqqq, Q))/5760.0;
      H_num = BCH(H2, Hnl, 4);
      std::cout << H_num;
      std::cout << (H_num-H);*/
  // Numerical search for the invariant
  /*  Id.identity(); Id[delta_] = tps(0.0, 0); Id[ct_] = tps(0.0, 0);
      H_num = H;
      for (i = 1; i <= 100; i++) {
      dH = H_num - H_num*Map; dH = dH*Inv(Map-Id);
      H_num = H_num + dH;
      }
      //  std::cout << H_num;
      std::cout << dH;*/

  // check the Lie generator
  //  danot_(no_tps); Map2 = LieExp(H_num, I());
  //  danot_(no_tps-1); TPSAEps(1e-9);
  //  std::cout << (Map*Inv(Map2)); TPSAEps(eps);
}

