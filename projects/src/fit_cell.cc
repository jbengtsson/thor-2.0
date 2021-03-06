#define ORDER 7

#include "thor_lib.h"

int no_tps = ORDER;


const int  max_Fnum = 50, max_ind = 25;

int       Fnum[max_Fnum], Fnum1[max_Fnum], ind[max_ind], n_prm, n_ind, n_iter;
double    alpha_rad_x, nu[2];


void get_emit(const long int i1, const long int i2)
{
  int           k;
  ss_vect<tps>  A;

  emittance_on = true; rad_on = true;

  get_COD(10, 1e-10, 0.0, true);

  danot_(1);
  // compute A and damping rate
  cout << endl;
  cout << "get_emit: computing damping rates" << endl;
  K = MapNorm(Map, g, A1, A0, Map_res, 1);

  // compute imaginary part of eigenvalues
  gettura_(nu_, rad_);

  // compute diffusion coeffs
  A = A1; A.propagate(i1, i2);

  for (k = 0; k <= 2; k++)
    eps_[k] = -D_[k]/(4.0*rad_[k]);

  cout << endl;
  cout << scientific << setprecision(3)
       << "eps_x = " << eps_[X_].cst() << endl;
}


void get_dnu(const ss_vect<tps> &A, double dnu[])
{   
  int  k;

  for (k = 0; k <= 1; k++) {
    dnu[k] = atan2(A[2*k][2*k+1], A[2*k][2*k])/(2.0*M_PI);
    if (dnu[k] < 0.0) dnu[k] += 1.0;
  }
} 


void get_ab(const ss_vect<tps> A1, tps ab[])
{
  /* parameter dependance for the periodic solution:


        -1       T           |  0  I |        T   | beta   -alpha |
       A   = -S A  S,    S = |       |,    A A  = |               |
                             | -I  0 |            | -alpha  gamma |

       alpha_x = -h_01000(A_Atp[x_]), alpha_y = -h_00010(A_Atp[y_]),
       beta_x  =  h_10000(A_Atp[x_]), beta_y  =  h_00100(A_Atp[y_])

       eta_x   =  h_00001(A[x_]),     eta_y   =  h_00001(A[y_]),
       eta'_x  =  h_00001(A[px_]),    eta'_y  =  h_00001(A[py_])

  */

  ss_vect<tps>  A1_A1tp;

  A1_A1tp = A1*tp_S(2, A1); ab[X_] = A1_A1tp[x_]; ab[Y_] = A1_A1tp[y_];
}


void calc_twiss(void)
{
  long int  j;
  int       k;
  double    dnu0[2], dnu[2];
  tps       ab[2];

  elem[0].Nu[X_] = 0.0; elem[0].Nu[Y_] = 0.0;
  for (j = 0; j < n_elem; j++) {
    for (k = 0; k <= 1; k++) {
      elem[j].Eta[k]  = h_ijklm(elem_tps[j].A1[2*k],  0, 0, 0, 0, 1);
      elem[j].Etap[k] = h_ijklm(elem_tps[j].A1[2*k+1], 0, 0, 0, 0, 1);
    }

    get_ab(elem_tps[j].A1, ab);

    elem[j].Alpha[X_] = -h_ijklm(ab[X_], 0, 1, 0, 0, 0);
    elem[j].Alpha[Y_] = -h_ijklm(ab[Y_], 0, 0, 0, 1, 0);
    elem[j].Beta[X_]  =  h_ijklm(ab[X_], 1, 0, 0, 0, 0);
    elem[j].Beta[Y_]  =  h_ijklm(ab[Y_], 0, 0, 1, 0, 0);

    if (j > 0) {
      get_dnu(elem_tps[j-1].A1, dnu0); get_dnu(elem_tps[j].A1, dnu);
      for (k = 0; k <= 1; k++) {
	dnu[k] -= dnu0[k];
	if ((dnu[k] < 0.0) && (elem[j].L > 0.0)) dnu[k] += 1.0;
	elem[j].Nu[k] = elem[j-1].Nu[k] + dnu[k];
      }

      elem[j].S = elem[j-1].S + elem[j].L;
    }
  }
}


void get_twiss(const long int i1, const long int i2)
{
  // computes linear dispersion for delta = 0
  ss_vect<tps>  A1_prm;

  emittance_on = true; rad_on = false;

  danot_(no_tps-1); Map.identity(); Map.propagate(i1, i2);

  danot_(no_tps); K = MapNorm(Map, g, A1, A0, Map_res, 1);

  nus_ = dHdJ(K); get_nu_ksi(nus_, nu_, ksi_);

  // compute A and diffusion coeffs
  A1_prm = LieExp(g, A0*A1); A1_prm.propagate(i1, i2);

  // radiation damping is constant
  eps_[X_] = -D_[X_]/(4.0*alpha_rad_x);
}


inline double le(const double x, const double x_max)
{

  return (fabs(x) <= fabs(x_max))? 0.0 : sqr(x-x_max);
}


void no_mpoles(void)
{
  int j, k;

  cout << endl;
  cout << "zeroing multipoles" << endl;
  for (j = 0; j < n_elem; j++)
    if (elem[j].kind == Mpole)
      for (k = Sext; k < mpole_max; k++)
        set_bn(elem[j].Fnum, elem[j].Knum, k, 0.0);
}


void get_ab_ms(const long int i1, const long int i2,
	       const double nu_x, const double nu_y,
	       double alpha[], double beta[])
{
  // compute alpha and nu from a mirror symmetric section for a given beta
  int     k;
  double  m11[2], m12[2], s, c, nu[2];

  danot_(1); Map.identity(); Map.propagate(i1, i2);

  m11[X_] = h_ijklm(Map[x_], 1, 0, 0, 0, 0);
  m12[X_] = h_ijklm(Map[x_], 0, 1, 0, 0, 0);
  m11[Y_] = h_ijklm(Map[y_], 0, 0, 1, 0, 0);
  m12[Y_] = h_ijklm(Map[y_], 0, 0, 0, 1, 0);

  nu[X_] = nu_x; nu[Y_] = nu_y;
  for (k = 0; k <= 1; k++) {
    s = sin(2.0*M_PI*nu[k]); c = sqrt(1.0-sqr(s));
    alpha[k] = (m11[k]-c)/s; beta[k] = m12[k]/s;
  }

  cout << fixed << setprecision(3)
       << "alpha = (" << setw(5) << alpha[X_]
       << ", " << setw(5) << alpha[Y_] << ")"
       << " beta = (" << setw(5) << beta[X_]
       << ", " << setw(5) << beta[Y_] << ")" << endl;
}


void prt_ab(const char *name, const tps ab[], const bool eoln)
{
  char  str[max_str];

  sprintf(str, "%-6s", name);

  cout << scientific << setprecision(1)
       << "alpha" << str << " = ("
       << setw(8) << -h_ijklm(ab[X_], 0, 1, 0, 0, 0) << ", "
       << setw(8) << -h_ijklm(ab[Y_], 0, 0, 0, 1, 0) << ")"
       << fixed << setprecision(3)
       << " beta" << str << " = ("
       << setw(6) << h_ijklm(ab[X_], 1, 0, 0, 0, 0) << ", "
       << setw(6) << h_ijklm(ab[Y_], 0, 0, 1, 0, 0) << ")";

  if (eoln) cout << endl;
}



void fit_prms(const bool fit_b2)
{
  n_prm = 0;

  if (fit_b2) {
    Fnum[n_prm++] = get_Fnum("qm1");
    Fnum[n_prm++] = get_Fnum("qm2");

    Fnum[n_prm++] = get_Fnum("ql1"); Fnum[n_prm++] = get_Fnum("ql2");
    Fnum[n_prm++] = get_Fnum("ql3");

    Fnum[n_prm++] = get_Fnum("qh1"); Fnum[n_prm++] = get_Fnum("qh2");
    Fnum[n_prm++] = get_Fnum("qh3");

//    Fnum[n_prm++] =  get_Fnum("q_dw");
  } //else {
    Fnum[n_prm] = -get_Fnum("du_qm1"); Fnum1[n_prm++] = -get_Fnum("dl_qm1");
    Fnum[n_prm] = -get_Fnum("du_qm2"); Fnum1[n_prm++] = -get_Fnum("dl_qm2");
//    Fnum[n_prm] = -get_Fnum("du_b1");  Fnum1[n_prm++] = -get_Fnum("dl_b1");

    Fnum[n_prm] = -get_Fnum("du_ql1"); Fnum1[n_prm++] = -get_Fnum("dl_ql1");
    Fnum[n_prm] = -get_Fnum("du_ql2"); Fnum1[n_prm++] = -get_Fnum("dl_ql2");
    Fnum[n_prm] = -get_Fnum("du_ql3"); Fnum1[n_prm++] = -get_Fnum("dl_ql3");

    Fnum[n_prm] = -get_Fnum("du_qh1"); Fnum1[n_prm++] = -get_Fnum("dl_qh1");
    Fnum[n_prm] = -get_Fnum("du_qh2"); Fnum1[n_prm++] = -get_Fnum("dl_qh2");
    Fnum[n_prm] = -get_Fnum("du_qh3"); Fnum1[n_prm++] = -get_Fnum("dl_qh3");
//  }

  if (n_prm > n_prm_max) {
    cout << "fit_prms: n_prm_max exceeded " << n_prm << "(" << n_prm_max
	 << ")" << endl;
    exit(1);
  }
}


void fit_cell(const double nu_x, const double nu_y, const bool fit_tune,
	      const double prm_tol, const bool fit_b2)
{

  int           i, j;
  long int      k;
  double        **A, *b, *prm_lim, *prm, *dprm, step, nu_fract[2], dprm_max;
  tps           ab[max_ind][2];
  ofstream      quad_out;

  const int     m          = 21;
  const double  s_cut      = 1e-3,  step0      = 0.5;
  const double  scl_eps_x  = 1e12;
  const double  scl_nu     = (fit_tune)? 1e6 : 0.0;
  const double  scl_ksi_x  = 1e3,   scl_ksi_y  = 1e3;
  const double  scl_beta   = 5e2,   scl_eta    = 1e6;
  const double  scl_beta_x = 0*1e0,   scl_beta_y = 0*1e0, scl_alpha = 1e5;

  const int     n_fold      = 8;
  const double  eps_x       = 2.0e-9;
  const double  ksi_max[]   = { -6.5, -2.4 };
  const double  beta_mp[]   = { 20.0, 10.0 };

  const double  beta_qm1[]  = {  5.0, 18.0 }; // beta_y is tunable

  const double  beta_ss[]   = {  1.0,  1.0 };
  const double  beta_ql2[]  = { 20.0, 10.0 };
  const double  beta_ql3[]  = { 10.0, 25.0 };

  const double  beta_ls[]   = { 15.0,  3.0 };
  const double  beta_qh2[]  = { 25.0, 10.0 };
  const double  beta_qh3[]  = { 10.0, 25.0 };

  b2_max = 2.3; ds_max = 0.4;

  nu_fract[X_] = fract(nu_x); nu_fract[Y_] = fract(nu_y);

  cout << endl;
  cout << fixed << setprecision(5)
       << "fit_cell: nu_x = " << nu_fract[X_] << ", nu_y = " << nu_fract[Y_]
       << endl;

  // => n_prm
  fit_prms(fit_b2);

  b = dvector(1, m); prm_lim = dvector(1, n_prm);
  prm = dvector(1, n_prm); dprm = dvector(1, n_prm);
  A = dmatrix(1, m, 1, n_prm);

  for (i = 1; i <= n_prm; i++) {
    if (Fnum[i-1] > 0) {
      prm_lim[i] = b2_max; prm[i] = get_bn(Fnum[i-1], 1, Quad);
    } else {
      prm_lim[i] = ds_max; k = get_loc(abs(Fnum[i-1]), 1) - 1;
      prm[i] = get_L(abs(Fnum[i-1]), 1);
    }
  }

  n_ind = 0;
  ind[n_ind++] = get_loc(get_Fnum("sm2"), 1) - 1;
  ind[n_ind++] = get_loc(get_Fnum("b1"), 1) - 2; // upstream
  ind[n_ind++] = get_loc(get_Fnum("b1"), 2) - 1;
  ind[n_ind++] = get_loc(get_Fnum("qm1"), 1) - 1;

  ind[n_ind++] = get_loc(get_Fnum("ss"), 1) - 1;
  ind[n_ind++] = get_loc(get_Fnum("ql2"), 1) - 2; // upstream
  ind[n_ind++] = get_loc(get_Fnum("ql3"), 1) - 2; // upstream

  ind[n_ind++] = get_loc(get_Fnum("ls"), 1) - 1;
  ind[n_ind++] = get_loc(get_Fnum("qh2"), 1) - 1;
  ind[n_ind++] = get_loc(get_Fnum("qh3"), 1) - 1;

  do {
    n_iter++; step = step0;
    for (i = 1; i <= n_prm; i++) {
      if (Fnum[i-1] > 0)
	set_bn_par(Fnum[i-1], Quad, 7);
      else
	set_s_pair(abs(Fnum[i-1]), abs(Fnum1[i-1]), 7);

      get_twiss(1, n_elem);

      for (j = 0; j < n_ind; j++)
	get_ab(elem_tps[ind[j]].A1, ab[j]);

      eps_[X_] = -D_[X_]/(4.0*alpha_rad_x);

      j = 0;
      // eps_x
      A[++j][i] = scl_eps_x*h_ijklm_p(eps_[X_], 0, 0, 0, 0, 0, 7);
      // nu
      A[++j][i] = scl_nu*h_ijklm_p(nus_[3], 0, 0, 0, 0, 0, 7);
      A[++j][i] = scl_nu*h_ijklm_p(nus_[4], 0, 0, 0, 0, 0, 7);
      // ksi
      A[++j][i] = scl_ksi_x*h_ijklm_p(nus_[3], 0, 0, 0, 0, 1, 7);
      A[++j][i] = scl_ksi_y*h_ijklm_p(nus_[4], 0, 0, 0, 0, 1, 7);

      A[++j][i] = -scl_alpha*h_ijklm_p(ab[0][X_], 0, 1, 0, 0, 0, 7);
      A[++j][i] = -scl_alpha*h_ijklm_p(ab[0][Y_], 0, 0, 0, 1, 0, 7);
      // beta_x at the center of the DBA
      A[++j][i] = scl_beta_x*h_ijklm_p(ab[0][X_], 1, 0, 0, 0, 0, 7);

      // dispersion at the DBA entrance
      A[++j][i]
	= scl_eta*h_ijklm_p(elem_tps[ind[1]].A1[x_],  0, 0, 0, 0, 1, 7);
      A[++j][i]
	= scl_eta*h_ijklm_p(elem_tps[ind[1]].A1[px_], 0, 0, 0, 0, 1, 7);
      // dispersion at the DBA exit (not automatically mirror symmetric)
      A[++j][i]
	= scl_eta*h_ijklm_p(elem_tps[ind[2]].A1[x_],  0, 0, 0, 0, 1, 7);
      A[++j][i]
	= scl_eta*h_ijklm_p(elem_tps[ind[2]].A1[px_], 0, 0, 0, 0, 1, 7);

      // beta_y at qm1
      A[++j][i] = scl_beta_y*h_ijklm_p(ab[3][Y_], 0, 0, 1, 0, 0, 7);

      // beta at the center of the short straigth
      A[++j][i] = scl_beta*h_ijklm_p(ab[4][X_], 1, 0, 0, 0, 0, 7);
      A[++j][i] = scl_beta*h_ijklm_p(ab[4][Y_], 0, 0, 1, 0, 0, 7);
      // beta in short matching section
      A[++j][i] = scl_beta_x*h_ijklm_p(ab[5][X_], 1, 0, 0, 0, 0, 7);
      A[++j][i] = scl_beta_y*h_ijklm_p(ab[6][Y_], 0, 0, 1, 0, 0, 7);

      // beta at the center of the long straigth
      A[++j][i] = scl_beta*h_ijklm_p(ab[7][X_], 1, 0, 0, 0, 0, 7);
      A[++j][i] = scl_beta*h_ijklm_p(ab[7][Y_], 0, 0, 1, 0, 0, 7);
      // beta in the long matching section
      A[++j][i] = scl_beta_x*h_ijklm_p(ab[8][X_], 1, 0, 0, 0, 0, 7);
      A[++j][i] = scl_beta_y*h_ijklm_p(ab[9][Y_], 0, 0, 1, 0, 0, 7);

      if (Fnum[i-1] > 0)
	clr_bn_par(Fnum[i-1], Quad);
      else
	clr_s_pair(abs(Fnum[i-1]), abs(Fnum1[i-1]));
    }

    j = 0;
    b[++j] = -scl_eps_x*(eps_[X_].cst()-eps_x);
    b[++j] = -scl_nu*(nus_[3].cst()-nu_fract[X_]);
    b[++j] = -scl_nu*(nus_[4].cst()-nu_fract[Y_]);
    b[++j] = -scl_ksi_x*(h_ijklm(nus_[3], 0, 0, 0, 0, 1)-ksi_max[X_]);
    b[++j] = -scl_ksi_y*(h_ijklm(nus_[4], 0, 0, 0, 0, 1)-ksi_max[Y_]);

    b[++j] = scl_alpha*h_ijklm(ab[0][X_], 0, 1, 0, 0, 0);
    b[++j] = scl_alpha*h_ijklm(ab[0][Y_], 0, 0, 0, 1, 0);
    b[++j] = -scl_beta_x*(h_ijklm(ab[0][X_], 1, 0, 0, 0, 0)-beta_mp[X_]);

    b[++j] = -scl_eta*h_ijklm(elem_tps[ind[1]].A1[x_],  0, 0, 0, 0, 1);
    b[++j] = -scl_eta*h_ijklm(elem_tps[ind[1]].A1[px_], 0, 0, 0, 0, 1);
    b[++j] = -scl_eta*h_ijklm(elem_tps[ind[2]].A1[x_],  0, 0, 0, 0, 1);
    b[++j] = -scl_eta*h_ijklm(elem_tps[ind[2]].A1[px_], 0, 0, 0, 0, 1);

    b[++j] = -scl_beta*(h_ijklm(ab[3][Y_], 0, 0, 1, 0, 0)-beta_qm1[Y_]);

    b[++j] = -scl_beta*(h_ijklm(ab[4][X_], 1, 0, 0, 0, 0)-beta_ss[X_]);
    b[++j] = -scl_beta*(h_ijklm(ab[4][Y_], 0, 0, 1, 0, 0)-beta_ss[Y_]);
 
    b[++j] = -scl_beta_x*(h_ijklm(ab[5][X_], 1, 0, 0, 0, 0)-beta_ql2[X_]);
    b[++j] = -scl_beta_y*(h_ijklm(ab[6][Y_], 0, 0, 1, 0, 0)-beta_ql3[Y_]);
 
    b[++j] = -scl_beta*(h_ijklm(ab[7][X_], 1, 0, 0, 0, 0)-beta_ls[X_]);
    b[++j] = -scl_beta*(h_ijklm(ab[7][Y_], 0, 0, 1, 0, 0)-beta_ls[Y_]);
 
    b[++j] = -scl_beta_x*(h_ijklm(ab[8][X_], 1, 0, 0, 0, 0)-beta_qh2[X_]);
    b[++j] = -scl_beta_y*(h_ijklm(ab[9][Y_], 0, 0, 1, 0, 0)-beta_qh3[Y_]);
 
    cout << endl;
    cout << " Ax = b:" << endl;
    for (j = 1; j <= n_prm; j++)
      if (j == 1)
	cout << setw(8) << j;
      else
	cout << setw(11) << j;
    cout << endl;
    for (i = 1; i <= m; i++) {
      cout << setw(2) << i;
      for (j = 1; j <= n_prm; j++)
	cout << scientific << setprecision(3) << setw(11) << A[i][j];
      cout << scientific << setprecision(3) << setw(11) << b[i] << endl;
    }

    SVD_lim(m, n_prm, A, b, prm_lim, s_cut, prm, dprm);

    cout << endl;
    dprm_max = 0.0;
    for (i = 1; i <= n_prm; i++) {
      if (Fnum[i-1] > 0) {
	set_dbn(Fnum[i-1], Quad, step*dprm[i]);
	prm_lim[i] = b2_max; prm[i] = get_bn(Fnum[i-1], 1, Quad);
      } else {
	set_dL(abs(Fnum[i-1]), step*dprm[i]);
	set_dL(abs(Fnum1[i-1]), -step*dprm[i]);
	prm[i] = get_L(abs(Fnum[i-1]), 1);
      }
      dprm_max = max(fabs(dprm[i]), dprm_max);
      cout << fixed << setprecision(5) << setw(11) << prm[i];
      if (i % n_fold == 0) cout << endl;
    }
    if (n_prm % n_fold != 0) cout << endl;

    get_twiss(1, n_elem);

    while (!stable) {
      // roll back
      for (i = 1; i <= n_prm; i++)
	if (Fnum[i-1] > 0)
	  set_dbn(Fnum[i-1], Quad, -step*dprm[i]);
	else {
	  set_dL(abs(Fnum[i-1]), -step*dprm[i]);
	  set_dL(abs(Fnum1[i-1]), step*dprm[i]);
	}

      step /= 2.0;
      cout << endl;
      cout << scientific << setprecision(3) << "step = " << step << endl;

      for (i = 1; i <= n_prm; i++)
	if (Fnum[i-1] > 0) {
	  set_dbn(Fnum[i-1], Quad, step*dprm[i]);
	  prm_lim[i] = b2_max; prm[i] = get_bn(Fnum[i-1], 1, Quad);
	} else {
	  set_dL(abs(Fnum[i-1]), step*dprm[i]);
	  set_dL(abs(Fnum1[i-1]), -step*dprm[i]);
	  prm[i] = get_L(abs(Fnum[i-1]), 1);
	}
	
      get_twiss(1, n_elem);
    }
      
    for (j = 0; j < n_ind; j++)
      get_ab(elem_tps[ind[j]].A1, ab[j]);

    eps_[X_] = -D_[X_]/(4.0*alpha_rad_x);

    cout << endl;
    cout << scientific << setprecision(3)
	 << "eps_x       = " << setw(9) << eps_[X_].cst()
	 << fixed << setprecision(5)
	 << " nu = (" << setw(7) << h_ijklm(nus_[3], 0, 0, 0, 0, 0)
	 << ", " << setw(7) << h_ijklm(nus_[4], 0, 0, 0, 0, 0) << ")"
    	 << " ksi = ("  << setw(6) << h_ijklm(nus_[3], 0, 0, 0, 0, 1)
	 << ", " << setw(6) << h_ijklm(nus_[4], 0, 0, 0, 0, 1) << ")" << endl;
    prt_ab("_mp", ab[0], false);
    cout << fixed << setprecision(3)
	 << " eta_x = "
	 << setw(5) << h_ijklm(elem_tps[ind[0]].A1[x_],  0, 0, 0, 0, 1)
	 << endl;
    cout << scientific << setprecision(1)
	 << "eta_x       =  "
	 << setw(8) << h_ijklm(elem_tps[ind[1]].A1[x_],  0, 0, 0, 0, 1)
	 << "            eta\'_x     =  "
	 << setw(8) << h_ijklm(elem_tps[ind[1]].A1[px_],  0, 0, 0, 0, 1)
	 << endl;
    cout << scientific << setprecision(1)
	 << "eta_x       =  "
	 << setw(8) << h_ijklm(elem_tps[ind[2]].A1[x_],  0, 0, 0, 0, 1)
	 << "            eta\'_x     =  "
	 << setw(8) << h_ijklm(elem_tps[ind[2]].A1[px_],  0, 0, 0, 0, 1)
	 << endl;
    prt_ab("_qm1", ab[3], true);
    prt_ab("_ql2", ab[5], true);
    prt_ab("_ql3", ab[6], true);
    prt_ab("_qh2", ab[8], true);
    prt_ab("_qh3", ab[9], true);
    prt_ab("_ss", ab[4], true);
    prt_ab("_ls", ab[7], true);

    calc_twiss(); prt_lat("linlat.out");
  } while (dprm_max > prm_tol);

  free_dvector(b, 1, m); free_dvector(prm_lim, 1, n_prm);
  free_dvector(prm, 1, n_prm); free_dvector(dprm, 1, n_prm);
  free_dmatrix(A, 1, m, 1, n_prm);

  prt_lat("linlat.out");
  prt_mfile("flat_file.dat");
}


void prt_lin_map(const int n_DOF, const ss_vect<tps> &map)
{
  int  i, j;

  cout << endl; 
  for (i = 0; i < 2*n_DOF; i++) {
    for (j = 0; j < 2*n_DOF; j++)
      if (true) 
        cout << scientific << setprecision(5) << setw(13) << map[i][j];
      else
        cout << scientific << setprecision(16) << setw(24) << map[i][j];
    cout << endl;
  }
}


ss_vect<tps> get_A(const double alpha[], const double beta[],
                   const double eta[], const double etap[])
{
  ss_vect<tps>  A, Id;
    
  Id.identity();
  
  A.identity();
  A[x_]  = sqrt(beta[X_])*Id[x_];
  A[px_] = -alpha[X_]/sqrt(beta[X_])*Id[x_] + 1.0/sqrt(beta[X_])*Id[px_];
  A[y_]  = sqrt(beta[Y_])*Id[y_];
  A[py_] = -alpha[Y_]/sqrt(beta[Y_])*Id[y_] + 1.0/sqrt(beta[Y_])*Id[py_];

  A[x_] += eta[X_]*Id[delta_]; A[px_] += etap[X_]*Id[delta_];

  return A;
} 


tps get_alpha_x(const tps ab, ss_vect<tps> Id)
{
  tps  alpha, ps;
 
  ps = tps(0.0, px_+1); Id[px_] = ps;
  alpha = -ab*Id; alpha -= h_ijklm(alpha, 0, 1, 0, 0, 0)*ps;
  Id[px_] = 1.0; alpha = alpha*Id;
  Id[px_] = 0.0;

  return alpha;
}


tps get_beta_x(const tps ab, ss_vect<tps> Id)
{
  tps  beta, ps;
 
  ps = tps(0.0, x_+1); Id[x_] = ps;
  beta = ab*Id; beta -= h_ijklm(beta, 1, 0, 0, 0, 0)*ps;
  Id[x_] = 1.0; beta = beta*Id;
  Id[x_] = 0.0;

  return beta;
}


void fit(const double nu_x, const double nu_y, const bool fit_tune)
{
  const int  m = 21;

  int           i, j, jj[ss_dim];
  double        nu_fract[2];
  tps           ab[max_ind][2], S[m], S_inv[m], b[m];
  ss_vect<tps>  Id;

  const double  scl_eps_x  = 1e12;
  const double  scl_nu     = (fit_tune)? 1e6 : 0.0;
  const double  scl_ksi_x  = 1e3,   scl_ksi_y  = 1e3;
  const double  scl_beta   = 5e2,   scl_eta    = 1e6;
  const double  scl_beta_x = 0*1e0,   scl_beta_y = 0*1e0, scl_alpha = 1e5;

  const double  eps_x       = 2.0e-9;
  const double  ksi_max[]   = { -6.5, -2.4 };
  const double  beta_mp[]   = { 20.0, 10.0 };

  const double  beta_qm1[]  = {  5.0, 18.0 }; // beta_y is tunable

  const double  beta_ss[]   = {  1.0,  1.0 };
  const double  beta_ql2[]  = { 20.0, 10.0 };
  const double  beta_ql3[]  = { 10.0, 25.0 };

  const double  beta_ls[]   = { 15.0,  3.0 };
  const double  beta_qh2[]  = { 25.0, 10.0 };
  const double  beta_qh3[]  = { 10.0, 25.0 };

  nu_fract[X_] = fract(nu_x); nu_fract[Y_] = fract(nu_y);

  cout << endl;
  cout << fixed << setprecision(5)
       << "fit: nu_x = " << nu_fract[X_] << ", nu_y = " << nu_fract[Y_]
       << endl;

  fit_prms(true);

  n_ind = 0;
  ind[n_ind++] = get_loc(get_Fnum("sm2"), 1) - 1;
  ind[n_ind++] = get_loc(get_Fnum("b1"), 1) - 2; // upstream
  ind[n_ind++] = get_loc(get_Fnum("b1"), 2) - 1;
  ind[n_ind++] = get_loc(get_Fnum("qm1"), 1) - 1;

  ind[n_ind++] = get_loc(get_Fnum("ss"), 1) - 1;
  ind[n_ind++] = get_loc(get_Fnum("ql2"), 1) - 2; // upstream
  ind[n_ind++] = get_loc(get_Fnum("ql3"), 1) - 2; // upstream

  ind[n_ind++] = get_loc(get_Fnum("ls"), 1) - 1;
  ind[n_ind++] = get_loc(get_Fnum("qh2"), 1) - 1;
  ind[n_ind++] = get_loc(get_Fnum("qh3"), 1) - 1;

//  for (i = 1; i <= n_prm; i++) {
  for (i = 1; i <= 4; i++) {
    if (Fnum[i-1] > 0)
      set_bn_par(Fnum[i-1], Quad, ps_dim+i);
    else
      set_s_pair(abs(Fnum[i-1]), abs(Fnum1[i-1]), ps_dim+i);
  }

  get_twiss(1, n_elem);

  for (j = 0; j < n_ind; j++)
    get_ab(elem_tps[ind[j]].A1, ab[j]);

  eps_[X_] = -D_[X_]/(4.0*alpha_rad_x);

  danot_(3);

  for (i = 0; i < ss_dim; i++)
    Id[i] = (i < ps_dim)? 0.0 : tps(0.0, i+1);

  for (i = 0; i < ps_dim; i++)
    S[i] = tps(0.0, i+1);

  j = 6;
  S[j++] = get_alpha_x(ab[0][X_], Id);
  S[j++] = get_beta_x(ab[0][X_], Id);

  cout << S[6] << S[7];

  for (i = j; i < ss_dim; i++)
    S[i] = tps(0.0, i+1);

  for (i = 0; i < ss_dim; i++)
    jj[i] = (i < ps_dim)? 0 : 1;

  danot_(2);

  PInv(ss_dim, S, ss_dim, S_inv, jj);

  cout << S_inv[6] << S_inv[7];

  exit(0);

  j = 0;
  b[++j] = -scl_eps_x*(eps_[X_].cst()-eps_x);
  b[++j] = -scl_nu*(nus_[3].cst()-nu_fract[X_]);
  b[++j] = -scl_nu*(nus_[4].cst()-nu_fract[Y_]);
  b[++j] = -scl_ksi_x*(h_ijklm(nus_[3], 0, 0, 0, 0, 1)-ksi_max[X_]);
  b[++j] = -scl_ksi_y*(h_ijklm(nus_[4], 0, 0, 0, 0, 1)-ksi_max[Y_]);

  b[++j] = scl_alpha*h_ijklm(ab[0][X_], 0, 1, 0, 0, 0);
  b[++j] = scl_alpha*h_ijklm(ab[0][Y_], 0, 0, 0, 1, 0);
  b[++j] = -scl_beta_x*(h_ijklm(ab[0][X_], 1, 0, 0, 0, 0)-beta_mp[X_]);

  b[++j] = -scl_eta*h_ijklm(elem_tps[ind[1]].A1[x_],  0, 0, 0, 0, 1);
  b[++j] = -scl_eta*h_ijklm(elem_tps[ind[1]].A1[px_], 0, 0, 0, 0, 1);
  b[++j] = -scl_eta*h_ijklm(elem_tps[ind[2]].A1[x_],  0, 0, 0, 0, 1);
  b[++j] = -scl_eta*h_ijklm(elem_tps[ind[2]].A1[px_], 0, 0, 0, 0, 1);

  b[++j] = -scl_beta*(h_ijklm(ab[3][Y_], 0, 0, 1, 0, 0)-beta_qm1[Y_]);

  b[++j] = -scl_beta*(h_ijklm(ab[4][X_], 1, 0, 0, 0, 0)-beta_ss[X_]);
  b[++j] = -scl_beta*(h_ijklm(ab[4][Y_], 0, 0, 1, 0, 0)-beta_ss[Y_]);
 
  b[++j] = -scl_beta_x*(h_ijklm(ab[5][X_], 1, 0, 0, 0, 0)-beta_ql2[X_]);
  b[++j] = -scl_beta_y*(h_ijklm(ab[6][Y_], 0, 0, 1, 0, 0)-beta_ql3[Y_]);
 
  b[++j] = -scl_beta*(h_ijklm(ab[7][X_], 1, 0, 0, 0, 0)-beta_ls[X_]);
  b[++j] = -scl_beta*(h_ijklm(ab[7][Y_], 0, 0, 1, 0, 0)-beta_ls[Y_]);
 
  b[++j] = -scl_beta_x*(h_ijklm(ab[8][X_], 1, 0, 0, 0, 0)-beta_qh2[X_]);
  b[++j] = -scl_beta_y*(h_ijklm(ab[9][Y_], 0, 0, 1, 0, 0)-beta_qh3[Y_]);
 
  cout << endl;
  cout << scientific << setprecision(3)
       << "eps_x       = " << setw(9) << eps_[X_].cst()
       << fixed << setprecision(5)
       << " nu = (" << setw(7) << h_ijklm(nus_[3], 0, 0, 0, 0, 0)
       << ", " << setw(7) << h_ijklm(nus_[4], 0, 0, 0, 0, 0) << ")"
       << " ksi = ("  << setw(6) << h_ijklm(nus_[3], 0, 0, 0, 0, 1)
       << ", " << setw(6) << h_ijklm(nus_[4], 0, 0, 0, 0, 1) << ")" << endl;
  prt_ab("_mp", ab[0], false);
  cout << fixed << setprecision(3)
       << " eta_x = "
       << setw(5) << h_ijklm(elem_tps[ind[0]].A1[x_],  0, 0, 0, 0, 1)
       << endl;
  cout << scientific << setprecision(1)
       << "eta_x       =  "
       << setw(8) << h_ijklm(elem_tps[ind[1]].A1[x_],  0, 0, 0, 0, 1)
       << "            eta\'_x     =  "
       << setw(8) << h_ijklm(elem_tps[ind[1]].A1[px_],  0, 0, 0, 0, 1)
       << endl;
  cout << scientific << setprecision(1)
       << "eta_x       =  "
       << setw(8) << h_ijklm(elem_tps[ind[2]].A1[x_],  0, 0, 0, 0, 1)
       << "            eta\'_x     =  "
       << setw(8) << h_ijklm(elem_tps[ind[2]].A1[px_],  0, 0, 0, 0, 1)
       << endl;
  prt_ab("_qm1", ab[3], true);
  prt_ab("_ql2", ab[5], true);
  prt_ab("_ql3", ab[6], true);
  prt_ab("_qh2", ab[8], true);
  prt_ab("_qh3", ab[9], true);
  prt_ab("_ss", ab[4], true);
  prt_ab("_ls", ab[7], true);
}


int main(int argc, const char *argv[])
{
  long int      i1, i2;
  int           k, n_step;
  double        alpha[2], beta[2], dL;
  ss_vect<tps>  dm, Id;

  rad_on    = false; H_exact        = false; totpath_on   = false;
  cavity_on = false; quad_fringe_on = false; emittance_on = false;
  IBS_on    = false;

  // disable from TPSALib- and LieLib log messages
  idprset(-1);

  pi = 4.0*atan(1.0);

  rd_mfile("flat_file.ref", elem); rd_mfile("flat_file.ref", elem_tps);

  // initialize the symplectic integrator after energy has been defined
  ini_si();

  fit(0.256, 0.143, false); exit(0);

  if (false) {
    i1 = get_loc(get_Fnum("b1"), 1)-1; i2 = get_loc(get_Fnum("b1"), 2);

    get_ab_ms(i1, i2, 0.804, 0.110, alpha, beta);

    danot_(no_tps-1); Map.identity(); Map.propagate(i1, i2);

    Id.identity(); dm.identity();
    dm[px_] += alpha[X_]/beta[X_]*Id[x_]; dm[py_] += alpha[Y_]/beta[Y_]*Id[y_];
  
    prt_lin_map(2, Map);

    exit(0);
  }

  // radiation damping is constant
  get_emit(1, n_elem); alpha_rad_x = rad_[X_];

  no_mpoles();

  get_twiss(1, n_elem); calc_twiss(); i1 = get_loc(get_Fnum("b1"), 2) - 1;

  cout << endl;
  cout << fixed << setprecision(5)
       << "alpha = (" << elem[i1].Alpha[X_]
       << ", " << elem[i1].Alpha[Y_] << ")"
       << " beta = (" << elem[i1].Beta[X_] << ", " << elem[i1].Beta[Y_] << ")"
       << endl;

//  fit_cell(0.256, 0.143, false, 1e-3, false); exit(0);

  do {
    n_iter = 0;

    fit_cell(0.256, 0.143, false, 1e-3, false);
    fit_cell(0.256, 0.143, false, 1e-3, true);
  } while (n_iter > 2);

  if (false) {
    n_step = 10;
    for (k = 0; k <= n_step; k++) {
      dL = k*0.5/n_step;
      set_L(get_Fnum("du_qm2"), -dL); set_L(get_Fnum("dl_qm2"), dL);
      fit_cell(2.2, 1.1, 1e-2, false, true);
    }
  }

  get_twiss(1, n_elem); calc_twiss(); prt_lat("linlat.out");
}
