#define NO 3

#include "thor_lib.h"

int no_tps   = NO,
    ndpt_tps = 5;


extern tps          K, g;
extern ss_vect<tps> Map, A0, A1, Map_res;


const int n_prm_max = 20;


void get_twiss(const double alpha[], const double beta[],
	       const double eta[], const double etap[])
{
  int    j;
  double dnu[2];

  // Include parameter dependence.
  danot_(2);

  A1 = get_A(alpha, beta, eta, etap);
  for (j = 1; j <= n_elem; j++) {
    A1.propagate(j, j);
    elem_tps[j-1].A1 = get_A_CS(2, A1, dnu);
  }
}


tps get_eps_x(void)
{
  // Compute hor. emittance with parameter dependence.
  int Fnum, j;
  tps eps_x, J_x, J_z;

  cavity_on = false; emittance_on = false;

  Fnum = get_Fnum("bh");
  for (j = 1; j <= get_n_Kids(Fnum); j++)
    set_bn_par(Fnum, j, Quad, 7);

  danot_(2);

  get_Map(); K = MapNorm(Map, g, A1, A0, Map_res, 1);

  // Include parameter dependence.
  A1 = LieExp(g, A1);

  // Include dispersion terms.
  A1[x_]  += A0[x_][delta_]*tps(0e0, delta_+1);
  A1[px_] += A0[px_][delta_]*tps(0e0, delta_+1);
  A1[ct_] += A0[ct_][x_]*tps(0e0, x_+1);
  A1[ct_] += A0[ct_][px_]*tps(0e0, px_+1);

  emittance_on = true;

  danot_(2);

  A1.propagate(1, n_elem);

  eps_x = 1470e0*sqr(E0)*I5/(I2-I4); J_x = 1e0 - I4/I2; J_z = 3e0 - J_x;

  printf("\neps_x = %5.3f pm.rad, J_x = %5.3f, J_z = %5.3f\n",
	 1e3*eps_x.cst(), J_x.cst(), J_z.cst());

  emittance_on = false;

  return eps_x;
}


void chk_lat(double nu[], double ksi[])
{
  double        alpha[2], beta[2];
  ss_vect<tps>  nus;

//  get_Map();
  get_COD(10, 1e-10, 0e0, true);
  cout << Map;
  K = MapNorm(Map, g, A1, A0, Map_res, 1);
  nus = dHdJ(K); get_nu_ksi(nus, nu, ksi); get_ab(alpha, beta, 0);
  cout << endl;
  printf("alpha_x = %5.3f, alpha_y = %5.3f, beta_x = %5.3f, beta_y  = %5.3f\n",
	 alpha[X_], alpha[Y_], beta[X_], beta[Y_]);
  prt_nu(nus);
}


void prt_system(const int m, const int n_b2, double **A, double *b)
{
  int i, j;

  printf("\n Ax = b:\n");
  for (j = 1; j <= n_b2; j++)
    printf("%8d", j);
  printf("\n");
  for (i = 1; i <= m; i++) {
    printf("%2d", i);
    for (j = 1; j <= n_b2; j++)
      printf("%11.3e", A[i][j]);
    printf("%11.3e\n", b[i]);
  }
}

	
void fit_emit(int &n_b2, int b2s[], const double nu_x, const double nu_y,
	      const double b2_max, const double db2_tol, const double s_cut,
	      const double step)
{
  // Optimize hor. emittance.
  const double scl_nu = 1e1;

  const int m_max = 3;

  int          i, j, m, loc;
  double       **A, *b, *b2_lim, *b2, *db2;
  double       db2_max;
  tps          eps_x;
  ss_vect<tps> nus;

  b = dvector(1, m_max); b2_lim = dvector(1, n_b2);
  b2 = dvector(1, n_b2); db2 = dvector(1, n_b2);
  A = dmatrix(1, m_max, 1, n_b2);

  for (i = 1; i <= n_b2; i++) {
    if (b2s[i-1] > 0) {
      b2_lim[i] = b2_max; b2[i] = get_bn(b2s[i-1], 1, Quad);
    } else {
      b2_lim[i] = ds_max;
      loc = get_loc(abs(b2s[i-1]), 1) - 1; b2[i] = get_L(elem[loc-1].Fnum, 1);
    }
  }

  danot_(2);

  do {
    for (i = 1; i <= n_b2; i++) {
      for (j = 1; j <= get_n_Kids(abs(b2s[i-1])); j++)
	if (b2s[i-1] > 0)
	  set_bn_par(b2s[i-1], j, Quad, 7);
	else
	  set_s_par(abs(b2s[i-1]), j, 7);

      get_Map(); K = MapNorm(Map, g, A1, A0, Map_res, 1); nus = dHdJ(K);

      eps_x = get_eps_x();
 
      m = 0;
      A[++m][i] = h_ijklm_p(eps_x, 0, 0, 0, 0, 0, 7);
      A[++m][i] = scl_nu*h_ijklm_p(nus[3], 0, 0, 0, 0, 0, 7);
      A[++m][i] = scl_nu*h_ijklm_p(nus[4], 0, 0, 0, 0, 0, 7);

      if (b2s[i-1] < 0)
	for (j = 1; j <= m; j++)
	  A[j][i] *= scl_ds;

      for (j = 1; j <= get_n_Kids(abs(b2s[i-1])); j++)
	if (b2s[i-1] > 0)
	  clr_bn_par(b2s[i-1], j, Quad);
	else
	  clr_s_par(abs(b2s[i-1]), j);
    }

    m = 0;
    b[++m] = -eps_x.cst();
    b[++m] = -scl_nu*(nus[3].cst()-nu_x);
    b[++m] = -scl_nu*(nus[4].cst()-nu_y);

    prt_system(m, n_b2, A, b);

    SVD_lim(m, n_b2, A, b, b2_lim, s_cut, b2, db2);

    printf("\n");
    db2_max = 0e0;
    for (i = 1; i <= n_b2; i++) {
      set_dbn_s(b2s[i-1], Quad, step*db2[i]);
      b2[i] = get_bn_s(b2s[i-1], 1, Quad);
      db2_max = max(fabs(step*db2[i]), db2_max);

      printf("%14.5e", b2[i]);
    }
    printf("\n");
  } while (db2_max > db2_tol);

  free_dvector(b, 1, m_max);
  free_dvector(b2_lim, 1, n_b2); free_dvector(b2, 1, n_b2);
  free_dvector(db2, 1, n_b2);
  free_dmatrix(A, 1, m_max, 1, n_b2);
}


void fit_match(int &n_b2, int b2s[], const double beta_x, const double beta_y,
	       const double eta_x, const double b2_max, const double db2_tol,
	       const double s_cut, const double step)
{
  // Match linear optics.

  const int m_max = 6;

  long int     loc;
  int          i, j, m;
  double       **A, *b, *b2_lim, *b2, *db2;
  double       db2_max;
  ss_vect<tps> AA_tp, A_disp;

  const double scl_eta = 1e1,
               alpha[] = {0e0,    0e0},
               beta[]  = {beta_x, beta_y},
               eta[]   = {eta_x,  0e0},
	       etap[]  = {0e0,    0e0};

  b = dvector(1, m_max); b2_lim = dvector(1, n_b2);
  b2 = dvector(1, n_b2); db2 = dvector(1, n_b2);
  A = dmatrix(1, m_max, 1, n_b2);

  for (i = 1; i <= n_b2; i++) {
    if (b2s[i-1] > 0) {
      b2_lim[i] = b2_max; b2[i] = get_bn(b2s[i-1], 1, Quad);
    } else {
      b2_lim[i] = ds_max;
      loc = get_loc(abs(b2s[i-1]), 1) - 1; b2[i] = get_L(elem[loc-1].Fnum, 1);
    }
  }

  danot_(3);

  loc = get_loc(get_Fnum("bm"), 1);

  do {
    for (i = 1; i <= n_b2; i++) {
      for (j = 1; j <= get_n_Kids(abs(b2s[i-1])); j++)
	if (b2s[i-1] > 0)
	  set_bn_par(b2s[i-1], j, Quad, 7);
	else
	  set_s_par(abs(b2s[i-1]), j, 7);

      get_twiss(alpha, beta, eta, etap);
      AA_tp = elem_tps[n_elem-1].A1*tp_S(2, elem_tps[n_elem-1].A1);
      A_disp = elem_tps[loc-1].A1;

      m = 0;
      A[++m][i] = h_ijklm_p(-AA_tp[x_], 0, 1, 0, 0, 0, 7);
      A[++m][i] = h_ijklm_p(-AA_tp[y_], 0, 0, 0, 1, 0, 7);
      A[++m][i] = h_ijklm_p(AA_tp[x_], 1, 0, 0, 0, 0, 7);
      A[++m][i] = h_ijklm_p(AA_tp[y_], 0, 0, 1, 0, 0, 7);
      A[++m][i] = scl_eta*h_ijklm_p(A_disp[x_],  0, 0, 0, 0, 1, 7);
      A[++m][i] = scl_eta*h_ijklm_p(A_disp[px_], 0, 0, 0, 0, 1, 7);

      if (b2s[i-1] < 0)
	for (j = 1; j <= m; j++)
	  A[j][i] *= scl_ds;

      for (j = 1; j <= get_n_Kids(abs(b2s[i-1])); j++)
	if (b2s[i-1] > 0)
	  clr_bn_par(b2s[i-1], j, Quad);
	else
	  clr_s_par(abs(b2s[i-1]), j);
    }

    m = 0;
    b[++m] = -h_ijklm(-AA_tp[x_], 0, 1, 0, 0, 0);
    b[++m] = -h_ijklm(-AA_tp[y_], 0, 0, 0, 1, 0);
    b[++m] = -(h_ijklm(AA_tp[x_], 1, 0, 0, 0, 0)-beta_x);
    b[++m] = -(h_ijklm(AA_tp[y_], 0, 0, 1, 0, 0)-beta_y);
    b[++m] = -scl_eta*h_ijklm(A_disp[x_],  0, 0, 0, 0, 1);
    b[++m] = -scl_eta*h_ijklm(A_disp[px_], 0, 0, 0, 0, 1);

    prt_system(m, n_b2, A, b);

    SVD_lim(m, n_b2, A, b, b2_lim, s_cut, b2, db2);

    printf("\n");
    db2_max = 0e0;
    for (i = 1; i <= n_b2; i++) {
      set_dbn_s(b2s[i-1], Quad, step*db2[i]);
      b2[i] = get_bn_s(b2s[i-1], 1, Quad);
      db2_max = max(fabs(step*db2[i]), db2_max);

      printf("%14.5e", b2[i]);
    }
    printf("\n");
  } while (db2_max > db2_tol);

  free_dvector(b, 1, m_max);
  free_dvector(b2_lim, 1, n_b2); free_dvector(b2, 1, n_b2);
  free_dvector(db2, 1, n_b2);
  free_dmatrix(A, 1, m_max, 1, n_b2);
}


void get_b2s(int &n_b2, int b2_Fams[])
{

  n_b2 = 0;
  if (false) {
    b2_Fams[n_b2++] = get_Fnum("bh");
    b2_Fams[n_b2++] = get_Fnum("qf");
  } else {
    b2_Fams[n_b2++] = get_Fnum("bm");
    b2_Fams[n_b2++] = get_Fnum("qfe");
    b2_Fams[n_b2++] = get_Fnum("qde");
    b2_Fams[n_b2++] = get_Fnum("qm");
    b2_Fams[n_b2++] = -get_Fnum("bm");
    b2_Fams[n_b2++] = -get_Fnum("qfe");
    b2_Fams[n_b2++] = -get_Fnum("qde");
    b2_Fams[n_b2++] = -get_Fnum("qm");
  }
}


void tst_eps_x()
{
  int    Fnum;
  double eps_x[2], b2, der;

  const double db2 = 1e-3;
  
  Fnum = get_Fnum("bh");
  b2 = get_bn(Fnum, 1, Quad);
  set_bn(Fnum, Quad, b2+db2); eps_x[1] = (get_eps_x()).cst();
  set_bn(Fnum, Quad, b2-db2); eps_x[0] = (get_eps_x()).cst();
  set_bn(Fnum, Quad, b2);

  der = (eps_x[1]-eps_x[0])/(2e0*db2);

  printf("\n%14.5e\n", der);
}


int main(int argc, char *argv[])
{
  int b2_Fams[n_prm_max], n_b2;
  tps eps_x;

  const double nu[] = {5.0/16.0, 1.0/16.0},
        beta[]      = {1.26669,  1.81291},
	  eta_x     = 0.00618198;

  rad_on    = false; H_exact        = false; totpath_on   = false;
  cavity_on = false; quad_fringe_on = false; emittance_on = false;
  IBS_on    = false;

  rd_mfile(argv[1], elem); rd_mfile(argv[1], elem_tps);
  
  // Initialize the symplectic integrator after energy has been defined.
  ini_si();

  // Disable log messages from TPSALib and LieLib.
  idprset(-1);

  danot_(1);

  // tst_eps_x();
  // eps_x = get_eps_x();
  // cout << eps_x << "\n"; 

  if (false) {
    get_b2s(n_b2, b2_Fams);
    fit_emit(n_b2, b2_Fams, nu[X_], nu[Y_], 100.0, 1e-4, 1e-4, 0.5);
  }

  if (true) {
    get_b2s(n_b2, b2_Fams);
    ds_max = 0.2, scl_ds = 0.01;
    fit_match(n_b2, b2_Fams, beta[X_], beta[Y_], eta_x, 100.0, 1e-5, 1e-4, 0.5);
  }
}
