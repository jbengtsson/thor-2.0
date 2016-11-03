#define NO 3

#include "thor_lib.h"

int  no_tps   = NO,
     ndpt_tps = 5;


extern tps          K, g;
extern ss_vect<tps> Map, A0, A1, Map_res;


const int n_prm_max = 20;


double get_eps_x(void)
{
  // Compute hor. emittance with parameter dependence.
  bool         cav, emit;
  double       eps_x;

  cav = cavity_on; emit = emittance_on;

  cavity_on = false; emittance_on = false;

  danot_(2);

  get_Map(); K = MapNorm(Map, g, A1, A0, Map_res, 1);

  // Include dispersion terms.
  A1[x_]  += A0[x_][delta_]*tps(0e0, delta_+1);
  A1[px_] += A0[px_][delta_]*tps(0e0, delta_+1);
  A1[ct_] += A0[ct_][x_]*tps(0e0, x_+1);
  A1[ct_] += A0[ct_][px_]*tps(0e0, px_+1);

  emittance_on = true;

  danot_(1);

  A1.propagate(1, n_elem);

  eps_x = 1470e0*pow(E0, 2)*I5/(I2-I4);

  printf("\neps_x = %5.3f pm.rad, J_x = %5.3f, J_z = %5.3f \n",
	 1e3*eps_x, 1e0-I4/I2, 2e0+I4/I2);

  cavity_on = cav; emittance_on = emit;

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
  cout << fixed << setprecision(3)
       << "alpha_x  = " << alpha[X_] << ", alpha_y = " << alpha[Y_]
       << ", beta_x = " << beta[X_] << ", beta_y  = " << beta[Y_] << endl;
  prt_nu(nus);
}

void prt_system(const int m, const int n_b2, double **A, double *b)
{
  int i, j;

  printf("\n Ax = b:\n\n");
  for (j = 1; j <= n_b2; j++)
    printf("%8d", j);
  printf("\n");
  for (i = 1; i <= m; i++) {
    printf("%2d", i);
    for (j = 1; j <= n_b2; j++)
      printf("%8.3f", A[i][j]);
    printf("%8.3f\n", b[i]);
  }
}

	
void fit_uc(int &n_b2, int b2s[], const double nu_x, const double nu_y,
	    const double b2_max, const double db2_tol, const double s_cut,
	    const double step)
{
  // Periodic unit cell: [nu_x, nu_y].

  const int m_max = 2;

  int          i, j, m, loc;
  double       **A, *b, *b2_lim, *b2, *db2;
  double       db2_max;
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

  danot_(3);

  get_Map(); K = MapNorm(Map, g, A1, A0, Map_res, 1); nus = dHdJ(K);

  do {
    for (i = 1; i <= n_b2; i++) {
      for (j = 1; j <= get_n_Kids(abs(b2s[i-1])); j++)
	if (b2s[i-1] > 0)
	  set_bn_par(b2s[i-1], j, Quad, 7);
	else
	  set_s_par(abs(b2s[i-1]), j, 7);

      get_Map(); K = MapNorm(Map, g, A1, A0, Map_res, 1); nus = dHdJ(K);
 
      m = 0;
      A[++m][i] = -h_ijklm_p(nus[3], 0, 0, 0, 0, 0, 7);
      A[++m][i] = -h_ijklm_p(nus[4], 0, 0, 0, 0, 0, 7);

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
    b[++m] = nus[3].cst() - nu_x;
    b[++m] = nus[4].cst() - nu_y;
    printf("\nnu  = [%8.5f, %8.5f]\n", nus[3].cst(), nus[4].cst());
    printf("dnu = [%8.5f, %8.5f]\n", b[1], b[2]);

    prt_system(m, n_b2, A, b);

    SVD_lim(m, n_b2, A, b, b2_lim, s_cut, b2, db2);

    db2_max = 0e0;
    for (i = 1; i <= n_b2; i++) {
      set_dbn_s(b2s[i-1], Quad, step*db2[i]);
      b2[i] = get_bn_s(b2s[i-1], 1, Quad);
      db2_max = max(fabs(step*db2[i]), db2_max);
    }
  } while (db2_max > db2_tol);

  free_dvector(b, 1, m_max);
  free_dvector(b2_lim, 1, n_b2); free_dvector(b2, 1, n_b2);
  free_dvector(db2, 1, n_b2);
  free_dmatrix(A, 1, m_max, 1, n_b2);
}


void get_b2s(int &n_b2, int b2_Fams[])
{

  n_b2 = 0;
  b2_Fams[n_b2++] = get_Fnum("bh");
  b2_Fams[n_b2++] = get_Fnum("qf");
}


int main(int argc, char *argv[])
{
  int b2_Fams[n_prm_max], n_b2 = 0;

  const double nu[] = {5.0/16.0, 1.0/16.0};

  rad_on    = false; H_exact        = false; totpath_on   = false;
  cavity_on = false; quad_fringe_on = false; emittance_on = false;
  IBS_on    = false;

  rd_mfile(argv[1], elem); rd_mfile(argv[1], elem_tps);
  
  // Initialize the symplectic integrator after energy has been defined.
  ini_si();

  // Disable log messages from TPSALib and LieLib.
  idprset(-1);

  danot_(1);

  get_eps_x();

  if (true) {
    get_b2s(n_b2, b2_Fams);
    fit_uc(n_b2, b2_Fams, nu[X_], nu[Y_], 100.0, 1e-4, 1e-4, 1.0);
  }
}
