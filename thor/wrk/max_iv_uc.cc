#define NO 3

#include "thor_lib.h"

int no_tps   = NO,
    ndpt_tps = 5;


extern tps          K, g;
extern ss_vect<tps> Map, A0, A1, Map_res;


struct param_type {
private:

public:
  int                 n_prm;
  double              bn_tol, svd_cut, step;
  std::vector<double> bn_max, bn_scl;
  std::vector<int>    Fnum, n;

  void add_prm(const std::string Fname, const int n,
	       const double bn_max, const double bn_scl);
  void ini_prm(double *bn, double *bn_lim);
  void set_prm(const int k) const;
  void clr_prm(const int k) const;
  void set_bn(double *db2, double *b2, double &db2_max) const;
};


void param_type::add_prm(const std::string Fname, const int n,
			 const double bn_max, const double bn_scl)
{
  Fnum.push_back(get_Fnum(Fname.c_str()));
  this->n.push_back(n);
  this->bn_max.push_back(bn_max);
  this->bn_scl.push_back(bn_scl);
  n_prm = Fnum.size();
}


void param_type::ini_prm(double *bn, double *bn_lim)
{
  int i;

  n_prm = Fnum.size();
  for (i = 1; i <= n_prm; i++) {
    bn_lim[i] = bn_max[i-1];
    if (n[i-1] > 0)
      bn[i] = get_bn(Fnum[i-1], 1, n[i-1]);
    else if (n[i-1] == -1)
      bn[i] = get_L(Fnum[i-1], 1);
    else if (n[i-1] == -2)
      bn[i] = get_bn_s(-Fnum[i-1], 1, n[i-1]);
  }
}


void param_type::set_prm(const int k) const
{
  int j;

  for (j = 1; j <= get_n_Kids(Fnum[k]); j++)
    if (n[k] > 0)
      set_bn_par(Fnum[k], j, n[k], 7);
    else if (n[k] == -1)
      set_L_par(Fnum[k], j, 7);
    else if (n[k] == -2)
      set_s_par(Fnum[k], j, 7);
}


void param_type::clr_prm(const int k) const
{
  int j;

  for (j = 1; j <= get_n_Kids(Fnum[k]); j++)
    if (n[k] > 0)
      clr_bn_par(Fnum[k], j, n[k]);
    else if (n[k] == -1)
      clr_L_par(Fnum[k], j);
    else if (n[k] == -2)
      clr_s_par(Fnum[k], j);
}


void param_type::set_bn(double *db2, double *b2, double &db2_max) const
{
  int i;

  printf("\n");
  db2_max = 0e0;
  for (i = 1; i <= n_prm; i++) {
    db2[i] *= step*bn_scl[i-1];
    if (n[i-1] > 0) {
      set_dbn(Fnum[i-1], n[i-1], db2[i]);
      b2[i] = get_bn(Fnum[i-1], 1, n[i-1]);
    } else if (n[i-1] == -1) {
      set_dL(Fnum[i-1], db2[i]);
      b2[i] = get_L(Fnum[i-1], 1);
    } else if (n[i-1] == -2) {
      set_dbn_s(-Fnum[i-1], n[i-1], db2[i]);
      b2[i] = get_bn_s(-Fnum[i-1], 1, n[i-1]);
    }
    db2_max = max(fabs(db2[i]), db2_max);
    printf("%10.5f", b2[i]);
  }
  printf("\n");
}


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


void tst_dL(const double beta_x, const double beta_y, const double eta_x)
{
  int          Fnum;
  double       der;
  ss_vect<tps> AA_tp, AA_tp0, AA_tp1;

  const double dL      = 1e-3,
               alpha[] = {0e0,    0e0},
               beta[]  = {beta_x, beta_y},
               eta[]   = {eta_x,  0e0},
	       etap[]  = {0e0,    0e0};
  
  const char elem_name[] = "bm";

  Fnum = get_Fnum(elem_name);
  
  set_s_par(Fnum, 1, 7);
  get_twiss(alpha, beta, eta, etap);
  clr_s_par(Fnum, 1);
  AA_tp = elem_tps[n_elem-1].A1*tp_S(2, elem_tps[n_elem-1].A1);

  // set_bn_s(-Fnum, Quad, L+dL);
  set_dbn_s(-Fnum, Quad, dL);

  get_twiss(alpha, beta, eta, etap);
  AA_tp1 = elem_tps[n_elem-1].A1*tp_S(2, elem_tps[n_elem-1].A1);

  // set_bn_s(-Fnum, Quad, L-dL);
  set_dbn_s(-Fnum, Quad, -2e0*dL);

  get_twiss(alpha, beta, eta, etap);
  AA_tp0 = elem_tps[n_elem-1].A1*tp_S(2, elem_tps[n_elem-1].A1);

  set_dbn_s(-Fnum, Quad, dL);

  der = (AA_tp1[x_][x_]-AA_tp0[x_][x_])/(2e0*dL);

  printf("\n%14.5e %14.5e\n", der, h_ijklm_p(AA_tp[x_], 1, 0, 0, 0, 0, 7));
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

  printf("eps_x = %5.3f pm.rad, J_x = %5.3f, J_z = %5.3f\n",
	 1e3*eps_x.cst(), J_x.cst(), J_z.cst());

  emittance_on = false;

  return eps_x;
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


void prt_system(const int m, const int n_b2, double **A, double *b)
{
  int i, j;

  printf("\n Ax = b:\n");
  for (j = 1; j <= n_b2; j++)
    if (j == 1)
      printf("%9d", j);
    else
      printf("%11d", j);
  printf("\n");
  for (i = 1; i <= m; i++) {
    printf("%2d", i);
    for (j = 1; j <= n_b2; j++)
      printf("%11.3e", A[i][j]);
    printf("%11.3e\n", b[i]);
  }
}

	
void fit_emit(param_type &b2_prms, const double eps_x,
	      const double nu_x, const double nu_y)
{
  // Optimize unit cell for hor. emittance.
  // 
  const double scl_nu = 1e1, scl_ksi = 1e-2, scl_eps = 1e1;

  const int m_max = 8;

  int          n_b2, i, j, m;
  double       **A, *b, *b2_lim, *b2, *db2;
  double       db2_max;
  tps          eps1_x;
  ss_vect<tps> nus;

  n_b2 = b2_prms.n_prm;
  printf("n_prm = %d\n", b2_prms.n_prm);

  b = dvector(1, m_max); b2_lim = dvector(1, n_b2); b2 = dvector(1, n_b2);
  db2 = dvector(1, n_b2); A = dmatrix(1, m_max, 1, n_b2);

  b2_prms.ini_prm(b2, b2_lim);

  danot_(2);

  do {
    printf("\n");
    for (i = 1; i <= n_b2; i++) {
      b2_prms.set_prm(i-1);

      get_Map(); K = MapNorm(Map, g, A1, A0, Map_res, 1); nus = dHdJ(K);

      eps1_x = get_eps_x();
 
      m = 0;
      A[++m][i] = scl_eps*h_ijklm_p(eps1_x, 0, 0, 0, 0, 0, 7);
      A[++m][i] = scl_nu*h_ijklm_p(nus[3], 0, 0, 0, 0, 0, 7);
      A[++m][i] = scl_nu*h_ijklm_p(nus[4], 0, 0, 0, 0, 0, 7);
      A[++m][i] = scl_ksi*h_ijklm_p(nus[3], 0, 0, 0, 0, 1, 7);
      A[++m][i] = scl_ksi*h_ijklm_p(nus[4], 0, 0, 0, 0, 1, 7);

      for (j = 1; j <= m; j++)
	A[j][i] *= b2_prms.bn_scl[i-1];

      b2_prms.clr_prm(i-1);
    }

    m = 0;
    b[++m] = -scl_eps*(eps1_x.cst()-eps_x);
    b[++m] = -scl_nu*(h_ijklm(nus[3], 0, 0, 0, 0, 0)-nu_x);
    b[++m] = -scl_nu*(h_ijklm(nus[4], 0, 0, 0, 0, 0)-nu_y);
    b[++m] = -scl_ksi*h_ijklm(nus[3], 0, 0, 0, 0, 1);
    b[++m] = -scl_ksi*h_ijklm(nus[4], 0, 0, 0, 0, 1);

    prt_system(m, n_b2, A, b);

    SVD_lim(m, n_b2, A, b, b2_lim, b2_prms.svd_cut, b2, db2);

    b2_prms.set_bn(db2, b2, db2_max);
  } while (db2_max > b2_prms.bn_tol);

  free_dvector(b, 1, m_max); free_dvector(b2_lim, 1, n_b2);
  free_dvector(b2, 1, n_b2); free_dvector(db2, 1, n_b2);
  free_dmatrix(A, 1, m_max, 1, n_b2);
}


void prt_match(const param_type &b2_prms, const double *b2)
{
  double l4, l5h, l6, l7h, l8;
  FILE *outf;

  std::string file_name = "match.out";

  outf = file_write(file_name.c_str());

  fprintf(outf, "bm:  bending, l = 0.166667, t = 0.5, k = %8.5f, t1 = 0.0"
	  ", t2 = 0.0,\n     gap = 0.00, N = Nbend, Method = Meth;\n", b2[1]);
  fprintf(outf, "qfe: quadrupole, l = 0.15, k = %9.5f, N = Nquad"
	  ", Method = Meth;\n", b2[2]);
  fprintf(outf, "qde: quadrupole, l = 0.1,  k = %9.5f, N = Nquad"
	  ", Method = Meth;\n", b2[3]);
  fprintf(outf, "qm:  quadrupole, l = 0.15, k = %9.5f, N = Nquad"
	  ", Method = Meth;\n", b2[4]);

  l4  = get_L(get_Fnum("l4"), 1);  l5h = get_L(get_Fnum("l5h"), 1);
  l6  = get_L(get_Fnum("l6"), 1);  l7h = get_L(get_Fnum("l7h"), 1);
  l8  = get_L(get_Fnum("l8"), 1);
  
  fprintf(outf, "\nl4:  drift, l = %7.5f;\n", l4+b2[6]);
  fprintf(outf, "l5h: drift, l = %7.5f;\n", l5h-b2[6]/2e0+b2[7]/2e0);
  fprintf(outf, "l6:  drift, l = %7.5f;\n", l6-b2[7]+b2[5]);
  fprintf(outf, "l7h: drift, l = %7.5f;\n", l7h-b2[5]/2e0+b2[8]/2e0);
  fprintf(outf, "l8:  drift, l = %7.5f;\n", l8-b2[8]);

  fclose(outf);
}


void fit_match(param_type &b2_prms,
	       const double beta0_x, const double beta0_y, const double eta0_x,
	       const double beta1_x, const double beta1_y)
{
  // Match linear optics to straight section.

  const int m_max = 8;

  long int     loc;
  int          n_b2, i, j, m;
  double       **A, *b, *b2_lim, *b2, *db2;
  double       db2_max;
  ss_vect<tps> AA_tp, A_disp;

  const double scl_eta  = 1e1,
               alpha0[] = {0e0,     0e0},
               beta0[]  = {beta0_x, beta0_y},
               eta0[]   = {eta0_x,  0e0},
	       etap0[]  = {0e0,     0e0};

  n_b2 = b2_prms.n_prm;

  b = dvector(1, m_max); b2_lim = dvector(1, n_b2); b2 = dvector(1, n_b2);
  db2 = dvector(1, n_b2); A = dmatrix(1, m_max, 1, n_b2);

  b2_prms.ini_prm(b2, b2_lim);

  danot_(3);

  loc = get_loc(get_Fnum("bm"), 1);

  do {
    for (i = 1; i <= n_b2; i++) {
      b2_prms.set_prm(i-1);

      get_twiss(alpha0, beta0, eta0, etap0);
      AA_tp = elem_tps[n_elem-1].A1*tp_S(2, elem_tps[n_elem-1].A1);
      A_disp = elem_tps[loc-1].A1;

      m = 0;
      A[++m][i] = h_ijklm_p(-AA_tp[x_], 0, 1, 0, 0, 0, 7);
      A[++m][i] = h_ijklm_p(-AA_tp[y_], 0, 0, 0, 1, 0, 7);
      A[++m][i] = scl_eta*h_ijklm_p(A_disp[x_],  0, 0, 0, 0, 1, 7);
      A[++m][i] = scl_eta*h_ijklm_p(A_disp[px_], 0, 0, 0, 0, 1, 7);
      A[++m][i] = h_ijklm_p(AA_tp[x_], 1, 0, 0, 0, 0, 7);
      A[++m][i] = h_ijklm_p(AA_tp[y_], 0, 0, 1, 0, 0, 7);

      for (j = 1; j <= m; j++)
	A[j][i] *= b2_prms.bn_scl[i-1];

      b2_prms.clr_prm(i-1);
    }

    m = 0;
    b[++m] = -h_ijklm(-AA_tp[x_], 0, 1, 0, 0, 0);
    b[++m] = -h_ijklm(-AA_tp[y_], 0, 0, 0, 1, 0);
    b[++m] = -scl_eta*h_ijklm(A_disp[x_],  0, 0, 0, 0, 1);
    b[++m] = -scl_eta*h_ijklm(A_disp[px_], 0, 0, 0, 0, 1);
    b[++m] = -(h_ijklm(AA_tp[x_], 1, 0, 0, 0, 0)-beta1_x);
    b[++m] = -(h_ijklm(AA_tp[y_], 0, 0, 1, 0, 0)-beta1_y);

    prt_system(m, n_b2, A, b);

    SVD_lim(m, n_b2, A, b, b2_lim, b2_prms.svd_cut, b2, db2);

    b2_prms.set_bn(db2, b2, db2_max);
  } while (db2_max > b2_prms.bn_tol);

  prt_match(b2_prms, b2);

  free_dvector(b, 1, m_max); free_dvector(b2_lim, 1, n_b2);
  free_dvector(b2, 1, n_b2); free_dvector(db2, 1, n_b2);
  free_dmatrix(A, 1, m_max, 1, n_b2);
}


void fit_3rd_achrom(std::vector<int> &bns,
		    const double bn_max, const double dbn_tol,
		    const double s_cut, const double step)
{
  // Third order achromat.

  const int m_max = 6;

  int          n_bn, i, j, m, loc;
  double       **A, *b, *bn_lim, *bn, *dbn;
  double       dbn_max;
  tps          K_re, K_im;
  ss_vect<tps> nus;

  n_bn = bns.size();

  b = dvector(1, m_max); bn_lim = dvector(1, n_bn);
  bn = dvector(1, n_bn); dbn = dvector(1, n_bn);
  A = dmatrix(1, m_max, 1, n_bn);

  for (i = 1; i <= n_bn; i++) {
    if (bns[i-1] > 0) {
      bn_lim[i] = bn_max; bn[i] = get_bn(bns[i-1], 1, Oct);
    } else {
      bn_lim[i] = ds_max;
      loc = get_loc(abs(bns[i-1]), 1) - 1; bn[i] = get_L(elem[loc-1].Fnum, 1);
    }
  }

  danot_(5);

  do {
    for (i = 1; i <= n_bn; i++) {
      for (j = 1; j <= get_n_Kids(abs(bns[i-1])); j++)
	if (bns[i-1] > 0)
	  set_bn_par(bns[i-1], j, Oct, 7);
	else
	  set_s_par(abs(bns[i-1]), j, 7);

      get_Map(); K = MapNorm(Map, g, A1, A0, Map_res, 1); CtoR(K, K_re, K_im);

      m = 0;
      A[++m][i] = h_ijklm_p(K_re, 2, 2, 0, 0, 0, 7);
      A[++m][i] = h_ijklm_p(K_re, 0, 0, 2, 2, 0, 7);
      A[++m][i] = h_ijklm_p(K_re, 1, 1, 1, 1, 0, 7);
      // A[++m][i] = h_ijklm_p(K_re, 1, 1, 0, 0, 2, 7);
      // A[++m][i] = h_ijklm_p(K_re, 0, 0, 1, 1, 2, 7);

      if (bns[i-1] < 0)
	for (j = 1; j <= m; j++)
	  A[j][i] *= scl_ds;

      for (j = 1; j <= get_n_Kids(abs(bns[i-1])); j++)
	if (bns[i-1] > 0)
	  clr_bn_par(bns[i-1], j, Oct);
	else
	  clr_s_par(abs(bns[i-1]), j);
    }

    m = 0;
    b[++m] = -h_ijklm(K_re, 2, 2, 0, 0, 0);
    b[++m] = -h_ijklm(K_re, 0, 0, 2, 2, 0);
    b[++m] = -h_ijklm(K_re, 1, 1, 1, 1, 0);
    // b[++m] = -h_ijklm(K_re, 1, 1, 0, 0, 2);
    // b[++m] = -h_ijklm(K_re, 0, 0, 1, 1, 2);

    prt_system(m, n_bn, A, b);

    SVD_lim(m, n_bn, A, b, bn_lim, s_cut, bn, dbn);

    printf("\n");
    dbn_max = 0e0;
    for (i = 1; i <= n_bn; i++) {
      set_dbn_s(bns[i-1], Oct, step*dbn[i]);
      bn[i] = get_bn_s(bns[i-1], 1, Oct);
      dbn_max = max(fabs(step*dbn[i]), dbn_max);

      printf("%14.5e", bn[i]);
    }
    printf("\n");
  } while (dbn_max > dbn_tol);

  free_dvector(b, 1, m_max);
  free_dvector(bn_lim, 1, n_bn); free_dvector(bn, 1, n_bn);
  free_dvector(dbn, 1, n_bn);
  free_dmatrix(A, 1, m_max, 1, n_bn);
}


void chk_lat(double nu[], double ksi[])
{
  double       alpha[2], beta[2];
  ss_vect<tps> nus;

//  get_Map();
  get_COD(10, 1e-10, 0e0, true);
  std::cout << Map;
  K = MapNorm(Map, g, A1, A0, Map_res, 1);
  nus = dHdJ(K); get_nu_ksi(nus, nu, ksi); get_ab(alpha, beta, 0);
   std::cout <<  std::endl;
  printf("alpha_x = %5.3f, alpha_y = %5.3f, beta_x = %5.3f, beta_y  = %5.3f\n",
	 alpha[X_], alpha[Y_], beta[X_], beta[Y_]);
  prt_nu(nus);
}


int main(int argc, char *argv[])
{
  double           Jx, Jy;
  tps              eps1_x, K_re, K_im, g_re, g_im, h_re, h_im, H_re, H_im;
  ss_vect<tps>     Id_scl;
  param_type       b2_prms, b3_prms;
  std::vector<int> b2_Fams, bn_Fams;
  std::ofstream    outf;
  
  const double max_Ax = 5e-3, max_Ay = 5e-3, max_delta = 3e-2;
  
  const double eps_x   = 15e-3,
               nu[]    = {4.0/15.0, 1.0/15.0},
               beta0[] = {0.38133, 5.00190},
	       eta0_x  = 0.00382435,
               beta1[] = {3.0, 3.0};
	    

  rad_on    = false; H_exact        = false; totpath_on   = false;
  cavity_on = false; quad_fringe_on = false; emittance_on = false;
  IBS_on    = false;

  rd_mfile(argv[1], elem); rd_mfile(argv[1], elem_tps);
  
  // Initialize the symplectic integrator after energy has been defined.
  ini_si();

  // Disable log messages from TPSALib and LieLib.
  idprset(-1);

  danot_(1);

  Jx = sqr(max_Ax)/(2.0*beta1[X_]); Jy = sqr(max_Ay)/(2.0*beta1[Y_]);

  Id_scl.identity();
  Id_scl[x_] *= sqrt(2.0*Jx); Id_scl[px_] *= sqrt(2.0*Jx);
  Id_scl[y_] *= sqrt(2.0*Jy); Id_scl[py_] *= sqrt(2.0*Jy);
  Id_scl[delta_] *= max_delta;

  // tst_eps_x();
  // eps1_x = get_eps_x();
  // cout << eps_x << "\n"; 

  // tst_dL(beta[X_], beta[Y_], eta_x);

  if (false) {
    b2_prms.add_prm("bh",  2, 10.0, 1.0);
    b2_prms.add_prm("qf",  2, 25.0, 1.0);
    b2_prms.add_prm("l2", -1,  0.5, 0.01);

    b2_prms.bn_tol = 1e-4; b2_prms.svd_cut = 1e-6;  b2_prms.step = 0.8;

    fit_emit(b2_prms, eps_x, nu[X_], nu[Y_]);

    prt_mfile("flat_file.dat");
  }

  if (true) {
    b2_prms.add_prm("bm",   2, 25.0, 1.0);
    b2_prms.add_prm("qfe",  2, 25.0, 1.0);
    b2_prms.add_prm("qde",  2, 50.0, 1.0);
    b2_prms.add_prm("qm",   2, 25.0, 1.0);

    b2_prms.add_prm("bm",  -2,  0.2, 0.01);
    b2_prms.add_prm("qfe", -2,  0.2, 0.01);
    b2_prms.add_prm("qde", -2,  0.2, 0.01);
    b2_prms.add_prm("qm",  -2,  0.2, 0.01);

    b2_prms.bn_tol = 1e-4; b2_prms.svd_cut = 1e-6;  b2_prms.step = 0.3;

    fit_match(b2_prms, beta0[X_], beta0[Y_], eta0_x, beta1[X_], beta1[Y_]);
  }

  if (false) {
    // b3_prms.add_prm("o1", 3);
    // b3_prms.add_prm("o2", 3);
    // b3_prms.add_prm("o3", 3);

    fit_3rd_achrom(bn_Fams, 1e5, 1e-5, 1e-4, 1.0);
  }

  if (false) {
    danot_(NO);

    get_Map(); K = MapNorm(Map, g, A1, A0, Map_res, no_tps);
    CtoR(K*Id_scl, K_re, K_im); CtoR(g*Id_scl, g_re, g_im);
    CtoR(get_h()*Id_scl, h_re, h_im);
    CtoR(get_H()*Id_scl, H_re, H_im);

    file_wr(outf, "map.dat"); outf << Map; outf.close();
    file_wr(outf, "K.dat"); outf << K_re*Id_scl; outf.close();
    file_wr(outf, "g.dat"); outf << g_im*Id_scl; outf.close();
    file_wr(outf, "h.dat"); outf << h_re*Id_scl; outf.close();
    file_wr(outf, "H.dat"); outf << H_re*Id_scl; outf.close();
  }
}
