#define NO 7

#include "thor_lib.h"

int no_tps   = NO,
    ndpt_tps = 5;


extern tps          K, g;
extern ss_vect<tps> Map, A0, A1, Map_res;


const int n_max = NO/2 - 1;

typedef double dnu1_type[n_max+1][n_max+1];
typedef double dnu2_type[2*n_max+1][2*n_max+1];

dnu2_type dnu_xy;

struct param_type {
private:

public:
  int                 m_constr, n_prm;
  double              bn_tol, svd_cut, step;
  double              *bn_lim, *bn, *dbn;
  std::vector<double> bn_max, bn_scl;
  std::vector<int>    Fnum, n;

  void add_prm(const std::string Fname, const int n,
	       const double bn_max, const double bn_scl);
  void ini_prm(void);
  void set_prm_dep(const int k) const;
  void clr_prm_dep(const int k) const;
  double set_dprm(void) const;
  void set_prm(void) const;
};


param_type   bn_prms;
int          n_iter, n_powell;
ss_vect<tps> Id_scl;

const int n_prm_max = 8;

const double scl_dnu = 1e5, scl_2d = 5e0;


void param_type::add_prm(const std::string Fname, const int n,
			 const double bn_max, const double bn_scl)
{
  Fnum.push_back(get_Fnum(Fname.c_str()));
  this->n.push_back(n);
  this->bn_max.push_back(bn_max);
  this->bn_scl.push_back(bn_scl);
  n_prm = Fnum.size();
}


void param_type::ini_prm(void)
{
  int i;

  n_prm = Fnum.size();

  bn_prms.bn_lim = dvector(1, n_prm); bn_prms.bn = dvector(1, n_prm);
  bn_prms.dbn = dvector(1, n_prm);

  printf("\nInitial bn (%d):\n", n_prm);
  for (i = 1; i <= n_prm; i++) {
    bn_lim[i] = bn_max[i-1];
    if (n[i-1] > 0)
      // Multipole.
      bn[i] = get_bn(Fnum[i-1], 1, n[i-1]);
    else if (n[i-1] == -1)
      // Drift.
      bn[i] = get_L(Fnum[i-1], 1);
    else if (n[i-1] == -2)
      // Location.
      bn[i] = get_bn_s(-Fnum[i-1], 1, n[i-1]);
    printf("%13.5e", bn[i]);
  }
  printf("\n");
}


void param_type::set_prm_dep(const int k) const
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


void param_type::clr_prm_dep(const int k) const
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


double param_type::set_dprm(void) const
{
  int    i;
  double dbn_max;

  printf("set_dprm:\n");
  dbn_max = 0e0;
  for (i = 1; i <= n_prm; i++) {
    dbn[i] *= step*bn_scl[i-1];
    if (n[i-1] > 0) {
      set_dbn(Fnum[i-1], n[i-1], dbn[i]);
      bn[i] = get_bn(Fnum[i-1], 1, n[i-1]);
    } else if (n[i-1] == -1) {
      set_dL(Fnum[i-1], dbn[i]);
      bn[i] = get_L(Fnum[i-1], 1);
    } else if (n[i-1] == -2) {
      set_dbn_s(-Fnum[i-1], n[i-1], dbn[i]);
      bn[i] = get_bn_s(-Fnum[i-1], 1, n[i-1]);
    }
    dbn_max = max(fabs(dbn[i]), dbn_max);
    printf(" %13.5e", bn[i]);
  }
  printf("\n");

  return dbn_max;
}


void param_type::set_prm(void) const
{
  int i;

  printf("set_prm:\n");
  for (i = 1; i <= n_prm; i++) {
    if (n[i-1] > 0)
      set_bn(Fnum[i-1], n[i-1], bn[i]);
    else if (n[i-1] == -1)
      set_L(Fnum[i-1], bn[i]);
    else if (n[i-1] == -2)
      set_bn_s(-Fnum[i-1], n[i-1], bn[i]);
    printf(" %13.5e", bn[i]);
  }
  printf("\n");
}


void no_mpoles(const int n)
{
  int j;

  printf("\nzeroing multipoles: %d\n", n);
  for (j = 0; j < n_elem; j++)
    if (elem[j].kind == Mpole)
	set_bn(elem[j].Fnum, elem[j].Knum, n, 0e0);
}


void get_S(void)
{
  int    j;
  double S;

  S = 0e0;
  for (j = 0; j < n_elem; j++) {
    S += elem[j].L;
    elem[j].S = S; elem_tps[j].S = S;
  }
}


void get_nu_ksi(void)
{
  ss_vect<tps> nus;

  danot_(2);
  get_Map();
  danot_(3);
  K = MapNorm(Map, g, A1, A0, Map_res, 1); nus = dHdJ(K);

  printf("\nnu  = [%8.5f, %8.5f]\n",
	 nus[0].cst(), nus[1].cst());
  printf("ksi = [%8.5f, %8.5f]\n",
	 h_ijklm(nus[0], 0, 0, 0, 0, 1),
	 h_ijklm(nus[1], 0, 0, 0, 0, 1));
}


void get_twiss(const ss_vect<tps> &A)
{
  int          j, k;
  double       alpha1[2], beta1[2], eta1[2], etap1[2], dnu1[2], dnu2[2];
  ss_vect<tps> A1;

  // Include parameter dependence.
  danot_(2);

  for (k = 0; k < 2; k++)
    dnu1[k] = 0e0;
  A1 = A;
  for (j = 1; j <= n_elem; j++) {
    A1.propagate(j, j);
    elem_tps[j-1].A1 = get_A_CS(2, A1, dnu2);

    // Store linear optics for convenience.
    get_ab(A1, alpha1, beta1, dnu2, eta1, etap1);
    for (k = 0; k < 2; k++) {
      elem[j-1].Alpha[k] = alpha1[k]; elem[j-1].Beta[k] = beta1[k];
      elem[j-1].Eta[k] = eta1[k]; elem[j-1].Etap[k] = etap1[k];
    }
    // Assumes dnu < 360 degrees.
    for (k = 0; k < 2; k++) {
      elem[j-1].Nu[k] = floor(elem[j-2].Nu[k]) + dnu2[k];
      if ((dnu2[k] < dnu1[k]) && (elem[j-1].L >= 0e0)) elem[j-1].Nu[k] += 1e0;
    }
    for (k = 0; k < 2; k++)
      dnu1[k] = dnu2[k];
  }
}


void get_twiss(const double alpha[], const double beta[],
	       const double eta[], const double etap[])
{
  ss_vect<tps> A;

  // Include parameter dependence.
  danot_(2);

  A = get_A(alpha, beta, eta, etap); get_twiss(A);
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

	
void prt_emit(const param_type &bn_prms)
{
  FILE *outf;

  std::string file_name = "emit.out";

  outf = file_write(file_name.c_str());
  fprintf(outf, "l2:  drift, l = %7.5f;\n\n", bn_prms.bn[3]);
  fprintf(outf, "bh:  bending, l = 0.166667, t = 0.5, k = %8.5f, t1 = 0.0"
	  ", t2 = 0.0,\n     gap = 0.0, N = Nbend, Method = Meth;\n",
	  bn_prms.bn[1]);
  fprintf(outf, "qf:  quadrupole, l = 0.08, k = %8.5f, N = Nquad"
	  ", Method = Meth;\n", bn_prms.bn[2]);
  fclose(outf);
}


void fit_emit(param_type &bn_prms, const double eps_x,
	      const double nu_x, const double nu_y)
{
  // Optimize unit cell for hor. emittance.
  // Lattice: unit cell.
  const double scl_eps = 1e2, scl_nu = 1e1, scl_ksi = 1e-1;

  const int m_max = 8;

  int          n_b2, i, j, m;
  double       **A, *b;
  double       db2_max;
  tps          eps1_x;
  ss_vect<tps> nus;

  n_b2 = bn_prms.n_prm;
  printf("n_prm = %d\n", bn_prms.n_prm);

  b = dvector(1, m_max); A = dmatrix(1, m_max, 1, n_b2);

  bn_prms.ini_prm();

  danot_(3);

  do {
    printf("\n");
    for (i = 1; i <= n_b2; i++) {
      bn_prms.set_prm_dep(i-1);

      get_Map(); K = MapNorm(Map, g, A1, A0, Map_res, 1); nus = dHdJ(K);

      eps1_x = get_eps_x();
 
      m = 0;
      A[++m][i] = scl_eps*h_ijklm_p(eps1_x, 0, 0, 0, 0, 0, 7);
      A[++m][i] = scl_nu*h_ijklm_p(nus[3], 0, 0, 0, 0, 0, 7);
      A[++m][i] = scl_nu*h_ijklm_p(nus[4], 0, 0, 0, 0, 0, 7);
      // A[++m][i] = scl_ksi*h_ijklm_p(nus[3], 0, 0, 0, 0, 1, 7);
      // A[++m][i] = scl_ksi*h_ijklm_p(nus[4], 0, 0, 0, 0, 1, 7);

      for (j = 1; j <= m; j++)
	A[j][i] *= bn_prms.bn_scl[i-1];

      bn_prms.clr_prm_dep(i-1);
    }

    m = 0;
    b[++m] = -scl_eps*(eps1_x.cst()-eps_x);
    b[++m] = -scl_nu*(h_ijklm(nus[3], 0, 0, 0, 0, 0)-nu_x);
    b[++m] = -scl_nu*(h_ijklm(nus[4], 0, 0, 0, 0, 0)-nu_y);
    // b[++m] = -scl_ksi*h_ijklm(nus[3], 0, 0, 0, 0, 1);
    // b[++m] = -scl_ksi*h_ijklm(nus[4], 0, 0, 0, 0, 1);

    prt_system(m, n_b2, A, b);

    SVD_lim(m, n_b2, A, b, bn_prms.bn_lim, bn_prms.svd_cut, bn_prms.bn,
	    bn_prms.dbn);

    db2_max = bn_prms.set_dprm();

    prt_mfile("flat_file.dat");
  } while (db2_max > bn_prms.bn_tol);

  prt_emit(bn_prms);

  free_dvector(b, 1, m_max); free_dmatrix(A, 1, m_max, 1, n_b2);
  free_dvector(bn_prms.bn_lim, 1, n_b2); free_dvector(bn_prms.bn, 1, n_b2);
  free_dvector(bn_prms.dbn, 1, n_b2);
}


void prt_match(const param_type &bn_prms)
{
  double l4, l5h, l6, l7h, l8;
  FILE *outf;

  std::string file_name = "match.out";

  outf = file_write(file_name.c_str());

  fprintf(outf, "bm:  bending, l = 0.166667, t = 0.5, k = %8.5f, t1 = 0.0"
	  ", t2 = 0.0,\n     gap = 0.00, N = Nbend, Method = Meth;\n",
	  bn_prms.bn[1]);
  fprintf(outf, "qfe: quadrupole, l = 0.15, k = %9.5f, N = Nquad"
	  ", Method = Meth;\n", bn_prms.bn[2]);
  fprintf(outf, "qde: quadrupole, l = 0.1,  k = %9.5f, N = Nquad"
	  ", Method = Meth;\n", bn_prms.bn[3]);
  fprintf(outf, "qm:  quadrupole, l = 0.15, k = %9.5f, N = Nquad"
	  ", Method = Meth;\n", bn_prms.bn[4]);

  l4  = get_L(get_Fnum("l4"), 1);  l5h = get_L(get_Fnum("l5h"), 1);
  l6  = get_L(get_Fnum("l6"), 1);  l7h = get_L(get_Fnum("l7h"), 1);
  l8  = get_L(get_Fnum("l8"), 1);
  
  fprintf(outf, "\nl4:  drift, l = %7.5f;\n", l4);
  fprintf(outf, "l5h: drift, l = %7.5f;\n", l5h+bn_prms.bn[10]/2e0);
  fprintf(outf, "l6:  drift, l = %7.5f;\n", l6-bn_prms.bn[10]+bn_prms.bn[9]);
  fprintf(outf, "l7h: drift, l = %7.5f;\n",
	  l7h-bn_prms.bn[9]/2e0+bn_prms.bn[11]/2e0);
  fprintf(outf, "l8:  drift, l = %7.5f;\n", l8-bn_prms.bn[11]);

  fprintf(outf, "\nL_bm  = %8.5f\n", bn_prms.bn[9]);
  fprintf(outf, "L_qfe = %8.5f\n", bn_prms.bn[10]);
  fprintf(outf, "L_qde = %8.5f\n", bn_prms.bn[11]);

  fclose(outf);
}


void fit_match(param_type &bn_prms,
	       const double beta0_x, const double beta0_y, const double eta0_x,
	       const double beta1_x, const double beta1_y)
{
  // Match linear optics to straight section.
  // Lattice: matching cell.

  const int m_max = 8;

  int          n_b2, i, j, m, loc;
  double       **A, *b, db2_max;
  ss_vect<tps> AA_tp, A_disp;

  const double scl_eta  = 1e1,
               alpha0[] = {0e0,     0e0},
               beta0[]  = {beta0_x, beta0_y},
               eta0[]   = {eta0_x,  0e0},
	       etap0[]  = {0e0,     0e0};

  n_b2 = bn_prms.n_prm;

  b = dvector(1, m_max); A = dmatrix(1, m_max, 1, n_b2);

  bn_prms.ini_prm();

  danot_(3);

  loc = get_loc(get_Fnum("bm"), 1);

  do {
    for (i = 1; i <= n_b2; i++) {
      bn_prms.set_prm_dep(i-1);

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
	A[j][i] *= bn_prms.bn_scl[i-1];

      bn_prms.clr_prm_dep(i-1);
    }

    m = 0;
    b[++m] = -h_ijklm(-AA_tp[x_], 0, 1, 0, 0, 0);
    b[++m] = -h_ijklm(-AA_tp[y_], 0, 0, 0, 1, 0);
    b[++m] = -scl_eta*h_ijklm(A_disp[x_],  0, 0, 0, 0, 1);
    b[++m] = -scl_eta*h_ijklm(A_disp[px_], 0, 0, 0, 0, 1);
    b[++m] = -(h_ijklm(AA_tp[x_], 1, 0, 0, 0, 0)-beta1_x);
    b[++m] = -(h_ijklm(AA_tp[y_], 0, 0, 1, 0, 0)-beta1_y);

    prt_system(m, n_b2, A, b);

    SVD_lim(m, n_b2, A, b, bn_prms.bn_lim, bn_prms.svd_cut, bn_prms.bn,
	    bn_prms.dbn);

    db2_max = bn_prms.set_dprm();
  } while (db2_max > bn_prms.bn_tol);

  prt_match(bn_prms);

  free_dvector(b, 1, m_max); free_dmatrix(A, 1, m_max, 1, n_b2);
  free_dvector(bn_prms.bn_lim, 1, n_b2); free_dvector(bn_prms.bn, 1, n_b2);
  free_dvector(bn_prms.dbn, 1, n_b2);
}


double f_match(double *b2)
{
  // Minimize linear linear chromaticity for super period.
  // Lattice: super period.

  static double chi2_ref = 1e30;

  int          i, loc1, loc2, loc3;
  double       eps1_x, tr[2], chi2;
  tps          K_re, K_im;
  ss_vect<tps> AA_tp1, AA_tp2, A_disp;

  const double eps_x   = 16e-3,
               scl_eps = 1e-1, scl_eta = 1e4, scl_ksi = 1e-1,
	       beta0[] = {3.0,     3.0},
               beta1[] = {1.35633, 1.91479};

  dvcopy(b2, bn_prms.n_prm, bn_prms.bn);
  bn_prms.set_prm();

  // Exit of 2nd BM.
  loc1 = get_loc(get_Fnum("bm"), 2);
  loc2 = get_loc(get_Fnum("sfh"), 1);
  loc3 = n_elem;

  danot_(1);
  eps1_x = get_eps_x().cst();

  danot_(2);
  get_Map();
  danot_(3);
  K = MapNorm(Map, g, A1, A0, Map_res, 1);
  CtoR(K, K_re, K_im);
  // Call to get_twiss sets no to 2.
  get_twiss(A1);
  A_disp = elem_tps[loc1-1].A1;
  AA_tp1 = elem_tps[loc2-1].A1*tp_S(2, elem_tps[loc2-1].A1);
  AA_tp2 = elem_tps[loc3-1].A1*tp_S(2, elem_tps[loc3-1].A1);

  tr[X_] = h_ijklm(Map[x_], 1, 0, 0, 0, 0)+h_ijklm(Map[px_], 0, 1, 0, 0, 0);
  tr[Y_] = h_ijklm(Map[y_], 0, 0, 1, 0, 0)+h_ijklm(Map[py_], 0, 0, 0, 1, 0);
  // printf("trace: %6.3f %6.3f\n", tr[X_], tr[Y_]);

  chi2 = 0e0;
  chi2 += sqr(scl_eps*(eps1_x-eps_x));
  chi2 += sqr(scl_eta*h_ijklm(A_disp[x_],  0, 0, 0, 0, 1));
  chi2 += sqr(scl_eta*h_ijklm(A_disp[px_], 0, 0, 0, 0, 1));
  chi2 += sqr(h_ijklm(AA_tp1[x_], 1, 0, 0, 0, 0)-beta1[X_]);
  chi2 += sqr(h_ijklm(AA_tp1[y_], 0, 0, 1, 0, 0)-beta1[Y_]);
  chi2 += sqr(h_ijklm(AA_tp2[x_], 1, 0, 0, 0, 0)-beta0[X_]);
  chi2 += sqr(h_ijklm(AA_tp2[y_], 0, 0, 1, 0, 0)-beta0[Y_]);
  chi2 += sqr(scl_ksi*(h_ijklm(K_re, 1, 1, 0, 0, 1)));
  chi2 += sqr(scl_ksi*(h_ijklm(K_re, 0, 0, 1, 1, 1)));
  if ((fabs(tr[X_]) > 2e0) || (fabs(tr[Y_]) > 2e0)) chi2 += 1e10;
  for (i = 1; i <= bn_prms.n_prm; i++)
    if (fabs(b2[i]) > bn_prms.bn_max[i-1]) chi2 += 1e10;

  if (chi2 < chi2_ref) {
    printf("\nchi2: %12.5e, %12.5e\n", chi2, chi2_ref);
    printf("b:    %10.3e %10.3e %10.3e %8.3f %8.3f %8.3f %8.3f %10.3e %10.3e\n",
	   eps1_x,
	   h_ijklm(A_disp[x_],  0, 0, 0, 0, 1),
	   h_ijklm(A_disp[px_], 0, 0, 0, 0, 1),
	   h_ijklm(AA_tp1[x_], 1, 0, 0, 0, 0),
	   h_ijklm(AA_tp1[y_], 0, 0, 1, 0, 0),
	   h_ijklm(AA_tp2[x_], 1, 0, 0, 0, 0),
	   h_ijklm(AA_tp2[y_], 0, 0, 1, 0, 0),
	   h_ijklm(K_re, 1, 1, 0, 0, 1),
	   h_ijklm(K_re, 0, 0, 1, 1, 1));
    printf("b2s: ");
    for (i = 1; i <= bn_prms.n_prm; i++)
      printf("%9.5f", b2[i]);
    printf("\n");

    prt_mfile("flat_file.dat");
    prt_match(bn_prms);

    get_S();
    prt_lat("linlat.out");
  }

  chi2_ref = min(chi2, chi2_ref);

  return chi2;
}


void opt_match(param_type &bn_prms)
{
  // Minimize linear linear chromaticity for super period.
  // Lattice: super period.

  int    n_b2, i, j, iter;
  double **xi, fret;

  n_b2 = bn_prms.n_prm;

  xi = dmatrix(1, n_b2, 1, n_b2);

  bn_prms.ini_prm();

  // Set initial directions (unit vectors).
  for (i = 1; i <= n_b2; i++)
    for (j = 1; j <= n_b2; j++)
      xi[i][j] = (i == j)? 1e0 : 0e0;

  dpowell(bn_prms.bn, xi, n_b2, bn_prms.bn_tol, &iter, &fret, f_match);

  free_dmatrix(xi, 1, n_b2, 1, n_b2);
}


void fit_3rd_achrom(param_type &bn_prms)
{
  // Third order achromat.

  int          n_b4, i, j, m;
  double       **A, *b;
  double       db4_max;
  tps          K_re, K_im;
  ss_vect<tps> nus;

  const int m_max = 6;

  n_b4 = bn_prms.n_prm;

  b = dvector(1, m_max); A = dmatrix(1, m_max, 1, n_b4);

  bn_prms.ini_prm();

  danot_(5);

  do {
    printf("\n");
    for (i = 1; i <= n_b4; i++) {
      bn_prms.set_prm_dep(i-1);

      get_Map(); K = MapNorm(Map, g, A1, A0, Map_res, 1); CtoR(K, K_re, K_im);

      m = 0;
      A[++m][i] = h_ijklm_p(K_re, 2, 2, 0, 0, 0, 7);
      A[++m][i] = h_ijklm_p(K_re, 0, 0, 2, 2, 0, 7);
      A[++m][i] = h_ijklm_p(K_re, 1, 1, 1, 1, 0, 7);
      // A[++m][i] = h_ijklm_p(K_re, 1, 1, 0, 0, 2, 7);
      // A[++m][i] = h_ijklm_p(K_re, 0, 0, 1, 1, 2, 7);

      for (j = 1; j <= m; j++)
	A[j][i] *= bn_prms.bn_scl[i-1];

      bn_prms.clr_prm_dep(i-1);
    }

    m = 0;
    b[++m] = -h_ijklm(K_re, 2, 2, 0, 0, 0);
    b[++m] = -h_ijklm(K_re, 0, 0, 2, 2, 0);
    b[++m] = -h_ijklm(K_re, 1, 1, 1, 1, 0);
    // b[++m] = -h_ijklm(K_re, 1, 1, 0, 0, 2);
    // b[++m] = -h_ijklm(K_re, 0, 0, 1, 1, 2);

    prt_system(m, n_b4, A, b);

    SVD_lim(m, n_b4, A, b, bn_prms.bn_lim, bn_prms.svd_cut, bn_prms.bn,
	    bn_prms.dbn);

    db4_max = bn_prms.set_dprm();
  } while (db4_max > bn_prms.bn_tol);

  free_dvector(b, 1, m_max); free_dmatrix(A, 1, m_max, 1, n_b4);
  free_dvector(bn_prms.bn_lim, 1, n_b4); free_dvector(bn_prms.bn, 1, n_b4);
  free_dvector(bn_prms.dbn, 1, n_b4);
}


void min_dnu_prt(const param_type &bn_prms)
{
  FILE *outf;

  std::string file_name = "dnu.out";

  outf = file_write(file_name.c_str());
  fprintf(outf, "o1: multipole, l = 0.0, N = Nsext, Method = Meth,"
	        "\n    HOM = (4, %12.5e, 0.0);\n", bn_prms.bn[1]);
  fprintf(outf, "o2: multipole, l = 0.0, N = Nsext, Method = Meth,"
	        "\n    HOM = (4, %12.5e, 0.0);\n", bn_prms.bn[2]);
  fprintf(outf, "o3: multipole, l = 0.0, N = Nsext, Method = Meth,"
	        "\n    HOM = (4, %12.5e, 0.0);\n", bn_prms.bn[3]);
  fprintf(outf, "o4: multipole, l = 0.0, N = Nsext, Method = Meth,"
	        "\n    HOM = (4, %12.5e, 0.0);\n", bn_prms.bn[4]);
  fclose(outf);
}


void min_dnu_prt2(const param_type &bn_prms)
{
  FILE *outf;

  std::string file_name = "dnu.out";

  outf = file_write(file_name.c_str());
  fprintf(outf, "o1: multipole, l = 0.0, N = Nsext, Method = Meth,"
	        "\n    HOM = (4, %12.5e, 0.0, 6, %12.5e, 0.0);\n",
	  bn_prms.bn[1], bn_prms.bn[5]);
  fprintf(outf, "o2: multipole, l = 0.0, N = Nsext, Method = Meth,"
	        "\n    HOM = (4, %12.5e, 0.0, 6, %12.5e, 0.0);\n",
	  bn_prms.bn[2], bn_prms.bn[6]);
  fprintf(outf, "o3: multipole, l = 0.0, N = Nsext, Method = Meth,"
	        "\n    HOM = (4, %12.5e, 0.0, 6, %12.5e, 0.0);\n",
	  bn_prms.bn[3], bn_prms.bn[7]);
  fprintf(outf, "o4: multipole, l = 0.0, N = Nsext, Method = Meth,"
	        "\n    HOM = (4, %12.5e, 0.0, 6, %12.5e, 0.0);\n",
	  bn_prms.bn[4], bn_prms.bn[8]);
  fclose(outf);
}


double get_chi2(const int n, const double data[])
{
  int    j;
  double chi2 = 0e0;

  for (j = 1; j <= n; j++)
    chi2 += sqr(data[j]);

  return chi2;
}


double min_dnu_f(double *b4s)
{
  int                 i;
  static double       chi2_ref = 1e30;
  double              chi2, twoJ[2];
  tps                 K_re, K_im;
  ss_vect<tps>        nus;
  std::vector<double> b;

  n_powell++;

  // Do not change parameters.
  for (i = 1; i <= bn_prms.n_prm; i++)
    set_bn(bn_prms.Fnum[i-1], bn_prms.n[i-1], b4s[i]);

  danot_(5);
  get_Map();
  danot_(6);
  K = MapNorm(Map, g, A1, A0, Map_res, 1);
  CtoR(K, K_re, K_im); nus = dHdJ(K);
  K_re = K_re*Id_scl; nus[3] = nus[3]*Id_scl; nus[4] = nus[4]*Id_scl;

  twoJ[X_] = 1e0; twoJ[Y_] = 1e0;
 
  b.push_back(scl_dnu*h_ijklm(K_re, 2, 2, 0, 0, 0));
  b.push_back(scl_dnu*h_ijklm(K_re, 0, 0, 2, 2, 0));
  b.push_back(scl_dnu*h_ijklm(K_re, 1, 1, 1, 1, 0));

  b.push_back(scl_dnu*h_ijklm(K_re, 3, 3, 0, 0, 0));
  b.push_back(scl_dnu*h_ijklm(K_re, 0, 0, 3, 3, 0));
  b.push_back(scl_dnu*h_ijklm(K_re, 2, 2, 1, 1, 0));
  b.push_back(scl_dnu*h_ijklm(K_re, 1, 1, 2, 2, 0));

  b.push_back(h_ijklm(nus[3], 1, 1, 0, 0, 0)
	      +2e0*h_ijklm(nus[3], 2, 2, 0, 0, 0)*twoJ[X_]);
  b.push_back(scl_2d*(h_ijklm(nus[3], 0, 0, 1, 1, 0)
		      +2e0*h_ijklm(nus[3], 0, 0, 2, 2, 0)*twoJ[Y_]));
  b.push_back(h_ijklm(nus[4], 0, 0, 1, 1, 0)
	      +2e0*h_ijklm(nus[4], 0, 0, 2, 2, 0)*twoJ[Y_]);
  b.push_back(scl_2d*(h_ijklm(nus[4], 1, 1, 0, 0, 0)
		      +2e0*h_ijklm(nus[4], 2, 2, 0, 0, 0)*twoJ[X_]));

  chi2 = 0e0;
  for (i = 0; i < (int)b.size(); i++)
    chi2 += sqr(b[i]);

  if (chi2 < chi2_ref) {
    printf("\n%3d chi2: %12.5e -> %12.5e\n", n_powell, chi2_ref, chi2);
    printf("b:\n");
    for (i = 0; i < (int)b.size(); i++)
      printf("%11.3e", b[i]);
    printf("\n");
    chi2 += sqr(b[i]);
    printf("bn:\n");
    for (i = 1; i <= bn_prms.n_prm; i++) 
      printf("%11.3e", b4s[i]);
    printf("\n");
  }

  chi2_ref = min(chi2, chi2_ref);

  return chi2;
}


void conj_grad(const int n_iter, double bn[], double dbn[],
	       double *g, double *h, double (*f)(double *))
{
  // Conjugate gradient method.
  int    n_bn, i;
  double fret, g2, gamma, dg2;

  n_bn = bn_prms.n_prm;

  if (n_iter == 1) {
    for (i = 1; i <= n_bn; i++) {
      g[i] = dbn[i]; h[i] = g[i];
    }
  } else {
    dg2 = g2 = 0e0;
    for (i = 1; i <= n_bn; i++) {
      g2 += sqr(g[i]); dg2 += (dbn[i]-g[i])*dbn[i];
    }
    if (g2 != 0e0) {
      gamma = dg2/g2;
      for (i = 1; i <= n_bn; i++) {
	g[i] = dbn[i]; dbn[i] = h[i] = g[i] + gamma*h[i];
      }
    } else {
      printf("g.g = 0\n");
      exit(0);
    }
  }

  printf("\n");
  d_linmin(bn, dbn, n_bn, &fret, f);
}


void min_dnu_grad(double &chi2, double &db4_max, double *g_, double *h_,
		  const bool cg_meth)
{
  int          n_b4, i, j, m;
  double       chi2_ref;
  double       **A, *b, *bn_ref;
  double       twoJ[2];
  tps          K_re, K_im;
  ss_vect<tps> nus;

  const int m_max = 15;

  n_b4 = bn_prms.n_prm;

  bn_ref = dvector(1, n_b4);
  b = dvector(1, m_max); A = dmatrix(1, m_max, 1, n_b4);

  twoJ[X_] = 1e0; twoJ[Y_] = 1e0;

  printf("\n");
  for (i = 1; i <= n_b4; i++) {
    bn_prms.set_prm_dep(i-1);

    danot_(6);
    get_Map();
    danot_(7);
    K = MapNorm(Map, g, A1, A0, Map_res, 1);
    CtoR(K, K_re, K_im); nus = dHdJ(K);
    K_re = K_re*Id_scl; nus[3] = nus[3]*Id_scl; nus[4] = nus[4]*Id_scl;

    m = 0;
    A[++m][i] = scl_dnu*h_ijklm_p(K_re, 2, 2, 0, 0, 0, 7);
    A[++m][i] = scl_dnu*h_ijklm_p(K_re, 0, 0, 2, 2, 0, 7);
    A[++m][i] = scl_dnu*h_ijklm_p(K_re, 1, 1, 1, 1, 0, 7);

    A[++m][i] = scl_dnu*h_ijklm_p(K_re, 3, 3, 0, 0, 0, 7);
    A[++m][i] = scl_dnu*h_ijklm_p(K_re, 0, 0, 3, 3, 0, 7);
    A[++m][i] = scl_dnu*h_ijklm_p(K_re, 2, 2, 1, 1, 0, 7);
    A[++m][i] = scl_dnu*h_ijklm_p(K_re, 1, 1, 2, 2, 0, 7);

    A[++m][i] =
      h_ijklm_p(nus[3], 1, 1, 0, 0, 0, 7)
      +2e0*h_ijklm_p(nus[3], 2, 2, 0, 0, 0, 7)*twoJ[X_];
    A[++m][i] =
      scl_2d*(h_ijklm_p(nus[3], 0, 0, 1, 1, 0, 7)
	      +2e0*h_ijklm_p(nus[3], 0, 0, 2, 2, 0, 7)*twoJ[Y_]);
    A[++m][i] =
       h_ijklm_p(nus[4], 0, 0, 1, 1, 0, 7)
      +2e0*h_ijklm_p(nus[4], 0, 0, 2, 2, 0, 7)*twoJ[Y_];
    A[++m][i] =
      scl_2d*(h_ijklm_p(nus[4], 1, 1, 0, 0, 0, 7)
	      +2e0*h_ijklm_p(nus[4], 2, 2, 0, 0, 0, 7)*twoJ[X_]);

    for (j = 1; j <= m; j++)
      A[j][i] *= bn_prms.bn_scl[i-1];

    bn_prms.clr_prm_dep(i-1);
  }

  m = 0;
  b[++m] = -scl_dnu*h_ijklm(K_re, 2, 2, 0, 0, 0);
  b[++m] = -scl_dnu*h_ijklm(K_re, 0, 0, 2, 2, 0);
  b[++m] = -scl_dnu*h_ijklm(K_re, 1, 1, 1, 1, 0);

  b[++m] = -scl_dnu*h_ijklm(K_re, 3, 3, 0, 0, 0);
  b[++m] = -scl_dnu*h_ijklm(K_re, 0, 0, 3, 3, 0);
  b[++m] = -scl_dnu*h_ijklm(K_re, 2, 2, 1, 1, 0);
  b[++m] = -scl_dnu*h_ijklm(K_re, 1, 1, 2, 2, 0);

  b[++m] =
    -(h_ijklm(nus[3], 1, 1, 0, 0, 0)
      +2e0*h_ijklm(nus[3], 2, 2, 0, 0, 0)*twoJ[X_]);
  b[++m] =
    -scl_2d*(h_ijklm(nus[3], 0, 0, 1, 1, 0)
	     +2e0*h_ijklm(nus[3], 0, 0, 2, 2, 0)*twoJ[Y_]);
  b[++m] =
    -(h_ijklm(nus[4], 0, 0, 1, 1, 0)
      +2e0*h_ijklm(nus[4], 0, 0, 2, 2, 0)*twoJ[Y_]);
  b[++m] =
    -scl_2d*(h_ijklm(nus[4], 1, 1, 0, 0, 0)
	     +2e0*h_ijklm(nus[4], 2, 2, 0, 0, 0)*twoJ[X_]);

  chi2_ref = chi2; chi2 = get_chi2(m, b);

  prt_system(m, n_b4, A, b);
  printf("\n%4d %12.5e -> %12.5e\n", n_iter, chi2_ref, chi2);

  SVD_lim(m, n_b4, A, b, bn_prms.bn_lim, bn_prms.svd_cut, bn_prms.bn,
	  bn_prms.dbn);

  dvcopy(bn_prms.bn, n_b4, bn_ref);
  if (cg_meth)
    conj_grad(n_iter, bn_prms.bn, bn_prms.dbn, g_, h_, min_dnu_f);
  else
    bn_prms.set_dprm();

  for (i = 1; i <= n_b4; i++)
    bn_prms.bn[i] = get_bn(bn_prms.Fnum[i-1], 1, bn_prms.n[i-1]);

  printf("\nbn & dbn:\n");
  for (i = 1; i <= n_b4; i++)
    printf("%13.5e", bn_prms.bn[i]);
  printf("\n");
  db4_max = 0e0;
  for (i = 1; i <= n_b4; i++) {
    db4_max = max(fabs((bn_prms.bn[i]-bn_ref[i])), db4_max);
    printf("%13.5e", bn_prms.bn[i]-bn_ref[i]);
  }
  printf("\n");

  free_dvector(bn_ref, 1, n_b4);
  free_dvector(b, 1, m_max); free_dmatrix(A, 1, m_max, 1, n_b4);
}


void min_dnu(const bool cg_meth)
{
  // Control tune foot print; conjugate gradient method.
  std::string str;
  int         n_b4;
  double      *g, *h, dbn_max, chi2;

  const int    n_iter_max = 1000;
  const double bn_tol     = 1e-5;

  n_b4 = bn_prms.n_prm;

  g = dvector(1, n_b4); h = dvector(1, n_b4);

  bn_prms.ini_prm();

  n_iter = 0; chi2 = 0e0;
  do {
    n_iter++;

    min_dnu_grad(chi2, dbn_max, g, h, cg_meth);

    prt_mfile("flat_file.fit");
    if (n_b4 <= 4)
      min_dnu_prt(bn_prms);
    else
      min_dnu_prt2(bn_prms);
  } while ((dbn_max > bn_tol) && (n_iter < n_iter_max));

  free_dvector(g, 1, n_b4); free_dvector(h, 1, n_b4);
  free_dvector(bn_prms.bn_lim, 1, n_b4); free_dvector(bn_prms.bn, 1, n_b4);
  free_dvector(bn_prms.dbn, 1, n_b4);
}


void min_dnu2(void)
{
  int    n_b4, i, j, iter;
  double **xi, fret;

  n_b4 = bn_prms.n_prm;

  xi = dmatrix(1, n_b4, 1, n_b4);

  bn_prms.ini_prm();

  // Set initial directions (unit vectors).
  for (i = 1; i <= n_b4; i++)
    for (j = 1; j <= n_b4; j++)
      xi[i][j] = (i == j)? 1e0 : 0e0;

  n_iter = 0;
  dpowell(bn_prms.bn, xi, n_b4, bn_prms.bn_tol, &iter, &fret, min_dnu_f);

  free_dmatrix(xi, 1, n_b4, 1, n_b4);
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
  int           j;
  double        twoJ[2];
  tps           eps1_x, K_re, K_im, g_re, g_im, h_re, h_im, H_re, H_im;
  std::ofstream outf;
  
  const double eps_x   = 15e-3,
               nu[]    = {4.0/15.0, 1.0/15.0},
               beta0[] = {0.38133, 5.00190},   // Center of dipole.
	       eta0_x  = 0.00382435, 
	       beta1[] = {1.35633, 1.91479},   // Center of SFh.
	       beta2[] = {3.0, 3.0},           // Center of straight.
               A_max[] = {1.0e-3, 1.0e-3}, delta = 3e-2;

  rad_on    = false; H_exact        = false; totpath_on   = false;
  cavity_on = false; quad_fringe_on = false; emittance_on = false;
  IBS_on    = false;

  rd_mfile(argv[1], elem); rd_mfile(argv[1], elem_tps);
  
  // Initialize the symplectic integrator after energy has been defined.
  ini_si();

  // Disable log messages from TPSALib and LieLib.
  idprset(-1);

  daeps_(1e-25);

  danot_(1);

  get_nu_ksi();

  // tst_eps_x();
  // eps1_x = get_eps_x();
  // cout << eps_x << "\n"; 

  if (false) {
    tst_dL(beta0[X_], beta0[Y_], eta0_x);
    exit(0);
  }

  if (false) {
    bn_prms.add_prm("bh",  2, 30.0, 1.0);
    bn_prms.add_prm("qfh",  2, 30.0, 1.0);
    bn_prms.add_prm("l2", -1,  0.5, 0.01);

    bn_prms.bn_tol = 1e-5; bn_prms.svd_cut = 1e-8;  bn_prms.step = 0.5;

    fit_emit(bn_prms, eps_x, nu[X_], nu[Y_]);
  }

  if (false) {
    bn_prms.add_prm("bm",   2, 25.0, 1.0);
    bn_prms.add_prm("qfe",  2, 25.0, 1.0);
    bn_prms.add_prm("qde",  2, 50.0, 1.0);
    bn_prms.add_prm("qm",   2, 25.0, 1.0);

    bn_prms.add_prm("bm",  -2,  0.2, 1.0);
    bn_prms.add_prm("qfe", -2,  0.2, 1.0);
    bn_prms.add_prm("qde", -2,  0.2, 1.0);
    bn_prms.add_prm("qm",  -2,  0.2, 1.0);

    bn_prms.bn_tol = 1e-6; bn_prms.svd_cut = 1e-6;  bn_prms.step = 0.5;

    fit_match(bn_prms, beta0[X_], beta0[Y_], eta0_x, beta1[X_], beta1[Y_]);
  }

  if (false) {
    bn_prms.add_prm("bm",   2, 25.0, 1.0);
    bn_prms.add_prm("qfe",  2, 25.0, 1.0);
    bn_prms.add_prm("qde",  2, 30.0, 1.0);
    bn_prms.add_prm("qm",   2, 25.0, 1.0);

    bn_prms.add_prm("l5h", -1, 0.2,  0.1);
    bn_prms.add_prm("l6",  -1, 0.5,  0.1);
    bn_prms.add_prm("l7h", -1, 0.5,  0.1);
    bn_prms.add_prm("l8",  -1, 0.5,  0.1);

    bn_prms.bn_tol = 1e-15; bn_prms.svd_cut = 0e0; bn_prms.step = 0.0;

    no_mpoles(Sext);
    opt_match(bn_prms);
  }

  if (false) {
    bn_prms.add_prm("o1", 4, 1e6, 1.0);
    bn_prms.add_prm("o2", 4, 1e6, 1.0);
    bn_prms.add_prm("o3", 4, 1e6, 1.0);

    bn_prms.bn_tol = 1e-5; bn_prms.svd_cut = 1e-3; bn_prms.step = 1.0;

    no_mpoles(Oct); no_mpoles(Dodec);

    fit_3rd_achrom(bn_prms);

    prt_mfile("flat_file.fit");
  }

  if (false) {
    for (j = 0; j < 2; j++)
      twoJ[j] = sqr(A_max[j])/beta2[j];

    Id_scl.identity();
    Id_scl[x_] *= sqrt(twoJ[X_]); Id_scl[px_] *= sqrt(twoJ[X_]);
    Id_scl[y_] *= sqrt(twoJ[Y_]); Id_scl[py_] *= sqrt(twoJ[Y_]);
    Id_scl[delta_] *= delta;

    bn_prms.add_prm("o1", 4, 1e6, 1.0);
    bn_prms.add_prm("o2", 4, 1e6, 1.0);
    bn_prms.add_prm("o3", 4, 1e6, 1.0);
    bn_prms.add_prm("o4", 4, 1e6, 1.0);

    bn_prms.add_prm("o1", 6, 1e9, 1.0);
    bn_prms.add_prm("o2", 6, 1e9, 1.0);
    bn_prms.add_prm("o3", 6, 1e9, 1.0);
    bn_prms.add_prm("o4", 6, 1e9, 1.0);

    bn_prms.bn_tol = 1e-4; bn_prms.svd_cut = 1e-10; bn_prms.step = 0.01;

    // no_mpoles(Oct); no_mpoles(Dodec);

    min_dnu(false);
  }

  if (true) {
    for (j = 0; j < 2; j++)
      twoJ[j] = sqr(A_max[j])/beta2[j];

    Id_scl.identity();
    Id_scl[x_] *= sqrt(twoJ[X_]); Id_scl[px_] *= sqrt(twoJ[X_]);
    Id_scl[y_] *= sqrt(twoJ[Y_]); Id_scl[py_] *= sqrt(twoJ[Y_]);
    Id_scl[delta_] *= delta;

    bn_prms.add_prm("o1", 4, 1e6, 1.0);
    bn_prms.add_prm("o2", 4, 1e6, 1.0);
    bn_prms.add_prm("o3", 4, 1e6, 1.0);
    bn_prms.add_prm("o4", 4, 1e6, 1.0);

    // bn_prms.add_prm("o1", 6, 1e9, 1.0);
    // bn_prms.add_prm("o2", 6, 1e9, 1.0);
    // bn_prms.add_prm("o3", 6, 1e9, 1.0);
    // bn_prms.add_prm("o4", 6, 1e9, 1.0);

    bn_prms.bn_tol = 1e-4; bn_prms.svd_cut = 1e-14; bn_prms.step = 1.0;

    no_mpoles(Oct); no_mpoles(Dodec);

    min_dnu2();
  }

  if (false) {
    danot_(NO);

    get_Map(); K = MapNorm(Map, g, A1, A0, Map_res, no_tps);
    CtoR(K*Id_scl, K_re, K_im); CtoR(g*Id_scl, g_re, g_im);
    CtoR(get_h()*Id_scl, h_re, h_im);
    CtoR(get_H()*Id_scl, H_re, H_im);

    file_wr(outf, "map.dat"); outf << Map; outf.close();
    file_wr(outf, "K.dat");   outf << K_re*Id_scl; outf.close();
    file_wr(outf, "g.dat");   outf << g_im*Id_scl; outf.close();
    file_wr(outf, "h.dat");   outf << h_re*Id_scl; outf.close();
    file_wr(outf, "H.dat");   outf << H_re*Id_scl; outf.close();
  }
}
