
#define NO 5

#include "thor_lib.h"

#include "param_type.h"

int
  no_tps   = NO,
  ndpt_tps = 5;


extern double       b2_max;
extern tps          K, g;
extern ss_vect<tps> Map, A0, A1, Map_res;


typedef struct {
  std::string
    label;
  double
    cst_scl,
    cst;
  std::vector<double>
    Jacobian;
} Lie_term;


const double
#if 0
  beta_inj[] = {2.7, 2.5},
#else
  beta_inj[] = {2.8, 2.6},
#endif
  A_max[]    = {3e-3, 1.5e-3},
  delta_max  = 2e-2,
  twoJ[]     = {sqr(A_max[X_])/beta_inj[X_], sqr(A_max[Y_])/beta_inj[Y_]};




tps get_h(ss_vect<tps> &map)
{
  // Parallel transport nonlinear kick to start of lattice.
  tps          h;
  ss_vect<tps> map1, R;

  if (true) {
    // Dragt-Finn factorization
    h = LieFact_DF(Inv(A0*A1)*map*A0*A1, R);
    R[6] = tps(0e0, 7); h = h*R;
    return h;
  } else {
    // single Lie exponent
    danot_(1); map1 = map; danot_(no_tps);
    return LieFact(Inv(A0*A1)*map*Inv(map1)*A0*A1);
  }
}


tps get_H(const tps &g)
{
  int           i;
  tps           H, gn;
  ss_vect<tps>  Id, Mn;

  const bool prt = false;

  // Compute the Lie generator.
  // K is in Dragt-Finn form but the Lie generators commute.
  Id.identity();
  H = K;
  for (i = no_tps; i >= 3; i--) {
    gn = Take(g, i); H = H*LieExp(-gn, Id);
  }

  if (prt) std::cout << (H-H*Inv(A0*A1)*Map*A0*A1);

  if (false) {
    // Normalize map (=> Map_res).
    for (i = 3; i <= no_tps; i++) {
      gn = Take(g, i); Mn = LieExp(gn, Id); Map = Inv(Mn)*Map*Mn;
    }
  }

  return H;
}


void chk_lat(ss_vect<tps> &map, ss_vect<tps> &map_res,
	     double nu[], double ksi[])
{
  double       alpha[2], beta[2];
  ss_vect<tps> nus;

  danot_(2);
  get_Map();
  danot_(3);
  K = MapNorm(map, g, A1, A0, map_res, 1);
  nus = dHdJ(K);
  get_nu_ksi(nus, nu, ksi);
  get_ab(alpha, beta, 0);
  printf("\n  alpha = [%6.3f, %6.3f]\n"
  	 "  beta  = [%6.3f, %6.3f]\n"
  	 "  nu    = [%6.3f, %6.3f]\n"
  	 "  ksi   = [%6.3f, %6.3f]\n",
  	 alpha[X_], alpha[Y_], beta[X_], beta[Y_], nu[X_], nu[Y_],
  	 ksi[X_], ksi[Y_]);
}


void prt_Lie_term(const Lie_term &k_ijklm)
{
  int k;

  printf(" %s %6.1e %10.3e",
	 k_ijklm.label.c_str(), k_ijklm.cst_scl, k_ijklm.cst);
  for (k = 0; k < (int)k_ijklm.Jacobian.size(); k++)
    printf(" %10.3e", k_ijklm.Jacobian[k]);
  printf("\n");
}


void prt_K_ijklm(const param_type &bns, const std::vector<Lie_term> &k_ijklm,
		 const bool tune_fp)
{
  int k;

  printf("\n           scl.      cst.");
  for (k = 0; k < bns.n_prm; k++)
    printf("      %-5s", bns.name[k].c_str());
  printf("\n                         ");
  for (k = 0; k < bns.n_prm; k++)
    printf("       %1d   ", bns.n[k]);
  printf("\n                         ");
  for (k = 0; k < bns.n_prm; k++)
    printf("    %7.1e", bns.bnL_scl[k]);

  printf("\nLinear chromaticity:\n");
  for (k = 0; k < 2; k++)
    prt_Lie_term(k_ijklm[k]);

  printf("\nChromatic terms:\n");
  for (k = 2; k < 5; k++)
    prt_Lie_term(k_ijklm[k]);

  printf("\nGeometric terms:\n");
  for (k = 5; k < 10; k++)
    prt_Lie_term(k_ijklm[k]);

  if (tune_fp) {
    printf("\nAnharmonic terms:\n");
    for (k = 10; k < 13; k++)
      prt_Lie_term(k_ijklm[k]);

    printf("\nAnharmonic terms:\n");
    for (k = 13; k < 17; k++)
      prt_Lie_term(k_ijklm[k]);

    printf("\n2nd order chromaticity:\n");
    for (k = 17; k < 19; k++)
      prt_Lie_term(k_ijklm[k]);

    printf("\n3rd order chromaticity:\n");
    for (k = 19; k < 21; k++)
      prt_Lie_term(k_ijklm[k]);
  }
}


void get_A(const int m, const int n, const std::vector<Lie_term> &k_ijklm,
	   double **A, double *w, double **U, double **V, double *b)
{
  int j, k;

  for (j = 0; j < m; j++) {
    b[j+1] = -k_ijklm[j].cst;
    for (k = 0; k < n; k++)
      A[j+1][k+1] = k_ijklm[j].Jacobian[k];
  }

  if (false) dmdump(stdout, (char *)"\nA:", A, m, n, (char *)" %10.3e");
}


void get_sing_val(const int n, double w[], const double svd_cut)
{
  int k;

  printf("\nsingular values:\n");
  for (k = 1; k <= n; k++) {
    printf("  %9.3e", w[k]);
    if (w[k] < svd_cut) {
      w[k] = 0e0;
      printf(" (zeroed)");
    }
    if (k % 6 == 0) printf("\n");
  }
  if (n % 6 != 0) printf("\n");
}


void set_bnL(const double *dbnL, param_type &bns)
{
  int    k;
  double bnL_ext;

  printf("\n      b_n       b_n*L\n");
  for (k = 0; k < bns.n_prm; k++) {
    bnL_ext =
      get_bnL(bns.Fnum[k], 1, bns.n[k]) + bns.bnL_scl[k]*dbnL[k+1];
    bns.bnL[k] = bnL_internal(bnL_ext, bns.bnL_min[k], bns.bnL_max[k]);
    if (false)
      printf("  %1d %11.3e %11.3e %11.3e\n",
	     k, bnL_ext, bns.bnL[k],
	     bnL_bounded(bns.bnL[k], bns.bnL_min[k], bns.bnL_max[k]));

    set_bnL(bns.Fnum[k], bns.n[k],
	    bnL_bounded(bns.bnL[k], bns.bnL_min[k], bns.bnL_max[k]));
    printf("  %10.3e %10.3e\n",
	   get_bn(bns.Fnum[k], 1, bns.n[k]), get_bnL(bns.Fnum[k], 1, bns.n[k]));
  }
}


void prt_bend(FILE *outf, const int loc, const int n)
{
  const elem_type<double> *elemp = &elem[loc];
  const int               Fnum = elemp->Fnum;

  fprintf(outf,
	  "%-8s: multipole, l = %7.5f, t = %7.5f, t1 = %7.5f, t2 = %7.5f,\n"
	  "          hom = (%d, %12.5e, 0e0,"
	  " %d, %12.5e, 0e0),\n"
	  "          n = nbend, Method = Meth;\n",
	  elemp->Name, elemp->L,
	  elemp->L*elemp->mpole->h_bend*180e0/M_PI,
	  elemp->mpole->edge1, elemp->mpole->edge2,
	  Quad, get_bn(Fnum, 1, Quad),
	  n, get_bn(Fnum, 1, n));
}


void prt_quad(FILE *outf, const int loc, const int n)
{
  const int Fnum = elem[loc].Fnum;

  fprintf(outf,
	  "%-8s: multipole, l = %7.5f,\n"
	  "          hom = (%d, %12.5e, 0e0,"
	  " %d, %12.5e, 0e0),\n"
	  "          n = 1, Method = Meth;\n",
	  elem[loc].Name, elem[loc].L, Quad, get_bn(Fnum, 1, Quad),
	  n, get_bn(Fnum, 1, n));
}


void prt_sext(FILE *outf, const int loc, const int n)
{
  const int Fnum = elem[loc].Fnum;

  if (n == Sext)
    fprintf(outf,
	    "%-8s: multipole, l = %7.5f,\n"
	    "          hom = (%d, %12.5e, 0e0),\n"
	    "          n = 1, Method = Meth;\n",
	    elem[loc].Name, elem[loc].L, Sext, get_bn(Fnum, 1, Sext));
  else
    fprintf(outf,
	    "%-8s: multipole, l = %7.5f,\n"
	    "          hom = (%d, %12.5e, 0e0,"
	    " %d, %12.5e, 0e0),\n"
	    "          n = 1, Method = Meth;\n",
	    elem[loc].Name, elem[loc].L, Sext, get_bn(Fnum, 1, Sext),
	    n, get_bn(Fnum, 1, n));
}


void prt_bn(const param_type &bns)
{
  long int loc;
  int      k;
  FILE     *outf;

  const std::string file_name = "b4.out";

  outf = file_write(file_name.c_str());

  fprintf(outf, "\n");
  for (k = 0; k < bns.n_prm; k++) {
    loc = get_loc(bns.Fnum[k], 1) - 1;
    if (elem[loc].mpole->n_design == Dip)
      prt_bend(outf, loc, bns.n[k]);
    else if (elem[loc].mpole->n_design == Quad)
      prt_quad(outf, loc, bns.n[k]);
    else if (elem[loc].mpole->n_design == Sext)
      prt_sext(outf, loc, bns.n[k]);
  }

  fclose(outf);
}


void correct(param_type &bns, const std::vector<Lie_term> &k_ijklm,
	     const double svd_cut, const double scl)
{
  int    k;
  double **A, **U, **V, *w, *b, *dbnL, *bnL_max, *bnL;

  const int
    m = k_ijklm.size(),
    n = bns.n_prm;

  printf("\nsvd:\n  m = %d n = %d\n", m, n);

  A = dmatrix(1, m, 1, n); U = dmatrix(1, m, 1, n); V = dmatrix(1, n, 1, n);
  w = dvector(1, n); b = dvector(1, m); dbnL = dvector(1, n);
  bnL_max = dvector(1, n); bnL = dvector(1, n);

  get_A(m, n, k_ijklm, A, w, U, V, b);

#if 0
  dmcopy(A, m, n, U);
  dsvdcmp(U, m, n, w, V);
  get_sing_val(n, w, svd_cut);

  dsvbksb(U, w, V, m, n, b, dbnL);
#else
  for (k = 0; k < n; k++) {
    bnL_max[k+1] = bns.bnL_max[k]*bns.L[k]/bns.bnL_scl[k];
    bnL[k+1] = get_bnL(bns.Fnum[k], 1, bns.n[k])/bns.bnL_scl[k];
  }

  SVD_lim(m, n, A, b, bnL_max, svd_cut, bnL, dbnL);

#endif

  set_bnL(dbnL, bns);
  bns.print();
  prt_bn(bns);

  free_dmatrix(A, 1, m, 1, n); free_dmatrix(U, 1, m, 1, n);
  free_dmatrix(V, 1, n, 1, n); free_dvector(w, 1, n); free_dvector(b, 1, m);
  free_dvector(dbnL, 1, n); free_dvector(bnL_max, 1, n);
  free_dvector(bnL, 1, n);
}


void get_h2_ijklm(const tps &h, const double scl, const int i, const int j,
		  const int k, const int l, const int m,
		  std::vector<Lie_term> &h2_ijklm)
{
  Lie_term           h2;
  std::ostringstream str;

  h2.cst_scl = scl;
  str << "k_" << i << j << k << l << m;
  h2.label = str.str();
  h2.cst = h2.cst_scl*h_ijklm(h, i, j, k, l, m);
  h2_ijklm.push_back(h2);
}


void get_constr(const ss_vect<tps> &Id_scl, std::vector<Lie_term> &k_ijklm,
		const bool tune_fp)
{
  double       nu[3], ksi[3];
  tps          g_re, g_im, K, K_re, K_im;
  ss_vect<tps> nus;
  Lie_term     h;

  const double
    scl_h     = 1e-1,
    scl_ksi[] = {0e0, 1e1, 1e0, 0e0},
    scl_a[]   = {1e0, 0e0};

  danot_(no_tps-1);
  get_Map();
  danot_(no_tps);
  K = MapNorm(Map, g, A1, A0, Map_res, 1);
  nus = dHdJ(K);
  get_nu_ksi(nus, nu, ksi);

  CtoR(g*Id_scl, g_re, g_im);
  CtoR(K*Id_scl, K_re, K_im);

  if (false)
    std::cout << std::scientific << std::setprecision(3)
	      << "\nK:\n" << K_re << "\ng:\n" << g_im;

  if (false)
    printf("\n  nu  = [%6.3f, %6.3f]\n"
	   "  ksi = [%6.3f, %6.3f]\n",
	   nu[X_], nu[Y_], ksi[X_], ksi[Y_]);

  k_ijklm.clear();

  get_h2_ijklm(K_re, scl_ksi[1], 1, 1, 0, 0, 1, k_ijklm);
  get_h2_ijklm(K_re, scl_ksi[1], 0, 0, 1, 1, 1, k_ijklm);

  get_h2_ijklm(g_im, scl_h, 1, 0, 0, 0, 2, k_ijklm);
  get_h2_ijklm(g_im, scl_h, 2, 0, 0, 0, 1, k_ijklm);
  get_h2_ijklm(g_im, scl_h, 0, 0, 2, 0, 1, k_ijklm);

  get_h2_ijklm(g_im, scl_h, 1, 0, 1, 1, 0, k_ijklm);
  get_h2_ijklm(g_im, scl_h, 2, 1, 0, 0, 0, k_ijklm);
  get_h2_ijklm(g_im, scl_h, 3, 0, 0, 0, 0, k_ijklm);
  get_h2_ijklm(g_im, scl_h, 1, 0, 0, 2, 0, k_ijklm);
  get_h2_ijklm(g_im, scl_h, 1, 0, 2, 0, 0, k_ijklm);

  if (tune_fp) {
    get_h2_ijklm(K_re, scl_a[0], 2, 2, 0, 0, 0, k_ijklm);
    get_h2_ijklm(K_re, scl_a[0], 1, 1, 1, 1, 0, k_ijklm);
    get_h2_ijklm(K_re, scl_a[0], 0, 0, 2, 2, 0, k_ijklm);

    get_h2_ijklm(K_re, scl_a[1], 3, 3, 0, 0, 0, k_ijklm);
    get_h2_ijklm(K_re, scl_a[1], 2, 2, 1, 1, 0, k_ijklm);
    get_h2_ijklm(K_re, scl_a[1], 1, 1, 2, 2, 0, k_ijklm);
    get_h2_ijklm(K_re, scl_a[1], 0, 0, 3, 3, 0, k_ijklm);

    get_h2_ijklm(K_re, scl_ksi[2], 1, 1, 0, 0, 2, k_ijklm);
    get_h2_ijklm(K_re, scl_ksi[2], 0, 0, 1, 1, 2, k_ijklm);

    get_h2_ijklm(K_re, scl_ksi[3], 1, 1, 0, 0, 3, k_ijklm);
    get_h2_ijklm(K_re, scl_ksi[3], 0, 0, 1, 1, 3, k_ijklm);
  }
}


void get_h1_ijklm(const tps &h1, const std::string &name, const int n,
		  const int i, const int j, const int k, const int l,
		  const int m, const double bnL_scl, Lie_term &k_ijklm)
{
  k_ijklm.Jacobian.push_back
    (bnL_scl*k_ijklm.cst_scl*h_ijklm_p(h1, i, j, k, l, m, 7));
}


void get_K_ijklm(const std::string &name, const int n,
		 const ss_vect<tps> &Id_scl, const double bnL_scl,
		 std::vector<Lie_term> &k_ijklm, const bool tune_fp)
{
  int k;
  tps g_re, g_im, K, K_re, K_im;

  const int Fnum = get_Fnum(name.c_str());

  set_bn_par(Fnum, n, 7);

  danot_(no_tps-1);
  get_Map();
  danot_(no_tps);
  K = MapNorm(Map, g, A1, A0, Map_res, 1);
  CtoR(g*Id_scl, g_re, g_im);
  CtoR(K*Id_scl, K_re, K_im);

  if (false)
    std::cout << std::scientific << std::setprecision(3)
	      << "\n" << std::setw(8) << name << std::setw(4) << Fnum << ":\n"
	      << std::scientific << std::setprecision(3) << K_re << K_im;

  clr_bn_par(Fnum, n);

  k = 0;
  get_h1_ijklm(K_re, name, n, 1, 1, 0, 0, 1, bnL_scl, k_ijklm[k++]);
  get_h1_ijklm(K_re, name, n, 0, 0, 1, 1, 1, bnL_scl, k_ijklm[k++]);

  get_h1_ijklm(g_im, name, n, 1, 0, 0, 0, 2, bnL_scl, k_ijklm[k++]);
  get_h1_ijklm(g_im, name, n, 2, 0, 0, 0, 1, bnL_scl, k_ijklm[k++]);
  get_h1_ijklm(g_im, name, n, 0, 0, 2, 0, 1, bnL_scl, k_ijklm[k++]);

  get_h1_ijklm(g_im, name, n, 1, 0, 1, 1, 0, bnL_scl, k_ijklm[k++]);
  get_h1_ijklm(g_im, name, n, 2, 1, 0, 0, 0, bnL_scl, k_ijklm[k++]);
  get_h1_ijklm(g_im, name, n, 3, 0, 0, 0, 0, bnL_scl, k_ijklm[k++]);
  get_h1_ijklm(g_im, name, n, 1, 0, 0, 2, 0, bnL_scl, k_ijklm[k++]);
  get_h1_ijklm(g_im, name, n, 1, 0, 2, 0, 0, bnL_scl, k_ijklm[k++]);

  if (tune_fp) {
    get_h1_ijklm(K_re, name, n, 2, 2, 0, 0, 0, bnL_scl, k_ijklm[k++]);
    get_h1_ijklm(K_re, name, n, 1, 1, 1, 1, 0, bnL_scl, k_ijklm[k++]);
    get_h1_ijklm(K_re, name, n, 0, 0, 2, 2, 0, bnL_scl, k_ijklm[k++]);

    get_h1_ijklm(K_re, name, n, 3, 3, 0, 0, 0, bnL_scl, k_ijklm[k++]);
    get_h1_ijklm(K_re, name, n, 2, 2, 1, 1, 0, bnL_scl, k_ijklm[k++]);
    get_h1_ijklm(K_re, name, n, 1, 1, 2, 2, 0, bnL_scl, k_ijklm[k++]);
    get_h1_ijklm(K_re, name, n, 0, 0, 3, 3, 0, bnL_scl, k_ijklm[k++]);

    get_h1_ijklm(K_re, name, n, 1, 1, 0, 0, 2, bnL_scl, k_ijklm[k++]);
    get_h1_ijklm(K_re, name, n, 0, 0, 1, 1, 2, bnL_scl, k_ijklm[k++]);

    get_h1_ijklm(K_re, name, n, 1, 1, 0, 0, 3, bnL_scl, k_ijklm[k++]);
    get_h1_ijklm(K_re, name, n, 0, 0, 1, 1, 3, bnL_scl, k_ijklm[k++]);
  }
}


void get_Jacobian(const ss_vect<tps> &Id_scl, const param_type &bns,
		std::vector<Lie_term> &k_ijklm, const bool tune_fp)
{
  int k;

  for (k = 0; k < bns.n_prm; k++)
    get_K_ijklm(bns.name[k], bns.n[k], Id_scl, bns.bnL_scl[k]/bns.L[k], k_ijklm,
		tune_fp);
}


void analyze(const ss_vect<tps> &Id_scl, const param_type &bns,
	     std::vector<Lie_term> &k_ijklm, const bool tune_fp)
{
  get_constr(Id_scl, k_ijklm, tune_fp);

  get_Jacobian(Id_scl, bns, k_ijklm, tune_fp);
  
  prt_K_ijklm(bns, k_ijklm, tune_fp);
}

void analyze_2(void)
{
  tps          h, h_re, h_im, dh;
  ss_vect<tps> Id;

  danot_(no_tps-1);
  get_Map();
  danot_(no_tps);
  K = MapNorm(Map, g, A1, A0, Map_res, 1);
  h = get_H(g);

  dh = (h-h*Inv(A0*A1)*Map*A0*A1);
  // Truncate to 4th order.
  danot_(3);
  h = 1e0*h;
  dh = 1e0*dh;
  // Remove momentum dependence.
  Id.identity();
  Id[delta_] = 0e0;
  h = h*Id;
  dh = dh*Id;
  
  CtoR(h, h_re, h_im);

  std::cout << std::scientific << std::setprecision(3)
	    << "\ndh:\n" << dh << "\nh:\n" << h
	    << "\nh_re:\n" << h_re << "\nh_im:\n" << h_im;
}


void get_bns(param_type &bns)
{
  const double
    bnL_scl[] = {0e0, 0e0, 0e0,  1e0,  1e2,  1e4},
    bnL_min[] = {0e0, 0e0, 0e0, -3e2, -1e4, -1e5},
    bnL_max[] = {0e0, 0e0, 0e0,  3e2,  1e4,  1e5};

  bns.add_prm("sf_h", Sext, bnL_min[Sext], bnL_max[Sext], bnL_scl[Sext]);
  bns.add_prm("sd_h", Sext, bnL_min[Sext], bnL_max[Sext], bnL_scl[Sext]);

  if (!false) {
    bns.add_prm("sf2_h", Sext, bnL_min[Sext], bnL_max[Sext], bnL_scl[Sext]);
    // bns.add_prm("sd2_h", Sext, bnL_min[Sext], bnL_max[Sext], bnL_scl[Sext]);
  }

  if (!false) {
    bns.add_prm("sf_h", Oct, bnL_min[Oct], bnL_max[Oct], bnL_scl[Oct]);
    bns.add_prm("sd_h", Oct, bnL_min[Oct], bnL_max[Oct], bnL_scl[Oct]);
  }
}


int main(int argc, char *argv[])
{
  int                   k;
  double                nu[3], ksi[3];
  ss_vect<tps>          Id_scl;
  param_type            bns;
  std::vector<Lie_term> k_ijklm;

  rad_on    = false; H_exact        = false; totpath_on   = false;
  cavity_on = false; quad_fringe_on = false; emittance_on = false;
  IBS_on    = false;

  const bool tune_fp = true;

  rd_mfile(argv[1], elem); rd_mfile(argv[1], elem_tps);

  // initialize the symplectic integrator after energy has been defined
  ini_si();

  // disable from TPSALib- and LieLib log messages
  idprset(-1);

  daeps_(1e-30);

  Id_scl.identity();
  for (k = 0; k < 4; k++)
    Id_scl[k] *= sqrt(twoJ[k/2]);
  Id_scl[delta_] *= delta_max;

  if (!false) chk_lat(Map, Map_res, nu, ksi);

  if (!false) {
    get_bns(bns);
    bns.ini_prm();
    bns.print();

    printf("\n");
    for (k = 1; k <= 30; k++) {
      printf("\nk = %d:", k);
      analyze(Id_scl, bns, k_ijklm, tune_fp);
      correct(bns, k_ijklm, 1e-10, 0.3);
      prt_mfile("flat_file.fit");
    }
    analyze(Id_scl, bns, k_ijklm, tune_fp);
  }

  if (false) analyze_2();
}
