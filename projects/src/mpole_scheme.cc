#define NO 5

#include "thor_lib.h"

int
  no_tps   = NO,
  ndpt_tps = 5;


extern double       b2_max;
extern tps          K, g;
extern ss_vect<tps> Map, A0, A1, Map_res;


typedef struct {
  std::string              label;
  double                   cst;
  std::vector<std::string> prms;
  std::vector<double>      Jacobian;
} Lie_term;


const double
  beta_inj[] = {2.8, 2.6},
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


void get_K_ijklm(const tps &K, std::vector<Lie_term> &k_ijklm)
{
  Lie_term h;

  h.label = "k_22000";
  h.cst = h_ijklm(K, 2, 2, 0, 0, 0);
  k_ijklm.push_back(h);
  h.label = "k_11110";
  h.cst = h_ijklm(K, 1, 1, 1, 1, 0);
  k_ijklm.push_back(h);
  h.label = "k_00220";
  h.cst = h_ijklm(K, 0, 0, 2, 2, 0);
  k_ijklm.push_back(h);

  h.label = "k_11002";
  h.cst = h_ijklm(K, 1, 1, 0, 0, 2);
  k_ijklm.push_back(h);
  h.label = "k_00112";
  h.cst = h_ijklm(K, 0, 0, 1, 1, 2);
  k_ijklm.push_back(h);
}


void get_K_ijklm_p(const std::string &name, const int n,
		   const ss_vect<tps> &Id_scl, std::vector<Lie_term> &k_ijklm)
{
  tps K, K_re, K_im;

  const int Fnum = get_Fnum(name.c_str());

  set_bn_par(Fnum, n, 7);

  danot_(no_tps-1);
  get_Map();
  danot_(no_tps);
  K = MapNorm(Map, g, A1, A0, Map_res, 1);
  CtoR(K*Id_scl, K_re, K_im);

  if (false)
    std::cout << std::scientific << std::setprecision(3)
	      << "\n" << std::setw(8) << name << std::setw(4) << Fnum << ":\n"
	      << std::scientific << std::setprecision(3) << K_re;

  clr_bn_par(Fnum, Oct);

  k_ijklm[0].prms.push_back(name);
  k_ijklm[0].Jacobian.push_back(h_ijklm_p(K_re, 2, 2, 0, 0, 0, 7));
  k_ijklm[1].prms.push_back(name);
  k_ijklm[1].Jacobian.push_back(h_ijklm_p(K_re, 1, 1, 1, 1, 0, 7));
  k_ijklm[2].prms.push_back(name);
  k_ijklm[2].Jacobian.push_back(h_ijklm_p(K_re, 0, 0, 2, 2, 0, 7));

  k_ijklm[3].prms.push_back(name);
  k_ijklm[3].Jacobian.push_back(h_ijklm_p(K_re, 1, 1, 0, 0, 2, 7));
  k_ijklm[4].prms.push_back(name);
  k_ijklm[4].Jacobian.push_back(h_ijklm_p(K_re, 0, 0, 1, 1, 2, 7));
}


void prt_Lie_term(const Lie_term &k_ijklm)
{
  int k;

  printf(" %s %10.3e", k_ijklm.label.c_str(), k_ijklm.cst);
  for (k = 0; k < (int)k_ijklm.Jacobian.size(); k++)
    printf(" %10.3e", k_ijklm.Jacobian[k]);
  printf("\n");
}


void prt_K_ijklm(const std::vector<Lie_term> &k_ijklm)
{
  int k;

  printf("\n             cst.");
  for (k = 0; k < k_ijklm[0].prms.size(); k++)
    printf(" %8s  ", k_ijklm[0].prms[k].c_str());
  printf("\nAnharmonic terms:\n");
  for (k = 0; k < 3; k++)
    prt_Lie_term(k_ijklm[k]);
  printf("2nd order chromaticity:\n");
  for (k = 3; k < 5; k++)
    prt_Lie_term(k_ijklm[k]);
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


void set_bn(const int n, const double *bn, const double scl,
	    const std::vector<Lie_term> &k_ijklm)
{
  // Integrated strengths: b_n*L.
  int k, Fnum;

  printf("\nb_n:\n  scaled by %3.1f\n ", scl);
  for (k = 0; k < n; k++) {
    Fnum = get_Fnum(k_ijklm[0].prms[k].c_str());
    set_bn(Fnum, Oct, scl*bn[k+1]);
    if (!true)
      printf(" %10.3e", get_bn(Fnum, 1, Oct));
    else
      printf(" %10.3e", get_bnL(Fnum, 1, Oct));
  }
  printf("\n");
}


void prt_bend(FILE *outf, const int loc)
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
	  Oct, get_bn(Fnum, 1, Oct));
}


void prt_quad(FILE *outf, const int loc)
{
  const int Fnum = elem[loc].Fnum;

  fprintf(outf,
	  "%-8s: multipole, l = %7.5f,\n"
	  "          hom = (%d, %12.5e, 0e0,"
	  " %d, %12.5e, 0e0),\n"
	  "          n = 1, Method = Meth;\n",
	  elem[loc].Name, elem[loc].L, Quad, get_bn(Fnum, 1, Quad),
	  Oct, get_bn(Fnum, 1, Oct));
}


void prt_sext(FILE *outf, const int loc)
{
  const int Fnum = elem[loc].Fnum;

  fprintf(outf,
	  "%-8s: multipole, l = %7.5f,\n"
	  "          hom = (%d, %12.5e, 0e0,"
	  " %d, %12.5e, 0e0),\n"
	  "          n = 1, Method = Meth;\n",
	  elem[loc].Name, elem[loc].L, Sext, get_bn(Fnum, 1, Sext),
	  Oct, get_bn(Fnum, 1, Oct));
}


void prt_bn(const std::vector<Lie_term> &k_ijklm)
{
  long int loc;
  int      k, Fnum;
  FILE     *outf;

  const std::string file_name = "b4.out";

  outf = file_write(file_name.c_str());

  fprintf(outf, "\n");
  for (k = 0; k < k_ijklm[0].prms.size(); k++) {
    Fnum = get_Fnum(k_ijklm[0].prms[k].c_str());
    loc = get_loc(Fnum, 1) - 1;
    if (elem[loc].mpole->n_design == Dip)
      prt_bend(outf, loc);
    else if (elem[loc].mpole->n_design == Quad)
      prt_quad(outf, loc);
    else if (elem[loc].mpole->n_design == Sext)
      prt_sext(outf, loc);
  }

  fclose(outf);
}


void correct(const std::vector<Lie_term> &k_ijklm, const double svd_cut,
	     const double scl)
{
  double **A, **U, **V, *w, *b, *bn;

  const int
    m = k_ijklm.size(),
    n = k_ijklm[0].Jacobian.size();

  printf("\nsvd:\n  m = %d n = %d\n", m, n);

  A = dmatrix(1, m, 1, n); U = dmatrix(1, m, 1, n); V = dmatrix(1, n, 1, n);
  w = dvector(1, n); b = dvector(1, m); bn = dvector(1, n);

  get_A(m, n, k_ijklm, A, w, U, V, b);

  dmcopy(A, m, n, U);
  dsvdcmp(U, m, n, w, V);
  get_sing_val(n, w, svd_cut);

  dsvbksb(U, w, V, m, n, b, bn);

  set_bn(n, bn, scl, k_ijklm);
  prt_bn(k_ijklm);

  free_dmatrix(A, 1, m, 1, n); free_dmatrix(U, 1, m, 1, n);
  free_dmatrix(V, 1, n, 1, n); free_dvector(w, 1, n); free_dvector(b, 1, m);
  free_dvector(bn, 1, n);
}


void analyze(const ss_vect<tps> &Id_scl, std::vector<Lie_term> &k_ijklm)
{
  double       nu[3], ksi[3];
  tps          g_re, g_im, K, K_re, K_im;
  ss_vect<tps> nus;

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

  get_K_ijklm(K_re, k_ijklm);

  if (!false) {
    get_K_ijklm_p("s1",  Oct, Id_scl, k_ijklm);
    get_K_ijklm_p("s2",  Oct, Id_scl, k_ijklm);
  }
  if (!false) {
    if (false) get_K_ijklm_p("q1",  Oct, Id_scl, k_ijklm);
    get_K_ijklm_p("uq1", Oct, Id_scl, k_ijklm);
    get_K_ijklm_p("uq2", Oct, Id_scl, k_ijklm);
    get_K_ijklm_p("uq3", Oct, Id_scl, k_ijklm);
  }
  if (false) {
    get_K_ijklm_p("b1",  Oct, Id_scl, k_ijklm);
    get_K_ijklm_p("b2",  Oct, Id_scl, k_ijklm);
    get_K_ijklm_p("mb1", Oct, Id_scl, k_ijklm);
    get_K_ijklm_p("mb2", Oct, Id_scl, k_ijklm);
  }

  prt_K_ijklm(k_ijklm);
}


int main(int argc, char *argv[])
{
  int                   k;
  double                nu[3], ksi[3];
  ss_vect<tps>          Id_scl;
  std::vector<Lie_term> k_ijklm;

  rad_on    = false; H_exact        = false; totpath_on   = false;
  cavity_on = false; quad_fringe_on = false; emittance_on = false;
  IBS_on    = false;

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

  if (!false) {
    chk_lat(Map, Map_res, nu, ksi);
  }

  analyze(Id_scl, k_ijklm);
  correct(k_ijklm, 1e-12, 0.9);
  prt_mfile("flat_file.fit");
  analyze(Id_scl, k_ijklm);
}
