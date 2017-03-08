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
double       twoJ[2];
ss_vect<tps> Id_scl;

const int n_prm_max = 8;

const double scl_dnu[] = {1e4, 1e-2};


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

  printf("\nInitial bn; (incl. scaling) (%d):\n", n_prm);
  for (i = 1; i <= n_prm; i++) {
    bn_lim[i] = bn_max[i-1];
    if (n[i-1] > 0)
      // Multipole.
      bn[i] = get_bn(Fnum[i-1], 1, n[i-1])/bn_scl[i-1];
    else if (n[i-1] == -1)
      // Drift.
      bn[i] = get_L(Fnum[i-1], 1)/bn_scl[i-1];
    else if (n[i-1] == -2)
      // Location.
      bn[i] = get_bn_s(-Fnum[i-1], 1, n[i-1])/bn_scl[i-1];
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
    dbn[i] *= bn_scl[i-1]*step;
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
    bn[i] /= bn_scl[i-1];
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
      set_bn(Fnum[i-1], n[i-1], bn_scl[i-1]*bn[i]);
    else if (n[i-1] == -1)
      set_L(Fnum[i-1], bn_scl[i-1]*bn[i]);
    else if (n[i-1] == -2)
      set_bn_s(-Fnum[i-1], n[i-1], bn_scl[i-1]*bn[i]);
    printf(" %13.5e", bn_scl[i-1]*bn[i]);
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

	
void min_dnu_prt(const param_type &bn_prms)
{
  FILE *outf;

  std::string file_name = "dnu.out";

  outf = file_write(file_name.c_str());
  fprintf(outf, "o1: multipole, l = 0.0, N = Nsext, Method = Meth,"
	  "\n    HOM = (4, %12.5e, 0.0);\n", bn_prms.bn_scl[0]*bn_prms.bn[1]);
  fprintf(outf, "o2: multipole, l = 0.0, N = Nsext, Method = Meth,"
	  "\n    HOM = (4, %12.5e, 0.0);\n", bn_prms.bn_scl[1]*bn_prms.bn[2]);
  fprintf(outf, "o3: multipole, l = 0.0, N = Nsext, Method = Meth,"
	  "\n    HOM = (4, %12.5e, 0.0);\n", bn_prms.bn_scl[2]*bn_prms.bn[3]);
  fprintf(outf, "o4: multipole, l = 0.0, N = Nsext, Method = Meth,"
	  "\n    HOM = (4, %12.5e, 0.0);\n", bn_prms.bn_scl[3]*bn_prms.bn[4]);
  fclose(outf);
}


void min_dnu_prt2(const param_type &bn_prms)
{
  FILE *outf;

  std::string file_name = "dnu.out";

  outf = file_write(file_name.c_str());
  fprintf(outf, "o1: multipole, l = 0.0, N = Nsext, Method = Meth,"
	  "\n    HOM = (4, %12.5e, 0.0, 6, %12.5e, 0.0);\n",
	  bn_prms.bn_scl[0]*bn_prms.bn[1], bn_prms.bn_scl[4]*bn_prms.bn[5]);
  fprintf(outf, "o2: multipole, l = 0.0, N = Nsext, Method = Meth,"
	  "\n    HOM = (4, %12.5e, 0.0, 6, %12.5e, 0.0);\n",
	  bn_prms.bn_scl[1]*bn_prms.bn[2], bn_prms.bn_scl[5]*bn_prms.bn[6]);
  fprintf(outf, "o3: multipole, l = 0.0, N = Nsext, Method = Meth,"
	  "\n    HOM = (4, %12.5e, 0.0, 6, %12.5e, 0.0);\n",
	  bn_prms.bn_scl[2]*bn_prms.bn[3], bn_prms.bn_scl[6]*bn_prms.bn[7]);
  fprintf(outf, "o4: multipole, l = 0.0, N = Nsext, Method = Meth,"
	  "\n    HOM = (4, %12.5e, 0.0, 6, %12.5e, 0.0);\n",
	  bn_prms.bn_scl[3]*bn_prms.bn[4], bn_prms.bn_scl[7]*bn_prms.bn[8]);
  fclose(outf);
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


double get_a(const double scale, const tps &t,
	     const int i, const int j, const int k, const int l, const int m)
{
  return scale*(h_ijklm_p(t, i, j, k, l, m, 7));
}


double get_b(const double scale, const tps &t,
	     const int i, const int j, const int k, const int l, const int m)
{
  return scale*(h_ijklm(t, i, j, k, l, m));
}


double min_dnu_f(double *b4s)
{
  int                 i;
  static double       chi2_ref = 1e30;
  double              chi2;
  tps                 K_re, K_im;
  ss_vect<tps>        nus;
  std::vector<double> b;

  const bool prt = false;

  n_powell++;

  // Do not change parameters.
  if (prt) printf("min_dnu_f (incl. scaling):\n");
  for (i = 1; i <= bn_prms.n_prm; i++) {
    set_bn(bn_prms.Fnum[i-1], bn_prms.n[i-1], bn_prms.bn_scl[i-1]*b4s[i]);
    if (prt) printf(" %13.5e", bn_prms.bn_scl[i-1]*b4s[i]);
  }
  if (prt) printf("\n");

  danot_(NO-1);
  get_Map();
  danot_(NO);
  K = MapNorm(Map, g, A1, A0, Map_res, 1);
  CtoR(K, K_re, K_im); nus = dHdJ(K);
  K_re = K_re*Id_scl;
  // nus_scl[3] = nus[3]*Id_scl; nus_scl[4] = nus[4]*Id_scl;

  b.push_back(get_b(scl_dnu[0], K_re, 2, 2, 0, 0, 0));
  b.push_back(get_b(scl_dnu[0], K_re, 0, 0, 2, 2, 0));
  b.push_back(get_b(scl_dnu[0], K_re, 1, 1, 1, 1, 0));

  b.push_back(get_b(scl_dnu[0], K_re, 2, 2, 1, 1, 0));
  b.push_back(get_b(scl_dnu[0], K_re, 1, 1, 2, 2, 0));

  // b.push_back(get_b(scl_dnu[0], K_re, 4, 4, 0, 0, 0));
  // b.push_back(get_b(scl_dnu[0], K_re, 0, 0, 4, 4, 0));
  // b.push_back(get_b(scl_dnu[0], K_re, 2, 2, 2, 2, 0));
  // b.push_back(get_b(scl_dnu[0], K_re, 3, 3, 1, 1, 0));
  // b.push_back(get_b(scl_dnu[0], K_re, 1, 1, 3, 3, 0));

  b.push_back(scl_dnu[1]*(get_b(1e0, K_re, 3, 3, 0, 0, 0)
			  +get_b(1e0, K_re, 2, 2, 0, 0, 0)/(3e0*twoJ[X_])));
  b.push_back(scl_dnu[1]*(get_b(1e0, K_re, 0, 0, 3, 3, 0)
			  +get_b(1e0, K_re, 0, 0, 2, 2, 0)/(3e0*twoJ[Y_])));

  b.push_back(scl_dnu[1]*(get_b(1e0, K_re, 1, 1, 2, 2, 0)
  			  +get_b(1e0, K_re, 1, 1, 1, 1, 0)/(2e0*twoJ[Y_])));
  b.push_back(scl_dnu[1]*(get_b(1e0, K_re, 2, 2, 1, 1, 0)
  			  +get_b(1e0, K_re, 1, 1, 1, 1, 0)/(2e0*twoJ[X_])));

  chi2 = 0e0;
  for (i = 0; i < (int)b.size(); i++)
    chi2 += sqr(b[i]);

  if (chi2 < chi2_ref) {
    printf("\n%3d chi2: %12.5e -> %12.5e\n", n_powell, chi2_ref, chi2);
    printf("b & bn:\n");
    for (i = 0; i < (int)b.size(); i++)
      printf("%11.3e", b[i]);
    printf("\n");
    for (i = 1; i <= bn_prms.n_prm; i++) 
      printf("%11.3e", b4s[i]);
    printf("\n");
  }

  chi2_ref = min(chi2, chi2_ref);

  return chi2;
}


void min_dnu_grad(double &chi2, double &db4_max, double *g_, double *h_,
		  const bool cg_meth)
{
  int           n_b4, i, j, m;
  double        chi2_ref, **A, *b, *bn_ref;
  tps           K_re, K_im;
  ss_vect<tps>  nus;
  std::ofstream outf;

  const int m_max = 20;

  n_b4 = bn_prms.n_prm;

  bn_ref = dvector(1, n_b4);
  b = dvector(1, m_max); A = dmatrix(1, m_max, 1, n_b4);

  printf("\n");
  for (i = 1; i <= n_b4; i++) {
    bn_prms.set_prm_dep(i-1);

    danot_(NO-1);
    get_Map();
    danot_(NO);
    K = MapNorm(Map, g, A1, A0, Map_res, 1);
    CtoR(K, K_re, K_im); nus = dHdJ(K);
    K_re = K_re*Id_scl;
    // nus[3] = nus[3]*Id_scl; nus[4] = nus[4]*Id_scl;

    // file_wr(outf, "nus.out");
    // outf << nus[3]*Id_scl << nus[4]*Id_scl;
    // outf.close();

    m = 0;
    A[++m][i] = get_a(scl_dnu[0], K_re, 2, 2, 0, 0, 0);
    A[++m][i] = get_a(scl_dnu[0], K_re, 0, 0, 2, 2, 0);

    A[++m][i] = get_a(scl_dnu[0], K_re, 1, 1, 1, 1, 0);
    A[++m][i] = get_a(scl_dnu[0], K_re, 2, 2, 1, 1, 0);
    A[++m][i] = get_a(scl_dnu[0], K_re, 1, 1, 2, 2, 0);

    // A[++m][i] = get_a(scl_dnu[0], K_re, 4, 4, 0, 0, 0);
    // A[++m][i] = get_a(scl_dnu[0], K_re, 0, 0, 4, 4, 0);
    // A[++m][i] = get_a(scl_dnu[0], K_re, 2, 2, 2, 2, 0);
    // A[++m][i] = get_a(scl_dnu[0], K_re, 3, 3, 1, 1, 0);
    // A[++m][i] = get_a(scl_dnu[0], K_re, 1, 1, 3, 3, 0);

    A[++m][i] = scl_dnu[1]*(get_a(1e0, K_re, 3, 3, 0, 0, 0)
			    +get_a(1e0/(3e0*twoJ[X_]), K_re, 2, 2, 0, 0, 0));
    A[++m][i] = scl_dnu[1]*(get_a(1e0, K_re, 0, 0, 3, 3, 0)
			    +get_a(1e0/(3e0*twoJ[Y_]), K_re, 0, 0, 2, 2, 0));

    A[++m][i] = scl_dnu[1]*(get_a(1e0, K_re, 1, 1, 2, 2, 0)
			    +get_a(1e0/(2e0*twoJ[Y_]), K_re, 1, 1, 1, 1, 0));
    A[++m][i] = scl_dnu[1]*(get_a(1e0, K_re, 2, 2, 1, 1, 0)
			    +get_a(1e0/(2e0*twoJ[X_]), K_re, 1, 1, 1, 1, 0));

    for (j = 1; j <= m; j++)
      A[j][i] *= bn_prms.bn_scl[i-1];

    bn_prms.clr_prm_dep(i-1);
  }

  m = 0;
  b[++m] = -get_b(scl_dnu[0], K_re, 2, 2, 0, 0, 0);
  b[++m] = -get_b(scl_dnu[0], K_re, 0, 0, 2, 2, 0);
  b[++m] = -get_b(scl_dnu[0], K_re, 1, 1, 1, 1, 0);

  b[++m] = -get_b(scl_dnu[0], K_re, 2, 2, 1, 1, 0);
  b[++m] = -get_b(scl_dnu[0], K_re, 1, 1, 2, 2, 0);

  // b[++m] = -get_b(scl_dnu[0], K_re, 4, 4, 0, 0, 0);
  // b[++m] = -get_b(scl_dnu[0], K_re, 0, 0, 4, 4, 0);
  // b[++m] = -get_b(scl_dnu[0], K_re, 2, 2, 2, 2, 0);
  // b[++m] = -get_b(scl_dnu[0], K_re, 3, 3, 1, 1, 0);
  // b[++m] = -get_b(scl_dnu[0], K_re, 1, 1, 3, 3, 0);

  b[++m] = -scl_dnu[1]*(get_b(1e0, K_re, 3, 3, 0, 0, 0)
			+get_b(1e0/(3e0*twoJ[X_]), K_re, 2, 2, 0, 0, 0));
  b[++m] = -scl_dnu[1]*(get_b(1e0, K_re, 0, 0, 3, 3, 0)
			+get_b(1e0/(3e0*twoJ[Y_]), K_re, 0, 0, 2, 2, 0));

  b[++m] = -scl_dnu[1]*(get_b(1e0, K_re, 1, 1, 2, 2, 0)
  			+get_b(1e0/(2e0*twoJ[Y_]), K_re, 1, 1, 1, 1, 0));
  b[++m] = -scl_dnu[1]*(get_b(1e0, K_re, 2, 2, 1, 1, 0)
  			+get_b(1e0/(2e0*twoJ[X_]), K_re, 1, 1, 1, 1, 0));

  chi2_ref = chi2;
  
  chi2 = 0e0;
  for (j = 1; j <= m; j++)
    chi2 += sqr(b[j]);

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
    bn_prms.bn[i] =
      get_bn(bn_prms.Fnum[i-1], 1, bn_prms.n[i-1])/bn_prms.bn_scl[i-1];

  printf("\nbn & dbn (incl. scaling):\n");
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


int main(int argc, char *argv[])
{
  int    j;
  tps    K_re, K_im;
  
  const double beta[]  = {3.0, 3.0},                     // Center of straight.
    A_max[] = {1.2e-3, 1.2e-3}, delta = 3e-2;

    rad_on    = false; H_exact        = false; totpath_on   = false;
    cavity_on = false; quad_fringe_on = false; emittance_on = false;
    IBS_on    = false;

    rd_mfile(argv[1], elem); rd_mfile(argv[1], elem_tps);
  
    // Initialize the symplectic integrator after energy has been defined.
    ini_si();

    // Disable log messages from TPSALib and LieLib.
    idprset(-1);

    daeps_(1e-30);

    danot_(1);

    get_nu_ksi();

    for (j = 0; j < 2; j++)
      twoJ[j] = sqr(A_max[j])/beta[j];

    Id_scl.identity();
    Id_scl[x_] *= sqrt(twoJ[X_]); Id_scl[px_] *= sqrt(twoJ[X_]);
    Id_scl[y_] *= sqrt(twoJ[Y_]); Id_scl[py_] *= sqrt(twoJ[Y_]);
    Id_scl[delta_] *= delta;

    bn_prms.add_prm("o1", 4, 5e5, 1.0);
    bn_prms.add_prm("o2", 4, 5e5, 1.0);
    bn_prms.add_prm("o3", 4, 5e5, 1.0);
    bn_prms.add_prm("o4", 4, 5e5, 1.0);

    bn_prms.add_prm("o1", 6, 1e10, 1.0);
    bn_prms.add_prm("o2", 6, 1e10, 1.0);
    bn_prms.add_prm("o3", 6, 1e10, 1.0);
    bn_prms.add_prm("o4", 6, 1e10, 1.0);

    bn_prms.bn_tol = 1e-1; bn_prms.svd_cut = 1e-16; bn_prms.step = 0.01;

    no_mpoles(Oct); no_mpoles(Dodec);

    min_dnu(true);

    // danot_(6);
    // get_Map();
    // danot_(7);
    // K = MapNorm(Map, g, A1, A0, Map_res, no_tps);
    // CtoR(K, K_re, K_im);

    // std::cout << K_re*Id_scl;
}
