#define NO 5

#include "thor_lib.h"

int no_tps   = NO,
    ndpt_tps = 5;


extern tps          K, g;
extern ss_vect<tps> Map, A0, A1, Map_res;

long int   rseed0, rseed;
double     normcut_, chi2 = 0e0;

const bool tune_conf = true;

const double scl_h[] = {1e0, 1e0}, scl_dnu[] = {1e10, 1e0, 0e0, 1e-9};

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


#define k 19
#define c 656329L
#define m 100000001

void iniranf(const long i)
{
  rseed0 = i; rseed = i;
}

void newseed(void)
{
  rseed0 = (k*rseed0+c) % m; rseed = (rseed0+54321) % m;
}

double ranf(void)
{
  // Random number [0, 1] with rectangular distribution.
  rseed = (k*rseed+c) % m; return (rseed/1e8);
}

#undef k
#undef c
#undef m  


void setrancut(const double cut)
{

  printf("\n");
  printf("setrancut: cut set to %3.1f\n", cut);

  normcut_ = cut;
}


double normranf(void)
{
  int     i, j;
  double  f, w;

  const int  maxiter = 100, n = 12;

  j = 0;
  do {
    j++;
    w = 0.0;
    for (i = 1; i <= n; i++)
      w += ranf();
    f = w - 6.0;
  } while (fabs(f) > fabs(normcut_) && j <= maxiter);

  if (j > maxiter)
    fprintf(stdout,"*** fatal error in normranf\n");
  return f;
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


void prt_map(void)
{
  tps           h, h_re, h_im, K_re, K_im;
  ss_vect<tps>  nus;
  std::ofstream outf;

  danot_(NO-1);
  get_Map();

  file_wr(outf, "map.out");
  outf << Map[ct_];
  outf.close();
}


void prt_h_K(void)
{
  tps           h, h_re, h_im, K_re, K_im;
  ss_vect<tps>  nus;
  std::ofstream outf;

  danot_(NO-1);
  get_Map();
  danot_(NO);
  K = MapNorm(Map, g, A1, A0, Map_res, no_tps); h = get_h();
  CtoR(K, K_re, K_im); CtoR(g, h_re, h_im); nus =  dHdJ(K);

  file_wr(outf, "h.out");
  outf << h_re*Id_scl << h_im*Id_scl;
  outf.close();
  file_wr(outf, "K.out");
  outf << K_re*Id_scl << K_im*Id_scl;
  outf.close();

  nus[3] = nus[3]*Id_scl; nus[4] = nus[4]*Id_scl;

  daeps_(1e-5);

  file_wr(outf, "nus.out");
  // Make a trivial evaluation for truncation.
  outf << 1e0*nus[3] << 1e0*nus[4];
  outf.close();

  daeps_(1e-30);
}


void prt_system(const int m, const int n_b2, double **A, double *b)
{
  int i, j;

  printf("\n Ax = b:\n");
  for (j = 1; j <= n_b2; j++)
    if (j == 1)
      printf("%11d", j);
    else
      printf("%11d", j);
  printf("\n");
  for (i = 1; i <= m; i++) {
    if (i == 1)
      printf("1st order chromatic\n");
    else if (i == 4)
      printf("1st order geometric\n");
    else if (i == 9)
      printf("2nd order geometric\n");
    else if (i == 17)
      printf("linear chromaticity\n");
    else if (i == 19)
      printf("ampl. dependant tune shift\n");
    else if (i == 22)
      printf("2nd order chromaticity\n");
    else if (i == 24)
      printf("cross terms\n");
    else if (i == 27)
      printf("3rd order chromaticity\n");
    else if (i == 29) {
      if (!tune_conf)
	printf("ampl. dependant tune shift\n");
      else
	printf("tune confinement\n");
    }

    printf("%4d", i);
    for (j = 1; j <= n_b2; j++)
      printf("%11.3e", A[i][j]);
    printf("%11.3e\n", b[i]);
  }
}

	
void prt_dnu(void)
{
  tps          K_re, K_im;
  ss_vect<tps> nus;

  danot_(NO-1);
  get_Map();
  danot_(NO);
  K = MapNorm(Map, g, A1, A0, Map_res, no_tps);
  CtoR(K, K_re, K_im); nus = dHdJ(K);
  nus[3] = nus[3]*Id_scl; nus[4] = nus[4]*Id_scl;

  printf("\ndnu:\n");
  printf(" %8.5f",   h_ijklm(nus[3], 1, 1, 0, 0, 0));
  printf(" %8.5f",   h_ijklm(nus[3], 0, 0, 1, 1, 0));
  printf(" %8.5f",   h_ijklm(nus[3], 2, 2, 0, 0, 0));
  printf(" %8.5f",   h_ijklm(nus[3], 1, 1, 1, 1, 0));
  printf(" %8.5f\n", h_ijklm(nus[3], 0, 0, 2, 2, 0));

  printf(" %8.5f",   h_ijklm(nus[4], 1, 1, 0, 0, 0));
  printf(" %8.5f",   h_ijklm(nus[4], 0, 0, 1, 1, 0));
  printf(" %8.5f",   h_ijklm(nus[4], 2, 2, 0, 0, 0));
  printf(" %8.5f",   h_ijklm(nus[4], 1, 1, 1, 1, 0));
  printf(" %8.5f\n", h_ijklm(nus[4], 0, 0, 2, 2, 0));

  printf("\n %8.5f", h_ijklm(nus[3], 0, 0, 0, 0, 2));
  printf(" %8.5f",   h_ijklm(nus[3], 1, 1, 0, 0, 2));
  printf(" %8.5f\n", h_ijklm(nus[3], 0, 0, 1, 1, 2));

  printf(" %8.5f",   h_ijklm(nus[4], 0, 0, 0, 0, 2));
  printf(" %8.5f",   h_ijklm(nus[4], 1, 1, 0, 0, 2));
  printf(" %8.5f\n", h_ijklm(nus[4], 0, 0, 1, 1, 2));

  printf("\nTune confinement:\n");
  printf(" %11.3e %11.3e\n",
	 h_ijklm(K_re, 3, 3, 0, 0, 0),
	 h_ijklm(K_re/(3e0*twoJ[X_]), 2, 2, 0, 0, 0));
  printf(" %11.3e %11.3e\n",
	 h_ijklm(K_re, 1, 1, 2, 2, 0),
	 h_ijklm(K_re/(2e0*twoJ[Y_]), 1, 1, 1, 1, 0));

  printf("\n %11.3e %11.3e\n",
	 h_ijklm(K_re, 0, 0, 3, 3, 0),
	 h_ijklm(K_re/(3e0*twoJ[Y_]), 0, 0, 2, 2, 0));
  printf(" %11.3e %11.3e\n",
	 h_ijklm(K_re, 2, 2, 1, 1, 0),
	 h_ijklm(K_re/(2e0*twoJ[X_]), 1, 1, 1, 1, 0));
 }


void prt_bn(const param_type &bn_prms)
{
  int  j, k;
  FILE *outf;

  const int n_sxt = 7;

  std::string file_name = "dnu.out";

  outf = file_write(file_name.c_str());

  k = 0;

  fprintf(outf, "sfh:  sextupole, l = 0.05, k = %12.5e, n = 4"
	  ", Method = Meth;\n", bn_prms.bn_scl[k]*bn_prms.bn[k+1]);
  k++;
  // if (bn_prms.n_prm == n_sxt)
  if (true)
    fprintf(outf, "sd:   sextupole, l = 0.10, k = %12.5e, n = 4"
	    ", Method = Meth;\n", bn_prms.bn_scl[k]*bn_prms.bn[k+1]);
  else {
    j = k;
  }
  k++;
  fprintf(outf, "sfmh: sextupole, l = 0.05, k = %12.5e, n = 4"
	  ", Method = Meth;\n", bn_prms.bn_scl[k]*bn_prms.bn[k+1]);
  k++;
  fprintf(outf, "sdm:  sextupole, l = 0.10, k = %12.5e, n = 4"
	  ", Method = Meth;\n", bn_prms.bn_scl[k]*bn_prms.bn[k+1]);

  k++;
  fprintf(outf, "sxxh: sextupole, l = 0.05, k = %12.5e, n = 4"
  	  ", Method = Meth;\n", bn_prms.bn_scl[k]*bn_prms.bn[k+1]);
  k++;
  fprintf(outf, "sxyh: sextupole, l = 0.05, k = %12.5e, n = 4"
  	  ", Method = Meth;\n", bn_prms.bn_scl[k]*bn_prms.bn[k+1]);
  k++;
  fprintf(outf, "syyh: sextupole, l = 0.05, k = %12.5e, n = 4"
  	  ", Method = Meth;\n", bn_prms.bn_scl[k]*bn_prms.bn[k+1]);

  if (bn_prms.n_prm > n_sxt) {
    // k++;
    // fprintf(outf, "ocx:  multipole, l = 0.0, n = 4, Method = Meth,"
    // 	    " HOM = (4, %12.5e, 0.0);\n", bn_prms.bn_scl[k]*bn_prms.bn[k+1]);
    // k++;
    // fprintf(outf, "ocxm: multipole, l = 0.0, n = 4, Method = Meth,"
    // 	    " HOM = (4, %12.5e, 0.0);\n", bn_prms.bn_scl[k]*bn_prms.bn[k+1]);
    // k++;
    // fprintf(outf, "sd:   multipole, l = 0.10, n = 4, Method = Meth,"
    // 	    "\n      HOM = (3, %12.5e, 0.0, 4, %12.5e, 0.0);\n",
    // 	    bn_prms.bn_scl[j]*bn_prms.bn[j+1],
    // 	    bn_prms.bn_scl[k]*bn_prms.bn[k+1]);
    // k++;
    // fprintf(outf, "ocxm: multipole, l = 0.0, n = 4, Method = Meth,"
    // 	    " HOM = (4, %12.5e, 0.0);\n", bn_prms.bn_scl[k]*bn_prms.bn[k+1]);

    k++;
    fprintf(outf, "\noxx:  multipole, l = 0.0, n = 4, Method = Meth,"
	    " HOM = (4, %12.5e, 0.0);\n", bn_prms.bn_scl[k]*bn_prms.bn[k+1]);
    k++;
    fprintf(outf, "oxy:  multipole, l = 0.0, n = 4, Method = Meth,"
	    " HOM = (4, %12.5e, 0.0);\n", bn_prms.bn_scl[k]*bn_prms.bn[k+1]);
    k++;
    fprintf(outf, "oyy:  multipole, l = 0.0, n = 4, Method = Meth,"
	    " HOM = (4, %12.5e, 0.0);\n", bn_prms.bn_scl[k]*bn_prms.bn[k+1]);
  }

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


void fit_ksi1(const double ksi_x, const double ksi_y)
{
  int          n_b4, i, m;
  double       **A, *b;
  ss_vect<tps> nus;

  const int m_max = 2;

  n_b4 = bn_prms.n_prm;

  b = dvector(1, m_max); A = dmatrix(1, m_max, 1, n_b4);

  no_mpoles(Sext); no_mpoles(Oct);

  printf("\n");
  for (i = 1; i <= n_b4; i++) {
    bn_prms.set_prm_dep(i-1);

    danot_(3);
    get_Map();
    danot_(4);
    K = MapNorm(Map, g, A1, A0, Map_res, 1); nus = dHdJ(K);

    m = 0;
    A[++m][i] = get_a(1e0, nus[3], 0, 0, 0, 0, 1);
    A[++m][i] = get_a(1e0, nus[4], 0, 0, 0, 0, 1);

    bn_prms.clr_prm_dep(i-1);
  }

  m = 0;
  b[++m] = -(get_b(1e0, nus[3], 0, 0, 0, 0, 1)-ksi_x);
  b[++m] = -(get_b(1e0, nus[4], 0, 0, 0, 0, 1)-ksi_y);

  prt_system(m, n_b4, A, b);

  SVD_lim(m, n_b4, A, b, bn_prms.bn_lim, bn_prms.svd_cut, bn_prms.bn,
	  bn_prms.dbn);

  bn_prms.set_dprm();

  for (i = 1; i <= n_b4; i++)
    bn_prms.bn[i] = get_bn(bn_prms.Fnum[i-1], 1, bn_prms.n[i-1]);

  printf("\nfit ksi:\n");
  for (i = 1; i <= n_b4; i++)
    printf("%13.5e", bn_prms.bn[i]);
  printf("\n");

  prt_mfile("flat_file.fit");

  free_dvector(b, 1, m_max); free_dmatrix(A, 1, m_max, 1, n_b4);
}


double min_dnu_f(double *b4s)
{
  int                 i;
  static double       chi2_ref = 1e30;
  double              chi2;
  tps                 h, h_re, h_im, K_re, K_im, K_re_scl;
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
  K = MapNorm(Map, g, A1, A0, Map_res, 1); h = get_h();
  CtoR(K, K_re, K_im); CtoR(g, h_re, h_im);
  K_re_scl = K_re*Id_scl; h_im = h_im*Id_scl;

  b.push_back(get_b(scl_h[0], h_im, 1, 0, 0, 0, 2));
  b.push_back(get_b(scl_h[0], h_im, 2, 0, 0, 0, 1));
  b.push_back(get_b(scl_h[0], h_im, 0, 0, 2, 0, 1));

  b.push_back(get_b(scl_h[0], h_im, 1, 0, 1, 1, 0));
  b.push_back(get_b(scl_h[0], h_im, 2, 1, 0, 0, 0));
  b.push_back(get_b(scl_h[0], h_im, 3, 0, 0, 0, 0));
  b.push_back(get_b(scl_h[0], h_im, 1, 0, 0, 2, 0));
  b.push_back(get_b(scl_h[0], h_im, 1, 0, 2, 0, 0));

  b.push_back(get_b(scl_h[1], h_im, 2, 0, 1, 1, 0));
  b.push_back(get_b(scl_h[1], h_im, 3, 1, 0, 0, 0));
  b.push_back(get_b(scl_h[1], h_im, 4, 0, 0, 0, 0));
  b.push_back(get_b(scl_h[1], h_im, 2, 0, 0, 2, 0));
  b.push_back(get_b(scl_h[1], h_im, 2, 0, 2, 0, 0));
  b.push_back(get_b(scl_h[1], h_im, 0, 0, 4, 0, 0));
  b.push_back(get_b(scl_h[1], h_im, 0, 0, 3, 1, 0));
  b.push_back(get_b(scl_h[1], h_im, 1, 1, 2, 0, 0));

  b.push_back(get_b(scl_dnu[0], K_re_scl, 1, 1, 0, 0, 1));
  b.push_back(get_b(scl_dnu[0], K_re_scl, 0, 0, 1, 1, 1));

  b.push_back(get_b(scl_dnu[1], K_re_scl, 2, 2, 0, 0, 0));
  b.push_back(get_b(scl_dnu[1], K_re_scl, 0, 0, 2, 2, 0));
  b.push_back(get_b(scl_dnu[1], K_re_scl, 1, 1, 1, 1, 0));

  b.push_back(get_b(scl_dnu[2], K_re_scl, 1, 1, 0, 0, 2));
  b.push_back(get_b(scl_dnu[2], K_re_scl, 0, 0, 1, 1, 2));

  if (NO >= 6) {
    b.push_back(get_b(scl_dnu[2], K_re_scl, 2, 2, 0, 0, 1));
    b.push_back(get_b(scl_dnu[2], K_re_scl, 0, 0, 2, 2, 1));
    b.push_back(get_b(scl_dnu[2], K_re_scl, 1, 1, 1, 1, 1));

    b.push_back(get_b(scl_dnu[2], K_re_scl, 1, 1, 0, 0, 3));
    b.push_back(get_b(scl_dnu[2], K_re_scl, 0, 0, 1, 1, 3));
  }

  if (NO >= 7) {
    if (!tune_conf) {
      b.push_back(get_b(scl_dnu[1], K_re_scl, 3, 3, 0, 0, 0));
      b.push_back(get_b(scl_dnu[1], K_re_scl, 2, 2, 1, 1, 0));
      b.push_back(get_b(scl_dnu[1], K_re_scl, 1, 1, 2, 2, 0));
      b.push_back(get_b(scl_dnu[1], K_re_scl, 0, 0, 3, 3, 0));
    } else {
      b.push_back(scl_dnu[3]
		  *(get_b(1e0, K_re, 3, 3, 0, 0, 0)
		    +get_b(1e0/(3e0*twoJ[X_]), K_re, 2, 2, 0, 0, 0)));
      b.push_back(scl_dnu[3]
		  *(get_b(1e0, K_re, 1, 1, 2, 2, 0)
		    +get_b(1e0/(2e0*twoJ[Y_]), K_re, 1, 1, 1, 1, 0)));

      b.push_back(scl_dnu[3]
		  *(get_b(1e0, K_re, 0, 0, 3, 3, 0)
		    +get_b(1e0/(3e0*twoJ[Y_]), K_re, 0, 0, 2, 2, 0)));
      b.push_back(scl_dnu[3]
		  *(get_b(1e0, K_re, 2, 2, 1, 1, 0)
		    +get_b(1e0/(2e0*twoJ[X_]), K_re, 1, 1, 1, 1, 0)));
    }
  }

  chi2 = 0e0;
  for (i = 0; i < (int)b.size(); i++)
    chi2 += sqr(b[i]);

  if (chi2 < chi2_ref) {
    prt_bn(bn_prms);
    prt_h_K();

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
  int    n_b4, i, j, m;
  double chi2_ref, **A, *b, *bn_ref;
  tps    h, h_re, h_im, K_re, K_im, K_re_scl;

  const int m_max = 50;

  n_b4 = bn_prms.n_prm;

  bn_ref = dvector(1, n_b4);
  b = dvector(1, m_max); A = dmatrix(1, m_max, 1, n_b4);

  printf("\n");
  for (i = 1; i <= n_b4; i++) {
    bn_prms.set_prm_dep(i-1);

    danot_(NO-1);
    get_Map();
    danot_(NO);
    K = MapNorm(Map, g, A1, A0, Map_res, 1); h = get_h();
    CtoR(K, K_re, K_im); CtoR(g, h_re, h_im);
    K_re_scl = K_re*Id_scl; h_im = h_im*Id_scl;

    m = 0;
    A[++m][i] = get_a(scl_h[0], h_im, 1, 0, 0, 0, 2);
    A[++m][i] = get_a(scl_h[0], h_im, 2, 0, 0, 0, 1);
    A[++m][i] = get_a(scl_h[0], h_im, 0, 0, 2, 0, 1);

    A[++m][i] = get_a(scl_h[0], h_im, 1, 0, 1, 1, 0);
    A[++m][i] = get_a(scl_h[0], h_im, 2, 1, 0, 0, 0);
    A[++m][i] = get_a(scl_h[0], h_im, 3, 0, 0, 0, 0);
    A[++m][i] = get_a(scl_h[0], h_im, 1, 0, 0, 2, 0);
    A[++m][i] = get_a(scl_h[0], h_im, 1, 0, 2, 0, 0);

    A[++m][i] = get_a(scl_h[1], h_im, 2, 0, 1, 1, 0);
    A[++m][i] = get_a(scl_h[1], h_im, 3, 1, 0, 0, 0);
    A[++m][i] = get_a(scl_h[1], h_im, 4, 0, 0, 0, 0);
    A[++m][i] = get_a(scl_h[1], h_im, 2, 0, 0, 2, 0);
    A[++m][i] = get_a(scl_h[1], h_im, 2, 0, 2, 0, 0);
    A[++m][i] = get_a(scl_h[1], h_im, 0, 0, 4, 0, 0);
    A[++m][i] = get_a(scl_h[1], h_im, 0, 0, 3, 1, 0);
    A[++m][i] = get_a(scl_h[1], h_im, 1, 1, 2, 0, 0);

    A[++m][i] = get_a(scl_dnu[0], K_re_scl, 1, 1, 0, 0, 1);
    A[++m][i] = get_a(scl_dnu[0], K_re_scl, 0, 0, 1, 1, 1);

    A[++m][i] = get_a(scl_dnu[1], K_re_scl, 2, 2, 0, 0, 0);
    A[++m][i] = get_a(scl_dnu[1], K_re_scl, 0, 0, 2, 2, 0);
    A[++m][i] = get_a(scl_dnu[1], K_re_scl, 1, 1, 1, 1, 0);

    A[++m][i] = get_a(scl_dnu[2], K_re_scl, 1, 1, 0, 0, 2);
    A[++m][i] = get_a(scl_dnu[2], K_re_scl, 0, 0, 1, 1, 2);

    if (NO >= 6) {
      A[++m][i] = get_a(scl_dnu[2], K_re_scl, 2, 2, 0, 0, 1);
      A[++m][i] = get_a(scl_dnu[2], K_re_scl, 0, 0, 2, 2, 1);
      A[++m][i] = get_a(scl_dnu[2], K_re_scl, 1, 1, 1, 1, 1);

      A[++m][i] = get_a(scl_dnu[2], K_re_scl, 1, 1, 0, 0, 3);
      A[++m][i] = get_a(scl_dnu[2], K_re_scl, 0, 0, 1, 1, 3);
    }

    if (NO >= 7) {
      if (!tune_conf) {
	A[++m][i] = get_a(scl_dnu[1], K_re_scl, 3, 3, 0, 0, 0);
	A[++m][i] = get_a(scl_dnu[1], K_re_scl, 2, 2, 1, 1, 0);
	A[++m][i] = get_a(scl_dnu[1], K_re_scl, 1, 1, 2, 2, 0);
	A[++m][i] = get_a(scl_dnu[1], K_re_scl, 0, 0, 3, 3, 0);
      } else {
	A[++m][i] =
	  scl_dnu[3]
	  *(get_a(1e0, K_re, 3, 3, 0, 0, 0)
	    +get_a(1e0/(3e0*twoJ[X_]), K_re, 2, 2, 0, 0, 0));
	A[++m][i] =
	  scl_dnu[3]
	  *(get_a(1e0, K_re, 1, 1, 2, 2, 0)
	    +get_a(1e0/(2e0*twoJ[Y_]), K_re, 1, 1, 1, 1, 0));

	A[++m][i] =
	  scl_dnu[3]
	  *(get_a(1e0, K_re, 0, 0, 3, 3, 0)
	    +get_a(1e0/(3e0*twoJ[Y_]), K_re, 0, 0, 2, 2, 0));
	A[++m][i] =
	  scl_dnu[3]
	  *(get_a(1e0, K_re, 2, 2, 1, 1, 0)
	    +get_a(1e0/(2e0*twoJ[X_]), K_re, 1, 1, 1, 1, 0));
       }
    }

    for (j = 1; j <= m; j++)
      A[j][i] *= bn_prms.bn_scl[i-1];

    bn_prms.clr_prm_dep(i-1);
  }

  m = 0;
  b[++m] = -get_b(scl_h[0], h_im, 1, 0, 0, 0, 2);
  b[++m] = -get_b(scl_h[0], h_im, 2, 0, 0, 0, 1);
  b[++m] = -get_b(scl_h[0], h_im, 0, 0, 2, 0, 1);

  b[++m] = -get_b(scl_h[0], h_im, 1, 0, 1, 1, 0);
  b[++m] = -get_b(scl_h[0], h_im, 2, 1, 0, 0, 0);
  b[++m] = -get_b(scl_h[0], h_im, 3, 0, 0, 0, 0);
  b[++m] = -get_b(scl_h[0], h_im, 1, 0, 0, 2, 0);
  b[++m] = -get_b(scl_h[0], h_im, 1, 0, 2, 0, 0);

  b[++m] = -get_b(scl_h[1], h_im, 2, 0, 1, 1, 0);
  b[++m] = -get_b(scl_h[1], h_im, 3, 1, 0, 0, 0);
  b[++m] = -get_b(scl_h[1], h_im, 4, 0, 0, 0, 0);
  b[++m] = -get_b(scl_h[1], h_im, 2, 0, 0, 2, 0);
  b[++m] = -get_b(scl_h[1], h_im, 2, 0, 2, 0, 0);
  b[++m] = -get_b(scl_h[1], h_im, 0, 0, 4, 0, 0);
  b[++m] = -get_b(scl_h[1], h_im, 0, 0, 3, 1, 0);
  b[++m] = -get_b(scl_h[1], h_im, 1, 1, 2, 0, 0);

  b[++m] = -get_b(scl_dnu[0], K_re_scl, 1, 1, 0, 0, 1);
  b[++m] = -get_b(scl_dnu[0], K_re_scl, 0, 0, 1, 1, 1);

  b[++m] = -get_b(scl_dnu[1], K_re_scl, 2, 2, 0, 0, 0);
  b[++m] = -get_b(scl_dnu[1], K_re_scl, 0, 0, 2, 2, 0);
  b[++m] = -get_b(scl_dnu[1], K_re_scl, 1, 1, 1, 1, 0);

  b[++m] = -get_b(scl_dnu[2], K_re_scl, 1, 1, 0, 0, 2);
  b[++m] = -get_b(scl_dnu[2], K_re_scl, 0, 0, 1, 1, 2);

  if (NO >= 6) {
    b[++m] = -get_b(scl_dnu[2], K_re_scl, 2, 2, 0, 0, 1);
    b[++m] = -get_b(scl_dnu[2], K_re_scl, 0, 0, 2, 2, 1);
    b[++m] = -get_b(scl_dnu[2], K_re_scl, 1, 1, 1, 1, 1);

    b[++m] = -get_b(scl_dnu[2], K_re_scl, 1, 1, 0, 0, 3);
    b[++m] = -get_b(scl_dnu[2], K_re_scl, 0, 0, 1, 1, 3);
  }

  if (NO >= 7) {
    if (!tune_conf) {
      b[++m] = -get_b(scl_dnu[1], K_re_scl, 3, 3, 0, 0, 0);
      b[++m] = -get_b(scl_dnu[1], K_re_scl, 2, 2, 1, 1, 0);
      b[++m] = -get_b(scl_dnu[1], K_re_scl, 1, 1, 2, 2, 0);
      b[++m] = -get_b(scl_dnu[1], K_re_scl, 0, 0, 3, 3, 0);
    } else {
      b[++m] =
	-scl_dnu[3]
	*(get_b(1e0, K_re, 3, 3, 0, 0, 0)
	  +get_b(1e0/(3e0*twoJ[X_]), K_re, 2, 2, 0, 0, 0));
      b[++m] =
	-scl_dnu[3]
	*(get_b(1e0, K_re, 1, 1, 2, 2, 0)
	  +get_b(1e0/(2e0*twoJ[Y_]), K_re, 1, 1, 1, 1, 0));

      b[++m] =
	-scl_dnu[3]
	*(get_b(1e0, K_re, 0, 0, 3, 3, 0)
	  +get_b(1e0/(3e0*twoJ[Y_]), K_re, 0, 0, 2, 2, 0));
      b[++m] =
	-scl_dnu[3]
	*(get_b(1e0, K_re, 2, 2, 1, 1, 0)
	  +get_b(1e0/(2e0*twoJ[X_]), K_re, 1, 1, 1, 1, 0));
    }
  }

  chi2_ref = chi2;
  
  chi2 = 0e0;
  for (j = 1; j <= m; j++)
    chi2 += sqr(b[j]);

  prt_system(m, n_b4, A, b);
  printf("\n%4d %12.5e -> %12.5e\n", n_iter, chi2_ref, chi2);

  prt_dnu();

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

  const int n_iter_max = 1000;

  n_b4 = bn_prms.n_prm;

  g = dvector(1, n_b4); h = dvector(1, n_b4);

  n_iter = 0; chi2 = 0e0;
  do {
    n_iter++;

    min_dnu_grad(chi2, dbn_max, g, h, cg_meth);

    prt_mfile("flat_file.fit");
    prt_bn(bn_prms);
    prt_h_K();
  } while ((dbn_max >  bn_prms.bn_tol) && (n_iter < n_iter_max));

  free_dvector(g, 1, n_b4); free_dvector(h, 1, n_b4);
  free_dvector(bn_prms.bn_lim, 1, n_b4); free_dvector(bn_prms.bn, 1, n_b4);
  free_dvector(bn_prms.dbn, 1, n_b4);
}


void min_dnu2(void)
{
  int    n_b3, i, j, iter;
  double **xi, fret;

  n_b3 = bn_prms.n_prm;

  xi = dmatrix(1, n_b3, 1, n_b3);

  // Set initial directions (unit vectors).
  for (i = 1; i <= n_b3; i++)
    for (j = 1; j <= n_b3; j++)
      xi[i][j] = (i == j)? 1e0 : 0e0;

  dpowell(bn_prms.bn, xi, n_b3, bn_prms.bn_tol, &iter, &fret, min_dnu_f);

  free_dmatrix(xi, 1, n_b3, 1, n_b3);
}


void prt_ps_long(const int n, const double twoJ[], const double delta)
{
  int              i, k;
  double           twoJ1[2], phi[2], ct_sum, ct_sum2, ct_mean, ct_sigma;
  double           delta_sum, delta_sum2, delta_mean, delta_sigma;
  ss_vect<double>  ps, ps_Fl;
  std::ofstream    outf;

  outf.open("ps_long.out");

  delta_sum = 0e0; delta_sum2 = 0e0;
  ct_sum = 0e0; ct_sum2 = 0e0;

  ps_Fl.zero();
  for (i = 1; i <= n; i++) {
    for (k = 0; k < 2; k++) {
      phi[k] = 2e0*M_PI*ranf();
      twoJ1[k] = twoJ[k]*normranf();
      while (twoJ1[k] < 0e0)
	twoJ1[k] = twoJ[k]*normranf();

      ps_Fl[2*k] = sqrt(twoJ1[k])*cos(phi[k]);
      ps_Fl[2*k+1] = -sqrt(twoJ1[k])*sin(phi[k]);
    }

    ps_Fl[delta_] = delta*normranf();

    delta_sum += ps_Fl[delta_]; delta_sum2 += sqr(ps_Fl[delta_]);

    ps = (A1*ps_Fl).cst();

    ps.propagate(1, n_elem);

    ct_sum += ps[ct_]; ct_sum2 += sqr(ps[ct_]);

    outf << std::scientific << std::setprecision(5)
	 << std::setw(4) << i << std::setw(13) << ps << std::endl;
  }

  delta_mean = delta_sum/n;
  delta_sigma = sqrt((n*delta_sum2-sqr(delta_sum))/(n*(n-1e0)));  
  ct_mean = ct_sum/n; ct_sigma = sqrt((n*ct_sum2-sqr(ct_sum))/(n*(n-1e0)));  

  std::cout << std::endl;
  std::cout << std::scientific << std::setprecision(3)
       << "delta = " << std::setw(11) << delta_mean
       << " +/-" << std::setw(10) << delta_sigma << std::endl;
  std::cout << std::scientific << std::setprecision(3)
       << "ct    = " << std::setw(11) << ct_mean
       << " +/-" << std::setw(10) << ct_sigma << std::endl;

  outf.close();
}


int main(int argc, char *argv[])
{
  int j;
 
  const long int  seed = 1121;
  // Center of straight.
  // const double beta[]  = {3.0, 3.0},
  //              A_max[] = {1.2e-3, 1.2e-3}, delta = 3e-2;
  const double beta[]  = {3.2, 2.1},
               A_max[] = {6e-3, 4e-3}, delta = 5e-2;

  rad_on    = false; H_exact        = false; totpath_on   = false;
  cavity_on = false; quad_fringe_on = false; emittance_on = false;
  IBS_on    = false;

  iniranf(seed); setrancut(3e0);

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

  if (false) {
    prt_H_long(10, M_PI, 10e-2, -544.7e3);
    exit(0);
  }

  if (false) {
    prt_ps_long(1000, twoJ, delta);
    exit(0);
  }

  if (true) {
    prt_map();
    exit(0);
  }

  if (false) {
    prt_h_K();
    exit(0);
  }

  if (false) {
    // MAX VI:
    bn_prms.add_prm("o1", 4, 5e5, 1.0);
    bn_prms.add_prm("o2", 4, 5e5, 1.0);
    bn_prms.add_prm("o3", 4, 5e5, 1.0);
    bn_prms.add_prm("o4", 4, 5e5, 1.0);

    bn_prms.add_prm("o1", 6, 5e10, 1.0);
    bn_prms.add_prm("o2", 6, 5e10, 1.0);
    bn_prms.add_prm("o3", 6, 5e10, 1.0);
    bn_prms.add_prm("o4", 6, 5e10, 1.0);
  } else {
    // SLS-2:
    bn_prms.add_prm("sfh",  3, 5e5, 1.0);
    bn_prms.add_prm("sd",   3, 5e5, 1.0);
    bn_prms.add_prm("sfmh", 3, 5e5, 1.0);
    bn_prms.add_prm("sdm",  3, 5e5, 1.0);

    bn_prms.add_prm("sxxh", 3, 5e5, 1.0);
    bn_prms.add_prm("sxyh", 3, 5e5, 1.0);
    bn_prms.add_prm("syyh", 3, 5e5, 1.0);

    // bn_prms.add_prm("ocx",  4, 5e10, 1.0);
    // bn_prms.add_prm("ocxm", 4, 5e10, 1.0);
    // bn_prms.add_prm("sd",   4, 5e10, 1.0);

    bn_prms.add_prm("oxx",  4, 5e10, 1.0);
    bn_prms.add_prm("oxy",  4, 5e10, 1.0);
    bn_prms.add_prm("oyy",  4, 5e10, 1.0);

    // bn_prms.add_prm("ocx",  6, 5e10, 1e5);
    // bn_prms.add_prm("ocxm", 6, 5e10, 1e5);

    // bn_prms.add_prm("oxx",  6, 5e10, 1e5);
    // bn_prms.add_prm("oxy",  6, 5e10, 1e5);
    // bn_prms.add_prm("oyy",  6, 5e10, 1e5);
  }

  // Step is 1.0 for conjugated gradient method.
  bn_prms.bn_tol = 1e-1; bn_prms.svd_cut = 1e-10; bn_prms.step = 1.0;

  no_mpoles(Sext);
  no_mpoles(Oct); no_mpoles(Dodec);

  bn_prms.ini_prm();

  fit_ksi1(0e0, 0e0);

  // min_dnu(true);
}
