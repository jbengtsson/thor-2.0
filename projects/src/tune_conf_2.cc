
#include <cfloat>

#define NO 7

#include "thor_lib.h"

int no_tps = NO,

#define DOF_3 0

#if !DOF_3
  ndpt_tps = 5;
#else
// Requires that cavity is turned on.
ndpt_tps = 0;
#endif


extern tps          K, g;
extern ss_vect<tps> Map, A0, A1, Map_res;

double       chi2 = 0e0, *f_lm, **A_lm;
tps          h_re, h_im, K_re, K_im;
ss_vect<tps> nus;

const bool   fit_ksi  = !true, symm  = true;
const int    n_cell   = 1,     n_cut = 0;
const double tpsa_eps = 1e-30;

// MAX-IV              1,
// SLS-2               2,
// DIAMOND             3,
// DIAMOND with VMX    4,
// DIAMOND-II 4-BA     5,
// DIAMOND-II 6-HMBA   6,
// DIAMOND-II 6-RB-BA  7,
// DIAMOND-II 8-BA     8,
// DIAMOND-II 8-HMBA   9,
// SLS-2              10.
const int lat_case = 10, n_prt = 8;

// Center of straight.
const double
  beta_inj[][2] =
    {{ 3.0, 3.0}, {3.4, 1.9}, { 9.9, 5.4}, { 9.9, 5.4}, {10.6, 8.6},
     {16.2, 4.6}, {5.3, 2.0}, {14.0, 4.5}, {10.5, 5.2}, { 3.4, 1.9}},
  A_max[][2] =
    {{1.2e-3, 1.2e-3}, {6e-3, 4e-3}, {15e-3, 8e-3}, {15e-3, 8e-3}, {5e-3, 3e-3},
     {6e-3, 4e-3},     {7e-3, 4e-3}, { 2e-3, 1e-3},  {2e-3, 1e-3}, {9e-3, 5e-3}
    },
  delta_max[] =
    {3e-2, 5e-2, 3e-2, 3e-2, 3e-2,
     3e-2, 3e-2, 3e-2, 3e-2, 3e-2};


#if true
// Sextupoles.
const bool   oct = false;
const double scl_h[]      = {1e0,  1e0,  1e-1},
             scl_dnu[]    = {1e-5, 1e-5, 1e-5, 1e-5},
             scl_ksi[]    = {1e5,  1e-1, 1e-5},
             scl_dnu_conf = 1e0;
#else
// Octupoles.
const bool   oct = true;
const double scl_h[]      = {1e0,  1e0,  1e0},
             scl_dnu[]    = {1e-1, 1e-5, 1e-5, 1e-5},
             scl_ksi[]    = {1e5,  1e-5, 1e-5},
             scl_dnu_conf = 1e-1;
#endif

struct param_type {
private:

public:
  int                 m_constr, n_prm, svd_n_cut;
  double              bn_tol, step;
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
  double L;

  const bool scale = true;
  const int  n_prt = 4;

  n_prm = Fnum.size();

  bn_prms.bn_lim = dvector(1, n_prm); bn_prms.bn = dvector(1, n_prm);
  bn_prms.dbn = dvector(1, n_prm);

  printf("\nInitial bn (scale factor in parenthesis):\n");
  printf("  No of Families: %1d\n", n_prm);
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

    if (scale) {
      L = get_L(Fnum[i-1], 1);
      if (L == 0e0) L = 1e0;
      bn_scl[i-1] = 1e0/sqrt(get_n_Kids(Fnum[i-1])*L);
    } else
      bn_scl[i-1] = 1e0;

    bn[i] /= bn_scl[i-1];
    printf(" %12.5e (%9.3e)", bn_scl[i-1]*bn[i], bn_scl[i-1]);
    if (i % n_prt == 0) printf("\n");
  }
  if (n_prm % n_prt != 0) printf("\n");
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
    printf(" %12.5e", bn[i]);
    if (i % n_prt == 0) printf("\n");
  }
  if (n_prm % n_prt != 0) printf("\n");

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
    printf(" %12.5e", bn_scl[i-1]*bn[i]);
    if (i % n_prt == 0) printf("\n");
  }
  if (n_prm % n_prt != 0) printf("\n");
}


// ---> Tune confinement.

static double xsav;
static tps (*func_save)(const double, const double);

tps gauss_quad(tps (*func)(const double), const double a, const double b)
{
  int    j;
  double xr, xm, dx;
  tps    s;

  static double x[] =
    {0.0, 0.1488743389, 0.4333953941, 0.6794095682, 0.8650633666, 0.9739065285};
  static double w[] =
    {0.0, 0.2955242247, 0.2692667193, 0.2190863625, 0.1494513491, 0.0666713443};

  xm = 0.5*(b+a);
  xr = 0.5*(b-a);
  s = 0;
  for (j = 1; j <= 5; j++) {
    dx = xr*x[j];
    s += w[j]*((*func)(xm+dx)+(*func)(xm-dx));
  }
  return s *= xr;
}


double gauss_quad_y0(const double x) { return 0e0; }


double gauss_quad_y1(const double x) { return twoJ[Y_]; }


tps gauss_quad_fy(const double y)  { return (*func_save)(xsav, y); }


tps gauss_quad_fx(const double x)
{
  xsav = x;
  return gauss_quad(gauss_quad_fy, gauss_quad_y0(x),
		    gauss_quad_y1(x));
}


tps gauss_quad_2d(tps (*func)(const double, const double),
		     const double x1, const double x2)
{
  func_save = func;
  return gauss_quad(gauss_quad_fx, x1, x2);
}


tps f_gauss_quad_2d_dnu(double x, double y)
{
  int          k;
  tps          dnu[2];
  ss_vect<tps> ps;

    ps.identity();
    ps[x_] = sqrt(x); ps[px_] = sqrt(x);
    ps[y_] = sqrt(y); ps[py_] = sqrt(y);
    ps[delta_] = 0*delta_max[lat_case-1];
    for (k = 0; k < 2; k++) {
      dnu[k] = (nus[k+3]-nus[k+3].cst())*ps;
      // Compute absolute value.
      if (dnu[k].cst() < 0e0) dnu[k] = -dnu[k];
    }

    return dnu[X_]*dnu[Y_]/(twoJ[X_]*twoJ[Y_]);
}



tps f_gauss_quad_2d(double x, double y)
{
  int          k, jj[ss_dim];
  tps          dK;
  ss_vect<tps> ps;

    ps.identity();
    ps[x_] = sqrt(x); ps[px_] = sqrt(x);
    ps[y_] = sqrt(y); ps[py_] = sqrt(y);
    ps[delta_] = 0*delta_max[lat_case-1];

    dK = K_re;
    for (k = 0; k < ss_dim; k++)
      jj[k] = 0;
    jj[x_] = 1; jj[px_] = 1;
    dK.pook(jj, 0e0);
    jj[x_] = 0; jj[px_] = 0; jj[y_] = 1; jj[py_] = 1;
    dK.pook(jj, 0e0);
    dK = dK*ps;
    // Compute absolute value.
    if (dK.cst() < 0e0) dK = -dK;
    // std::cout << std::scientific << std::setprecision(3)
    // 	      << "\n |dK| = " << dK << "\n";

    return dK/(twoJ[X_]*twoJ[Y_]);
}

// <--- Tune confinement.


void no_mpoles(const int n)
{
  int j;

  printf("\nzeroing multipoles: %d\n", n);
  for (j = 0; j < n_elem; j++)
    if (elem[j].kind == Mpole)
      set_bn(elem[j].Fnum, elem[j].Knum, n, 0e0);
}


void get_map_n(const int n)
{
  int          k;
  ss_vect<tps> map_n;

  get_Map(); map_n.identity();
  for (k = 0; k < n; k++)
    map_n = map_n*Map;
  Map = map_n;
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
 std::ofstream outf;

  danot_(NO-1);
  get_map_n(n_cell);

  file_wr(outf, "map.out");
  outf << Map[ct_];
  outf.close();
}


void prt_h_K(void)
{
  std::ofstream outf;

  file_wr(outf, "h.out");
  outf << h_re*Id_scl << h_im*Id_scl;
  outf.close();

  file_wr(outf, "K.out");
  outf << K_re*Id_scl << K_im*Id_scl;
  outf.close();

  file_wr(outf, "nus.out");
  nus = dHdJ(K);
  // Remove numeric noise.
  daeps_(1e-5);
  outf << nus[3]*Id_scl << nus[4]*Id_scl;
  daeps_(tpsa_eps);
  outf.close();
}


void prt_h_K_sym(void)
{
  int          k;
  double       mu[2];
  ss_vect<tps> R_half;

  danot_(NO-1);
  get_map_n(n_cell);
  danot_(NO);
  K = MapNorm(Map, g, A1, A0, Map_res, no_tps);
  CtoR(K, K_re, K_im); nus = dHdJ(K);

#if 0
  // Compute h at symmetry point of lattice; for a lattice with two super
  // periods.
  printf("\nprt_h_K: nu/2 = [%7.5f, %7.5f]\n",
	 nus[0].cst()/2e0, nus[1].cst()/2e0);
  R_half.identity();
  mu[0] = 2e0*M_PI*nus[0].cst()/2e0; mu[1] = 2e0*M_PI*nus[1].cst()/2e0;
  for (k = 0; k < 2; k++) {
    R_half[2*k]   =  cos(mu[k])*tps(0e0, 2*k+1) + sin(mu[k])*tps(0e0, 2*k+2);
    R_half[2*k+1] = -sin(mu[k])*tps(0e0, 2*k+1) + cos(mu[k])*tps(0e0, 2*k+2);
  }
  CtoR(get_h()*Inv(R_half), h_re, h_im);
#else
  CtoR(get_h(), h_re, h_im);
#endif

  prt_h_K();
}


void prt_dnu(tps &K)
{
  tps          K_re, K_im, nu_scl[2];
  ss_vect<tps> nus;

  CtoR(K, K_re, K_im); nus = dHdJ(K);
  nu_scl[0] = nus[3]*Id_scl; nu_scl[1] = nus[4]*Id_scl;

  printf("\ndnu:\n");
  printf(" %8.5f",   h_ijklm(nu_scl[X_], 1, 1, 0, 0, 0));
  printf(" %8.5f,",  h_ijklm(nu_scl[X_], 0, 0, 1, 1, 0));

  printf(" %8.5f",   h_ijklm(nu_scl[X_], 2, 2, 0, 0, 0));
  printf(" %8.5f",   h_ijklm(nu_scl[X_], 1, 1, 1, 1, 0));
  printf(" %8.5f,",  h_ijklm(nu_scl[X_], 0, 0, 2, 2, 0));

  printf(" %8.5f",   h_ijklm(nu_scl[X_], 3, 3, 0, 0, 0));
  printf(" %8.5f",   h_ijklm(nu_scl[X_], 2, 2, 1, 1, 0));
  printf(" %8.5f",   h_ijklm(nu_scl[X_], 1, 1, 2, 2, 0));
  printf(" %8.5f\n", h_ijklm(nu_scl[X_], 0, 0, 3, 3, 0));

  printf(" %8.5f",   h_ijklm(nu_scl[Y_], 1, 1, 0, 0, 0));
  printf(" %8.5f,",  h_ijklm(nu_scl[Y_], 0, 0, 1, 1, 0));

  printf(" %8.5f",   h_ijklm(nu_scl[Y_], 2, 2, 0, 0, 0));
  printf(" %8.5f",   h_ijklm(nu_scl[Y_], 1, 1, 1, 1, 0));
  printf(" %8.5f,",  h_ijklm(nu_scl[Y_], 0, 0, 2, 2, 0));

  printf(" %8.5f",   h_ijklm(nu_scl[Y_], 3, 3, 0, 0, 0));
  printf(" %8.5f",   h_ijklm(nu_scl[Y_], 2, 2, 1, 1, 0));
  printf(" %8.5f",   h_ijklm(nu_scl[Y_], 1, 1, 2, 2, 0));
  printf(" %8.5f\n", h_ijklm(nu_scl[Y_], 0, 0, 3, 3, 0));

  printf("\n %8.5f", h_ijklm(nu_scl[X_], 0, 0, 0, 0, 2));
  printf(" %8.5f",   h_ijklm(nu_scl[X_], 1, 1, 0, 0, 2));
  printf(" %8.5f\n", h_ijklm(nu_scl[X_], 0, 0, 1, 1, 2));

  printf(" %8.5f",   h_ijklm(nu_scl[Y_], 0, 0, 0, 0, 2));
  printf(" %8.5f",   h_ijklm(nu_scl[Y_], 1, 1, 0, 0, 2));
  printf(" %8.5f\n", h_ijklm(nu_scl[Y_], 0, 0, 1, 1, 2));

  printf("\nTune confinement:\n");
  printf(" %11.3e %11.3e\n",
	 h_ijklm(K_re/(3e0*twoJ[X_]), 2, 2, 0, 0, 0),
	 h_ijklm(K_re, 3, 3, 0, 0, 0));
  printf(" %11.3e %11.3e\n",
	 h_ijklm(K_re/(3e0*twoJ[Y_]), 0, 0, 2, 2, 0),
	 h_ijklm(K_re, 0, 0, 3, 3, 0));
  printf(" %11.3e %11.3e %11.3e\n",
	 h_ijklm(K_re*Id_scl, 1, 1, 1, 1, 0),
	 h_ijklm(K_re*Id_scl, 2, 2, 1, 1, 0),
	 h_ijklm(K_re*Id_scl, 1, 1, 2, 2, 0));
}


void prt_system(const int m, const int n_b2, double **A, double *b)
{
  int i, j, d, n_h;

  printf("\n Ax = b:\n");
  for (j = 1; j <= n_b2; j++)
    printf("%11d", j);
  printf("\n");
  d = (symm)? 0 : 16;
  for (i = 1; i <= m; i++) {
    if (i-1 == 0)
      printf("1st order chromatic\n");
    else if (i-1 == 3)
      printf("1st order geometric\n");
    else if ((NO >= 5) && (i-1 == 3+5))
      printf("2nd order geometric\n");
    else if ((NO >= 6) && (i-1 == 3+5+8))
      printf("3rd order geometric\n");

    n_h = 0;
    if (NO >= 3+1) n_h += 3 + 5;   // 8.
    if (NO >= 4+1) n_h += 8;       // 16.
    if (NO >= 5+1) n_h += 14;      // 30.
    if (!symm)     n_h += 16 + 14;

    if (i-1 == n_h)
      printf("linear chromaticity\n");
    else if (i-1 == n_h+2)
      printf("ampl. dependant tune shift\n");
    else if (i-1 == n_h+2+3)
      printf("2nd order chromaticity\n");
    else if (i-1 == n_h+2+3+2)
      printf("|dnu|\n");
    else if (i-1 == n_h+2+3+2+1)
      printf("cross terms\n");
    else if (i-1 == n_h+2+3+2+1+3)
      printf("3rd order chromaticity\n");
    else if (i-1 == n_h+2+3+2+1+3+2) {
      printf("ampl. dependant tune shift\n");
    }

    printf("%4d", i);
    for (j = 1; j <= n_b2; j++)
      printf("%11.3e", A[i][j]);
    printf("%11.3e\n", b[i]);
  }

  prt_dnu(K);
}


void prt_b2(const std::vector<int> &b2_Fnum)
{
  long int loc;
  int      k;
  FILE     *outf;

  const std::string file_name = "dnu.out";

  outf = file_write(file_name.c_str());

  fprintf(outf, "\n");
  for (k = 0; k < (int)b2_Fnum.size(); k++) {
    loc = get_loc(b2_Fnum[k], 1) - 1;
    fprintf(outf,
	    "%-8s: quadrupole, l = %7.5f"
	    ", k = %12.5e, n = nquad, Method = Meth;\n",
	    elem[loc].Name, elem[loc].L, get_bn(b2_Fnum[k], 1, Quad));
  }

  fclose(outf);
}


void prt_bn(const param_type &bn_prms)
{
  long int loc;
  int      k;
  FILE     *outf;

  const std::string file_name = "dnu.out";

  outf = file_write(file_name.c_str());

  fprintf(outf, "\n");
  for (k = 0; k < bn_prms.n_prm; k++) {
    loc = get_loc(bn_prms.Fnum[k], 1) - 1;
    if (!oct)
      fprintf(outf,
	      "%-8s: sextupole, l = %7.5f"
	      ", k = %12.5e, n = nsext, Method = Meth;\n",
	      elem[loc].Name, elem[loc].L, bn_prms.bn_scl[k]*bn_prms.bn[k+1]);
      else
	fprintf(outf,
		"%-8s: multipole, l = %7.5f,"
		"\n          hom = (3, %12.5e, 0e0, 4, %12.5e, 0e0),"
		"\n          n = nsext, Method = Meth;\n",
		elem[loc].Name, elem[loc].L, elem[loc].mpole->bn[Sext-1],
		bn_prms.bn_scl[k]*bn_prms.bn[k+1]);
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
  int          n_bn, i, j, m;
  double       **A, *b;
  ss_vect<tps> nus;

  const int    m_max = 2;
  const double s_cut = 1e-10;

  n_bn = bn_prms.n_prm;

  b = dvector(1, m_max); A = dmatrix(1, m_max, 1, n_bn);

  no_mpoles(Sext); no_mpoles(Oct);

  printf("\n");
  for (i = 1; i <= n_bn; i++) {
    bn_prms.set_prm_dep(i-1);

    danot_(3);
    get_map_n(n_cell);
    danot_(4);
    K = MapNorm(Map, g, A1, A0, Map_res, 1); nus = dHdJ(K);

    m = 0;
    A[++m][i] = get_a(1e0, nus[3], 0, 0, 0, 0, 1);
    A[++m][i] = get_a(1e0, nus[4], 0, 0, 0, 0, 1);

    for (j = 1; j <= m; j++)
      A[j][i] *= bn_prms.bn_scl[i-1];

    bn_prms.clr_prm_dep(i-1);
  }

  m = 0;
  b[++m] = -(get_b(1e0, nus[3], 0, 0, 0, 0, 1)-ksi_x);
  b[++m] = -(get_b(1e0, nus[4], 0, 0, 0, 0, 1)-ksi_y);

  prt_system(m, n_bn, A, b);

  SVD_lim(m, n_bn, A, b, bn_prms.bn_lim, s_cut, bn_prms.bn, bn_prms.dbn);

  bn_prms.set_dprm();

  printf("\nfit ksi:\n");
  for (i = 1; i <= n_bn; i++) {
    printf(" %12.5e", bn_prms.bn_scl[i-1]*bn_prms.bn[i]);
    if (i % n_prt == 0) printf("\n");
  }
  if (n_bn % n_prt != 0) printf("\n");

  prt_mfile("flat_file.fit");
  prt_bn(bn_prms);

  danot_(3);
  get_map_n(n_cell);
  danot_(4);
  K = MapNorm(Map, g, A1, A0, Map_res, 1);
  CtoR(get_h(), h_re, h_im); nus = dHdJ(K);
  prt_h_K();

  free_dvector(b, 1, m_max); free_dmatrix(A, 1, m_max, 1, n_bn);
}


double get_f(double *bns)
{
  int                 i;
  static double       chi2_ref = 1e30;
  double              chi2;
  tps                 K_re_scl, h_re_scl, h_im_scl, dnu;
  std::vector<double> b;

  const bool prt = false;

  n_powell++;

  // Do not change parameters.
  if (prt) printf("get_f:\n");
  for (i = 1; i <= bn_prms.n_prm; i++) {
    set_bn(bn_prms.Fnum[i-1], bn_prms.n[i-1], bn_prms.bn_scl[i-1]*bns[i]);
    if (prt) printf(" %12.5e", bn_prms.bn_scl[i-1]*bns[i]);
  }
  if (prt) printf("\n");

  danot_(NO-1);
  get_map_n(n_cell);
  danot_(NO);
  K = MapNorm(Map, g, A1, A0, Map_res, 1); nus = dHdJ(K);
  CtoR(K, K_re, K_im); K_re_scl = K_re*Id_scl;
  CtoR(get_h(), h_re, h_im); h_re_scl = h_re*Id_scl; h_im_scl = h_im*Id_scl;

  dnu = gauss_quad_2d(f_gauss_quad_2d, 0e0, twoJ[X_]);
  printf("\n|dnu| = %9.3e\n", dnu.cst());

  b.push_back(get_b(scl_h[0], h_re_scl, 1, 0, 0, 0, 2));
  b.push_back(get_b(scl_h[0], h_re_scl, 2, 0, 0, 0, 1));
  b.push_back(get_b(scl_h[0], h_re_scl, 0, 0, 2, 0, 1));

  b.push_back(get_b(scl_h[0], h_re_scl, 1, 0, 1, 1, 0));
  b.push_back(get_b(scl_h[0], h_re_scl, 2, 1, 0, 0, 0));
  b.push_back(get_b(scl_h[0], h_re_scl, 3, 0, 0, 0, 0));
  b.push_back(get_b(scl_h[0], h_re_scl, 1, 0, 0, 2, 0));
  b.push_back(get_b(scl_h[0], h_re_scl, 1, 0, 2, 0, 0));

  if (NO >= 5) {
    b.push_back(get_b(scl_h[1], h_re_scl, 2, 0, 1, 1, 0));
    b.push_back(get_b(scl_h[1], h_re_scl, 3, 1, 0, 0, 0));
    b.push_back(get_b(scl_h[1], h_re_scl, 4, 0, 0, 0, 0));
    b.push_back(get_b(scl_h[1], h_re_scl, 2, 0, 0, 2, 0));
    b.push_back(get_b(scl_h[1], h_re_scl, 2, 0, 2, 0, 0));
    b.push_back(get_b(scl_h[1], h_re_scl, 0, 0, 4, 0, 0));
    b.push_back(get_b(scl_h[1], h_re_scl, 0, 0, 3, 1, 0));
    b.push_back(get_b(scl_h[1], h_re_scl, 1, 1, 2, 0, 0));
  }

  if (NO >= 6) {
    b.push_back(get_b(scl_h[2], h_re_scl, 5, 0, 0, 0, 0));
    b.push_back(get_b(scl_h[2], h_re_scl, 4, 1, 0, 0, 0));
    b.push_back(get_b(scl_h[2], h_re_scl, 3, 2, 0, 0, 0));
    b.push_back(get_b(scl_h[2], h_re_scl, 3, 0, 2, 0, 0));
    b.push_back(get_b(scl_h[2], h_re_scl, 2, 1, 2, 0, 0));
    b.push_back(get_b(scl_h[2], h_re_scl, 2, 1, 0, 2, 0));
    b.push_back(get_b(scl_h[2], h_re_scl, 3, 0, 0, 2, 0));
    b.push_back(get_b(scl_h[2], h_re_scl, 1, 0, 4, 0, 0));
    b.push_back(get_b(scl_h[2], h_re_scl, 1, 0, 0, 4, 0));
    b.push_back(get_b(scl_h[2], h_re_scl, 3, 0, 1, 1, 0));
    b.push_back(get_b(scl_h[2], h_re_scl, 2, 1, 1, 1, 0));
    b.push_back(get_b(scl_h[2], h_re_scl, 1, 0, 3, 1, 0));
    b.push_back(get_b(scl_h[2], h_re_scl, 1, 0, 1, 3, 0));
    b.push_back(get_b(scl_h[2], h_re_scl, 1, 0, 2, 2, 0));
  }

  if (!symm) {
    b.push_back(get_b(scl_h[0], h_im_scl, 1, 0, 0, 0, 2));
    b.push_back(get_b(scl_h[0], h_im_scl, 2, 0, 0, 0, 1));
    b.push_back(get_b(scl_h[0], h_im_scl, 0, 0, 2, 0, 1));

    b.push_back(get_b(scl_h[0], h_im_scl, 1, 0, 1, 1, 0));
    b.push_back(get_b(scl_h[0], h_im_scl, 2, 1, 0, 0, 0));
    b.push_back(get_b(scl_h[0], h_im_scl, 3, 0, 0, 0, 0));
    b.push_back(get_b(scl_h[0], h_im_scl, 1, 0, 0, 2, 0));
    b.push_back(get_b(scl_h[0], h_im_scl, 1, 0, 2, 0, 0));

    b.push_back(get_b(scl_h[1], h_im_scl, 2, 0, 1, 1, 0));
    b.push_back(get_b(scl_h[1], h_im_scl, 3, 1, 0, 0, 0));
    b.push_back(get_b(scl_h[1], h_im_scl, 4, 0, 0, 0, 0));
    b.push_back(get_b(scl_h[1], h_im_scl, 2, 0, 0, 2, 0));
    b.push_back(get_b(scl_h[1], h_im_scl, 2, 0, 2, 0, 0));
    b.push_back(get_b(scl_h[1], h_im_scl, 0, 0, 4, 0, 0));
    b.push_back(get_b(scl_h[1], h_im_scl, 0, 0, 3, 1, 0));
    b.push_back(get_b(scl_h[1], h_im_scl, 1, 1, 2, 0, 0));

    if (NO >= 6) {
      b.push_back(get_b(scl_h[2], h_im_scl, 5, 0, 0, 0, 0));
      b.push_back(get_b(scl_h[2], h_im_scl, 4, 1, 0, 0, 0));
      b.push_back(get_b(scl_h[2], h_im_scl, 3, 2, 0, 0, 0));
      b.push_back(get_b(scl_h[2], h_im_scl, 3, 0, 2, 0, 0));
      b.push_back(get_b(scl_h[2], h_im_scl, 2, 1, 2, 0, 0));
      b.push_back(get_b(scl_h[2], h_im_scl, 2, 1, 0, 2, 0));
      b.push_back(get_b(scl_h[2], h_im_scl, 3, 0, 0, 2, 0));
      b.push_back(get_b(scl_h[2], h_im_scl, 1, 0, 4, 0, 0));
      b.push_back(get_b(scl_h[2], h_im_scl, 1, 0, 0, 4, 0));
      b.push_back(get_b(scl_h[2], h_im_scl, 3, 0, 1, 1, 0));
      b.push_back(get_b(scl_h[2], h_im_scl, 2, 1, 1, 1, 0));
      b.push_back(get_b(scl_h[2], h_im_scl, 1, 0, 3, 1, 0));
      b.push_back(get_b(scl_h[2], h_im_scl, 1, 0, 1, 3, 0));
      b.push_back(get_b(scl_h[2], h_im_scl, 1, 0, 2, 2, 0));
    }
  }

  b.push_back(get_b(scl_ksi[0], K_re_scl, 1, 1, 0, 0, 1));
  b.push_back(get_b(scl_ksi[0], K_re_scl, 0, 0, 1, 1, 1));

  if (NO >= 5) {
    b.push_back(get_b(scl_dnu[0],   K_re_scl, 2, 2, 0, 0, 0));
    b.push_back(get_b(scl_dnu[0],   K_re_scl, 0, 0, 2, 2, 0));
    b.push_back(get_b(scl_dnu[0],   K_re_scl, 1, 1, 1, 1, 0));

    b.push_back(get_b(scl_ksi[1],   K_re_scl, 1, 1, 0, 0, 2));
    b.push_back(get_b(scl_ksi[1],   K_re_scl, 0, 0, 1, 1, 2));

    b.push_back(get_b(scl_dnu_conf, dnu,      0, 0, 0, 0, 0));
  }

  if (NO >= 6) {
    b.push_back(get_b(scl_dnu[1], K_re_scl, 2, 2, 0, 0, 1));
    b.push_back(get_b(scl_dnu[1], K_re_scl, 0, 0, 2, 2, 1));
    b.push_back(get_b(scl_dnu[1], K_re_scl, 1, 1, 1, 1, 1));

    b.push_back(get_b(scl_ksi[2], K_re_scl, 1, 1, 0, 0, 3));
    b.push_back(get_b(scl_ksi[2], K_re_scl, 0, 0, 1, 1, 3));
  }

  if (NO >= 7) {
    b.push_back(get_b(scl_dnu[2], K_re_scl, 3, 3, 0, 0, 0));
    b.push_back(get_b(scl_dnu[2], K_re_scl, 2, 2, 1, 1, 0));
    b.push_back(get_b(scl_dnu[2], K_re_scl, 1, 1, 2, 2, 0));
    b.push_back(get_b(scl_dnu[2], K_re_scl, 0, 0, 3, 3, 0));
  }

  if (NO >= 9) {
    b.push_back(get_b(scl_dnu[3], K_re_scl, 4, 4, 0, 0, 0));
    b.push_back(get_b(scl_dnu[3], K_re_scl, 3, 3, 1, 1, 0));
    b.push_back(get_b(scl_dnu[3], K_re_scl, 2, 2, 2, 2, 0));
    b.push_back(get_b(scl_dnu[3], K_re_scl, 1, 1, 3, 3, 0));
    b.push_back(get_b(scl_dnu[3], K_re_scl, 0, 0, 4, 4, 0));
  }

  chi2 = 0e0;
  for (i = 0; i < (int)b.size(); i++)
    chi2 += sqr(b[i]);

  if (chi2 < chi2_ref) {
    prt_bn(bn_prms);

    printf("\n%3d %12.5e -> %12.5e\n", n_powell, chi2_ref, chi2);
    printf("b & bn:\n");
    for (i = 0; i < (int)b.size(); i++)
      printf("%11.3e", b[i]);
    printf("\n");
    for (i = 1; i <= bn_prms.n_prm; i++) 
      printf("%11.3e", bns[i]);
    printf("\n");
  }

  chi2_ref = min(chi2, chi2_ref);

  return chi2;
}


void get_f_grad(const int n_bn, double *f, double **A, double &chi2, int &m)
{
  int    i, j;
  tps    K_re_scl, h_re_scl, h_im_scl, dnu;

  // printf("\n");
  for (i = 1; i <= n_bn; i++) {
    bn_prms.set_prm_dep(i-1);

    danot_(NO-1);
    get_map_n(n_cell);
    danot_(NO);
    K = MapNorm(Map, g, A1, A0, Map_res, 1); nus = dHdJ(K);
    CtoR(K, K_re, K_im); K_re_scl = K_re*Id_scl;
    CtoR(get_h(), h_re, h_im); h_re_scl = h_re*Id_scl; h_im_scl = h_im*Id_scl;

    dnu = gauss_quad_2d(f_gauss_quad_2d, 0e0, twoJ[X_]);
    // std::cout << std::scientific << std::setprecision(3)
    // 	      << "\n |dnu| = " << dnu << "\n";

    m = 0;
    A[++m][i] = get_a(scl_h[0], h_re_scl, 1, 0, 0, 0, 2);
    A[++m][i] = get_a(scl_h[0], h_re_scl, 2, 0, 0, 0, 1);
    A[++m][i] = get_a(scl_h[0], h_re_scl, 0, 0, 2, 0, 1);

    A[++m][i] = get_a(scl_h[0], h_re_scl, 1, 0, 1, 1, 0);
    A[++m][i] = get_a(scl_h[0], h_re_scl, 2, 1, 0, 0, 0);
    A[++m][i] = get_a(scl_h[0], h_re_scl, 3, 0, 0, 0, 0);
    A[++m][i] = get_a(scl_h[0], h_re_scl, 1, 0, 0, 2, 0);
    A[++m][i] = get_a(scl_h[0], h_re_scl, 1, 0, 2, 0, 0);

    if (NO >= 5) {
      A[++m][i] = get_a(scl_h[1], h_re_scl, 2, 0, 1, 1, 0);
      A[++m][i] = get_a(scl_h[1], h_re_scl, 3, 1, 0, 0, 0);
      A[++m][i] = get_a(scl_h[1], h_re_scl, 4, 0, 0, 0, 0);
      A[++m][i] = get_a(scl_h[1], h_re_scl, 2, 0, 0, 2, 0);
      A[++m][i] = get_a(scl_h[1], h_re_scl, 2, 0, 2, 0, 0);
      A[++m][i] = get_a(scl_h[1], h_re_scl, 0, 0, 4, 0, 0);
      A[++m][i] = get_a(scl_h[1], h_re_scl, 0, 0, 3, 1, 0);
      A[++m][i] = get_a(scl_h[1], h_re_scl, 1, 1, 2, 0, 0);
    }

    if (NO >= 6) {
      A[++m][i] = get_a(scl_h[2], h_re_scl, 5, 0, 0, 0, 0);
      A[++m][i] = get_a(scl_h[2], h_re_scl, 4, 1, 0, 0, 0);
      A[++m][i] = get_a(scl_h[2], h_re_scl, 3, 2, 0, 0, 0);
      A[++m][i] = get_a(scl_h[2], h_re_scl, 3, 0, 2, 0, 0);
      A[++m][i] = get_a(scl_h[2], h_re_scl, 2, 1, 2, 0, 0);
      A[++m][i] = get_a(scl_h[2], h_re_scl, 2, 1, 0, 2, 0);
      A[++m][i] = get_a(scl_h[2], h_re_scl, 3, 0, 0, 2, 0);
      A[++m][i] = get_a(scl_h[2], h_re_scl, 1, 0, 4, 0, 0);
      A[++m][i] = get_a(scl_h[2], h_re_scl, 1, 0, 0, 4, 0);
      A[++m][i] = get_a(scl_h[2], h_re_scl, 3, 0, 1, 1, 0);
      A[++m][i] = get_a(scl_h[2], h_re_scl, 2, 1, 1, 1, 0);
      A[++m][i] = get_a(scl_h[2], h_re_scl, 1, 0, 3, 1, 0);
      A[++m][i] = get_a(scl_h[2], h_re_scl, 1, 0, 1, 3, 0);
      A[++m][i] = get_a(scl_h[2], h_re_scl, 1, 0, 2, 2, 0);
    }

    if (!symm) {
      A[++m][i] = get_a(scl_h[0], h_im_scl, 1, 0, 0, 0, 2);
      A[++m][i] = get_a(scl_h[0], h_im_scl, 2, 0, 0, 0, 1);
      A[++m][i] = get_a(scl_h[0], h_im_scl, 0, 0, 2, 0, 1);

      A[++m][i] = get_a(scl_h[0], h_im_scl, 1, 0, 1, 1, 0);
      A[++m][i] = get_a(scl_h[0], h_im_scl, 2, 1, 0, 0, 0);
      A[++m][i] = get_a(scl_h[0], h_im_scl, 3, 0, 0, 0, 0);
      A[++m][i] = get_a(scl_h[0], h_im_scl, 1, 0, 0, 2, 0);
      A[++m][i] = get_a(scl_h[0], h_im_scl, 1, 0, 2, 0, 0);

      A[++m][i] = get_a(scl_h[1], h_im_scl, 2, 0, 1, 1, 0);
      A[++m][i] = get_a(scl_h[1], h_im_scl, 3, 1, 0, 0, 0);
      A[++m][i] = get_a(scl_h[1], h_im_scl, 4, 0, 0, 0, 0);
      A[++m][i] = get_a(scl_h[1], h_im_scl, 2, 0, 0, 2, 0);
      A[++m][i] = get_a(scl_h[1], h_im_scl, 2, 0, 2, 0, 0);
      A[++m][i] = get_a(scl_h[1], h_im_scl, 0, 0, 4, 0, 0);
      A[++m][i] = get_a(scl_h[1], h_im_scl, 0, 0, 3, 1, 0);
      A[++m][i] = get_a(scl_h[1], h_im_scl, 1, 1, 2, 0, 0);

      if (NO >= 6) {
	A[++m][i] = get_a(scl_h[2], h_im_scl, 5, 0, 0, 0, 0);
	A[++m][i] = get_a(scl_h[2], h_im_scl, 4, 1, 0, 0, 0);
	A[++m][i] = get_a(scl_h[2], h_im_scl, 3, 2, 0, 0, 0);
	A[++m][i] = get_a(scl_h[2], h_im_scl, 3, 0, 2, 0, 0);
	A[++m][i] = get_a(scl_h[2], h_im_scl, 2, 1, 2, 0, 0);
	A[++m][i] = get_a(scl_h[2], h_im_scl, 2, 1, 0, 2, 0);
	A[++m][i] = get_a(scl_h[2], h_im_scl, 3, 0, 0, 2, 0);
	A[++m][i] = get_a(scl_h[2], h_im_scl, 1, 0, 4, 0, 0);
	A[++m][i] = get_a(scl_h[2], h_im_scl, 1, 0, 0, 4, 0);
	A[++m][i] = get_a(scl_h[2], h_im_scl, 3, 0, 1, 1, 0);
	A[++m][i] = get_a(scl_h[2], h_im_scl, 2, 1, 1, 1, 0);
	A[++m][i] = get_a(scl_h[2], h_im_scl, 1, 0, 3, 1, 0);
	A[++m][i] = get_a(scl_h[2], h_im_scl, 1, 0, 1, 3, 0);
	A[++m][i] = get_a(scl_h[2], h_im_scl, 1, 0, 2, 2, 0);
      }
    }

    A[++m][i] = get_a(scl_ksi[0], K_re_scl, 1, 1, 0, 0, 1);
    A[++m][i] = get_a(scl_ksi[0], K_re_scl, 0, 0, 1, 1, 1);

    if (NO >= 5) {
      A[++m][i] = get_a(scl_dnu[0],   K_re_scl, 2, 2, 0, 0, 0);
      A[++m][i] = get_a(scl_dnu[0],   K_re_scl, 0, 0, 2, 2, 0);
      A[++m][i] = get_a(scl_dnu[0],   K_re_scl, 1, 1, 1, 1, 0);

      A[++m][i] = get_a(scl_ksi[1],   K_re_scl, 1, 1, 0, 0, 2);
      A[++m][i] = get_a(scl_ksi[1],   K_re_scl, 0, 0, 1, 1, 2);

      A[++m][i] = get_a(scl_dnu_conf, dnu,      0, 0, 0, 0, 0);
    }

    if (NO >= 6) {
      A[++m][i] = get_a(scl_dnu[1], K_re_scl, 2, 2, 0, 0, 1);
      A[++m][i] = get_a(scl_dnu[1], K_re_scl, 0, 0, 2, 2, 1);
      A[++m][i] = get_a(scl_dnu[1], K_re_scl, 1, 1, 1, 1, 1);

      A[++m][i] = get_a(scl_ksi[2], K_re_scl, 1, 1, 0, 0, 3);
      A[++m][i] = get_a(scl_ksi[2], K_re_scl, 0, 0, 1, 1, 3);
    }

    if (NO >= 7) {
      A[++m][i] = get_a(scl_dnu[2], K_re_scl, 3, 3, 0, 0, 0);
      A[++m][i] = get_a(scl_dnu[2], K_re_scl, 2, 2, 1, 1, 0);
      A[++m][i] = get_a(scl_dnu[2], K_re_scl, 1, 1, 2, 2, 0);
      A[++m][i] = get_a(scl_dnu[2], K_re_scl, 0, 0, 3, 3, 0);
    }

    if (NO >= 9) {
      A[++m][i] = get_a(scl_dnu[3], K_re_scl, 4, 4, 0, 0, 0);
      A[++m][i] = get_a(scl_dnu[3], K_re_scl, 3, 3, 1, 1, 0);
      A[++m][i] = get_a(scl_dnu[3], K_re_scl, 2, 2, 2, 2, 0);
      A[++m][i] = get_a(scl_dnu[3], K_re_scl, 1, 1, 3, 3, 0);
      A[++m][i] = get_a(scl_dnu[3], K_re_scl, 0, 0, 4, 4, 0);
    }

    for (j = 1; j <= m; j++)
      A[j][i] *= bn_prms.bn_scl[i-1];

    bn_prms.clr_prm_dep(i-1);
  }

  printf("\n dnu = %10.3e\n", dnu.cst());

  m = 0;
  f[++m] = get_b(scl_h[0], h_re_scl, 1, 0, 0, 0, 2);
  f[++m] = get_b(scl_h[0], h_re_scl, 2, 0, 0, 0, 1);
  f[++m] = get_b(scl_h[0], h_re_scl, 0, 0, 2, 0, 1);

  f[++m] = get_b(scl_h[0], h_re_scl, 1, 0, 1, 1, 0);
  f[++m] = get_b(scl_h[0], h_re_scl, 2, 1, 0, 0, 0);
  f[++m] = get_b(scl_h[0], h_re_scl, 3, 0, 0, 0, 0);
  f[++m] = get_b(scl_h[0], h_re_scl, 1, 0, 0, 2, 0);
  f[++m] = get_b(scl_h[0], h_re_scl, 1, 0, 2, 0, 0);

  if (NO >= 5) {
    f[++m] = get_b(scl_h[1], h_re_scl, 2, 0, 1, 1, 0);
    f[++m] = get_b(scl_h[1], h_re_scl, 3, 1, 0, 0, 0);
    f[++m] = get_b(scl_h[1], h_re_scl, 4, 0, 0, 0, 0);
    f[++m] = get_b(scl_h[1], h_re_scl, 2, 0, 0, 2, 0);
    f[++m] = get_b(scl_h[1], h_re_scl, 2, 0, 2, 0, 0);
    f[++m] = get_b(scl_h[1], h_re_scl, 0, 0, 4, 0, 0);
    f[++m] = get_b(scl_h[1], h_re_scl, 0, 0, 3, 1, 0);
    f[++m] = get_b(scl_h[1], h_re_scl, 1, 1, 2, 0, 0);
  }

  if (NO >= 6) {
    f[++m] = get_b(scl_h[2], h_re_scl, 5, 0, 0, 0, 0);
    f[++m] = get_b(scl_h[2], h_re_scl, 4, 1, 0, 0, 0);
    f[++m] = get_b(scl_h[2], h_re_scl, 3, 2, 0, 0, 0);
    f[++m] = get_b(scl_h[2], h_re_scl, 3, 0, 2, 0, 0);
    f[++m] = get_b(scl_h[2], h_re_scl, 2, 1, 2, 0, 0);
    f[++m] = get_b(scl_h[2], h_re_scl, 2, 1, 0, 2, 0);
    f[++m] = get_b(scl_h[2], h_re_scl, 3, 0, 0, 2, 0);
    f[++m] = get_b(scl_h[2], h_re_scl, 1, 0, 4, 0, 0);
    f[++m] = get_b(scl_h[2], h_re_scl, 1, 0, 0, 4, 0);
    f[++m] = get_b(scl_h[2], h_re_scl, 3, 0, 1, 1, 0);
    f[++m] = get_b(scl_h[2], h_re_scl, 2, 1, 1, 1, 0);
    f[++m] = get_b(scl_h[2], h_re_scl, 1, 0, 3, 1, 0);
    f[++m] = get_b(scl_h[2], h_re_scl, 1, 0, 1, 3, 0);
    f[++m] = get_b(scl_h[2], h_re_scl, 1, 0, 2, 2, 0);
  }

  if (!symm) {
    f[++m] = get_b(scl_h[0], h_im_scl, 1, 0, 0, 0, 2);
    f[++m] = get_b(scl_h[0], h_im_scl, 2, 0, 0, 0, 1);
    f[++m] = get_b(scl_h[0], h_im_scl, 0, 0, 2, 0, 1);

    f[++m] = get_b(scl_h[0], h_im_scl, 1, 0, 1, 1, 0);
    f[++m] = get_b(scl_h[0], h_im_scl, 2, 1, 0, 0, 0);
    f[++m] = get_b(scl_h[0], h_im_scl, 3, 0, 0, 0, 0);
    f[++m] = get_b(scl_h[0], h_im_scl, 1, 0, 0, 2, 0);
    f[++m] = get_b(scl_h[0], h_im_scl, 1, 0, 2, 0, 0);

    f[++m] = get_b(scl_h[1], h_im_scl, 2, 0, 1, 1, 0);
    f[++m] = get_b(scl_h[1], h_im_scl, 3, 1, 0, 0, 0);
    f[++m] = get_b(scl_h[1], h_im_scl, 4, 0, 0, 0, 0);
    f[++m] = get_b(scl_h[1], h_im_scl, 2, 0, 0, 2, 0);
    f[++m] = get_b(scl_h[1], h_im_scl, 2, 0, 2, 0, 0);
    f[++m] = get_b(scl_h[1], h_im_scl, 0, 0, 4, 0, 0);
    f[++m] = get_b(scl_h[1], h_im_scl, 0, 0, 3, 1, 0);
    f[++m] = get_b(scl_h[1], h_im_scl, 1, 1, 2, 0, 0);

    if (NO >= 6) {
      f[++m] = get_b(scl_h[2], h_im_scl, 5, 0, 0, 0, 0);
      f[++m] = get_b(scl_h[2], h_im_scl, 4, 1, 0, 0, 0);
      f[++m] = get_b(scl_h[2], h_im_scl, 3, 2, 0, 0, 0);
      f[++m] = get_b(scl_h[2], h_im_scl, 3, 0, 2, 0, 0);
      f[++m] = get_b(scl_h[2], h_im_scl, 2, 1, 2, 0, 0);
      f[++m] = get_b(scl_h[2], h_im_scl, 2, 1, 0, 2, 0);
      f[++m] = get_b(scl_h[2], h_im_scl, 3, 0, 0, 2, 0);
      f[++m] = get_b(scl_h[2], h_im_scl, 1, 0, 4, 0, 0);
      f[++m] = get_b(scl_h[2], h_im_scl, 1, 0, 0, 4, 0);
      f[++m] = get_b(scl_h[2], h_im_scl, 3, 0, 1, 1, 0);
      f[++m] = get_b(scl_h[2], h_im_scl, 2, 1, 1, 1, 0);
      f[++m] = get_b(scl_h[2], h_im_scl, 1, 0, 3, 1, 0);
      f[++m] = get_b(scl_h[2], h_im_scl, 1, 0, 1, 3, 0);
      f[++m] = get_b(scl_h[2], h_im_scl, 1, 0, 2, 2, 0);
    }
  }

  f[++m] = get_b(scl_ksi[0], K_re_scl, 1, 1, 0, 0, 1);
  f[++m] = get_b(scl_ksi[0], K_re_scl, 0, 0, 1, 1, 1);

  if (NO >= 5) {
    f[++m] = get_b(scl_dnu[0],   K_re_scl, 2, 2, 0, 0, 0);
    f[++m] = get_b(scl_dnu[0],   K_re_scl, 0, 0, 2, 2, 0);
    f[++m] = get_b(scl_dnu[0],   K_re_scl, 1, 1, 1, 1, 0);

    f[++m] = get_b(scl_ksi[1],   K_re_scl, 1, 1, 0, 0, 2);
    f[++m] = get_b(scl_ksi[1],   K_re_scl, 0, 0, 1, 1, 2);

    f[++m] = get_b(scl_dnu_conf, dnu,      0, 0, 0, 0, 0);
  }

  if (NO >= 6) {
    f[++m] = get_b(scl_dnu[1], K_re_scl, 2, 2, 0, 0, 1);
    f[++m] = get_b(scl_dnu[1], K_re_scl, 0, 0, 2, 2, 1);
    f[++m] = get_b(scl_dnu[1], K_re_scl, 1, 1, 1, 1, 1);

    f[++m] = get_b(scl_ksi[2], K_re_scl, 1, 1, 0, 0, 3);
    f[++m] = get_b(scl_ksi[2], K_re_scl, 0, 0, 1, 1, 3);
  }

  if (NO >= 7) {
    f[++m] = get_b(scl_dnu[2], K_re_scl, 3, 3, 0, 0, 0);
    f[++m] = get_b(scl_dnu[2], K_re_scl, 2, 2, 1, 1, 0);
    f[++m] = get_b(scl_dnu[2], K_re_scl, 1, 1, 2, 2, 0);
    f[++m] = get_b(scl_dnu[2], K_re_scl, 0, 0, 3, 3, 0);
  }

  if (NO >= 9) {
    f[++m] = get_b(scl_dnu[3], K_re_scl, 4, 4, 0, 0, 0);
    f[++m] = get_b(scl_dnu[3], K_re_scl, 3, 3, 1, 1, 0);
    f[++m] = get_b(scl_dnu[3], K_re_scl, 2, 2, 2, 2, 0);
    f[++m] = get_b(scl_dnu[3], K_re_scl, 1, 1, 3, 3, 0);
    f[++m] = get_b(scl_dnu[3], K_re_scl, 0, 0, 4, 4, 0);
  }

  prt_h_K();

  chi2 = 0e0;
  for (j = 1; j <= m; j++)
    chi2 += sqr(f[j]);
}


void min_conj_grad(double &chi2, double &dbn_max, double *g_, double *h_,
		   const bool cg_meth)
{
  int    n_bn, i, m;
  double chi2_ref, **A, *f, *b, *bn_ref;

  const int m_max = 100;

  n_bn = bn_prms.n_prm;

  bn_ref = dvector(1, n_bn); f = dvector(1, m_max); b = dvector(1, m_max);
  A = dmatrix(1, m_max, 1, n_bn);

  chi2_ref = chi2;

  get_f_grad(n_bn, f, A, chi2, m);
  prt_system(m, n_bn, A, f);

  printf("\n%4d chi2: %12.5e -> %12.5e\n", n_iter, chi2_ref, chi2);

  for (i = 1; i <= m; i++)
    b[i] = -f[i];

  SVD_lim(m, n_bn, A, b, bn_prms.bn_lim, bn_prms.svd_n_cut, bn_prms.bn,
	  bn_prms.dbn);

  dvcopy(bn_prms.bn, n_bn, bn_ref);
  if (cg_meth)
    conj_grad(n_iter, bn_prms.bn, bn_prms.dbn, g_, h_, get_f);
  else
    bn_prms.set_dprm();

  printf("\nbn & dbn:\n");
  for (i = 1; i <= n_bn; i++)
    printf(" %12.5e", bn_prms.bn_scl[i-1]*bn_prms.bn[i]);
  printf("\n");
  dbn_max = 0e0;
  for (i = 1; i <= n_bn; i++) {
    dbn_max = max(fabs((bn_prms.bn[i]-bn_ref[i])), dbn_max);
    printf(" %12.5e", bn_prms.bn[i]-bn_ref[i]);
  }
  printf("\n");

  free_dvector(bn_ref, 1, n_bn); free_dvector(f, 1, m_max);
  free_dvector(b, 1, m_max); free_dmatrix(A, 1, m_max, 1, n_bn);
}


void min_conj_grad(const bool cg_meth)
{
  // Control tune foot print; conjugate gradient method.
  std::string str;
  int         n_bn;
  double      *g, *h, dbn_max, chi2;

  const int n_iter_max = 1000;

  n_bn = bn_prms.n_prm;

  g = dvector(1, n_bn); h = dvector(1, n_bn);

  n_iter = 0; chi2 = 0e0;
  do {
    n_iter++;

    min_conj_grad(chi2, dbn_max, g, h, cg_meth);

    prt_mfile("flat_file.fit");
    prt_bn(bn_prms);
  } while ((dbn_max >  bn_prms.bn_tol) && (n_iter < n_iter_max));
}


void min_powell(void)
{
  int    n_bn, i, j, iter;
  double **xi, fret;

  n_bn = bn_prms.n_prm;

  xi = dmatrix(1, n_bn, 1, n_bn);

  // Set initial directions (unit vectors).
  for (i = 1; i <= n_bn; i++)
    for (j = 1; j <= n_bn; j++)
      xi[i][j] = (i == j)? 1e0 : 0e0;

  dpowell(bn_prms.bn, xi, n_bn, bn_prms.bn_tol, &iter, &fret, get_f);

  free_dmatrix(xi, 1, n_bn, 1, n_bn);
}


void prt_lev_marq(const int m, const int n)
{
  int i;

  prt_system(m, n, A_lm, f_lm);

  prt_bn(bn_prms);

  printf("\n%d bn:\n", n_powell);
  for (i = 1; i <= bn_prms.n_prm; i++) {
    bn_prms.bn[i] =
      get_bn(bn_prms.Fnum[i-1], 1, bn_prms.n[i-1]/bn_prms.bn_scl[i-1]);
    printf("%11.3e", bn_prms.bn_scl[i-1]*bn_prms.bn[i]);
    if (i % n_prt == 0) printf("\n");
  }
  if (bn_prms.n_prm % n_prt != 0) printf("\n");
  prt_mfile("flat_file.fit");
}


void get_f_der(double x, double *bn, double *yfit, double *dyda, int n)
{
  int        i, m1;
  static int m;
  double     chi2;

  const bool prt = false;

  m1 = (int)(x+0.5);
  if (prt) printf(" %d", m1);

  if (m1 == 1) {
    n_powell++;
    for (i = 1; i <= n; i++)
      set_bn(bn_prms.Fnum[i-1], bn_prms.n[i-1], bn_prms.bn_scl[i-1]*bn[i]);
    get_f_grad(n, f_lm, A_lm, chi2, m);
  }

  if (prt && (m1 == m)) {
    printf("\n");
    for (i = 1; i <= n; i++) {
      printf(" %12.5e", bn_prms.bn_scl[i-1]*bn[i]);
      if (i % n_prt == 0) printf("\n");
    }
    if (n % n_prt != 0) printf("\n");
  }

  *yfit = f_lm[m1];

  for (i = 1; i <= n; i++)
    dyda[i] = A_lm[m1][i];
}


void min_lev_marq(void)
{
  int    n_data, n_bn, i, n, *ia;
  double *x, *y, *sigma, **covar, **alpha, chisq, alambda, alambda0;

  n_data = 0;
  if (NO >= 3+1) n_data += 3 + 5 + 2;     // 10.
  if (NO >= 4+1) n_data += 8 + 3 + 2 + 1; // 23 + 1.
  if (NO >= 5+1) n_data += 14 + 3 + 2;    // 42.
  if (NO >= 6+1) n_data += 4;             // 46.
  if (NO >= 8+1) n_data += 5;             // 51.
  if (!symm)     n_data += 16 + 14;

  n_bn = bn_prms.n_prm;

  ia = ivector(1, n_bn);
  x = dvector(1, n_data); y = dvector(1, n_data); sigma = dvector(1, n_data);
  covar = dmatrix(1, n_bn, 1, n_bn); alpha = dmatrix(1, n_bn, 1, n_bn);
  f_lm = dvector(1, n_data); A_lm = dmatrix(1, n_data, 1, n_bn);

  get_f_grad(n_bn, f_lm, A_lm, chi2, n_data);
  prt_system(n_data, n_bn, A_lm, f_lm);

  for (i = 1; i <= n_bn; i++)
    ia[i] = 1;

  for (i = 1; i <= n_data; i++) {
    sigma[i] = 1e0; x[i] = i; y[i] = 0e0;
  }

  alambda = -1e0; alambda0 = 1e-3;
  dmrqmin(x, y, sigma, n_data, bn_prms.bn, ia, n_bn, covar, alpha, &chisq,
	  get_f_der, &alambda);
  printf("\nalambda = %7.1e, chi2 = %9.3e\n", alambda, chisq);
  if (alambda < alambda0) prt_lev_marq(n_data, n_bn);
  alambda0 = alambda;

  n = 0;
  do {
    n++;
    dmrqmin(x, y, sigma, n_data, bn_prms.bn, ia, n_bn,  covar, alpha, &chisq,
	    get_f_der, &alambda);
    printf("\nalambda = %7.1e, chi2 = %9.3e\n", alambda, chisq);
    if (alambda < alambda0) prt_lev_marq(n_data, n_bn);
    alambda0 = alambda;
  } while (n < 25);

  alambda = 0e0;
  dmrqmin(x, y, sigma, n_data, bn_prms.bn, ia, n_bn,  covar, alpha, &chisq,
	  get_f_der, &alambda);

  free_dvector(f_lm, 1, n_data); free_dmatrix(A_lm, 1, n_data, 1, n_bn);
  free_ivector(ia, 1, n_bn);
  free_dvector(x, 1, n_data); free_dvector(y, 1, n_data);
  free_dvector(sigma, 1, n_data);
  free_dmatrix(covar, 1, n_bn, 1, n_bn); free_dmatrix(alpha, 1, n_bn, 1, n_bn);
}


void prt_ct(const int n, const double delta)
{
  int             k;
  double          delta1;
  ss_vect<double> ps;
  FILE            *outf;

  outf = file_write("ct.out");

  // cavity_on = true;

  for (k = -n; k <= n; k++) {
    delta1 = (double)k/(double)n*delta;
    ps.zero();
    ps[x_] = 0*2.6e-3; ps[delta_] = delta1;
    ps.propagate(1, n_elem);

    fprintf(outf, "%3d %12.5e %12.5e\n", k, delta1, ps[ct_]);
  }

  fclose(outf);
}


void fit_tune(const double nu_x, const double nu_y,
	      const std::vector<int> &b2_Fam,
	      const double eps, const bool prt)
{
  // Periodic solution: [nu_x, nu_y, beta_x, beta_y]

  int           i, j, n_b2;
  double        **A, *b, *b2_lim, *b2, *db2, step;
  double        nu_fract[2], dnu[2];
  ss_vect<tps>  nus, dnus;
  std::ofstream quad_out;

  const bool    debug = true;
  const int     m     = 2, n_cut = 0;
  const double  step0 = 0.1;

  b2_max = 10e0;

  n_b2 = b2_Fam.size();

  b = dvector(1, m); b2_lim = dvector(1, n_b2);
  b2 = dvector(1, n_b2); db2 = dvector(1, n_b2);
  A = dmatrix(1, m, 1, n_b2);

  nu_fract[X_] = fract(nu_x); nu_fract[Y_] = fract(nu_y);
  printf("\nfit_tune nu = [%7.5f, %7.5f]\n", nu_fract[X_], nu_fract[Y_]);
  printf("\ninitial b2 (%d):\n", n_b2);
  for (i = 1; i <= n_b2; i++) {
    b2_lim[i] = b2_max; b2[i] = get_bn(b2_Fam[i-1], 1, Quad);
    printf(" %9.5f", b2[i]);
  }
  printf("\n");

  danot_(3);

  get_Map(); K = MapNorm(Map, g, A1, A0, Map_res, 1); nus = dHdJ(K);

  dnu[X_] = nus[3].cst() - nu_fract[X_];
  dnu[Y_] = nus[4].cst() - nu_fract[Y_];
  
  printf("\ndnu = [%7.5f, %7.5f]\n", dnu[X_], dnu[Y_]);
  while ((fabs(dnu[X_]) > eps) || (fabs(dnu[Y_]) > eps)) {
    step = step0;
    for (i = 1; i <= n_b2; i++) {
      for (j = 1; j <= get_n_Kids(b2_Fam[i-1]); j++)
	set_bn_par(b2_Fam[i-1], j, Quad, 7);

      get_Map(); K = MapNorm(Map, g, A1, A0, Map_res, 1); nus = dHdJ(K);

      A[1][i] = h_ijklm_p(nus[3], 0, 0, 0, 0, 0, 7);
      A[2][i] = h_ijklm_p(nus[4], 0, 0, 0, 0, 0, 7);

      for (j = 1; j <= get_n_Kids(b2_Fam[i-1]); j++)
	clr_bn_par(b2_Fam[i-1], j, Quad);
    }

    b[1] = -dnu[X_]; b[2] = -dnu[Y_];

    SVD_lim(m, n_b2, A, b, b2_lim, n_cut, b2, db2);

    for (i = 1; i <= n_b2; i++) {
      set_dbn(b2_Fam[i-1], Quad, step*db2[i]);
      b2[i] = get_bn(b2_Fam[i-1], 1, Quad);
    }

    if (debug) {
      printf("\n Ax = b:\n");
      for (i = 1; i <= m; i++) {
	for (j = 1; j <= n_b2; j++)
	  printf("%11.3e", A[i][j]);
        printf("%11.3e\n", b[i]);
      }
    }

    // Evaluate if stable.
    get_Map(); K = MapNorm(Map, g, A1, A0, Map_res, 1); nus = dHdJ(K);
    while (!stable) {
      // Roll back.
      for (i = 1; i <= n_b2; i++) {
	set_dbn(b2_Fam[i-1], Quad, -step*db2[i]);
	b2[i] = get_bn_s(b2_Fam[i-1], 1, Quad);
      }

      step /= 2.0;
      printf("\nstep = %5.3f\n", step);

      for (i = 1; i <= n_b2; i++) {
	set_dbn_s(b2_Fam[i-1], Quad, step*db2[i]);
	b2[i] = get_bn_s(b2_Fam[i-1], 1, Quad);
      }
	
      get_Map(); K = MapNorm(Map, g, A1, A0, Map_res, 1); nus = dHdJ(K);
    }
      
    dnu[X_] = nus[3].cst() - nu_fract[X_];
    dnu[Y_] = nus[4].cst() - nu_fract[Y_];
    
    if (debug) {
      printf("\ndnu = [%8.5f, %8.5f]\n", dnu[X_], dnu[Y_]);
    }
    printf("nu  = [%8.5f, %8.5f]\n", nus[3].cst(), nus[4].cst());
  }

  if (prt) {
    quad_out.open("fit_tune.dat", std::ios::out);
    quad_out << "\nn = 1:" << "\n";
    for (i = 1; i <= n_b2; i++)
      for (j = 1; j <= get_n_Kids(b2_Fam[i-1]); j++)
	if (b2_Fam[i-1] > 0)
	  quad_out << std::fixed << std::setprecision(7) 
		   << std::setw(6) << get_Name(b2_Fam[i-1]) << "(" << j
		   << ") = " << std::setw(11) << get_bn(b2_Fam[i-1], j, Quad)
		   << std::setw(2) << Quad << "\n";
    quad_out.close();
  }

  free_dvector(b, 1, m);
  free_dvector(b2_lim, 1, n_b2); free_dvector(b2, 1, n_b2);
  free_dvector(db2, 1, n_b2);
  free_dmatrix(A, 1, m, 1, n_b2);
}


void fit_tune(const double nu_x, const double nu_y)
{
  std::vector<int> b2_Fnum;

  const double eps = 1e-5;

  // DIAMOND-II 8-BA by Hossein.
  // Zero dispersion.
  b2_Fnum.push_back(get_Fnum("q1"));
  b2_Fnum.push_back(get_Fnum("q2"));
  b2_Fnum.push_back(get_Fnum("q3"));

  // b2_Fnum.push_back(get_Fnum("q4"));
  // b2_Fnum.push_back(get_Fnum("q5"));
  // b2_Fnum.push_back(get_Fnum("q6"));
  // b2_Fnum.push_back(get_Fnum("q7"));
  // b2_Fnum.push_back(get_Fnum("q8"));
  // b2_Fnum.push_back(get_Fnum("q9"));
  // b2_Fnum.push_back(get_Fnum("q10"));
  // b2_Fnum.push_back(get_Fnum("qu1"));
  // b2_Fnum.push_back(get_Fnum("qu2"));
  // b2_Fnum.push_back(get_Fnum("qu3"));
  // b2_Fnum.push_back(get_Fnum("qu4"));
  // b2_Fnum.push_back(get_Fnum("qu5"));
  // b2_Fnum.push_back(get_Fnum("qu6"));
  // b2_Fnum.push_back(get_Fnum("qu7"));
  // b2_Fnum.push_back(get_Fnum("qu8"));
  // b2_Fnum.push_back(get_Fnum("qu9"));
 
  fit_tune(nu_x, nu_y, b2_Fnum, eps, true);

  prt_b2(b2_Fnum);
}


int main(int argc, char *argv[])
{
  int j;

    rad_on    = false; H_exact        = false; totpath_on   = false;
    cavity_on = false; quad_fringe_on = false; emittance_on = false;
    IBS_on    = false;

    rd_mfile(argv[1], elem); rd_mfile(argv[1], elem_tps);

    // Initialize the symplectic integrator after energy has been defined.
    ini_si();

    // Disable log messages from TPSALib and LieLib.
#if !DOF_3
    idprset(-1);
#else
    idprset(1);
    cavity_on = true;
#endif

    // printf("\nDBL_EPSILON = %9.2e\n", DBL_EPSILON);
    // daeps_(DBL_EPSILON);
    daeps_(tpsa_eps);

    danot_(1);

    printf("\nscl_h:          %7.1e, %7.1e, %7.1e\n",
	   scl_h[0], scl_h[1], scl_h[2]);
    printf("scl_dnu:        %7.1e, %7.1e, %7.1e, %7.1e\n",
	   scl_dnu[0], scl_dnu[1], scl_dnu[2], scl_dnu[3]);
    printf("scl_ksi:        %7.1e, %7.1e, %7.1e\n",
	   scl_ksi[0], scl_ksi[1], scl_ksi[2]);
    printf("scl_dnu_conf:   %7.1e\n", scl_dnu_conf);
    printf("n_cut:          %d\n", n_cut);
    printf("symmetric:      %d\n", symm);
    printf("\nA_max:          %7.1e, %7.1e\n",
	   A_max[lat_case-1][X_], A_max[lat_case-1][Y_]);
    printf("delta_max:      %7.1e\n", delta_max[lat_case-1]);
    printf("beta_inj:       %7.1e, %7.1e\n",
	   beta_inj[lat_case-1][X_], beta_inj[lat_case-1][Y_]);

    get_nu_ksi();

    if (false) {
      // fit_tune(57.15/6.0, 22.25/6.0);
      fit_tune(0.530831725+1e-4, 0.685735574-0*1e-4);
      get_nu_ksi();
      exit(0);
    }

    for (j = 0; j < 2; j++)
      twoJ[j] =	sqr(A_max[lat_case-1][j])/beta_inj[lat_case-1][j];

    Id_scl.identity();
    Id_scl[x_] *= sqrt(twoJ[X_]); Id_scl[px_] *= sqrt(twoJ[X_]);
    Id_scl[y_] *= sqrt(twoJ[Y_]); Id_scl[py_] *= sqrt(twoJ[Y_]);
    Id_scl[delta_] *= delta_max[lat_case-1];

    if (false) {
      danot_(NO-1);
      cavity_on = true; rad_on = true;
      get_map_n(n_cell);
      // MAX-VI:
      // prt_H_long(10, M_PI, 10e-2, -405.6e3, false);
      // SLS-2:
      prt_H_long(10, M_PI, 10e-2, -544.7e3, true);
      prt_alphac();
      exit(0);
    }

    if (false) {
      if (false) {
	no_mpoles(Sext);
	no_mpoles(Oct);
      }

      danot_(NO-1);
      get_map_n(n_cell);
      prt_alphac();

      prt_ct(100, 8e-2);
      exit(0);
    }

    if (false) {
      prt_map();
      exit(0);
    }

    if (false) {
      prt_h_K();
      exit(0);
    }

    switch (lat_case) {
    case 1:
      // MAX VI:
      bn_prms.add_prm("o1", 4, 5e5, 1.0);
      bn_prms.add_prm("o2", 4, 5e5, 1.0);
      bn_prms.add_prm("o3", 4, 5e5, 1.0);
      bn_prms.add_prm("o4", 4, 5e5, 1.0);

      // bn_prms.add_prm("o1", 6, 5e10, 1.0);
      // bn_prms.add_prm("o2", 6, 5e10, 1.0);
      // bn_prms.add_prm("o3", 6, 5e10, 1.0);
      // bn_prms.add_prm("o4", 6, 5e10, 1.0);
      break;
    case 2:
      // SLS-2:
      if (true) {
	bn_prms.add_prm("sfh",  3, 5e5, 1.0);
	bn_prms.add_prm("sdh",  3, 5e5, 1.0);
	bn_prms.add_prm("sfmh", 3, 5e5, 1.0);
	bn_prms.add_prm("sdmh", 3, 5e5, 1.0);

	bn_prms.add_prm("sxxh", 3, 5e5, 1.0);
	bn_prms.add_prm("sxyh", 3, 5e5, 1.0);
	bn_prms.add_prm("syyh", 3, 5e5, 1.0);
      } else {
	bn_prms.add_prm("ocx",  4, 5e10, 1.0);
	bn_prms.add_prm("ocxm", 4, 5e10, 1.0);
	bn_prms.add_prm("ocy",  4, 5e10, 1.0);
	bn_prms.add_prm("ocym", 4, 5e10, 1.0);

	bn_prms.add_prm("oxx",  4, 5e10, 1.0);
	bn_prms.add_prm("oxy",  4, 5e10, 1.0);
	bn_prms.add_prm("oyy",  4, 5e10, 1.0);
      }
      break;
    case 3 ... 4:
      // DIAMOND:
      if (true) {
	bn_prms.add_prm("ts1a",  3, 5e5, 1.0);
	bn_prms.add_prm("ts1ab", 3, 5e5, 1.0);
	bn_prms.add_prm("ts2a",  3, 5e5, 1.0);
	bn_prms.add_prm("ts2ab", 3, 5e5, 1.0);
	bn_prms.add_prm("ts1b",  3, 5e5, 1.0);
	bn_prms.add_prm("ts2b",  3, 5e5, 1.0);
	bn_prms.add_prm("ts1c",  3, 5e5, 1.0);
	bn_prms.add_prm("ts2c",  3, 5e5, 1.0);
	// bn_prms.add_prm("ts1d",  3, 5e5, 1.0);
	bn_prms.add_prm("ts2d",  3, 5e5, 1.0);
	bn_prms.add_prm("ts1e",  3, 5e5, 1.0);
	bn_prms.add_prm("ts2e",  3, 5e5, 1.0);

	if (lat_case == 4) {
	  // VMX.
	  bn_prms.add_prm("s1", 3, 5e5, 1.0);
	  bn_prms.add_prm("s2", 3, 5e5, 1.0);
	  bn_prms.add_prm("s3", 3, 5e5, 1.0);
	  bn_prms.add_prm("s4", 3, 5e5, 1.0);
	  bn_prms.add_prm("s5", 3, 5e5, 1.0);
	}
      } else {
	bn_prms.add_prm("ts1a",  3, 5e5, 1.0);
	bn_prms.add_prm("ts2a",  3, 5e5, 1.0);
	bn_prms.add_prm("ts1b",  3, 5e5, 1.0);
	bn_prms.add_prm("ts2b",  3, 5e5, 1.0);
	bn_prms.add_prm("ts1c",  3, 5e5, 1.0);
	bn_prms.add_prm("ts2c",  3, 5e5, 1.0);

	// bn_prms.add_prm("s1",    3, 5e5, 1.0);
	// bn_prms.add_prm("s2",    3, 5e5, 1.0);
	// bn_prms.add_prm("s3",    3, 5e5, 1.0);
	// bn_prms.add_prm("s4",    3, 5e5, 1.0);
	// bn_prms.add_prm("s5",    3, 5e5, 1.0);
      }
      break;
    case 5:
      // DIAMOND-II, 4-BA:
      bn_prms.add_prm("s1b", 3, 5e5, 1.0);
      bn_prms.add_prm("s1d", 3, 5e5, 1.0);
      bn_prms.add_prm("s2b", 3, 5e5, 1.0);
      bn_prms.add_prm("s2d", 3, 5e5, 1.0);
      bn_prms.add_prm("sx1", 3, 5e5, 1.0);
      bn_prms.add_prm("sy1", 3, 5e5, 1.0);

      bn_prms.add_prm("s3",  3, 5e5, 1.0);
      break;
    case 6:
      // DIAMOND-II, 6-HMBA:
      // bn_prms.add_prm("sd1",  3, 5e5, 1.0);
      // bn_prms.add_prm("sd2",  3, 5e5, 1.0);
      // bn_prms.add_prm("sd3",  3, 5e5, 1.0);
      // bn_prms.add_prm("sf21", 3, 5e5, 1.0);
      // bn_prms.add_prm("sd31", 3, 5e5, 1.0);
      // bn_prms.add_prm("sf1",  3, 5e5, 1.0);
      // bn_prms.add_prm("sh1a", 3, 5e5, 1.0);
      // bn_prms.add_prm("sh1e", 3, 5e5, 1.0);
      // break;
      bn_prms.add_prm("sfa",  3, 5e5, 1.0);
      bn_prms.add_prm("sfb",  3, 5e5, 1.0);
      bn_prms.add_prm("sda",  3, 5e5, 1.0);
      bn_prms.add_prm("sdb",  3, 5e5, 1.0);

      if (!fit_ksi) {
	bn_prms.add_prm("s1a", 3, 5e5, 1.0);
	bn_prms.add_prm("s1b", 3, 5e5, 1.0);
	bn_prms.add_prm("s2a", 3, 5e5, 1.0);
	bn_prms.add_prm("s2b", 3, 5e5, 1.0);
	bn_prms.add_prm("s3",  3, 5e5, 1.0);
	bn_prms.add_prm("s4",  3, 5e5, 1.0);
	bn_prms.add_prm("s5",  3, 5e5, 1.0);
      }
      break;
    case 7:
      // DIAMOND-II, 6-RB-BA:
      bn_prms.add_prm("sd",   3, 5e5, 1.0);
      bn_prms.add_prm("sfm",  3, 5e5, 1.0);
      bn_prms.add_prm("sdm",  3, 5e5, 1.0);
      bn_prms.add_prm("sxx",  3, 5e5, 1.0);
      bn_prms.add_prm("sxya", 3, 5e5, 1.0);
      bn_prms.add_prm("sxyb", 3, 5e5, 1.0);
      bn_prms.add_prm("syy",  3, 5e5, 1.0);
      break;
    case 8:
      // DIAMOND-II, 8-BA:
      bn_prms.add_prm("sfh",  3, 5e5, 1.0);
      bn_prms.add_prm("sdh",  3, 5e5, 1.0);
      bn_prms.add_prm("sfmh", 3, 5e5, 1.0);
      bn_prms.add_prm("sdmh", 3, 5e5, 1.0);
      if (true) {
	bn_prms.add_prm("sxxh", 3, 5e5, 1.0);
	bn_prms.add_prm("sxyh", 3, 5e5, 1.0);
	bn_prms.add_prm("syyh", 3, 5e5, 1.0);
      }
      break;
    case 9:
      // DIAMOND-II, 8-BA by Hossein:
      if (!oct) {
	bn_prms.add_prm("s1",  3, 5e5, 1.0);
	bn_prms.add_prm("s2",  3, 5e5, 1.0);
	bn_prms.add_prm("s3",  3, 5e5, 1.0);
	bn_prms.add_prm("s4",  3, 5e5, 1.0);
	bn_prms.add_prm("s5",  3, 5e5, 1.0);
	bn_prms.add_prm("s6",  3, 5e5, 1.0);
	bn_prms.add_prm("s7",  3, 5e5, 1.0);
	bn_prms.add_prm("s8",  3, 5e5, 1.0);
	bn_prms.add_prm("s9",  3, 5e5, 1.0);
	bn_prms.add_prm("s10", 3, 5e5, 1.0);
      } else {
	bn_prms.add_prm("s1",  4, 5e5, 1.0);
	bn_prms.add_prm("s2",  4, 5e5, 1.0);
	bn_prms.add_prm("s3",  4, 5e5, 1.0);
	bn_prms.add_prm("s4",  4, 5e5, 1.0);
	bn_prms.add_prm("s5",  4, 5e5, 1.0);
	bn_prms.add_prm("s6",  4, 5e5, 1.0);
	bn_prms.add_prm("s7",  4, 5e5, 1.0);
	bn_prms.add_prm("s8",  4, 5e5, 1.0);
	bn_prms.add_prm("s9",  4, 5e5, 1.0);
	bn_prms.add_prm("s10", 4, 5e5, 1.0);
      }
      break;
    case 10:
      // SLS-2:
      if (!oct) {
	bn_prms.add_prm("sdmh", 3, 5e5, 1.0);
	bn_prms.add_prm("sfmh", 3, 5e5, 1.0);
	bn_prms.add_prm("sdh",  3, 5e5, 1.0);
	bn_prms.add_prm("sfh",  3, 5e5, 1.0);
	if (!fit_ksi) {
	  bn_prms.add_prm("sxxh", 3, 5e5, 1.0);
	  bn_prms.add_prm("sxyh", 3, 5e5, 1.0);
	  bn_prms.add_prm("syyh", 3, 5e5, 1.0);
	}
      } else {
	bn_prms.add_prm("oxx",  4, 5e5, 1.0);
	bn_prms.add_prm("oxy",  4, 5e5, 1.0);
	bn_prms.add_prm("oyy",  4, 5e5, 1.0);
	bn_prms.add_prm("ocxm", 4, 5e5, 1.0);
	bn_prms.add_prm("ocx1", 4, 5e5, 1.0);
	bn_prms.add_prm("ocx2", 4, 5e5, 1.0);
       }
      break;
    }

    // Step is 1.0 for conjugated gradient method.
    bn_prms.bn_tol = 1e-1; bn_prms.svd_n_cut = 0; bn_prms.step = 1.0;

    if (fit_ksi) {
      no_mpoles(Sext); no_mpoles(Oct); no_mpoles(Dodec);
    }

    bn_prms.ini_prm();

    prt_bn(bn_prms);

    if (fit_ksi) {
      bn_prms.svd_n_cut = 0;
      fit_ksi1(0e0, 0e0);
      exit(0);
    }

    no_mpoles(Oct); no_mpoles(Dodec);

    if (true) {
      bn_prms.svd_n_cut = n_cut;
      min_conj_grad(true);
    } else
      min_lev_marq();

    prt_h_K();
}
