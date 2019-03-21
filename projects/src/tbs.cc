#include <cfloat>

#define NO 7

#include "thor_lib.h"

int no_tps = NO,

#define DOF_3 0

#if !DOF_3
  ndpt_tps = 5;
#else
  // The cavity must be turned on.
  ndpt_tps = 0;
#endif

#define DNU 0


extern tps          K, g;
extern ss_vect<tps> Map, A0, A1, Map_res;

const double tpsa_eps = 1e-30;

std::vector<std::string> drv_term;

int          n_iter;
double       chi2 = 1e30;
tps          h_re, h_im, h_re_scl, h_im_scl, K_re, K_im, K_re_scl;
tps          K_re_delta_scl;
ss_vect<tps> nus, nus_scl, Id_scl, Id_delta_scl;

const bool scale = false;

const int n_prt  = 8;

// Center of straight.
const double
  beta_inj[]   = {8.7, 2.1},
  A_max[]      = {6e-3, 2e-3},
  twoJ[]       = {sqr(A_max[X_])/beta_inj[X_], sqr(A_max[Y_])/beta_inj[Y_]},
  twoJ_delta[] = {sqr(0.5e-3)/beta_inj[X_], sqr(0.1e-3)/beta_inj[Y_]},
  delta_max    = 3e-2;

const double
  scl_h[]            = {0e0, 0e0, 0e0},
  scl_dnu[]          = {1e0, 1e0, 1e0, 0e0, 0e0},
  scl_ksi[]          = {0e0, 1e8, 0e0, 0e0, 0e0, 0e0}, // 1st not used.
  delta_scl          = 0e0,
  scl_dnu_conf       = 1e5,
  scl_dnu_delta_conf = 0e0;


double bn_internal(const double bn_bounded,
		   const double bn_min, const double bn_max);
double bn_bounded(const double bn_internal,
		  const double bn_min, const double bn_max);

struct param_type {
private:

public:
  int                 n_bn;
  std::vector<double> bn_min, bn_max, bn_scl;
  std::vector<int>    Fnum, n;

  void add_prm(const std::string Fname, const int n,
	       const double bn_min, const double bn_max, const double bn_scl);
  void ini_prm(double *bn);
  void set_prm_dep(const int k) const;
  void clr_prm_dep(const int k) const;
  void set_prm(double *bn) const;
  void set_dparam(const int k, const double eps) const;
  void prt_bn(double *bn) const;
  void prt_bn_lat(void) const;
};


param_type bn_prms;


void param_type::add_prm(const std::string Fname, const int n,
			 const double bn_min, const double bn_max,
			 const double bn_scl)
{
  Fnum.push_back(get_Fnum(Fname.c_str()));
  this->n.push_back(n);
  this->bn_min.push_back(bn_min);
  this->bn_max.push_back(bn_max);
  this->bn_scl.push_back(bn_scl);
  n_bn = Fnum.size();
}


void param_type::ini_prm(double *bn)
{
  int i, loc;

  printf("\nInitial bn (scale factor in parenthesis):\n");
  printf("  No of Families: %1d\n", n_bn);
  for (i = 1; i <= n_bn; i++) {
    bn[i] = get_bn(Fnum[i-1], 1, n[i-1]);

    // Bounded.
    if ((bn_min[i-1] <= bn[i]) && (bn[i] <= bn_max[i-1]))
      bn[i] = bn_internal(bn[i], bn_min[i-1], bn_max[i-1]);
    else {
      loc = get_loc(Fnum[i-1], 1);
      printf("\nini_prm:\n outside range ");
      printf(" %s %10.3e [%10.3e, %10.3e]\n",
	     elem[loc].Name, bn[i], bn_min[i-1], bn_max[i-1]);
      exit(1);
    }
  }

  prt_bn(bn);
}


void param_type::set_prm_dep(const int k) const
{
  int j;

  for (j = 1; j <= get_n_Kids(Fnum[k]); j++)
    set_bn_par(Fnum[k], j, n[k], 7);
}


void param_type::clr_prm_dep(const int k) const
{
  int j;

  for (j = 1; j <= get_n_Kids(Fnum[k]); j++)
    clr_bn_par(Fnum[k], j, n[k]);
}


void param_type::set_prm(double *bn) const
{
  int    i;
  double bn_ext;

  const bool prt = false;

  for (i = 1; i <= n_bn; i++) {
    // Bounded.
    bn_ext = bn_bounded(bn[i], bn_min[i-1], bn_max[i-1]);
    set_bn(Fnum[i-1], n[i-1], bn_ext);
  }

  if (prt) {
    printf("set_prm:\n");
    prt_bn(bn);
  }
}


void param_type::set_dparam(const int k, double eps) const
{
  const bool prt = false;

  if (prt) printf("set_dparam: %12.5e\n", eps);
  set_dbn(Fnum[k-1], n[k-1], eps);
}


void param_type::prt_bn(double *bn) const
{
  int    i;
  double bn_ext;

  for (i = 1; i <= n_bn; i++) {
    // Bounded.
    bn_ext = bn_bounded(bn[i], bn_min[i-1], bn_max[i-1]);
    printf(" %12.5e", bn_ext);
  }
  printf("\n");
}


void param_type::prt_bn_lat(void) const
{
  bool     first = true;
  long int loc;
  int      j, k, n_ord;
  double   bn;
  FILE     *outf;

  const std::string file_name = "dnu.out";

  outf = file_write(file_name.c_str());

  fprintf(outf, "\n");
  for (k = 0; k < n_bn; k++) {
    loc = get_loc(Fnum[k], 1) - 1;
    bn = get_bn(Fnum[k], 1, n[k]);
    if (n[k] == Sext)
      fprintf(outf,
	      "%-8s: sextupole, l = %7.5f"
	      ", k = %12.5e, n = nsext, Method = Meth;\n",
	      elem[loc].Name, elem[loc].L, bn);
    else {
      n_ord = elem[loc].mpole->order;
      for (j = Sext; j <= n_ord; j++) {
	bn = get_bn(Fnum[k], 1, j);
	if (first && (bn != 0e0)) {
	  first = false;
	  fprintf(outf,
		  "%-8s: multipole, l = %7.5f,"
		  "\n          hom = (%1d, %12.5e, 0e0",
		  elem[loc].Name, elem[loc].L, j, bn);
	} else if (bn != 0e0)
	  fprintf(outf,
		  ",\n"
		  "                 %1d, %12.5e, 0e0", j, bn);
	if (j == n_ord) {
	  fprintf(outf,
		  "),"
		  "\n          n = nsext, Method = Meth;\n");
	  first = true;
	}
      }
    }
  }

  fclose(outf);
}


double bn_internal(const double bn_bounded,
		   const double bn_min, const double bn_max)
{
  return asin((2e0*(bn_bounded-bn_min))/(bn_max-bn_min)-1e0);
}


double bn_bounded(const double bn_internal,
		  const double bn_min, const double bn_max)
{
  return bn_min + (sin(bn_internal)+1e0)*(bn_max-bn_min)/2e0;
}


// ---> Tune confinement.

static double xsav_2D;
static tps (*func_save_2D)(const double, const double);

tps gauss_quad_2D(tps (*func)(const double), const double a, const double b)
{
  int    j;
  double xr, xm, dx;
  tps    s;

  static double
    x[] =
    {0.0, 0.1488743389, 0.4333953941, 0.6794095682, 0.8650633666, 0.9739065285},
    w[] =
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

double gauss_quad_2D_y0(const double x) { return 0e0; }

double gauss_quad_2D_y1(const double x) { return twoJ[Y_]; }

tps gauss_quad_2D_fy(const double y)  { return (*func_save_2D)(xsav_2D, y); }

tps gauss_quad_2D_fx(const double x)
{
  xsav_2D = x;
  return gauss_quad_2D(gauss_quad_2D_fy, gauss_quad_2D_y0(x),
		       gauss_quad_2D_y1(x));
}

tps gauss_quad_2D(tps (*func)(const double, const double),
		  const double x0, const double x1)
{
  func_save_2D = func;
  return gauss_quad_2D(gauss_quad_2D_fx, x0, x1);
}

tps f_gauss_quad_2D(double x, double y)
{
  int          k;
  long int     jj[ss_dim];
  tps          dK, dnu[2], K_re_no_m_one, K_re_no;
  ss_vect<tps> ps, Id;

  // printf("f_gauss_quad_2D\n");
  ps.identity();
  ps[x_] = ps[px_] = sqrt(x);
  ps[y_] = ps[py_] = sqrt(y);
  ps[delta_] = 0e0;

#if DNU
  for (k = 0; k < 2; k++) {
    dnu[k] = (nus[k+3]-nus[k+3].cst())*ps;
    // Compute absolute value.
    if (dnu[k].cst() < 0e0) dnu[k] = -dnu[k];
  }

  return dnu[X_]*dnu[Y_]/(twoJ[X_]*twoJ[Y_]);
#else
  Id.identity(); Id[6] = 0e0;
  danot_(NO-1);
  K_re_no_m_one = K_re*Id;
  danot_(NO);
  K_re_no = K_re*Id - K_re_no_m_one;
  dK = K_re - K_re_no;
  Id.identity(); Id[delta_] = 0e0;
  dK = dK*Id;
  for (k = 0; k < ss_dim; k++)
    jj[k] = 0;
  jj[x_] = 1; jj[px_] = 1;
  dK.pook(jj, 0e0);
  jj[x_] = 0; jj[px_] = 0; jj[y_] = 1; jj[py_] = 1;
  dK.pook(jj, 0e0);
  // std::cout << std::scientific << std::setprecision(3)
  // 	    << std::setw(11) << dK << "\n";
  dK = dK*ps;
  // Compute absolute value.
  if (dK.cst() < 0e0) dK = -dK;
  // std::cout << std::scientific << std::setprecision(3)
  // 	      << "\n |dK| = " << dK << "\n";

  return dK/(twoJ[X_]*twoJ[Y_]);
#endif
}

static double xsav_3D, ysav_3D;
static tps (*func_save_3D)(const double, const double, const double);

tps gauss_quad_3D(tps (*func)(const double), const double a, const double b)
{
  int    j;
  double xr, xm, dx;
  tps    s;

  static double
    x[] =
    {0.0, 0.1488743389, 0.4333953941, 0.6794095682, 0.8650633666, 0.9739065285},
    w[] =
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

double gauss_quad_3D_y0(const double x) { return 0e0; }

double gauss_quad_3D_y1(const double x) { return twoJ[Y_]; }

double gauss_quad_3D_z0(const double x, const double y)
{
  return -delta_max;
}

double gauss_quad_3D_z1(const double x, const double y)
{
  return delta_max;
}

tps gauss_quad_3D_fz(const double z)
{
  return (*func_save_3D)(xsav_3D, ysav_3D, z);
}

tps gauss_quad_3D_fy(const double y)
{
  ysav_3D = y;
  return gauss_quad_3D(gauss_quad_3D_fz,
		       gauss_quad_3D_z0(xsav_3D, y),
		       gauss_quad_3D_z1(xsav_3D, y));
}

tps gauss_quad_3D_fx(const double x)
{
  xsav_3D = x;
  return gauss_quad_3D(gauss_quad_3D_fy, gauss_quad_3D_y0(x),
		       gauss_quad_3D_y1(x));
}

tps gauss_quad_3D(tps (*func)(const double, const double, const double),
		  const double x0, const double x1)
{
  func_save_3D = func;
  return gauss_quad_3D(gauss_quad_3D_fx, x0, x1);
}

tps f_gauss_quad_3D(double x, double y, double z)
{
  int          k, jj[ss_dim];
  tps          dK, dnu[2], K_re_no_m_one, K_re_no, K_re_no_J;
  ss_vect<tps> Id, ps;

  // printf("f_gauss_quad_3D\n");

  ps.identity();
  ps[x_] = ps[px_] = sqrt(x);
  ps[y_] = ps[py_] = sqrt(y);
  ps[delta_] = z;

#if DNU
  for (k = 0; k < 2; k++) {
    dnu[k] = (nus[k+3]-nus[k+3].cst())*ps;
    // Compute absolute value.
    if (dnu[k].cst() < 0e0) dnu[k] = -dnu[k];
  }

  return dnu[X_]*dnu[Y_]/(twoJ[X_]*twoJ[Y_]*2e0*delta_max);
#else
  Id.identity(); Id[6] = 0e0;
  danot_(NO-1);
  K_re_no_m_one = K_re*Id;
  danot_(NO);
  K_re_no = K_re*Id - K_re_no_m_one;
  dK = K_re - K_re_no;
  Id.identity(); Id[delta_] = 0e0;
  dK = dK - dK*Id;
  Id.identity(); Id[x_] = Id[px_] = Id[y_] = Id[py_] = 0e0;
  dK = dK - dK*Id;
  // std::cout << std::scientific << std::setprecision(3)
  // 	    << "\n |dK| = " << dK << "\n";
  dK = dK*ps;
  // Compute absolute value.
  if (dK.cst() < 0e0) dK = -dK;
  // std::cout << std::scientific << std::setprecision(3)
  // 	    << "\n |dK| = " << dK << "\n";

  return dK/(twoJ[X_]*twoJ[Y_]*2e0*delta_max);
#endif
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


void prt_name(FILE *outf, const char *name, const std::string &str,
	      const int len)
{
  int j, k;

  j = 0;
  do {
    fprintf(outf, "%c", name[j]);
    j++;
  } while (name[j] != ' ');
  fprintf(outf, str.c_str());
  for (k = j; k < len; k++)
    fprintf(outf, " ");
}


double get_a(const double scale, const tps &t,
	     const int i, const int j, const int k, const int l, const int m)
{
  return scale*(h_ijklm_p(t, i, j, k, l, m, 7));
}


double get_b(const std::string &label, const double scale, const tps &t,
	     const int i, const int j, const int k, const int l, const int m)
{
  drv_term.push_back(label);
  return scale*(h_ijklm(t, i, j, k, l, m));
}


double get_chi2(void)
{
  bool                prt_ln = false;
  int                 k;
  double              chi2_1;
  tps                 dnu, dnu_delta;
  std::vector<double> dnu_vec;

  const bool chrom = false;

  danot_(NO-1);
  get_Map();
  danot_(NO);
  K = MapNorm(Map, g, A1, A0, Map_res, 1);
  nus = dHdJ(K); nus_scl = nus*Id_scl;
  CtoR(K, K_re, K_im);
  K_re_scl = K_re*Id_scl; K_re_delta_scl = K_re*Id_delta_scl;
  CtoR(get_h(), h_re, h_im); h_re_scl = h_re*Id_scl; h_im_scl = h_im*Id_scl;

  dnu = gauss_quad_2D(f_gauss_quad_2D, 0e0, twoJ[X_]);
  dnu_delta = gauss_quad_3D(f_gauss_quad_3D, 0e0, twoJ[X_]);

  dnu_vec.push_back(h_ijklm(K_re_scl, 1, 1, 0, 0, 1));
  dnu_vec.push_back(h_ijklm(K_re_scl, 0, 0, 1, 1, 1));

  if (!true) {
    dnu_vec.push_back(h_ijklm(K_re_scl, 2, 2, 0, 0, 0));
    dnu_vec.push_back(h_ijklm(K_re_scl, 1, 1, 1, 1, 0));
    dnu_vec.push_back(h_ijklm(K_re_scl, 0, 0, 2, 2, 0));

    dnu_vec.push_back(h_ijklm(K_re_scl, 3, 3, 0, 0, 0));
    dnu_vec.push_back(h_ijklm(K_re_scl, 2, 2, 1, 1, 0));
    dnu_vec.push_back(h_ijklm(K_re_scl, 1, 1, 2, 2, 0));
    dnu_vec.push_back(h_ijklm(K_re_scl, 0, 0, 3, 3, 0));

    dnu_vec.push_back(h_ijklm(K_re_scl, 4, 4, 0, 0, 0));
    dnu_vec.push_back(h_ijklm(K_re_scl, 3, 3, 1, 1, 0));
    dnu_vec.push_back(h_ijklm(K_re_scl, 2, 2, 2, 2, 0));
    dnu_vec.push_back(h_ijklm(K_re_scl, 1, 1, 3, 3, 0));
    dnu_vec.push_back(h_ijklm(K_re_scl, 0, 0, 4, 4, 0));
  }

  if (chrom) {
    dnu_vec.push_back(h_ijklm(K_re_scl, 1, 1, 0, 0, 2));
    dnu_vec.push_back(h_ijklm(K_re_scl, 0, 0, 1, 1, 2));

    dnu_vec.push_back(h_ijklm(K_re_scl, 1, 1, 0, 0, 3));
    dnu_vec.push_back(h_ijklm(K_re_scl, 0, 0, 1, 1, 3));

    dnu_vec.push_back(h_ijklm(K_re_scl, 1, 1, 0, 0, 4));
    dnu_vec.push_back(h_ijklm(K_re_scl, 0, 0, 1, 1, 4));
  }

  dnu_vec.push_back(dnu.cst());
  dnu_vec.push_back(dnu_delta.cst());

  chi2_1 = 0e0; k = 0;

  chi2_1 += scl_ksi[1]*sqr(dnu_vec[k]); k++;
  chi2_1 += scl_ksi[1]*sqr(dnu_vec[k]); k++;

  chi2_1 += scl_dnu_conf*sqr(dnu_vec[k]); k++;
  chi2_1 += scl_dnu_delta_conf*sqr(dnu_vec[k]); k++;

  if (!true) {
    chi2_1 += scl_dnu[0]*sqr(dnu_vec[k]); k++;
    chi2_1 += scl_dnu[0]*sqr(dnu_vec[k]); k++;
    chi2_1 += scl_dnu[0]*sqr(dnu_vec[k]); k++;

    chi2_1 += scl_dnu[1]*sqr(dnu_vec[k]); k++;
    chi2_1 += scl_dnu[1]*sqr(dnu_vec[k]); k++;
    chi2_1 += scl_dnu[1]*sqr(dnu_vec[k]); k++;
    chi2_1 += scl_dnu[1]*sqr(dnu_vec[k]); k++;

    chi2_1 += scl_dnu[2]*sqr(dnu_vec[k]); k++;
    chi2_1 += scl_dnu[2]*sqr(dnu_vec[k]); k++;
    chi2_1 += scl_dnu[2]*sqr(dnu_vec[k]); k++;
    chi2_1 += scl_dnu[2]*sqr(dnu_vec[k]); k++;
    chi2_1 += scl_dnu[2]*sqr(dnu_vec[k]); k++;
  }

  if (chrom) {
    chi2_1 += scl_ksi[2]*sqr(dnu_vec[k]); k++;
    chi2_1 += scl_ksi[2]*sqr(dnu_vec[k]); k++;

    chi2_1 += scl_ksi[3]*sqr(dnu_vec[k]); k++;
    chi2_1 += scl_ksi[3]*sqr(dnu_vec[k]); k++;

    chi2_1 += scl_ksi[4]*sqr(dnu_vec[k]); k++;
    chi2_1 += scl_ksi[4]*sqr(dnu_vec[k]); k++;
  }

  if (chi2_1 < chi2) {
    printf("\nchi2: %21.15e -> %21.15e\n", chi2, chi2_1);
    for (k = 1; k <= (int)dnu_vec.size(); k++) {
      printf(" %10.3e", dnu_vec[k-1]);
      switch (k) {
      case 2:
      case 4:
      case 9:
      case 14:
      case 16:
      case 18: 
      case 20: 
	printf("\n");
	prt_ln = true;
	break;
      }
    }
    if (!prt_ln) printf("\n");
    printf("\n|dnu|       = %9.3e\n", dnu.cst());
    printf("|dnu_delta| = %9.3e\n", dnu_delta.cst());
  }

  return chi2_1;
}


void df_nl(double *bn, double *df)
{
  int    k;
  double eps;

  const bool prt = !false;

  bn_prms.set_prm(bn);
  for (k = 1; k <= bn_prms.n_bn; k++) {
    // eps = (k <= 0)? 1e-2 : 1e0;
    eps = (k <= 3)? 1e0 : 1e2;
    bn_prms.set_dparam(k, eps);
    df[k] = get_chi2();
    bn_prms.set_dparam(k, -2e0*eps);
    df[k] -= get_chi2();
    df[k] /= 2e0*eps;

    bn_prms.set_dparam(k, eps);
  }

  if (prt)
    dvdump(stdout, (char *)"\ndf_nl:", df, bn_prms.n_bn, (char *)" %12.5e");
}


double f_nl(double bn[])
{
  double chi2_1;

  n_iter++;
  bn_prms.set_prm(bn);

  chi2_1 = get_chi2();
  if (chi2_1 < chi2) {
    printf("bn:\n");
    bn_prms.prt_bn(bn);
    chi2 = chi2_1;

    prt_mfile("flat_file.fit");
    bn_prms.prt_bn_lat();
  }

  return chi2_1;
}


void conj_grad(param_type &bn_prms, double (*f)(double *),
		   void df(double *, double *))
{
  int          iter;
  double       *bn, fret;
  ss_vect<tps> A;

  const double ftol = 1e-15;

  bn = dvector(1, bn_prms.n_bn);

  bn_prms.ini_prm(bn);

  dfrprmn(bn, bn_prms.n_bn, ftol, &iter, &fret, f, df);

  prt_mfile("flat_file.fit");
  bn_prms.prt_bn_lat();

 free_dvector(bn, 1, bn_prms.n_bn);
}


void powell(param_type &bn_prms, double (*f)(double *))
{
  int    n_bn, i, j, iter;
  double *bn, **xi, fret;

  double const bn_tol = 1e-10;

  n_bn = bn_prms.n_bn;

  bn = dvector(1, n_bn); xi = dmatrix(1, n_bn, 1, n_bn);

  bn_prms.ini_prm(bn);

  // Set initial directions (unit vectors).
  for (i = 1; i <= n_bn; i++)
    for (j = 1; j <= n_bn; j++)
      xi[i][j] = (i == j)? 1e0 : 0e0;

  dpowell(bn, xi, n_bn, bn_tol, &iter, &fret, f);

  prt_mfile("flat_file.fit");
  bn_prms.prt_bn_lat();

  free_dvector(bn, 1, bn_prms.n_bn); free_dmatrix(xi, 1, n_bn, 1, n_bn);
}


void lat_select(void)
{

  if (!true) {
    bn_prms.add_prm("sf1", 3, -5e3, 5e3, 1e0);
    bn_prms.add_prm("sd1", 3, -5e3, 5e3, 1e0);
    bn_prms.add_prm("sd2", 3, -5e3, 5e3, 1e0);
  }

  if (!false) {
    // Sextupole Length is 0.1 m.

    // bn_prms.add_prm("sh1a", 4, -5e3, 5e3, 1.0);
    // bn_prms.add_prm("sh1b", 4, -5e3, 5e3, 1.0);

    bn_prms.add_prm("sh2",  4, -1e3/0.1, 1e3/0.1, 1e0);
    bn_prms.add_prm("s",    4, -1e3/0.1, 1e3/0.1, 1e0);
    bn_prms.add_prm("of1",  4, -1e3,     1e3,     1e0);

    // bn_prms.add_prm("sh2",  6, -1e6/0.1, 1e6/0.1, 1e0);
    // bn_prms.add_prm("s",    6, -1e6/0.1, 1e6/0.1, 1e0);
    // bn_prms.add_prm("of1",  6, -1e6,     1e6,     1e0);
  }

  if (false) {
    bn_prms.add_prm("sf1", 4, -1e3, 1e3, 1e0);
    bn_prms.add_prm("sd1", 4, -1e3, 1e3, 1e0);
    bn_prms.add_prm("sd2", 4, -1e3, 1e3, 1e0);
  }

  if (false) {
    bn_prms.add_prm("sf1", 5, -1e4, 1e4, 1e0);
    bn_prms.add_prm("sd1", 5, -1e4, 1e4, 1e0);
    bn_prms.add_prm("sd2", 5, -1e4, 1e4, 1e0);
  }

  // bn_prms.add_prm("sf1", 6, 1e4, 1e0);
  // bn_prms.add_prm("sd1", 6, 1e4, 1e0);
  // bn_prms.add_prm("sd2", 6, 1e4, 1e0);
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

  daeps_(tpsa_eps);

  Id_scl.identity();
  for (j = 0; j < 4; j++)
    Id_scl[j] *= sqrt(twoJ[j/2]);
  Id_scl[delta_] *= delta_max;

  Id_delta_scl.identity();
  for (j = 0; j < 4; j++)
    Id_delta_scl[j] *= sqrt(twoJ_delta[j/2]);
  Id_delta_scl[delta_] *= delta_max;

  lat_select();

  if (true)
    conj_grad(bn_prms, f_nl, df_nl);
  else
    powell(bn_prms, f_nl);
}
