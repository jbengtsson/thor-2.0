#include <cfloat>

#define NO 4

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
tps          h_re, h_im, h_re_scl, h_im_scl, K_re, K_im, K_re_scl;
tps          K_re_delta_scl;
ss_vect<tps> nus, nus_scl, Id_scl, Id_delta_scl;

const bool
  fit_ksi = false,
  scale   = false;

const int n_prt  = 8;

// Center of straight.
const double
  beta_inj[]   = {8.7, 2.1},
  twoJ[]       = {sqr(8e-3)/beta_inj[X_], sqr(2e-3)/beta_inj[Y_]},
  twoJ_delta[] = {sqr(0.5e-3)/beta_inj[X_], sqr(0.1e-3)/beta_inj[Y_]},
  delta_max    = 3e-2;

const double
  scl_h[]            = {0e-1,  0e-2, 0e-2},
  scl_dnu[]          = {1e-2, 1e-2, 1e-2, 1e-2},
  scl_ksi[]          = {1e5,  1e-2, 1e-2, 1e-2, 1e-2},
  delta_scl          = 0e0;


struct param_type {
private:

public:
  int                 m_constr, n_bn;
  double              bn_tol;
  std::vector<double> bn_max, bn_scl;
  std::vector<int>    Fnum, n, svd_list;

  void add_prm(const std::string Fname, const int n,
	       const double bn_max, const double bn_scl);
  void ini_prm(double *bn);
  void set_prm_dep(const int k) const;
  void clr_prm_dep(const int k) const;
  void set_prm(double *bn) const;
  void set_dparam(const int k, const double eps) const;
  void prt_bn(double *bn) const;
};


param_type bn_prms;


void param_type::add_prm(const std::string Fname, const int n,
			 const double bn_max, const double bn_scl)
{
  Fnum.push_back(get_Fnum(Fname.c_str()));
  this->n.push_back(n);
  this->bn_max.push_back(bn_max);
  this->bn_scl.push_back(bn_scl);
  n_bn = Fnum.size();
}


void param_type::ini_prm(double *bn)
{
  int    i;
  double L;

  const int n_prt = 4;

  printf("\nInitial bn (scale factor in parenthesis):\n");
  printf("  No of Families: %1d\n", n_bn);
  for (i = 1; i <= n_bn; i++) {
    bn[i] = get_bn(Fnum[i-1], 1, n[i-1]);

    if (scale) {
      L = get_L(Fnum[i-1], 1);
      if (L == 0e0) L = 1e0;
      bn_scl[i-1] = 1e0/sqrt(get_n_Kids(Fnum[i-1])*L);
    }

    bn[i] /= bn_scl[i-1];
    printf(" %12.5e (%9.3e)", bn_scl[i-1]*bn[i], bn_scl[i-1]);
    if (i % n_prt == 0) printf("\n");
  }
  if (n_bn % n_prt != 0) printf("\n");
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
  int i;

  const bool prt = false;

  if (prt) printf("set_prm:\n");
  for (i = 1; i <= n_bn; i++) {
    set_bn(Fnum[i-1], n[i-1], bn_scl[i-1]*bn[i]);
    if (prt) {
      printf(" %12.5e", bn_scl[i-1]*bn[i]);
      if (i % n_prt == 0) printf("\n");
    }
  }
  if (prt && (n_bn % n_prt != 0)) printf("\n");
}


void param_type::set_dparam(const int k, double eps) const
{
  const bool prt = false;

  if (prt) printf("set_dparam: %12.5e\n", eps);
  set_dbn(Fnum[k-1], n[k-1], eps);
}


void param_type::prt_bn(double *bn) const
{
  int k;

  for (k = 1; k <= n_bn; k++)
    printf(" %10.3e", bn[k]);
  printf("\n");
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


void prt_dnu(void)
{

  printf("\ndnu(2J=[%10.3e, %10.3e]):\n", twoJ[X_], twoJ[Y_]);
  printf("   k_22000    k_11110    k_00220\n");
  printf("   k_33000    k_22110    k_11220    k_00330\n");
  printf("   k_44000    k_33110    k_22220    k_11330    k_00440\n");
  printf("   k_55000    k_44110    k_33220    k_22330    k_11440    k_00550\n");
  printf(" %10.3e %10.3e %10.3e\n",
	 h_ijklm(K_re_scl, 2, 2, 0, 0, 0), h_ijklm(K_re_scl, 1, 1, 1, 1, 0),
	 h_ijklm(K_re_scl, 0, 0, 2, 2, 0));
  printf(" %10.3e %10.3e %10.3e %10.3e\n",
	 h_ijklm(K_re_scl, 3, 3, 0, 0, 0), h_ijklm(K_re_scl, 2, 2, 1, 1, 0),
	 h_ijklm(K_re_scl, 1, 1, 2, 2, 0), h_ijklm(K_re_scl, 0, 0, 3, 3, 0));
  printf(" %10.3e %10.3e %10.3e %10.3e %10.3e\n",
	 h_ijklm(K_re_scl, 4, 4, 0, 0, 0), h_ijklm(K_re_scl, 3, 3, 1, 1, 0),
	 h_ijklm(K_re_scl, 2, 2, 2, 2, 0), h_ijklm(K_re_scl, 1, 1, 3, 3, 0),
	 h_ijklm(K_re_scl, 0, 0, 4, 4, 0));
  printf(" %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e\n",
	 h_ijklm(K_re_scl, 5, 5, 0, 0, 0), h_ijklm(K_re_scl, 4, 4, 1, 1, 0),
	 h_ijklm(K_re_scl, 3, 3, 2, 2, 0), h_ijklm(K_re_scl, 2, 2, 3, 3, 0),
	 h_ijklm(K_re_scl, 1, 1, 4, 4, 0), h_ijklm(K_re_scl, 0, 0, 5, 5, 0));
  printf("dnu_x:\n %8.5f %8.5f\n",
	 h_ijklm(nus_scl[3], 1, 1, 0, 0, 0),
	 h_ijklm(nus_scl[3], 0, 0, 1, 1, 0));
  printf(" %8.5f %8.5f %8.5f\n",
	 h_ijklm(nus_scl[3], 2, 2, 0, 0, 0),
	 h_ijklm(nus_scl[3], 1, 1, 1, 1, 0),
	 h_ijklm(nus_scl[3], 0, 0, 2, 2, 0));
  printf(" %8.5f %8.5f %8.5f %8.5f\n",
	 h_ijklm(nus_scl[3], 3, 3, 0, 0, 0),
	 h_ijklm(nus_scl[3], 2, 2, 1, 1, 0),
	 h_ijklm(nus_scl[3], 1, 1, 2, 2, 0),
	 h_ijklm(nus_scl[3], 0, 0, 3, 3, 0));
  printf(" %8.5f %8.5f %8.5f %8.5f %8.5f\n",
	 h_ijklm(nus_scl[3], 4, 4, 0, 0, 0),
	 h_ijklm(nus_scl[3], 3, 3, 1, 1, 0),
	 h_ijklm(nus_scl[3], 2, 2, 2, 2, 0),
	 h_ijklm(nus_scl[3], 1, 1, 3, 3, 0),
	 h_ijklm(nus_scl[3], 0, 0, 4, 4, 0));
  printf("dnu_y:\n %8.5f %8.5f\n",
	 h_ijklm(nus_scl[4], 1, 1, 0, 0, 0),
	 h_ijklm(nus_scl[4], 0, 0, 1, 1, 0));
  printf(" %8.5f %8.5f %8.5f\n",
	 h_ijklm(nus_scl[4], 2, 2, 0, 0, 0),
	 h_ijklm(nus_scl[4], 1, 1, 1, 1, 0),
	 h_ijklm(nus_scl[4], 0, 0, 2, 2, 0));
  printf(" %8.5f %8.5f %8.5f %8.5f\n",
	 h_ijklm(nus_scl[4], 3, 3, 0, 0, 0),
	 h_ijklm(nus_scl[4], 2, 2, 1, 1, 0),
	 h_ijklm(nus_scl[4], 1, 1, 2, 2, 0),
	 h_ijklm(nus_scl[4], 0, 0, 3, 3, 0));
  printf(" %8.5f %8.5f %8.5f %8.5f %8.5f\n",
	 h_ijklm(nus_scl[4], 4, 4, 0, 0, 0),
	 h_ijklm(nus_scl[4], 3, 3, 1, 1, 0),
	 h_ijklm(nus_scl[4], 2, 2, 2, 2, 0),
	 h_ijklm(nus_scl[4], 1, 1, 3, 3, 0),
	 h_ijklm(nus_scl[4], 0, 0, 4, 4, 0));

  printf("\nksi(delta=%3.1f%%):\n", 1e2*delta_max);
  printf("   k_11001    k_00111    k_22001    k_11111    k_00221\n");
  printf("   k_11002    k_00112    k_22002    k_11112    k_00222"
	 "   k_33002    k_22112    k_11222    k_00332\n");
  printf("   k_11003    k_00113    k_22003    k_11113    k_00223"
	 "   k_33003    k_22113    k_11223    k_00333\n");
  printf("   k_11004    k_00114    k_22004    k_11114    k_00224"
	 "   k_33004    k_22114    k_11224    k_00334\n");
  printf("   k_11005    k_00115    k_22005    k_11115    k_00225"
	 "   k_33005    k_22115    k_11225    k_00335\n");
  printf(" %10.3e %10.3e %10.3e %10.3e %10.3e\n",
	 h_ijklm(K_re_scl, 1, 1, 0, 0, 1), h_ijklm(K_re_scl, 0, 0, 1, 1, 1),
	 h_ijklm(K_re_scl, 2, 2, 0, 0, 1), h_ijklm(K_re_scl, 1, 1, 1, 1, 1),
	 h_ijklm(K_re_scl, 0, 0, 2, 2, 1));
  printf(" %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e\n",
	 h_ijklm(K_re_scl, 1, 1, 0, 0, 2), h_ijklm(K_re_scl, 0, 0, 1, 1, 2),
	 h_ijklm(K_re_scl, 2, 2, 0, 0, 2), h_ijklm(K_re_scl, 1, 1, 1, 1, 2),
	 h_ijklm(K_re_scl, 0, 0, 2, 2, 2), h_ijklm(K_re_scl, 3, 3, 0, 0, 2),
	 h_ijklm(K_re_scl, 2, 2, 1, 1, 2), h_ijklm(K_re_scl, 1, 1, 2, 2, 2),
	 h_ijklm(K_re_scl, 0, 0, 3, 3, 2));
  printf(" %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e\n",
	 h_ijklm(K_re_scl, 1, 1, 0, 0, 3), h_ijklm(K_re_scl, 0, 0, 1, 1, 3),
	 h_ijklm(K_re_scl, 2, 2, 0, 0, 3), h_ijklm(K_re_scl, 1, 1, 1, 1, 3),
	 h_ijklm(K_re_scl, 0, 0, 2, 2, 3), h_ijklm(K_re_scl, 3, 3, 0, 0, 3),
	 h_ijklm(K_re_scl, 2, 2, 1, 1, 3), h_ijklm(K_re_scl, 1, 1, 2, 2, 3),
	 h_ijklm(K_re_scl, 0, 0, 3, 3, 3));
  printf(" %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e\n",
	 h_ijklm(K_re_scl, 1, 1, 0, 0, 4), h_ijklm(K_re_scl, 0, 0, 1, 1, 4),
	 h_ijklm(K_re_scl, 2, 2, 0, 0, 4), h_ijklm(K_re_scl, 1, 1, 1, 1, 4),
	 h_ijklm(K_re_scl, 0, 0, 2, 2, 4), h_ijklm(K_re_scl, 3, 3, 0, 0, 4),
	 h_ijklm(K_re_scl, 2, 2, 1, 1, 4), h_ijklm(K_re_scl, 1, 1, 2, 2, 4),
	 h_ijklm(K_re_scl, 0, 0, 3, 3, 4));
  printf(" %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e\n",
	 h_ijklm(K_re_scl, 1, 1, 0, 0, 5), h_ijklm(K_re_scl, 0, 0, 1, 1, 5),
	 h_ijklm(K_re_scl, 2, 2, 0, 0, 5), h_ijklm(K_re_scl, 1, 1, 1, 1, 5),
	 h_ijklm(K_re_scl, 0, 0, 2, 2, 5), h_ijklm(K_re_scl, 3, 3, 0, 0, 5),
	 h_ijklm(K_re_scl, 2, 2, 1, 1, 5), h_ijklm(K_re_scl, 1, 1, 2, 2, 5),
	 h_ijklm(K_re_scl, 0, 0, 3, 3, 5));
  printf("ksi_x:\n %8.5f %8.5f %8.5f\n",
	 h_ijklm(nus_scl[3], 0, 0, 0, 0, 1),
	 h_ijklm(nus_scl[3], 1, 1, 0, 0, 1),
	 h_ijklm(nus_scl[3], 0, 0, 1, 1, 1));
  printf(" %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f\n",
	 h_ijklm(nus_scl[3], 0, 0, 0, 0, 2),
	 h_ijklm(nus_scl[3], 1, 1, 0, 0, 2),
	 h_ijklm(nus_scl[3], 0, 0, 1, 1, 2),
	 h_ijklm(nus_scl[3], 2, 2, 0, 0, 2),
	 h_ijklm(nus_scl[3], 1, 1, 1, 1, 2),
	 h_ijklm(nus_scl[3], 0, 0, 2, 2, 2));
  printf(" %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f\n",
	 h_ijklm(nus_scl[3], 0, 0, 0, 0, 3),
	 h_ijklm(nus_scl[3], 1, 1, 0, 0, 3),
	 h_ijklm(nus_scl[3], 0, 0, 1, 1, 3),
	 h_ijklm(nus_scl[3], 2, 2, 0, 0, 3),
	 h_ijklm(nus_scl[3], 1, 1, 1, 1, 3),
	 h_ijklm(nus_scl[3], 0, 0, 2, 2, 3));
  printf(" %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f\n",
	 h_ijklm(nus_scl[3], 0, 0, 0, 0, 4),
	 h_ijklm(nus_scl[3], 1, 1, 0, 0, 4),
	 h_ijklm(nus_scl[3], 0, 0, 1, 1, 4),
	 h_ijklm(nus_scl[3], 2, 2, 0, 0, 4),
	 h_ijklm(nus_scl[3], 1, 1, 1, 1, 4),
	 h_ijklm(nus_scl[3], 0, 0, 2, 2, 4));
  printf(" %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f\n",
	 h_ijklm(nus_scl[3], 0, 0, 0, 0, 5),
	 h_ijklm(nus_scl[3], 1, 1, 0, 0, 5),
	 h_ijklm(nus_scl[3], 0, 0, 1, 1, 5),
	 h_ijklm(nus_scl[3], 2, 2, 0, 0, 5),
	 h_ijklm(nus_scl[3], 1, 1, 1, 1, 5),
	 h_ijklm(nus_scl[3], 0, 0, 2, 2, 5));
  printf("ksi_y:\n %8.5f %8.5f %8.5f\n",
	 h_ijklm(nus_scl[4], 0, 0, 0, 0, 1),
	 h_ijklm(nus_scl[4], 1, 1, 0, 0, 1),
	 h_ijklm(nus_scl[4], 0, 0, 1, 1, 1));
  printf(" %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f\n",
	 h_ijklm(nus_scl[4], 0, 0, 0, 0, 2),
	 h_ijklm(nus_scl[4], 1, 1, 0, 0, 2),
	 h_ijklm(nus_scl[4], 0, 0, 1, 1, 2),
	 h_ijklm(nus_scl[4], 2, 2, 0, 0, 2),
	 h_ijklm(nus_scl[4], 1, 1, 1, 1, 2),
	 h_ijklm(nus_scl[4], 0, 0, 2, 2, 2));
  printf(" %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f\n",
	 h_ijklm(nus_scl[4], 0, 0, 0, 0, 3),
	 h_ijklm(nus_scl[4], 1, 1, 0, 0, 3),
	 h_ijklm(nus_scl[4], 0, 0, 1, 1, 3),
	 h_ijklm(nus_scl[4], 2, 2, 0, 0, 3),
	 h_ijklm(nus_scl[4], 1, 1, 1, 1, 3),
	 h_ijklm(nus_scl[4], 0, 0, 2, 2, 3));
  printf(" %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f\n",
	 h_ijklm(nus_scl[4], 0, 0, 0, 0, 4),
	 h_ijklm(nus_scl[4], 1, 1, 0, 0, 4),
	 h_ijklm(nus_scl[4], 0, 0, 1, 1, 4),
	 h_ijklm(nus_scl[4], 2, 2, 0, 0, 4),
	 h_ijklm(nus_scl[4], 1, 1, 1, 1, 4),
	 h_ijklm(nus_scl[4], 0, 0, 2, 2, 4));
  printf(" %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f\n",
	 h_ijklm(nus_scl[4], 0, 0, 0, 0, 5),
	 h_ijklm(nus_scl[4], 1, 1, 0, 0, 5),
	 h_ijklm(nus_scl[4], 0, 0, 1, 1, 5),
	 h_ijklm(nus_scl[4], 2, 2, 0, 0, 5),
	 h_ijklm(nus_scl[4], 1, 1, 1, 1, 5),
	 h_ijklm(nus_scl[4], 0, 0, 2, 2, 5));

  printf("\nTune confinement:\n");
  // printf(" %11.3e %11.3e\n",
  // 	 h_ijklm(K_re/(3e0*twoJ[X_]), 2, 2, 0, 0, 0),
  // 	 h_ijklm(K_re, 3, 3, 0, 0, 0));
  // printf(" %11.3e %11.3e\n",
  // 	 h_ijklm(K_re/(3e0*twoJ[Y_]), 0, 0, 2, 2, 0),
  // 	 h_ijklm(K_re, 0, 0, 3, 3, 0));
  printf(" %11.3e %11.3e\n",
	 h_ijklm(nus_scl[3], 1, 1, 0, 0, 0),
	 h_ijklm(nus_scl[3], 2, 2, 0, 0, 0));
  printf(" %11.3e %11.3e\n",
	 h_ijklm(nus_scl[3], 0, 0, 1, 1, 0),
	 h_ijklm(nus_scl[3], 0, 0, 2, 2, 0));
  printf(" %11.3e %11.3e\n",
	 h_ijklm(nus_scl[4], 0, 0, 1, 1, 0),
	 h_ijklm(nus_scl[4], 0, 0, 2, 2, 0));
  printf(" %11.3e %11.3e\n",
	 h_ijklm(nus_scl[4], 1, 1, 0, 0, 0),
	 h_ijklm(nus_scl[4], 2, 2, 0, 0, 0));

  printf("\n %11.3e %11.3e %11.3e\n",
	 h_ijklm(K_re*Id_scl, 1, 1, 1, 1, 0),
	 h_ijklm(K_re*Id_scl, 2, 2, 1, 1, 0),
	 h_ijklm(K_re*Id_scl, 1, 1, 2, 2, 0));
}


void prt_system(const int m, const int n_b2, double **A, double *b)
{
  int i, j;

  printf("\n Ax = b:\n          ");
  for (j = 1; j <= n_b2; j++)
    printf("%11d", j);
  printf("\n");
  for (i = 1; i <= m; i++) {
    printf("%4d %10s", i, drv_term[i-1].c_str());
    for (j = 1; j <= n_b2; j++)
      printf("%11.3e", A[i][j]);
    printf("%11.3e\n", b[i]);
  }

  prt_dnu();
}


void prt_bn(const param_type &bn_prms)
{
  bool     first = true;
  long int loc;
  int      j, k, n, n_bn;
  double   bn;
  FILE     *outf;

  const std::string file_name = "dnu.out";

  outf = file_write(file_name.c_str());

  n_bn = bn_prms.n_bn;
  fprintf(outf, "\n");
  for (k = 0; k < n_bn; k++) {
    loc = get_loc(bn_prms.Fnum[k], 1) - 1;
    bn = get_bn(bn_prms.Fnum[k], 1, bn_prms.n[k]);
    if (bn_prms.n[k] == Sext)
      fprintf(outf,
	      "%-8s: sextupole, l = %7.5f"
	      ", k = %12.5e, n = nsext, Method = Meth;\n",
	      elem[loc].Name, elem[loc].L, bn);
    else {
      n = elem[loc].mpole->order;
      for (j = Sext; j <= n; j++) {
	bn = get_bn(bn_prms.Fnum[k], 1, j);
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
	if (j == n) {
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


double get_chi2(void)
{
  int                 k;
  double              chi2;
  std::vector<double> dnu;

  const bool prt = false;

  danot_(NO-1);
  get_Map();
  danot_(NO);
  K = MapNorm(Map, g, A1, A0, Map_res, 1);
  nus = dHdJ(K); nus_scl = nus*Id_scl;
  CtoR(K, K_re, K_im);
  K_re_scl = K_re*Id_scl; K_re_delta_scl = K_re*Id_delta_scl;
  CtoR(get_h(), h_re, h_im); h_re_scl = h_re*Id_scl; h_im_scl = h_im*Id_scl;

  dnu.push_back(h_ijklm(K_re_scl, 1, 1, 0, 0, 1));
  dnu.push_back(h_ijklm(K_re_scl, 0, 0, 1, 1, 1));

  dnu.push_back(h_ijklm(K_re_scl, 2, 2, 0, 0, 0));
  dnu.push_back(h_ijklm(K_re_scl, 1, 1, 1, 1, 0));
  dnu.push_back(h_ijklm(K_re_scl, 0, 0, 2, 2, 0));

  dnu.push_back(h_ijklm(K_re_scl, 1, 1, 0, 0, 2));
  dnu.push_back(h_ijklm(K_re_scl, 0, 0, 1, 1, 2));

  chi2 = 0e0;
  chi2 += scl_ksi[0]*sqr(dnu[0]);
  chi2 += scl_ksi[0]*sqr(dnu[1]);

  chi2 += scl_dnu[0]*sqr(dnu[2]);
  chi2 += scl_dnu[0]*sqr(dnu[3]);
  chi2 += scl_dnu[0]*sqr(dnu[4]);

  chi2 += scl_ksi[2]*sqr(dnu[5]);
  chi2 += scl_ksi[2]*sqr(dnu[6]);

  if (prt) {
    printf("\nget_chi2:\n ");
    for (k = 0; k < dnu.size(); k++)
      printf(" %10.3e", dnu[k]);
    printf("\n");
    printf("  chi2: %9.3e\n", chi2);
  }

  return chi2;
}


void df_nl(double *bn, double *df)
{
  int    k;
  double eps;

  const bool   prt = !false;

  bn_prms.set_prm(bn);
  for (k = 1; k <= bn_prms.n_bn; k++) {
    eps = (k <= 3)? 1e-2 : 1e1;
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
  int           k;
  double        chi2;
  static double chi20 = 1e30;

  const bool prt = !false;

  n_iter++;
  bn_prms.set_prm(bn);

  chi2 = get_chi2();
  if (chi2 < chi20) {
    printf("\nf_nl:\n");
    bn_prms.prt_bn(bn);
    printf("  chi2 = %9.3e -> %9.3e\n", chi20, chi2);
  }
  chi20 = min(chi2, chi20);

  return chi2;
}


void fit_conj_grad(param_type &bn_prms, double (*f)(double *),
		   void df(double *, double *))
{
  int          iter, k;
  double       *bn, fret;
  ss_vect<tps> A;

  const double ftol = 1e-8;

  bn = dvector(1, bn_prms.n_bn);

  bn_prms.ini_prm(bn);
  f(bn);

  dfrprmn(bn, bn_prms.n_bn, ftol, &iter, &fret, f, df);

  // printf("\n  iter = %d fret = %12.5e\n", iter, fret);
  // printf("bns:\n");
  // bn_prms.set_prm(bn);
  // f_prt(bn);

  prt_mfile("flat_file.fit");

 free_dvector(bn, 1, bn_prms.n_bn);
}


void lat_select(void)
{
  bn_prms.add_prm("sf1", 3, 1e4, 1.0);
  bn_prms.add_prm("sd1", 3, 1e4, 1.0);
  bn_prms.add_prm("sd2", 3, 1e4, 1.0);

  bn_prms.add_prm("sf1", 4, 1e4, 1.0);
  bn_prms.add_prm("sd1", 4, 1e4, 1.0);
  bn_prms.add_prm("sd2", 4, 1e4, 1.0);

  // bn_prms.add_prm("sh1a", 4, 1e4, 1.0);
  // bn_prms.add_prm("sh1b", 4, 1e4, 1.0);

  bn_prms.add_prm("sh2",  4, 1e4, 1.0);
  bn_prms.add_prm("s",    4, 1e4, 1.0);
  bn_prms.add_prm("of1",  4, 1e4, 1.0);
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

  Id_scl.identity();
  for (j = 0; j < 4; j++)
    Id_scl[j] *= sqrt(twoJ[j/2]);
  Id_scl[delta_] *= delta_max;

  Id_delta_scl.identity();
  for (j = 0; j < 4; j++)
    Id_delta_scl[j] *= sqrt(twoJ_delta[j/2]);
  Id_delta_scl[delta_] *= delta_max;

  lat_select();
  fit_conj_grad(bn_prms, f_nl, df_nl);
}
