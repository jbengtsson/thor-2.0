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

#define DNU 1


extern tps          K, g;
extern ss_vect<tps> Map, A0, A1, Map_res;

const double tpsa_eps = 1e-30;

std::vector<std::string> drv_term;

int          n_iter;
double       chi2 = 1e30;
tps          h_re, h_im, h_re_scl, h_im_scl, K_re, K_im, K_re_scl;
tps          K_re_delta_scl;
ss_vect<tps> nus, nus_scl, Id_scl, Id_delta_scl;

const int n_prt = 8;

// Center of straight.
const double
  beta_inj[]   = {8.7, 2.1},
  A_max[]      = {3e-3, 1e-3},
  twoJ[]       = {sqr(A_max[X_])/beta_inj[X_], sqr(A_max[Y_])/beta_inj[Y_]},
  twoJ_delta[] = {sqr(0.5e-3)/beta_inj[X_], sqr(0.1e-3)/beta_inj[Y_]},
  delta_max    = 2.5e-2;

const double
  scl_h[]            = {0e0, 0e0, 0e0},
  scl_dnu[]          = {1e0, 1e0, 1e0, 0e0, 0e0},
  scl_ksi[]          = {0e0, 0e5, 0e0, 0e0, 0e0, 0e0}, // 1st not used.
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
  std::vector<double> bn_min, bn_max, dbn;
  std::vector<int>    Fnum, n;

  void add_prm(const std::string Fname, const int n,
	       const double bn_min, const double bn_max, const double dbn);
  void ini_prm(double *bn);
  void set_prm(double *bn) const;
  void set_dparam(const int k, const double eps) const;
  void prt_bn(double *bn) const;
  void prt_bn_lat(void) const;
};


param_type bn_prms;


void param_type::add_prm(const std::string Fname, const int n,
			 const double bn_min, const double bn_max,
			 const double dbn)
{
  Fnum.push_back(get_Fnum(Fname.c_str()));
  this->n.push_back(n);
  this->bn_min.push_back(bn_min);
  this->bn_max.push_back(bn_max);
  this->dbn.push_back(dbn);
  n_bn = Fnum.size();
}


void param_type::ini_prm(double *bn)
{
  int i, loc;

  printf("\nInitial bn (scale factor in parenthesis):\n");
  printf("  No of Families: %1d\n", n_bn);
  for (i = 0; i < n_bn; i++) {
    bn[i+1] = get_bn(Fnum[i], 1, n[i]);

    // Bounded.
    if ((bn_min[i] <= bn[i+1]) && (bn[i+1] <= bn_max[i]))
      bn[i+1] = bn_internal(bn[i+1], bn_min[i], bn_max[i]);
    else {
      loc = get_loc(Fnum[i], 1);
      printf("\nini_prm:\n outside range ");
      printf(" %s %10.3e [%10.3e, %10.3e]\n",
	     elem[loc].Name, bn[i+1], bn_min[i], bn_max[i]);
      exit(1);
    }
  }

  prt_bn(bn);
}


void param_type::set_prm(double *bn) const
{
  int    i;
  double bn_ext;

  const bool prt = false;

  for (i = 0; i < n_bn; i++) {
    // Bounded.
    bn_ext = bn_bounded(bn[i+1], bn_min[i], bn_max[i]);
    set_bn(Fnum[i], n[i], bn_ext);
  }

  if (prt) {
    printf("set_prm:\n");
    prt_bn(bn);
  }
}


void param_type::set_dparam(const int k, double eps) const
{
  char ch;

  const bool prt = false;

  if (prt) {
    ch = (eps >= 0e0)? '+' : '-';
    printf("\nset_dparam:\n  %12.5e %c %11.5e",
	   get_bn(Fnum[k], 1, n[k]), ch, fabs(eps));
  }
  set_dbn(Fnum[k], n[k], eps);
  if (prt) printf(" = %12.5e\n", get_bn(Fnum[k], 1, n[k]));
}


void param_type::prt_bn(double *bn) const
{
  int    i;
  double bn_ext;

  for (i = 0; i < n_bn; i++) {
    // Bounded.
    bn_ext = bn_bounded(bn[i+1], bn_min[i], bn_max[i]);
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
static double (*func_save_2D)(const double, const double);

double gauss_quad_2D(double (*func)(const double),
		     const double a, const double b)
{
  int    j;
  double xr, xm, dx, s;

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

double gauss_quad_2D_fy(const double y)  { return (*func_save_2D)(xsav_2D, y); }

double gauss_quad_2D_fx(const double x)
{
  xsav_2D = x;
  return gauss_quad_2D(gauss_quad_2D_fy, gauss_quad_2D_y0(x),
		       gauss_quad_2D_y1(x));
}

double gauss_quad_2D(double (*func)(const double, const double),
		     const double x0, const double x1)
{
  func_save_2D = func;
  return gauss_quad_2D(gauss_quad_2D_fx, x0, x1);
}

double f_gauss_quad_2D(double x, double y)
{
  int             k;
  long int        jj[ss_dim];
  double          dnu[2];
  tps             dK;
  ss_vect<double> ps;
  ss_vect<tps>    Id;

  // printf("f_gauss_quad_2D\n");
  ps.zero();
  ps[x_] = ps[px_] = sqrt(x);
  ps[y_] = ps[py_] = sqrt(y);
  ps[delta_] = 0e0;

#if DNU
  for (k = 0; k < 2; k++) {
    dnu[k] = ((nus[k+3]-nus[k+3].cst())*ps).cst();
    // Compute absolute value.
    if (dnu[k] < 0e0) dnu[k] = -dnu[k];
  }

  return dnu[X_]*dnu[Y_]/(twoJ[X_]*twoJ[Y_]);
#else
  dK = K_re;
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

  return (dK/(twoJ[X_]*twoJ[Y_])).cst();
#endif
}

static double xsav_3D, ysav_3D;
static double (*func_save_3D)(const double, const double, const double);

double gauss_quad_3D(double (*func)(const double),
		     const double a, const double b)
{
  int    j;
  double xr, xm, dx, s;

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

double gauss_quad_3D_fz(const double z)
{
  return (*func_save_3D)(xsav_3D, ysav_3D, z);
}

double gauss_quad_3D_fy(const double y)
{
  ysav_3D = y;
  return gauss_quad_3D(gauss_quad_3D_fz,
		       gauss_quad_3D_z0(xsav_3D, y),
		       gauss_quad_3D_z1(xsav_3D, y));
}

double gauss_quad_3D_fx(const double x)
{
  xsav_3D = x;
  return gauss_quad_3D(gauss_quad_3D_fy, gauss_quad_3D_y0(x),
		       gauss_quad_3D_y1(x));
}

double gauss_quad_3D(double (*func)(const double, const double, const double),
		     const double x0, const double x1)
{
  func_save_3D = func;
  return gauss_quad_3D(gauss_quad_3D_fx, x0, x1);
}

double f_gauss_quad_3D(double x, double y, double z)
{
  int             k;
  long int        jj[ss_dim];
  double          dnu[2];
  tps             dK;
  ss_vect<double> ps;
  ss_vect<tps>    Id;

  // printf("f_gauss_quad_3D\n");

  ps.zero();
  ps[x_] = ps[px_] = sqrt(x);
  ps[y_] = ps[py_] = sqrt(y);
  ps[delta_] = z;

#if DNU
  for (k = 0; k < 2; k++) {
    dnu[k] = ((nus[k+3]-nus[k+3].cst())*ps).cst();
    // Compute absolute value.
    if (dnu[k] < 0e0) dnu[k] = -dnu[k];
  }

  return dnu[X_]*dnu[Y_]/(twoJ[X_]*twoJ[Y_]*2e0*delta_max);
#else
  dK = K_re;
  for (k = 0; k < ss_dim; k++)
    jj[k] = 0;
  jj[x_] = 1; jj[px_] = 1;
  dK.pook(jj, 0e0);
  jj[x_] = 0; jj[px_] = 0; jj[y_] = 1; jj[py_] = 1;
  dK.pook(jj, 0e0);
  // std::cout << std::scientific << std::setprecision(3)
  // 	    << "\n |dK| = " << dK << "\n";
  dK = dK*ps;
  // Compute absolute value.
  if (dK.cst() < 0e0) dK = -dK;
  // std::cout << std::scientific << std::setprecision(3)
  // 	    << "\n |dK| = " << dK << "\n";

  return (dK/(twoJ[X_]*twoJ[Y_]*2e0*delta_max)).cst();
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


void get_dK(std::vector<double> &dK)
{

  dK.push_back(h_ijklm(K_re, 1, 1, 0, 0, 1));
  dK.push_back(h_ijklm(K_re, 0, 0, 1, 1, 1));

  dK.push_back(h_ijklm(K_re_scl, 2, 2, 0, 0, 0));
  dK.push_back(h_ijklm(K_re_scl, 1, 1, 1, 1, 0));
  dK.push_back(h_ijklm(K_re_scl, 0, 0, 2, 2, 0));

  dK.push_back(h_ijklm(K_re_scl, 3, 3, 0, 0, 0));
  dK.push_back(h_ijklm(K_re_scl, 2, 2, 1, 1, 0));
  dK.push_back(h_ijklm(K_re_scl, 1, 1, 2, 2, 0));
  dK.push_back(h_ijklm(K_re_scl, 0, 0, 3, 3, 0));

  dK.push_back(h_ijklm(K_re_scl, 4, 4, 0, 0, 0));
  dK.push_back(h_ijklm(K_re_scl, 3, 3, 1, 1, 0));
  dK.push_back(h_ijklm(K_re_scl, 2, 2, 2, 2, 0));
  dK.push_back(h_ijklm(K_re_scl, 1, 1, 3, 3, 0));
  dK.push_back(h_ijklm(K_re_scl, 0, 0, 4, 4, 0));

  dK.push_back(h_ijklm(K_re_scl, 1, 1, 0, 0, 2));
  dK.push_back(h_ijklm(K_re_scl, 0, 0, 1, 1, 2));

  dK.push_back(h_ijklm(K_re_scl, 1, 1, 0, 0, 3));
  dK.push_back(h_ijklm(K_re_scl, 0, 0, 1, 1, 3));

  dK.push_back(h_ijklm(K_re_scl, 1, 1, 0, 0, 4));
  dK.push_back(h_ijklm(K_re_scl, 0, 0, 1, 1, 4));
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


double get_chi2(const bool prt)
{
  int                 k;
  double              chi2_1, dnu, dnu_delta;
  std::vector<double> dK;

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

  get_dK(dK);

  chi2_1 = 0e0; k = 0;

  chi2_1 += scl_ksi[1]*sqr(dK[k]); k++;
  chi2_1 += scl_ksi[1]*sqr(dK[k]); k++;

  chi2_1 += scl_dnu_conf*sqr(dnu);
  chi2_1 += scl_dnu_delta_conf*sqr(dnu_delta);

  if (!true) {
    chi2_1 += scl_dnu[0]*sqr(dK[k]); k++;
    chi2_1 += scl_dnu[0]*sqr(dK[k]); k++;
    chi2_1 += scl_dnu[0]*sqr(dK[k]); k++;

    chi2_1 += scl_dnu[1]*sqr(dK[k]); k++;
    chi2_1 += scl_dnu[1]*sqr(dK[k]); k++;
    chi2_1 += scl_dnu[1]*sqr(dK[k]); k++;
    chi2_1 += scl_dnu[1]*sqr(dK[k]); k++;

    chi2_1 += scl_dnu[2]*sqr(dK[k]); k++;
    chi2_1 += scl_dnu[2]*sqr(dK[k]); k++;
    chi2_1 += scl_dnu[2]*sqr(dK[k]); k++;
    chi2_1 += scl_dnu[2]*sqr(dK[k]); k++;
    chi2_1 += scl_dnu[2]*sqr(dK[k]); k++;
  }

  if (chrom) {
    chi2_1 += scl_ksi[2]*sqr(dK[k]); k++;
    chi2_1 += scl_ksi[2]*sqr(dK[k]); k++;

    chi2_1 += scl_ksi[3]*sqr(dK[k]); k++;
    chi2_1 += scl_ksi[3]*sqr(dK[k]); k++;

    chi2_1 += scl_ksi[4]*sqr(dK[k]); k++;
    chi2_1 += scl_ksi[4]*sqr(dK[k]); k++;
  }

  if (prt && (chi2_1 < chi2)) {
    printf("\nchi2: %21.15e -> %21.15e\n", chi2, chi2_1);
    prt_dnu();
  }

  return chi2_1;
}


void df_nl(double *bn, double *df)
{
  int    k;
  double eps;

  const bool prt = !false;

  bn_prms.set_prm(bn);
  for (k = 0; k < bn_prms.n_bn; k++) {
    eps = bn_prms.dbn[k];
    bn_prms.set_dparam(k, eps);
    df[k+1] = get_chi2(false);
    bn_prms.set_dparam(k, -2e0*eps);
    df[k+1] -= get_chi2(false);
    df[k+1] /= 2e0*eps;

    bn_prms.set_dparam(k, eps);
  }

  if (prt)
    dvdump(stdout, (char *)"\ndf_nl:", df, bn_prms.n_bn, (char *)" %12.5e");
}


double f_nl(double *bn)
{
  double chi2_1;

  n_iter++;
  bn_prms.set_prm(bn);

  chi2_1 = get_chi2(true);
  if (chi2_1 < chi2) {
    printf("\nbn:\n");
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
    bn_prms.add_prm("sf1", 3, -5e3, 5e3, 1e-2);
    bn_prms.add_prm("sd1", 3, -5e3, 5e3, 1e-2);
    bn_prms.add_prm("sd2", 3, -5e3, 5e3, 1e-2);
  }

  if (!false) {
    // Sextupole Length is 0.1 m.

    // bn_prms.add_prm("sh1a", 4, -5e3/0.1, 5e3/0.1, 1e1);
    // bn_prms.add_prm("sh1b", 4, -5e3/0.1, 5e3/0.1, 1e1);

    bn_prms.add_prm("sh2",  4, -1e3/0.1, 1e3/0.1, 1e1);
    bn_prms.add_prm("s",    4, -1e3/0.1, 1e3/0.1, 1e1);
    bn_prms.add_prm("of1",  4, -1e3,     1e3,     1e1);

    // bn_prms.add_prm("sh2",  6, -1e6/0.1, 1e6/0.1, 1e2);
    // bn_prms.add_prm("s",    6, -1e6/0.1, 1e6/0.1, 1e2);
    // bn_prms.add_prm("of1",  6, -1e6,     1e6,     1e2);
  }

  if (false) {
    bn_prms.add_prm("sf1", 4, -1e3, 1e3, 1e0);
    bn_prms.add_prm("sd1", 4, -1e3, 1e3, 1e0);
    bn_prms.add_prm("sd2", 4, -1e3, 1e3, 1e0);
  }

  if (false) {
    bn_prms.add_prm("sf1", 5, -1e4, 1e4, 1e1);
    bn_prms.add_prm("sd1", 5, -1e4, 1e4, 1e1);
    bn_prms.add_prm("sd2", 5, -1e4, 1e4, 1e1);
  }

  // bn_prms.add_prm("sf1", 6, 1e4, 1e2);
  // bn_prms.add_prm("sd1", 6, 1e4, 1e2);
  // bn_prms.add_prm("sd2", 6, 1e4, 1e2);
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
