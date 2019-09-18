#include <cfloat>

#define NO 4

#include "thor_lib.h"

#include "Powell/src/newuoa.h"

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


typedef std::vector< std::vector<double> > perf_vec;


const double tpsa_eps = 1e-30;

int          n_iter;
double       chi2_ref = 1e30;
tps          h_re, h_im, h_re_scl, h_im_scl, K_re, K_im, K_re_scl;
tps          K_re_delta_scl;
ss_vect<tps> nus, nus_scl, Id_scl, Id_delta_scl;

#define LAT_CASE 1

// Center of straight.
const double
#if LAT_CASE == 1
// M-H6BA-17E-79pm-01.03-01.
  beta_inj[]     = {7.9, 3.1},
#elif LAT_CASE == 2
// M-H6BA-17E-69pm-04.02-01
  beta_inj[]     = {11.1, 5.5},
#endif
  A_max[]        = {3e-3, 1.5e-3},
  delta_max      = 2e-2,
  twoJ[]         = {sqr(A_max[X_])/beta_inj[X_], sqr(A_max[Y_])/beta_inj[Y_]},
  twoJ_delta[]   = {sqr(0.5e-3)/beta_inj[X_], sqr(0.1e-3)/beta_inj[Y_]};

const double
  scl_h[]        = {0e0, 0e0, 0e0},
  scl_dnu[]      = {0e-2, 0e-2, 0e-2},
  scl_ksi[]      = {0e0, 1e0, 0e0, 0e0, 0e0, 0e0}, // 1st not used.
  delta_scl      = 0e0,
  dx_dJ_scl      = 1e4,
  // Negative: minimize,
  // Positive: maintain opposite signs;
  // increase weight on remaining until opposite signs are obtained.
#define CASE_DNU 4
#if CASE_DNU == 1
  scl_dnu_conf[] = {-1e1, -1e1, -1e1, -1e1, -1e1, -1e1,
                     0e1,  0e1},
#elif CASE_DNU == 2
  scl_dnu_conf[] = {1e1, 1e1, 1e1, 1e1, 1e1, 1e1,
                    0e1, 0e1},
#elif CASE_DNU == 3
  scl_dnu_conf[] = {1e1, 1e1, 1e1, 1e1, 0e1, 0e1,
                    0e1, 0e1},
#elif CASE_DNU == 4
  scl_dnu_conf[] = {-1e1, -1e1, -1e1, -1e1, -1e1, -1e1,
                    0e-1, 0e-1},
#endif
#if DNU
  scl_dnu_2d     = 1e6,
#else
  scl_dnu_2d     = 0e14,
#endif
  scl_dnu_3d     = 0e0;


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
  void set_prm_dep(const int k) const;
  void clr_prm_dep(const int k) const;
  void set_dprm(const int k, const double eps) const;
  void prt_bn(double *bn) const;
  void prt_bn_lat(const std::string &fnam) const;
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


void param_type::set_dprm(const int k, double eps) const
{
  char ch;

  const bool prt = false;

  if (prt) {
    ch = (eps >= 0e0)? '+' : '-';
    printf("\nset_dprm:\n  %12.5e %c %11.5e",
	   get_bn(Fnum[k], 1, n[k]), ch, fabs(eps));
  }
  set_dbn(Fnum[k], n[k], eps);
  if (prt) printf(" = %12.5e\n", get_bn(Fnum[k], 1, n[k]));
}


void param_type::prt_bn(double *bn) const
{
  int    i;
  double bn_ext;

  const int n_prt = 6;

  printf("  ");
  for (i = 0; i < n_bn; i++) {
    // Bounded.
    bn_ext = bn_bounded(bn[i+1], bn_min[i], bn_max[i]);
    printf(" %12.5e", bn_ext);
    if ((i+1) % n_prt == 0) printf("\n  ");
  }
  if (n_bn % n_prt != 0) printf("\n");
}


void param_type::prt_bn_lat(const std::string &fnam) const
{
  bool     first = true;
  long int loc;
  int      j, k, n_ord;
  double   bn;
  FILE     *outf;

  outf = file_write(fnam.c_str());

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

double gauss_quad_2D_fy(const double y) { return (*func_save_2D)(xsav_2D, y); }

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
  double          dnu_xy, dK_abs;
  tps             dK;
  ss_vect<double> ps;

  const bool prt = false;

  if (prt) printf("f_gauss_quad_2D\n");
  ps.zero();
  ps[x_] = ps[px_] = sqrt(x); ps[y_] = ps[py_] = sqrt(y); ps[delta_] = 0e0;

#if DNU
  dnu_xy = abs(((nus[3]-nus[3].cst())*(nus[4]-nus[4].cst()))*ps);

  return dnu_xy/(twoJ[X_]*twoJ[Y_]);
#else
  dK = K_re;
  for (k = 0; k < ss_dim; k++)
    jj[k] = 0;
  jj[x_] = 1; jj[px_] = 1;
  dK.pook(jj, 0e0);
  jj[x_] = 0; jj[px_] = 0; jj[y_] = 1; jj[py_] = 1;
  dK.pook(jj, 0e0);
  if (prt)
    std::cout << std::scientific << std::setprecision(3)
	      << std::setw(11) << dK << "\n";
  dK_abs = abs(dK*ps);
  if (prt)
    std::cout << std::scientific << std::setprecision(3)
  	      << "\n |dK| = " << dK_abs << "\n";

  return dK_abs/(twoJ[X_]*twoJ[Y_]);
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

double gauss_quad_3D_y1(const double x) { return twoJ_delta[Y_]; }

double gauss_quad_3D_z0(const double x, const double y) { return -delta_max; }

double gauss_quad_3D_z1(const double x, const double y) { return delta_max; }

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
  double          dnu_xy, dK_abs;
  tps             dK;
  ss_vect<double> ps;

  const bool prt = false;

  if (prt) printf("f_gauss_quad_3D\n");

  ps.zero();
  ps[x_] = ps[px_] = sqrt(x); ps[y_] = ps[py_] = sqrt(y); ps[delta_] = z;
  if (prt)
    std::cout << std::scientific << std::setprecision(3)
	      << "\nps:\n" << std::setw(11) << ps << "\n";

#if DNU
  dnu_xy = abs(((nus[3]-nus[3].cst())*(nus[4]-nus[4].cst()))*ps);

  return dnu_xy/(twoJ[X_]*twoJ[Y_]*2e0*delta_max);
#else
  dK = K_re;
  for (k = 0; k < ss_dim; k++)
    jj[k] = 0;
  jj[x_] = 1; jj[px_] = 1;
  dK.pook(jj, 0e0);
  jj[x_] = 0; jj[px_] = 0; jj[y_] = 1; jj[py_] = 1;
  dK.pook(jj, 0e0);
  if (prt)
    std::cout << std::scientific << std::setprecision(3)
	      << "\n |dK| = " << dK << "\n";
  dK_abs = abs(dK*ps);
  if (prt)
    std::cout << std::scientific << std::setprecision(3)
	      << "\n |dK| = " << dK_abs << "\n";

  return dK_abs/(twoJ[X_]*twoJ[Y_]*2e0*delta_max);
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


tps get_a(const tps &t,
	  const int i, const int j, const int k, const int l, const int m)
{
  return h_ijklm(t, i, j, k, l, m) + h_ijklm_p(t, i, j, k, l, m, 7)*tps(0e0, 7);
}


tps dacfu1(const tps &a, double (*func)(const long int []))
{
  char    name[11];
  int     j, n;
  long int jj[ss_dim], ibuf1[bufsize], ibuf2[bufsize];
  double  rbuf[bufsize];
  tps     b;

  a.exprt(rbuf, ibuf1, ibuf2, name); n = (int)rbuf[0];

  for (j = 0; j < n; j++) {
    dehash_(no_tps, ss_dim, ibuf1[j], ibuf2[j], jj);

    rbuf[j+1] *= (*func)(jj);
  }

  b.imprt(n, rbuf, ibuf1, ibuf2);

  // Remove zeroes.
  return 1e0*b;
}


double f_kernel(const long int jj[])
{

  return
    (((jj[x_] != 0) || (jj[y_] != 0)) &&
     (jj[x_] == jj[px_]) && (jj[y_] == jj[py_]))?
    1e0 : 0e0;
}


void get_dx_dJ(double dx2[], const bool prt)
{
  int           j, k;
  double        dx[2];
  ss_vect<tps>  Id, Id_scl1, dx_fl, dx_re, dx_im, M;
  std::ofstream outf;

  const int no_b3 = 3;

  if (prt) outf.open("dx_dJ.out", std::ios::out);

  Id.identity();

  Id_scl1.identity();
  for (j = 0; j < 4; j++)
    Id_scl1[j] *= sqrt(twoJ[j/2]);
  Id_scl1[delta_] = 0e0;

  danot_(no_b3-1);
  Map.identity(); Map.propagate(1, n_elem);
  danot_(no_b3);

  M.identity();
  for (k = 0; k < 2; k++)
    dx2[k] = 0e0;
  for (j = 1; j <= n_elem; j++) {
    M.propagate(j, j);
    if ((elem[j-1].kind == Mpole) && (elem[j-1].mpole->bn[Sext-1] != 0e0)) {
      K = MapNorm(M*Map*Inv(M), g, A1, A0, Map_res, 1);
#if 0
      dx_fl = LieExp(g, Id);
#else
      for (k = 0; k < 4; k++)
      	dx_fl[k] = PB(g, Id[k]);
#endif
      for (k = 0; k < 4; k++)
	CtoR(dx_fl[k], dx_re[k], dx_im[k]);
      dx_re = A1*dx_re; dx_im = A1*dx_im;
      dx[X_] = h_ijklm(dx_re[x_]*Id_scl, 1, 1, 0, 0, 0);
      dx[Y_] = h_ijklm(dx_re[x_]*Id_scl, 0, 0, 1, 1, 0);
      for (k = 0; k < 2; k++)
	// dx2[k] += sqr(dx[k]);
	dx2[k] +=
	  elem[j-1].mpole->bn[Sext-1]*elem[j-1].L*elem[j-1].Beta[k]*dx[k];
      if (prt) {
	outf << std::setw(4) << j << std::fixed << std::setprecision(3)
	     << std::setw(8) << elem[j-1].S
	     << " " << std::setw(8) << elem[j-1].Name;
	for (k = 0; k < 2; k++)
	  outf << std::scientific << std::setprecision(5)
	       << std::setw(13)
	       << elem[j-1].mpole->bn[Sext-1]*elem[j-1].L*dx[k];
	outf << "\n";
      }
    }
  }

  for (k = 0; k < 2; k++)
    dx2[k] = sqr(dx2[k]);

  printf("  dx(J) = [%9.3e, %9.3e]\n", dx2[X_], dx2[Y_]);
  if (prt) outf.close();
}


void prt_dnu(void)
{

  const bool chrom = !true;

  printf("\ndnu(2J=[%9.3e, %9.3e]):\n", twoJ[X_], twoJ[Y_]);
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

  if (chrom) {
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
  }

  printf("\nTune confinement:\n");
  printf(" %11.3e %11.3e\n",
	 h_ijklm(nus_scl[3], 1, 1, 0, 0, 0),
	 2e0*h_ijklm(nus_scl[3], 2, 2, 0, 0, 0));
  printf(" %11.3e %11.3e\n",
	 h_ijklm(nus_scl[3], 0, 0, 1, 1, 0),
	 2e0*h_ijklm(nus_scl[3], 0, 0, 2, 2, 0));
  printf(" %11.3e %11.3e\n",
	 h_ijklm(nus_scl[4], 0, 0, 1, 1, 0),
	 2e0*h_ijklm(nus_scl[4], 0, 0, 2, 2, 0));
  printf(" %11.3e %11.3e\n",
	 h_ijklm(nus_scl[4], 1, 1, 0, 0, 0),
	 2e0*h_ijklm(nus_scl[4], 2, 2, 0, 0, 0));

  printf("\n %11.3e %11.3e\n",
	 h_ijklm(nus_scl[3], 0, 0, 0, 0, 2),
	 2e0*h_ijklm(nus_scl[3], 0, 0, 0, 0, 4));
  printf(" %11.3e %11.3e\n",
	 h_ijklm(nus_scl[4], 0, 0, 0, 0, 2),
	 2e0*h_ijklm(nus_scl[4], 0, 0, 0, 0, 4));

  printf("\n %11.3e %11.3e\n",
	 h_ijklm(nus_scl[3], 0, 0, 0, 0, 3),
	 5e0/3e0*h_ijklm(nus_scl[3], 0, 0, 0, 0, 5));
  printf(" %11.3e %11.3e\n",
	 h_ijklm(nus_scl[4], 0, 0, 0, 0, 3),
	 5e0/3e0*h_ijklm(nus_scl[4], 0, 0, 0, 0, 5));
}


void get_dK(std::vector<tps> &dK)
{
  double dnu, dnu_delta;

  const double K_scl = 1e6;

  danot_(NO-1);
  get_Map();
  danot_(NO);
  K = MapNorm(Map, g, A1, A0, Map_res, 1);
  nus = dHdJ(K); nus_scl = nus*Id_scl;
  CtoR(K, K_re, K_im);
  K_re_scl = K_scl*K_re*Id_scl; K_re_delta_scl = K_scl*K_re*Id_delta_scl;
  CtoR(get_h(), h_re, h_im); h_re_scl = h_re*Id_scl; h_im_scl = h_im*Id_scl;

  dnu = gauss_quad_2D(f_gauss_quad_2D, 0e0, twoJ[X_]);
  dnu_delta = gauss_quad_3D(f_gauss_quad_3D, 0e0, twoJ_delta[X_]);

  dK.clear();

  dK.push_back(get_a(nus[3], 0, 0, 0, 0, 1));
  dK.push_back(get_a(nus[4], 0, 0, 0, 0, 1));

  dK.push_back(get_a(nus_scl[3], 1, 1, 0, 0, 0));
  dK.push_back(get_a(nus_scl[3], 2, 2, 0, 0, 0));
  dK.push_back(get_a(nus_scl[3], 0, 0, 1, 1, 0));
  dK.push_back(get_a(nus_scl[3], 0, 0, 2, 2, 0));

  dK.push_back(get_a(nus_scl[4], 0, 0, 1, 1, 0));
  dK.push_back(get_a(nus_scl[4], 0, 0, 2, 2, 0));
  dK.push_back(get_a(nus_scl[4], 1, 1, 0, 0, 0));
  dK.push_back(get_a(nus_scl[4], 2, 2, 0, 0, 0));

  dK.push_back(get_a(nus_scl[3], 0, 0, 0, 0, 2));
  dK.push_back(get_a(nus_scl[3], 0, 0, 0, 0, 4));

  dK.push_back(get_a(nus_scl[4], 0, 0, 0, 0, 2));
  dK.push_back(get_a(nus_scl[4], 0, 0, 0, 0, 4));

  dK.push_back(get_a(nus_scl[3], 0, 0, 0, 0, 3));
  dK.push_back(get_a(nus_scl[4], 0, 0, 0, 0, 3));

  dK.push_back(dnu);
  dK.push_back(dnu_delta);

  dK.push_back(get_a(K_re_scl, 2, 2, 0, 0, 0));
  dK.push_back(get_a(K_re_scl, 1, 1, 1, 1, 0));
  dK.push_back(get_a(K_re_scl, 0, 0, 2, 2, 0));

  dK.push_back(get_a(K_re_scl, 3, 3, 0, 0, 0));
  dK.push_back(get_a(K_re_scl, 2, 2, 1, 1, 0));
  dK.push_back(get_a(K_re_scl, 1, 1, 2, 2, 0));
  dK.push_back(get_a(K_re_scl, 0, 0, 3, 3, 0));

  dK.push_back(get_a(K_re_scl, 4, 4, 0, 0, 0));
  dK.push_back(get_a(K_re_scl, 3, 3, 1, 1, 0));
  dK.push_back(get_a(K_re_scl, 2, 2, 2, 2, 0));
  dK.push_back(get_a(K_re_scl, 1, 1, 3, 3, 0));
  dK.push_back(get_a(K_re_scl, 0, 0, 4, 4, 0));

  dK.push_back(get_a(K_re_scl, 1, 1, 0, 0, 2));
  dK.push_back(get_a(K_re_scl, 0, 0, 1, 1, 2));

  dK.push_back(get_a(K_re_scl, 1, 1, 0, 0, 3));
  dK.push_back(get_a(K_re_scl, 0, 0, 1, 1, 3));

  dK.push_back(get_a(K_re_scl, 1, 1, 0, 0, 4));
  dK.push_back(get_a(K_re_scl, 0, 0, 1, 1, 4));
}


tps tps_abs(const tps &a) { return (a.cst() > 0e0)? a : -a; }


template<typename T>
void dK_shift(const double scl, const T dnu1, const T dnu2, std::vector<T> &b)
{
  // scl > 0: maintain tune confinement; 
  T val;

  if ((sgn(dnu1.cst()) != sgn(dnu2.cst())) || (scl < 0e0))
    val = fabs(scl)*sqr(dnu1+2e0*dnu2);
  else
    val = scl*1e30;

  b.push_back(val);
}


template<typename T>
void get_b(std::vector<T> &dK, std::vector<T> &b)
{
  int k;

  const bool chrom = false;

  k = 0;
  b.clear();

  b.push_back(scl_ksi[1]*sqr(dK[k])); k++;
  b.push_back(scl_ksi[1]*sqr(dK[k])); k++;

  dK_shift(scl_dnu_conf[0], dK[k], dK[k+1], b); k += 2;
  dK_shift(scl_dnu_conf[1], dK[k], dK[k+1], b); k += 2;
  dK_shift(scl_dnu_conf[2], dK[k], dK[k+1], b); k += 2;
  dK_shift(scl_dnu_conf[3], dK[k], dK[k+1], b); k += 2;

  dK_shift(scl_dnu_conf[4], dK[k], dK[k+1], b); k += 2;
  dK_shift(scl_dnu_conf[5], dK[k], dK[k+1], b); k += 2;

  b.push_back(scl_dnu_conf[6]*sqr(dK[k])); k++;
  b.push_back(scl_dnu_conf[7]*sqr(dK[k])); k++;

  b.push_back(scl_dnu_2d*sqr(dK[k])); k++;
  b.push_back(scl_dnu_3d*sqr(dK[k])); k++;

  if (true) {
    b.push_back(scl_dnu[0]*sqr(dK[k])); k++;
    b.push_back(scl_dnu[0]*sqr(dK[k])); k++;
    b.push_back(scl_dnu[0]*sqr(dK[k])); k++;

    b.push_back(scl_dnu[1]*sqr(dK[k])); k++;
    b.push_back(scl_dnu[1]*sqr(dK[k])); k++;
    b.push_back(scl_dnu[1]*sqr(dK[k])); k++;
    b.push_back(scl_dnu[1]*sqr(dK[k])); k++;

    b.push_back(scl_dnu[2]*sqr(dK[k])); k++;
    b.push_back(scl_dnu[2]*sqr(dK[k])); k++;
    b.push_back(scl_dnu[2]*sqr(dK[k])); k++;
    b.push_back(scl_dnu[2]*sqr(dK[k])); k++;
    b.push_back(scl_dnu[2]*sqr(dK[k])); k++;
  }

  if (chrom) {
    b.push_back(scl_ksi[2]*sqr(dK[k])); k++;
    b.push_back(scl_ksi[2]*sqr(dK[k])); k++;

    b.push_back(scl_ksi[3]*sqr(dK[k])); k++;
    b.push_back(scl_ksi[3]*sqr(dK[k])); k++;

    b.push_back(scl_ksi[4]*sqr(dK[k])); k++;
    b.push_back(scl_ksi[4]*sqr(dK[k])); k++;
  }
}


double get_chi2(const bool prt)
{
  int              n, j, k, n_extra;
  double           chi2, bn, dx2[2];
  std::vector<int> Fnum_extra;
  std::vector<tps> dK, b, b_extra;
  static bool      first = true;

  const int n_prt = 4;

#define CASE_SCL 6

  // First minimize, then balance.
#if CASE_SCL == 1
  // Equalize.
  const bool   chi2_extra = true;
  const double scl        = 1e3;
#elif CASE_SCL == 2
  const bool   chi2_extra = true;
  const double scl        = 1e2;
#elif CASE_SCL == 3
  const bool   chi2_extra = true;
  const double scl        = 1e1;
#elif CASE_SCL == 4
  const bool   chi2_extra = true;
  const double scl        = 1e0;
#elif CASE_SCL == 5
  const bool   chi2_extra = true;
  const double scl        = 1e-2;
#elif CASE_SCL == 6
  const bool   chi2_extra = false;
  const double scl        = 0e0;
#endif

  get_dK(dK);
  get_b(dK, b);
  get_dx_dJ(dx2, first);

  chi2 = 0e0;
  n = (int)b.size();
  for (k = 0; k < n; k++)
    chi2 += b[k].cst();

  for (k = 0; k < 2; k++)
    chi2 += dx_dJ_scl*dx2[k];


  if (chi2_extra) {
    b_extra.clear();

    if (!false) {
      b_extra.push_back(scl*sqr(h_ijklm(nus_scl[3], 1, 1, 0, 0, 0)));
      b_extra.push_back(scl*sqr(2e0*h_ijklm(nus_scl[3], 2, 2, 0, 0, 0)));
      b_extra.push_back(scl*sqr(h_ijklm(nus_scl[3], 0, 0, 1, 1, 0)));
      b_extra.push_back(scl*sqr(2e0*h_ijklm(nus_scl[3], 0, 0, 2, 2, 0)));

      b_extra.push_back(scl*sqr(h_ijklm(nus_scl[4], 0, 0, 1, 1, 0)));
      b_extra.push_back(scl*sqr(2e0*h_ijklm(nus_scl[4], 0, 0, 2, 2, 0)));
      b_extra.push_back(scl*sqr(h_ijklm(nus_scl[4], 1, 1, 0, 0, 0)));
      b_extra.push_back(scl*sqr(2e0*h_ijklm(nus_scl[4], 2, 2, 0, 0, 0)));

      b_extra.push_back(scl*sqr(h_ijklm(nus_scl[3], 0, 0, 0, 0, 2)));
      b_extra.push_back(scl*sqr(h_ijklm(nus_scl[3], 0, 0, 0, 0, 3)));
      b_extra.push_back(scl*sqr(h_ijklm(nus_scl[3], 0, 0, 0, 0, 4)));

      b_extra.push_back(scl*sqr(h_ijklm(nus_scl[4], 0, 0, 0, 0, 2)));
      b_extra.push_back(scl*sqr(h_ijklm(nus_scl[4], 0, 0, 0, 0, 3)));
      b_extra.push_back(scl*sqr(h_ijklm(nus_scl[4], 0, 0, 0, 0, 4)));
    }

    for (k = 0; k < (int)b_extra.size(); k++)
      chi2 += b_extra[k].cst();
   }

  if (prt && (first || (chi2 < chi2_ref))) {
    first = false;
    printf("\nget_chi2(%1d): scl = %9.3e\n", n, scl);

    get_dx_dJ(dx2, true);
    prt_dnu();

    k = 0;
    printf("\n  ksi1        = [%7.5f, %7.5f] ([%9.3e, %9.3e])\n",
	   h_ijklm(nus[3], 0, 0, 0, 0, 1), h_ijklm(nus[4], 0, 0, 0, 0, 1),
	   b[k].cst(), b[k+1].cst());
    k += 2;
    printf("\n  Tune Conf.:   %10.3e %10.3e %10.3e %10.3e\n",
    	   b[k].cst(), b[k+1].cst(), b[k+2].cst(), b[k+3].cst());
    printf("                %10.3e %10.3e\n", b[k+4].cst(), b[k+5].cst());
    printf("                %10.3e %10.3e\n", b[k+6].cst(), b[k+7].cst());
    k += 8;
    if (chi2_extra) {
      printf("\n  b_extra     =");
      n_extra = b_extra.size();
      for (j = 0; j < n_extra; j++) {
	printf(" %10.3e", b_extra[j].cst());
	if ((j != n_extra-1) && ((j+1) % n_prt == 0))
	  printf("\n               ");
      }
      printf("\n");
    }
    printf("\n  |dnu|       = %15.8e\n", b[k].cst());
    printf("  |dnu_delta| = %15.8e\n", b[k+1].cst());
    k += 2;
    printf("\n  K:            %10.3e %10.3e %10.3e\n",
    	   b[k].cst(), b[k+1].cst(), b[k+2].cst());
    k += 3;
    printf("                %10.3e %10.3e %10.3e %10.3e\n",
    	   b[k].cst(), b[k+1].cst(), b[k+2].cst(), b[k+3].cst());
    k += 4;
    printf("                %10.3e %10.3e %10.3e %10.3e %10.3e\n",
    	   b[k].cst(), b[k+1].cst(), b[k+2].cst(), b[k+3].cst(), b[k+4].cst());
    k += 5;

    printf("\n  %-4d chi2: %21.15e -> %21.15e\n", n_iter, chi2_ref, chi2);
  }

  return chi2;
}


void df_nl2(double *bn, double *df)
{
  int              i, k;
  tps              chi2;
  std::vector<tps> dK, b, c;

  const bool prt = false;

  bn_prms.set_prm(bn);
  for (i = 0; i < bn_prms.n_bn; i++) {
    bn_prms.set_prm_dep(i);

    get_dK(dK);
    get_b(dK, b);

    chi2 = 0e0;
    for (k = 0; k < (int)b.size(); k++)
      chi2 += b[k];
    df[i+1] = h_ijklm_p(chi2, 0, 0, 0, 0, 0, 7);
 
    bn_prms.clr_prm_dep(i);
  }

  if (prt)
    dvdump(stdout, (char *)"\ndf_nl2:", df, bn_prms.n_bn, (char *)" %12.5e");
}


void df_nl(double *bn, double *df)
{
  int    k;
  double eps;

  const bool prt = !false;

  bn_prms.set_prm(bn);
  for (k = 0; k < bn_prms.n_bn; k++) {
    eps = bn_prms.dbn[k];
    bn_prms.set_dprm(k, eps);
    df[k+1] = get_chi2(false);
    bn_prms.set_dprm(k, -2e0*eps);
    df[k+1] -= get_chi2(false);
    df[k+1] /= 2e0*eps;
    bn_prms.set_dprm(k, eps);
  }

  df_nl2(bn, df);

  if (prt)
    dvdump(stdout, (char *)"\ndf_nl:", df, bn_prms.n_bn, (char *)" %12.5e");
}


double f_nl(double *bn)
 {
   double chi2;

   bn_prms.set_prm(bn);

   chi2 = get_chi2(true);
   if (chi2 < chi2_ref) {
    n_iter++;

    printf("\n  bn:\n");
    bn_prms.prt_bn(bn);
    fflush(stdout);
    bn_prms.prt_bn_lat("dnu.out");
    prt_mfile("flat_file.fit");

    chi2_ref = chi2;
  }

   return chi2;
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

  n_iter = 1;

  dfrprmn(bn, bn_prms.n_bn, ftol, &iter, &fret, f, df);

  prt_mfile("flat_file.fit");
  bn_prms.prt_bn_lat("dnu.out");

 free_dvector(bn, 1, bn_prms.n_bn);
}


void powell(param_type &bn_prms, double (*f)(double *))
{
  const int
    n_bn    = bn_prms.n_bn,
    n_pt    = n_bn + 2,
    n_w     = (n_pt+13)*(n_pt+n_bn)+3*n_bn*(n_bn+3)/2,
    n_prt   = 2,
    max_fun = 5000;

  const double
    rho_beg = bn_prms.dbn[0],
    rho_end = 1e-6;

  int    i, j, iter;
  double *bn, **xi, fret, w[n_w], x[n_bn];

  double const bn_tol = 1e-10;

  bn = dvector(1, n_bn); xi = dmatrix(1, n_bn, 1, n_bn);

  bn_prms.ini_prm(bn);

  n_iter = 1;

  if (true) {
    // Set initial directions (unit vectors).
    for (i = 1; i <= n_bn; i++)
      for (j = 1; j <= n_bn; j++)
	xi[i][j] = (i == j)? 1e0 : 0e0;

    dpowell(bn, xi, n_bn, bn_tol, &iter, &fret, f);
  } else {
    for (i = 1; i <= n_bn; i++)
      x[i-1] = bn[i];
    newuoa_(n_bn, n_pt, x, rho_beg, rho_end, n_prt, max_fun, w);
  }

  prt_mfile("flat_file.fit");
  bn_prms.prt_bn_lat("dnu.out");

  free_dvector(bn, 1, bn_prms.n_bn); free_dmatrix(xi, 1, n_bn, 1, n_bn);
}


void fit_ksi1(const double ksi_x, const double ksi_y,
	      const std::vector<int> &Fnum)
{
  int    i, n_svd;
  double **A, **U, **V, *w, *b, *dbn;

  const bool prt = false;

  const int
    m    = 2,
    n_b3 = Fnum.size();

  const double svd_cut = 1e-10;

  A = dmatrix(1, m, 1, n_b3); U = dmatrix(1, m, 1, n_b3);
  V = dmatrix(1, n_b3, 1, n_b3);
  w = dvector(1, n_b3); b = dvector(1, m); dbn = dvector(1, n_b3);

  for (i = 1; i <= n_b3; i++) {
    set_bn_par(Fnum[i-1], Sext, 7);

    danot_(3);
    get_Map();
    danot_(4);
    K = MapNorm(Map, g, A1, A0, Map_res, 1); nus = dHdJ(K);

    A[1][i] = h_ijklm_p(nus[3], 0, 0, 0, 0, 1, 7);
    A[2][i] = h_ijklm_p(nus[4], 0, 0, 0, 0, 1, 7);

    clr_bn_par(Fnum[i-1], Sext);
  }

  b[1] = -(h_ijklm(nus[3], 0, 0, 0, 0, 1)-ksi_x);
  b[2] = -(h_ijklm(nus[4], 0, 0, 0, 0, 1)-ksi_y);

  dmcopy(A, m, n_b3, U); dsvdcmp(U, m, n_b3, w, V);

  if (prt) printf("  singular values:\n   ");
  n_svd = 0;
  for (i = 1; i <= n_b3; i++) {
    if (prt) printf("%10.3e", w[i]);
    if (w[i] < svd_cut) {
      w[i] = 0e0;
      if (prt) printf(" (zeroed)");
    } else {
      if (n_svd > 2) {
	if (prt) printf("fit_ksi1: more than 2 non-zero singular values");
	exit(1);
      }
    }
  }
  if (prt) printf("\n");

  dsvbksb(U, w, V, m, n_b3, b, dbn);

  for (i = 1; i <= n_b3; i++)
    set_dbn(Fnum[i-1], Sext, dbn[i]);

  if (prt) {
    printf("  b3:\n  ");
    for (i = 1; i <= n_b3; i++)
      printf(" %12.5e (%12.5e)", get_bn(Fnum[i-1], 1, Sext), dbn[i]);
    printf("\n");
  }

  free_dmatrix(A, 1, m, 1, n_b3); free_dmatrix(U, 1, m, 1, n_b3);
  free_dmatrix(V, 1, n_b3, 1, n_b3);
  free_dvector(w, 1, n_b3); free_dvector(b, 1, m); free_dvector(dbn, 1, n_b3);
}


void Bubble_Sort2(const int ind, perf_vec &w)
{
  bool   swapped;
  int    j, k, n;
  double val;

  do {
    swapped = false;
    n = w[0].size();
    for (k = 0; k < (int)w.size()-1; k++) {
      if (fabs(w[k][ind-1]) < fabs(w[k+1][ind-1])) {
	for (j = 0; j < n; j++) {
	  val = w[k][j]; w[k][j] = w[k+1][j]; w[k+1][j] = val;
	}
	swapped = true;
      }
    }
  } while (swapped);
}


double rnd(const double x_min, const double x_max)
{
  return (x_max-x_min)*(double)rand()/(double)RAND_MAX + x_min;
}


void get_perf(int &n_good, double &chi2)
{
    danot_(NO-1);
    get_Map();
    danot_(NO);
    K = MapNorm(Map, g, A1, A0, Map_res, 1);
    nus = dHdJ(K); nus_scl = nus*Id_scl;

    n_good = 0;
    if (sgn(h_ijklm(nus_scl[3], 1, 1, 0, 0, 0))
	!= sgn(h_ijklm(nus_scl[3], 2, 2, 0, 0, 0)))
      n_good++;
    if (sgn(h_ijklm(nus_scl[3], 0, 0, 1, 1, 0))
	!= sgn(h_ijklm(nus_scl[3], 0, 0, 2, 2, 0)))
      n_good++;
    if (sgn(h_ijklm(nus_scl[4], 0, 0, 1, 1, 0))
	!= sgn(h_ijklm(nus_scl[4], 0, 0, 2, 2, 0)))
      n_good++;
    if (sgn(h_ijklm(nus_scl[4], 1, 1, 0, 0, 0))
	!= sgn(h_ijklm(nus_scl[4], 2, 2, 0, 0, 0)))
      n_good++;

    if (sgn(h_ijklm(nus_scl[3], 0, 0, 0, 0, 2))
	!= sgn(h_ijklm(nus_scl[3], 0, 0, 0, 0, 4)))
      n_good++;
    if (sgn(h_ijklm(nus_scl[4], 0, 0, 0, 0, 2))
	!= sgn(h_ijklm(nus_scl[4], 0, 0, 0, 0, 4)))
      n_good++;

    chi2 = get_chi2(false);
}


void prt_perf(std::vector<double> p)
{
  int k, n;

  const int ind = 16;

  printf("\n");
  n = p.size();
  for (k = 0; k < n; k++) {
    if (k == n-ind)
      printf(" %10d", (int)(p[k]+0.5));
    else
      printf(" %10.3e", p[k]);
    if ((k == n-ind+1) || (k == n-3)) printf("\n");
  }
  printf("\n");
  fflush(stdout);
}


void bn_mc(const int n_stats, const int ind, const int n_ksi)
{
  int                 j, k, n_good;
  double              chi2, r = 0e0;
  std::vector<int>    Fnum_ksi1, sgns;
  std::vector<double> p;
  perf_vec            perf;

 for (k = 0; k < n_ksi; k++)
    Fnum_ksi1.push_back(bn_prms.Fnum[k]);

  printf("\nbn_mc1:\n");
  for (j = 1; j <= n_stats; j++) {
    sgns.clear();
    for (k = n_ksi; k < bn_prms.n_bn; k++) {
      r = rnd(bn_prms.bn_min[k], bn_prms.bn_max[k]);
      set_bn(bn_prms.Fnum[k], bn_prms.n[k], r);
    }
    if (n_ksi != 0) fit_ksi1(0e0, 0e0, Fnum_ksi1);

    get_perf(n_good, chi2);

    p.clear();
    for (k = 0; k < bn_prms.n_bn; k++)
      p.push_back(get_bn(bn_prms.Fnum[k], 1, bn_prms.n[k]));

    p.push_back(n_good);

    p.push_back(chi2);

    p.push_back(h_ijklm(nus_scl[3], 1, 1, 0, 0, 0));
    p.push_back(h_ijklm(2e0*nus_scl[3], 2, 2, 0, 0, 0));
    p.push_back(h_ijklm(nus_scl[3], 0, 0, 1, 1, 0));
    p.push_back(h_ijklm(2e0*nus_scl[3], 0, 0, 2, 2, 0));

    p.push_back(h_ijklm(nus_scl[4], 0, 0, 1, 1, 0));
    p.push_back(h_ijklm(2e0*nus_scl[4], 0, 0, 2, 2, 0));
    p.push_back(h_ijklm(nus_scl[4], 1, 1, 0, 0, 0));
    p.push_back(h_ijklm(2e0*nus_scl[4], 2, 2, 0, 0, 0));

    p.push_back(h_ijklm(nus_scl[3], 0, 0, 0, 0, 2));
    p.push_back(h_ijklm(2e0*nus_scl[3], 0, 0, 0, 0, 4));
    p.push_back(h_ijklm(nus_scl[4], 0, 0, 0, 0, 2));
    p.push_back(h_ijklm(2e0*nus_scl[4], 0, 0, 0, 0, 4));

    p.push_back(h_ijklm(nus_scl[3], 0, 0, 0, 0, 3));
    p.push_back(h_ijklm(nus_scl[4], 0, 0, 0, 0, 3));

    perf.push_back(p);

    prt_perf(p);
  }

  Bubble_Sort2(ind, perf);

  printf("\n");
  for (j = 0; j < (int)perf.size(); j++)
    prt_perf(perf[j]);
}


void m_c(const int n)
{

  const int rand_seed = 100001;

  srand(rand_seed);

  switch (1) {
  case 1:
    bn_prms.add_prm("sf1", 3, -2e2,   2e2, 1e-2);
    bn_prms.add_prm("sd1", 3, -2.5e2, 2.5e2, 1e-2);
    bn_prms.add_prm("sd2", 3, -2e2,   0e2,   1e-2);

    bn_prms.add_prm("s",   3, -1.5e2, 1.5e2, 1e-2);
    bn_prms.add_prm("sh2", 3, -1.5e2, 1.5e2, 1e-2);
 
    bn_mc(n, bn_prms.n_bn+2, 2);
    break;
  case 2:
    bn_prms.add_prm("sf1", 3, -4.5e2,  4.5e2, 1e-2);
    bn_prms.add_prm("sd1", 3, -4.5e2,  4.5e2, 1e-2);
    bn_prms.add_prm("sd2", 3, -3e2,    0e2,   1e-2);

    bn_prms.add_prm("s",   3, -2e2, 2e2, 1e-2);
    bn_prms.add_prm("sh2", 3, -2e2, 2e2, 1e-2);

    bn_prms.add_prm("sf1", 4, -1e3, 1e3, 1e-2);
    bn_prms.add_prm("sd1", 4, -1e3, 1e3, 1e-2);
    bn_prms.add_prm("sd2", 4, -1e3, 1e3, 1e-2);

    bn_mc(n, bn_prms.n_bn+2, 2);
    break;
  case 3:
    bn_prms.add_prm("sf1", 3, -4.5e2,  4.5e2, 1e-2);
    bn_prms.add_prm("sd1", 3, -4.5e2,  4.5e2, 1e-2);
    bn_prms.add_prm("sd2", 3, -3e2,    0e2,   1e-2);

    bn_prms.add_prm("s",   3, -2e2, 2e2, 1e-2);
    bn_prms.add_prm("sh2", 3, -2e2, 2e2, 1e-2);

    bn_prms.add_prm("sf1", 4, -1e3, 1e3, 1e-2);
    bn_prms.add_prm("sd1", 4, -1e3, 1e3, 1e-2);
    bn_prms.add_prm("sd2", 4, -1e3, 1e3, 1e-2);

    bn_prms.add_prm("of1", 4, -1e3, 1e3, 1e-2);

    bn_mc(n, bn_prms.n_bn+2, 2);
    break;
  default:
    printf("\nm_c: unknown case\n");
    break;
  }
}


void lat_select(void)
{

  const double
    //                         b_3   b_4  b_5  b_6
    bn_max[] = {0e0, 0e0, 0e0, 2e3,  1e6, 5e7, 1e9},
    dbn[]    = {0e0, 0e0, 0e0, 1e-2, 1e0, 1e1, 1e0};

  switch (6) {
  case 1:
    // First minimize magnitude of tune footprint.
    // 3+0 b_3.
    bn_prms.add_prm("sf1", 3, -bn_max[3], bn_max[3], dbn[3]);
    bn_prms.add_prm("sd1", 3, -bn_max[3], bn_max[3], dbn[3]);
    bn_prms.add_prm("sd2", 3, -bn_max[3], bn_max[3], dbn[3]);
    break;
  case 2:
    // 3+0 b_3, 3+0 b_4.
    bn_prms.add_prm("sf1", 4, -bn_max[4], bn_max[4], dbn[4]);
    bn_prms.add_prm("sd1", 4, -bn_max[4], bn_max[4], dbn[4]);
    bn_prms.add_prm("sd2", 4, -bn_max[4], bn_max[4], dbn[4]);

    bn_prms.add_prm("sf1", 3, -bn_max[3], bn_max[3], dbn[3]);
    bn_prms.add_prm("sd1", 3, -bn_max[3], bn_max[3], dbn[3]);
    bn_prms.add_prm("sd2", 3, -bn_max[3], bn_max[3], dbn[3]);
    break;
  case 3:
    // Then balance terms.
    // 3+2 b_3+0.
    bn_prms.add_prm("s",   3, -bn_max[3], bn_max[3], dbn[3]);
    bn_prms.add_prm("sh2", 3, -bn_max[3], bn_max[3], dbn[3]);

    bn_prms.add_prm("sf1", 3, -bn_max[3], bn_max[3], dbn[3]);
    bn_prms.add_prm("sd1", 3, -bn_max[3], bn_max[3], dbn[3]);
    bn_prms.add_prm("sd2", 3, -bn_max[3], bn_max[3], dbn[3]);
    break;
  case 4:
    // 3+2 b_3, 3+0 b_4.
    bn_prms.add_prm("sf1", 4, -bn_max[4], bn_max[4], dbn[4]);
    bn_prms.add_prm("sd1", 4, -bn_max[4], bn_max[4], dbn[4]);
    bn_prms.add_prm("sd2", 4, -bn_max[4], bn_max[4], dbn[4]);

    bn_prms.add_prm("s",   3, -bn_max[3], bn_max[3], dbn[3]);
    bn_prms.add_prm("sh2", 3, -bn_max[3], bn_max[3], dbn[3]);

    bn_prms.add_prm("sf1", 3, -bn_max[3], bn_max[3], dbn[3]);
    bn_prms.add_prm("sd1", 3, -bn_max[3], bn_max[3], dbn[3]);
    bn_prms.add_prm("sd2", 3, -bn_max[3], bn_max[3], dbn[3]);
    break;
  case 5:
    // 3+2 b_3, 3+2 b_4.
    bn_prms.add_prm("sf1", 4, -bn_max[4], bn_max[4], dbn[4]);
    bn_prms.add_prm("sd1", 4, -bn_max[4], bn_max[4], dbn[4]);
    bn_prms.add_prm("sd2", 4, -bn_max[4], bn_max[4], dbn[4]);

    bn_prms.add_prm("s",   4, -bn_max[4], bn_max[4], dbn[4]);
    bn_prms.add_prm("sh2", 4, -bn_max[4], bn_max[4], dbn[4]);

    bn_prms.add_prm("s",   3, -bn_max[3], bn_max[3], dbn[3]);
    bn_prms.add_prm("sh2", 3, -bn_max[3], bn_max[3], dbn[3]);

    bn_prms.add_prm("sf1", 3, -bn_max[3], bn_max[3], dbn[3]);
    bn_prms.add_prm("sd1", 3, -bn_max[3], bn_max[3], dbn[3]);
    bn_prms.add_prm("sd2", 3, -bn_max[3], bn_max[3], dbn[3]);
    break;
  case 6:
    // 3+4 b_3.
    bn_prms.add_prm("sh1a", 3, -bn_max[3], bn_max[3], dbn[3]);
    bn_prms.add_prm("sh1b", 3, -bn_max[3], bn_max[3], dbn[3]);
    bn_prms.add_prm("s",    3, -bn_max[3], bn_max[3], dbn[3]);
    bn_prms.add_prm("sh2",  3, -bn_max[3], bn_max[3], dbn[3]);

    bn_prms.add_prm("sf1",  3, -bn_max[3], bn_max[3], dbn[3]);
    bn_prms.add_prm("sd1",  3, -bn_max[3], bn_max[3], dbn[3]);
    bn_prms.add_prm("sd2",  3, -bn_max[3], bn_max[3], dbn[3]);
    break;
  case 7:
    // 3+4 b_3, 3+0 b_4.
    bn_prms.add_prm("sf1",  4, -bn_max[4], bn_max[4], dbn[4]);
    bn_prms.add_prm("sd1",  4, -bn_max[4], bn_max[4], dbn[4]);
    bn_prms.add_prm("sd2",  4, -bn_max[4], bn_max[4], dbn[4]);

    bn_prms.add_prm("sh1a", 3, -bn_max[3], bn_max[3], dbn[3]);
    bn_prms.add_prm("sh1b", 3, -bn_max[3], bn_max[3], dbn[3]);
    bn_prms.add_prm("s",    3, -bn_max[3], bn_max[3], dbn[3]);
    bn_prms.add_prm("sh2",  3, -bn_max[3], bn_max[3], dbn[3]);

    bn_prms.add_prm("sf1",  3, -bn_max[3], bn_max[3], dbn[3]);
    bn_prms.add_prm("sd1",  3, -bn_max[3], bn_max[3], dbn[3]);
    bn_prms.add_prm("sd2",  3, -bn_max[3], bn_max[3], dbn[3]);
    break;
  case 8:
    bn_prms.add_prm("sh1a", 4, -bn_max[4], bn_max[4], dbn[4]);
    bn_prms.add_prm("sh1b", 4, -bn_max[4], bn_max[4], dbn[4]);

    bn_prms.add_prm("s",    4, -bn_max[4], bn_max[4], dbn[4]);
    bn_prms.add_prm("sh2",  4, -bn_max[4], bn_max[4], dbn[4]);

    bn_prms.add_prm("sf1",  3, -bn_max[3], bn_max[3], dbn[3]);
    bn_prms.add_prm("sd1",  3, -bn_max[3], bn_max[3], dbn[3]);
    bn_prms.add_prm("sd2",  3, -bn_max[3], bn_max[3], dbn[3]);
    break;
  case 9:
    bn_prms.add_prm("sf1", 4, -bn_max[4], bn_max[4], dbn[4]);
    bn_prms.add_prm("sd1", 4, -bn_max[4], bn_max[4], dbn[4]);
    bn_prms.add_prm("sd2", 4, -bn_max[4], bn_max[4], dbn[4]);

    bn_prms.add_prm("sf1", 5, -bn_max[5], bn_max[5], dbn[5]);
    bn_prms.add_prm("sd1", 5, -bn_max[5], bn_max[5], dbn[5]);
    bn_prms.add_prm("sd2", 5, -bn_max[5], bn_max[5], dbn[5]);

    bn_prms.add_prm("s",   3, -bn_max[3], bn_max[3], dbn[3]);
    bn_prms.add_prm("sh2", 3, -bn_max[3], bn_max[3], dbn[3]);

    bn_prms.add_prm("sf1", 3, -bn_max[3], bn_max[3], dbn[3]);
    bn_prms.add_prm("sd1", 3, -bn_max[3], bn_max[3], dbn[3]);
    bn_prms.add_prm("sd2", 3, -bn_max[3], bn_max[3], dbn[3]);
    break;
  case 10:
    bn_prms.add_prm("s",   4, -bn_max[4], bn_max[4], dbn[4]);
    bn_prms.add_prm("sh2", 4, -bn_max[4], bn_max[4], dbn[4]);
    bn_prms.add_prm("of1", 4, -bn_max[4], bn_max[4], dbn[4]);

    bn_prms.add_prm("sf1", 4, -bn_max[4], bn_max[4], dbn[4]);
    bn_prms.add_prm("sd1", 4, -bn_max[4], bn_max[4], dbn[4]);
    bn_prms.add_prm("sd2", 4, -bn_max[4], bn_max[4], dbn[4]);

    bn_prms.add_prm("s",   3, -bn_max[3], bn_max[3], dbn[3]);
    bn_prms.add_prm("sh2", 3, -bn_max[3], bn_max[3], dbn[3]);

    bn_prms.add_prm("sf1", 3, -bn_max[3], bn_max[3], dbn[3]);
    bn_prms.add_prm("sd1", 3, -bn_max[3], bn_max[3], dbn[3]);
    bn_prms.add_prm("sd2", 3, -bn_max[3], bn_max[3], dbn[3]);
    break;
  case 11:
    // ALS-U.
    // bn_prms.add_prm("sf", 3, -bn_max[3], bn_max[3], dbn[3]);
    // bn_prms.add_prm("sd", 3, -bn_max[3], bn_max[3], dbn[3]);

    bn_prms.add_prm("sh1", 3, -bn_max[3], bn_max[3], dbn[3]);
    bn_prms.add_prm("sh2", 3, -bn_max[3], bn_max[3], dbn[3]);
    bn_prms.add_prm("sh3", 3, -bn_max[3], bn_max[3], dbn[3]);
    bn_prms.add_prm("sh4", 3, -bn_max[3], bn_max[3], dbn[3]);
    break;
  default:
    printf("\nlat_select: unknown case\n");
    break;
  }
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

  if (false) {
    no_mpoles(3);
    no_mpoles(4);
    m_c(1000);
    exit(0);
  }

  if (false) {
    no_mpoles(3);
    bn_prms.add_prm("sf1", 3, -5e3, 5e3, 1e-2);
    bn_prms.add_prm("sd1", 3, -5e3, 5e3, 1e-2);
    bn_prms.add_prm("sd2", 3, -5e3, 5e3, 1e-2);
    fit_ksi1(0.0, 0.0, bn_prms.Fnum);
    bn_prms.prt_bn_lat("ksi1.out");
    exit(0);
  }

  lat_select();

  if (false) {
    const double bn[] =
      {2.023e+01, -8.931e+01, 3.101e+02, -2.322e+02, -2.930e+02};

    for (j = 0; j < 5; j++)
      set_bn(bn_prms.Fnum[j], bn_prms.n[j], bn[j]);

    bn_prms.prt_bn_lat("bn_ini.out");
    // exit(0);
  }

  if (false) no_mpoles(3);

  if (!true)
    conj_grad(bn_prms, f_nl, df_nl);
  else
    powell(bn_prms, f_nl);
}
