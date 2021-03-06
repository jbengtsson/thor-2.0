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

std::vector<std::string> drv_term;
int                      n_cell, n_cut;
double                   chi2 = 0e0, *f_lm, **A_lm;
tps                      h_re, h_im, h_re_scl, h_im_scl, K_re, K_im;
tps                      K_re_scl, K_re_delta_scl;
ss_vect<tps>             nus, nus_scl;

const bool
  fit_ksi = false,
  symm    = true,
  scale   = false,
  c_g     = true;

const double tpsa_eps = 1e-30;

// MAX-V         1,
// SLS-2         2,
// DIAMOND       3,
// DIAMOND & VMX 4,
// M-H6BAi       5,
// M-H6BA        6,
// ALS-U         7.

const int
  lat_case = 6,
  n_prt    = 8;

// Center of straight.
const double
  beta_inj[]   = {8.1, 2.3},
  twoJ[]       = {sqr(8e-3)/beta_inj[X_], sqr(2e-3)/beta_inj[Y_]},
  twoJ_delta[] = {sqr(0.5e-3)/beta_inj[X_], sqr(0.1e-3)/beta_inj[Y_]},
  delta_max    = 2e-2;


const double
  scl_h[]            = {0e-1,  0e-2, 0e-2},
  scl_dnu[]          = {1e-2, 1e-2, 1e-2, 1e-2},
  scl_ksi[]          = {1e5,  1e-2, 1e-2, 1e-2, 1e-2},
  delta_scl          = 0e0,
  scl_dnu_conf       = 0e0,
  scl_dnu_delta_conf = 0e0;


struct param_type {
private:

public:
  int                 m_constr, n_prm, svd_n_cut;
  double              bn_tol, step;
  double              *bn_lim, *bn, *dbn;
  std::vector<double> bn_max, bn_scl;
  std::vector<int>    Fnum, n, svd_list;

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
ss_vect<tps> Id_scl, Id_delta_scl;

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
    }

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

static double xsav_2D;
static tps (*func_save_2D)(const double, const double);

tps gauss_quad_2D(tps (*func)(const double), const double a, const double b)
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
		    gauss_quad_3D_z0(xsav_3D, y), gauss_quad_3D_z1(xsav_3D, y));
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


void get_twiss(void)
{
  int      j, k;
  double   alpha1[2], beta1[2], eta1[2], etap1[2], dnu1[2], dnu2[2];
  std::ofstream outf;

  const std::string file_name = "linlat.out";

  file_wr(outf, file_name.c_str());

  danot_(1); get_Map(); danot_(2);
  K = MapNorm(Map, g, A1, A0, Map_res, 1);
  get_ab(A1, alpha1, beta1, dnu1, eta1, etap1);

  // Crucial; to only store linear part of A.
  danot_(1);

  for (k = 0; k < 2; k++)
    dnu1[k] = 0e0;

  for (j = 1; j <= n_elem; j++) {
    A1.propagate(j, j);
    elem_tps[j-1].A1 = get_A_CS(2, A1, dnu2);
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

    outf << std::fixed << std::setprecision(4)
	 << std::setw(7) << elem[j-1].S
	 << std::setprecision(4)
	 << std::setw(9) << elem[j-1].Alpha[X_]
	      << std::setw(8) << elem[j-1].Beta[X_]
	 << std::setw(8) << elem[j-1].Nu[X_]
	 << std::setw(9) << elem[j-1].Eta[X_]
	      << std::setw(9) << elem[j-1].Etap[X_]
	 << std::setw(9) << elem[j-1].Alpha[Y_]
	      << std::setw(8) << elem[j-1].Beta[Y_]
	 << std::setw(8) << elem[j-1].Nu[Y_]
	 << std::setw(9) << elem[j-1].Eta[Y_]
	      << std::setw(9) << elem[j-1].Etap[Y_]
	 << std::endl;
  }

  outf.close();
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

  danot_(2); get_Map(); danot_(3);
  K = MapNorm(Map, g, A1, A0, Map_res, 1); nus = dHdJ(K);
  CtoR(K, K_re, K_im); CtoR(get_h(), h_re, h_im);

  printf("\nnu  = [%8.5f, %8.5f]\n", nus[0].cst(), nus[1].cst());
  printf("ksi = [%8.5f, %8.5f]\n",
	 h_ijklm(nus[0], 0, 0, 0, 0, 1),
	 h_ijklm(nus[1], 0, 0, 0, 0, 1));
  // printf("\nnu  = [%8.5f, %8.5f]\n",
  // 	 -h_ijklm(K_re, 1, 1, 0, 0, 0)/M_PI,
  // 	 -h_ijklm(K_re, 0, 0, 1, 1, 0)/M_PI);
  // printf("ksi = [%8.5f, %8.5f]\n",
  // 	 -h_ijklm(h_re, 1, 1, 0, 0, 1)/M_PI,
  // 	 -h_ijklm(h_re, 0, 0, 1, 1, 1)/M_PI);
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

  // Remove numeric noise.
  daeps_(1e-20);

  file_wr(outf, "h.out");
  outf << h_re*Id_scl;
  outf.close();

  file_wr(outf, "K.out");
  outf << K_re*Id_scl;
  outf.close();

  file_wr(outf, "nus.out");
  nus = dHdJ(K);
  daeps_(1e-5);
  outf << nus[3]*Id_scl << nus[4]*Id_scl;
  outf.close();

  daeps_(tpsa_eps);
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
  bool     first = true;
  long int loc;
  int      j, k, n, n_prm;
  double   bn;
  FILE     *outf;

  const std::string file_name = "dnu.out";

  outf = file_write(file_name.c_str());

  n_prm = bn_prms.n_prm;
  fprintf(outf, "\n");
  for (k = 0; k < n_prm; k++) {
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


double get_b(const std::string &label, const double scale, const tps &t,
	     const int i, const int j, const int k, const int l, const int m)
{
  drv_term.push_back(label);
  return scale*(h_ijklm(t, i, j, k, l, m));
}


void fit_ksi1(const double ksi_x, const double ksi_y)
{
  int    n_bn, i, j, m;
  double **A, *b, L;

  const int    m_max = 2;
  const double s_cut = 1e-10;

  n_bn = bn_prms.n_prm;

  b = dvector(1, m_max); A = dmatrix(1, m_max, 1, n_bn);

  no_mpoles(Sext); no_mpoles(Oct);

  printf("\n");
  for (i = 1; i <= n_bn; i++) {
    bn_prms.set_prm_dep(i-1);

    danot_(3);
    get_Map();
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
  b[++m] = -(get_b("ksi1_x", 1e0, nus[3], 0, 0, 0, 0, 1)-ksi_x);
  b[++m] = -(get_b("ksi1_y", 1e0, nus[4], 0, 0, 0, 0, 1)-ksi_y);

  prt_system(m, n_bn, A, b);

  SVD_lim(m, n_bn, A, b, bn_prms.bn_lim, s_cut, bn_prms.bn, bn_prms.dbn);

  bn_prms.set_dprm();

  printf("\nfit ksi (integrates strength in parenthesis):\n");
  for (i = 1; i <= n_bn; i++) {
    L = get_L(bn_prms.Fnum[i-1], 1);
    printf(" %12.5e", bn_prms.bn_scl[i-1]*bn_prms.bn[i]);
    if (i % n_prt == 0) printf("\n");
  }
  if (n_bn % n_prt != 0) printf("\n");
  printf("\nItegrated strenghts:\n");
  for (i = 1; i <= n_bn; i++) {
    L = get_L(bn_prms.Fnum[i-1], 1);
    printf(" %9.5f", bn_prms.bn_scl[i-1]*bn_prms.bn[i]*L);
    if (i % n_prt == 0) printf("\n");
  }
  if (n_bn % n_prt != 0) printf("\n");

  prt_mfile("flat_file.fit");
  prt_bn(bn_prms);

  danot_(3);
  get_Map();
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
  tps                 dnu, dnu_delta;
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
  K = MapNorm(Map, g, A1, A0, Map_res, 1);
  nus = dHdJ(K); nus_scl = nus*Id_scl;
  CtoR(K, K_re, K_im);
  K_re_scl = K_re*Id_scl; K_re_delta_scl = K_re*Id_delta_scl;
  CtoR(get_h(), h_re, h_im); h_re_scl = h_re*Id_scl; h_im_scl = h_im*Id_scl;

  dnu = gauss_quad_2D(f_gauss_quad_2D, 0e0, twoJ[X_]);
  dnu_delta = gauss_quad_3D(f_gauss_quad_3D, 0e0, twoJ[X_]);
  printf("\n|dnu|       = %9.3e\n", dnu.cst());
  printf("|dnu_delta| = %9.3e\n", dnu_delta.cst());

  if (scl_h[0] != 0e0) {
    b.push_back(get_b("", scl_h[0], h_re_scl, 1, 0, 0, 0, 2));
    b.push_back(get_b("", scl_h[0], h_re_scl, 2, 0, 0, 0, 1));
    b.push_back(get_b("", scl_h[0], h_re_scl, 0, 0, 2, 0, 1));

    if (!symm) {
      b.push_back(get_b("", scl_h[0], h_im_scl, 1, 0, 0, 0, 2));
      b.push_back(get_b("", scl_h[0], h_im_scl, 2, 0, 0, 0, 1));
      b.push_back(get_b("", scl_h[0], h_im_scl, 0, 0, 2, 0, 1));
    }

    b.push_back(get_b("", scl_h[0], h_re_scl, 1, 0, 1, 1, 0));
    b.push_back(get_b("", scl_h[0], h_re_scl, 2, 1, 0, 0, 0));
    b.push_back(get_b("", scl_h[0], h_re_scl, 3, 0, 0, 0, 0));
    b.push_back(get_b("", scl_h[0], h_re_scl, 1, 0, 0, 2, 0));
    b.push_back(get_b("", scl_h[0], h_re_scl, 1, 0, 2, 0, 0));

    if (!symm) {
      b.push_back(get_b("", scl_h[0], h_im_scl, 1, 0, 1, 1, 0));
      b.push_back(get_b("", scl_h[0], h_im_scl, 2, 1, 0, 0, 0));
      b.push_back(get_b("", scl_h[0], h_im_scl, 3, 0, 0, 0, 0));
      b.push_back(get_b("", scl_h[0], h_im_scl, 1, 0, 0, 2, 0));
      b.push_back(get_b("", scl_h[0], h_im_scl, 1, 0, 2, 0, 0));
    }
  }

  if (scl_h[1] != 0e0) {
    if (NO >= 5) {
      b.push_back(get_b("", scl_h[1], h_re_scl, 2, 0, 1, 1, 0));
      b.push_back(get_b("", scl_h[1], h_re_scl, 3, 1, 0, 0, 0));
      b.push_back(get_b("", scl_h[1], h_re_scl, 4, 0, 0, 0, 0));
      b.push_back(get_b("", scl_h[1], h_re_scl, 2, 0, 0, 2, 0));
      b.push_back(get_b("", scl_h[1], h_re_scl, 2, 0, 2, 0, 0));
      b.push_back(get_b("", scl_h[1], h_re_scl, 0, 0, 4, 0, 0));
      b.push_back(get_b("", scl_h[1], h_re_scl, 0, 0, 3, 1, 0));
      b.push_back(get_b("", scl_h[1], h_re_scl, 1, 1, 2, 0, 0));

      if (!symm) {
	b.push_back(get_b("", scl_h[1], h_im_scl, 2, 0, 1, 1, 0));
	b.push_back(get_b("", scl_h[1], h_im_scl, 3, 1, 0, 0, 0));
	b.push_back(get_b("", scl_h[1], h_im_scl, 4, 0, 0, 0, 0));
	b.push_back(get_b("", scl_h[1], h_im_scl, 2, 0, 0, 2, 0));
	b.push_back(get_b("", scl_h[1], h_im_scl, 2, 0, 2, 0, 0));
	b.push_back(get_b("", scl_h[1], h_im_scl, 0, 0, 4, 0, 0));
	b.push_back(get_b("", scl_h[1], h_im_scl, 0, 0, 3, 1, 0));
	b.push_back(get_b("", scl_h[1], h_im_scl, 1, 1, 2, 0, 0));
      }
    }
  }

  if (scl_h[2] != 0e0) {
    if (NO >= 6) {
      b.push_back(get_b("", scl_h[2], h_re_scl, 5, 0, 0, 0, 0));
      b.push_back(get_b("", scl_h[2], h_re_scl, 4, 1, 0, 0, 0));
      b.push_back(get_b("", scl_h[2], h_re_scl, 3, 2, 0, 0, 0));
      b.push_back(get_b("", scl_h[2], h_re_scl, 3, 0, 2, 0, 0));
      b.push_back(get_b("", scl_h[2], h_re_scl, 2, 1, 2, 0, 0));
      b.push_back(get_b("", scl_h[2], h_re_scl, 2, 1, 0, 2, 0));
      b.push_back(get_b("", scl_h[2], h_re_scl, 3, 0, 0, 2, 0));
      b.push_back(get_b("", scl_h[2], h_re_scl, 1, 0, 4, 0, 0));
      b.push_back(get_b("", scl_h[2], h_re_scl, 1, 0, 0, 4, 0));
      b.push_back(get_b("", scl_h[2], h_re_scl, 3, 0, 1, 1, 0));
      b.push_back(get_b("", scl_h[2], h_re_scl, 2, 1, 1, 1, 0));
      b.push_back(get_b("", scl_h[2], h_re_scl, 1, 0, 3, 1, 0));
      b.push_back(get_b("", scl_h[2], h_re_scl, 1, 0, 1, 3, 0));
      b.push_back(get_b("", scl_h[2], h_re_scl, 1, 0, 2, 2, 0));

      if (!symm) {
	b.push_back(get_b("", scl_h[2], h_im_scl, 5, 0, 0, 0, 0));
	b.push_back(get_b("", scl_h[2], h_im_scl, 4, 1, 0, 0, 0));
	b.push_back(get_b("", scl_h[2], h_im_scl, 3, 2, 0, 0, 0));
	b.push_back(get_b("", scl_h[2], h_im_scl, 3, 0, 2, 0, 0));
	b.push_back(get_b("", scl_h[2], h_im_scl, 2, 1, 2, 0, 0));
	b.push_back(get_b("", scl_h[2], h_im_scl, 2, 1, 0, 2, 0));
	b.push_back(get_b("", scl_h[2], h_im_scl, 3, 0, 0, 2, 0));
	b.push_back(get_b("", scl_h[2], h_im_scl, 1, 0, 4, 0, 0));
	b.push_back(get_b("", scl_h[2], h_im_scl, 1, 0, 0, 4, 0));
	b.push_back(get_b("", scl_h[2], h_im_scl, 3, 0, 1, 1, 0));
	b.push_back(get_b("", scl_h[2], h_im_scl, 2, 1, 1, 1, 0));
	b.push_back(get_b("", scl_h[2], h_im_scl, 1, 0, 3, 1, 0));
	b.push_back(get_b("", scl_h[2], h_im_scl, 1, 0, 1, 3, 0));
	b.push_back(get_b("", scl_h[2], h_im_scl, 1, 0, 2, 2, 0));
      }
    }
  }

  if (scl_ksi[0] != 0e0) {
    b.push_back(get_b("", scl_ksi[0], K_re_scl, 1, 1, 0, 0, 1));
    b.push_back(get_b("", scl_ksi[0], K_re_scl, 0, 0, 1, 1, 1));
  }

  if (NO >= 5) {
    if (scl_dnu[0] != 0e0) {
      b.push_back(get_b("", scl_dnu[0], K_re_scl, 2, 2, 0, 0, 0));
      b.push_back(get_b("", scl_dnu[0], K_re_scl, 1, 1, 1, 1, 0));
      b.push_back(get_b("", scl_dnu[0], K_re_scl, 0, 0, 2, 2, 0));
    }

    if (scl_ksi[2] != 0e0) {
      b.push_back(get_b("", scl_ksi[2], K_re_scl, 1, 1, 0, 0, 2));
      b.push_back(get_b("", scl_ksi[2], K_re_scl, 0, 0, 1, 1, 2));
    }

    if (scl_dnu_conf != 0e0)
      b.push_back(get_b("", scl_dnu_conf, dnu, 0, 0, 0, 0, 0));
    if (scl_dnu_delta_conf != 0e0)
      b.push_back(get_b("", scl_dnu_delta_conf, dnu_delta, 0, 0, 0, 0, 0));
  }

  if (NO >= 6) {
    if (delta_scl*scl_ksi[1] != 0e0) {
      b.push_back(get_b("", scl_ksi[1], K_re_delta_scl, 2, 2, 0, 0, 1));
      b.push_back(get_b("", scl_ksi[1], K_re_delta_scl, 1, 1, 1, 1, 1));
      b.push_back(get_b("", scl_ksi[1], K_re_delta_scl, 0, 0, 2, 2, 1));
    }

    if (scl_ksi[3] != 0e0) {
      b.push_back(get_b("", scl_ksi[3], K_re_scl, 1, 1, 0, 0, 3));
      b.push_back(get_b("", scl_ksi[3], K_re_scl, 0, 0, 1, 1, 3));
    }
  }

  if (NO >= 7) {
    if (scl_dnu[1] != 0e0) {
      b.push_back(get_b("", scl_dnu[1], K_re_scl, 3, 3, 0, 0, 0));
      b.push_back(get_b("", scl_dnu[1], K_re_scl, 2, 2, 1, 1, 0));
      b.push_back(get_b("", scl_dnu[1], K_re_scl, 1, 1, 2, 2, 0));
      b.push_back(get_b("", scl_dnu[1], K_re_scl, 0, 0, 3, 3, 0));
    }

    if (delta_scl*scl_ksi[2] != 0e0) {
      b.push_back(get_b("", scl_ksi[2], K_re_delta_scl, 2, 2, 0, 0, 2));
      b.push_back(get_b("", scl_ksi[2], K_re_delta_scl, 1, 1, 1, 1, 2));
      b.push_back(get_b("", scl_ksi[2], K_re_delta_scl, 0, 0, 2, 2, 2));
    }

    if (scl_ksi[4] != 0e0) {
      b.push_back(get_b("", scl_ksi[4], K_re_scl, 1, 1, 0, 0, 4));
      b.push_back(get_b("", scl_ksi[4], K_re_scl, 0, 0, 1, 1, 4));
    }
  }

  if (NO >= 8) {
    if (delta_scl*scl_ksi[1] != 0e0) {
      b.push_back(get_b("", scl_ksi[1], K_re_delta_scl, 3, 3, 0, 0, 1));
      b.push_back(get_b("", scl_ksi[1], K_re_delta_scl, 2, 2, 1, 1, 1));
      b.push_back(get_b("", scl_ksi[1], K_re_delta_scl, 1, 1, 2, 2, 1));
      b.push_back(get_b("", scl_ksi[1], K_re_delta_scl, 0, 0, 3, 3, 1));
    }

    if (delta_scl*scl_ksi[3] != 0e0) {
      b.push_back(get_b("", scl_ksi[3], K_re_delta_scl, 2, 2, 0, 0, 3));
      b.push_back(get_b("", scl_ksi[3], K_re_delta_scl, 1, 1, 1, 1, 3));
      b.push_back(get_b("", scl_ksi[3], K_re_delta_scl, 0, 0, 2, 2, 3));
    }
  }

  if (NO >= 9) {
    if (scl_dnu[2] != 0e0) {
      b.push_back(get_b("", scl_dnu[2], K_re_scl, 4, 4, 0, 0, 0));
      b.push_back(get_b("", scl_dnu[2], K_re_scl, 3, 3, 1, 1, 0));
      b.push_back(get_b("", scl_dnu[2], K_re_scl, 2, 2, 2, 2, 0));
      b.push_back(get_b("", scl_dnu[2], K_re_scl, 1, 1, 3, 3, 0));
      b.push_back(get_b("", scl_dnu[2], K_re_scl, 0, 0, 4, 4, 0));
    }

    if (delta_scl*scl_ksi[2] != 0e0) {
      b.push_back(get_b("", scl_ksi[2], K_re_delta_scl, 3, 3, 0, 0, 2));
      b.push_back(get_b("", scl_ksi[2], K_re_delta_scl, 2, 2, 1, 1, 2));
      b.push_back(get_b("", scl_ksi[2], K_re_delta_scl, 1, 1, 2, 2, 2));
      b.push_back(get_b("", scl_ksi[2], K_re_delta_scl, 0, 0, 3, 3, 2));
    }
  }

  if (NO >= 11) {
    if (scl_dnu[3] != 0e0) {
      b.push_back(get_b("", scl_dnu[3], K_re_scl, 5, 5, 0, 0, 0));
      b.push_back(get_b("", scl_dnu[3], K_re_scl, 4, 4, 1, 1, 0));
      b.push_back(get_b("", scl_dnu[3], K_re_scl, 3, 3, 2, 2, 0));
      b.push_back(get_b("", scl_dnu[3], K_re_scl, 2, 2, 3, 3, 0));
      b.push_back(get_b("", scl_dnu[3], K_re_scl, 1, 1, 4, 4, 0));
      b.push_back(get_b("", scl_dnu[3], K_re_scl, 0, 0, 5, 5, 0));
    }
  }

  chi2 = 0e0;
  for (i = 0; i < (int)b.size(); i++)
    chi2 += sqr(b[i]);

  if (chi2 < chi2_ref) {
    prt_bn(bn_prms);

    printf("\n%3d %12.5e -> %12.5e\n", n_powell, chi2_ref, chi2);
    for (i = 1; i <= bn_prms.n_prm; i++) {
      printf("%11.3e", bns[i]);
      if (i % n_prt == 0) printf("\n");
    }
    if (bn_prms.n_prm % n_prt != 0) printf("\n");
  }

  chi2_ref = min(chi2, chi2_ref);

  return chi2;
}


void get_f_grad(const int n_bn, double *f, double **A, double &chi2, int &m)
{
  int    i, j;
  tps    dnu, dnu_delta;

  // printf("\n");
  for (i = 1; i <= n_bn; i++) {
    bn_prms.set_prm_dep(i-1);

    danot_(NO-1);
    get_map_n(n_cell);
    danot_(NO);
    K = MapNorm(Map, g, A1, A0, Map_res, 1);
    nus = dHdJ(K); nus_scl = nus*Id_scl;
    CtoR(K, K_re, K_im);
    K_re_scl = K_re*Id_scl; K_re_delta_scl = K_re*Id_delta_scl;
    CtoR(get_h(), h_re, h_im); h_re_scl = h_re*Id_scl; h_im_scl = h_im*Id_scl;

    dnu = gauss_quad_2D(f_gauss_quad_2D, 0e0, twoJ[X_]);
    dnu_delta = gauss_quad_3D(f_gauss_quad_3D, 0e0, twoJ[X_]);
    // std::cout << std::scientific << std::setprecision(3)
    // 	      << "\n |dnu| = " << dnu << "\n";

    m = 0;
    if (scl_h[0] != 0e0) {
      A[++m][i] = get_a(scl_h[0], h_re_scl, 1, 0, 0, 0, 2);
      A[++m][i] = get_a(scl_h[0], h_re_scl, 2, 0, 0, 0, 1);
      A[++m][i] = get_a(scl_h[0], h_re_scl, 0, 0, 2, 0, 1);

      if (!symm) {
	A[++m][i] = get_a(scl_h[0], h_im_scl, 1, 0, 0, 0, 2);
	A[++m][i] = get_a(scl_h[0], h_im_scl, 2, 0, 0, 0, 1);
	A[++m][i] = get_a(scl_h[0], h_im_scl, 0, 0, 2, 0, 1);
      }

      A[++m][i] = get_a(scl_h[0], h_re_scl, 1, 0, 1, 1, 0);
      A[++m][i] = get_a(scl_h[0], h_re_scl, 2, 1, 0, 0, 0);
      A[++m][i] = get_a(scl_h[0], h_re_scl, 3, 0, 0, 0, 0);
      A[++m][i] = get_a(scl_h[0], h_re_scl, 1, 0, 0, 2, 0);
      A[++m][i] = get_a(scl_h[0], h_re_scl, 1, 0, 2, 0, 0);

      if (!symm) {
	A[++m][i] = get_a(scl_h[0], h_im_scl, 1, 0, 1, 1, 0);
	A[++m][i] = get_a(scl_h[0], h_im_scl, 2, 1, 0, 0, 0);
	A[++m][i] = get_a(scl_h[0], h_im_scl, 3, 0, 0, 0, 0);
	A[++m][i] = get_a(scl_h[0], h_im_scl, 1, 0, 0, 2, 0);
	A[++m][i] = get_a(scl_h[0], h_im_scl, 1, 0, 2, 0, 0);
      }
    }

    if (scl_h[1] != 0e0) {
      if (NO >= 5) {
	A[++m][i] = get_a(scl_h[1], h_re_scl, 2, 0, 1, 1, 0);
	A[++m][i] = get_a(scl_h[1], h_re_scl, 3, 1, 0, 0, 0);
	A[++m][i] = get_a(scl_h[1], h_re_scl, 4, 0, 0, 0, 0);
	A[++m][i] = get_a(scl_h[1], h_re_scl, 2, 0, 0, 2, 0);
	A[++m][i] = get_a(scl_h[1], h_re_scl, 2, 0, 2, 0, 0);
	A[++m][i] = get_a(scl_h[1], h_re_scl, 0, 0, 4, 0, 0);
	A[++m][i] = get_a(scl_h[1], h_re_scl, 0, 0, 3, 1, 0);
	A[++m][i] = get_a(scl_h[1], h_re_scl, 1, 1, 2, 0, 0);

	if (!symm) {
	  A[++m][i] = get_a(scl_h[1], h_im_scl, 2, 0, 1, 1, 0);
	  A[++m][i] = get_a(scl_h[1], h_im_scl, 3, 1, 0, 0, 0);
	  A[++m][i] = get_a(scl_h[1], h_im_scl, 4, 0, 0, 0, 0);
	  A[++m][i] = get_a(scl_h[1], h_im_scl, 2, 0, 0, 2, 0);
	  A[++m][i] = get_a(scl_h[1], h_im_scl, 2, 0, 2, 0, 0);
	  A[++m][i] = get_a(scl_h[1], h_im_scl, 0, 0, 4, 0, 0);
	  A[++m][i] = get_a(scl_h[1], h_im_scl, 0, 0, 3, 1, 0);
	  A[++m][i] = get_a(scl_h[1], h_im_scl, 1, 1, 2, 0, 0);
	}
      }
    }

    if (scl_h[2] != 0e0) {
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

	if (!symm) {
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
    }

    if (scl_ksi[0] != 0e0) {
      A[++m][i] = get_a(scl_ksi[0], K_re_scl, 1, 1, 0, 0, 1);
      A[++m][i] = get_a(scl_ksi[0], K_re_scl, 0, 0, 1, 1, 1);
    }

    if (NO >= 5) {
      if (scl_dnu[0] != 0e0) {
	A[++m][i] = get_a(scl_dnu[0], K_re_scl, 2, 2, 0, 0, 0);
	A[++m][i] = get_a(scl_dnu[0], K_re_scl, 1, 1, 1, 1, 0);
	A[++m][i] = get_a(scl_dnu[0], K_re_scl, 0, 0, 2, 2, 0);
      }

      if (scl_ksi[2] != 0e0) {
	A[++m][i] = get_a(scl_ksi[2], K_re_scl, 1, 1, 0, 0, 2);
	A[++m][i] = get_a(scl_ksi[2], K_re_scl, 0, 0, 1, 1, 2);
      }

      if (scl_dnu_conf != 0e0)
	A[++m][i] = get_a(scl_dnu_conf, dnu, 0, 0, 0, 0, 0);
      if (scl_dnu_delta_conf != 0e0)
	A[++m][i] = get_a(scl_dnu_delta_conf, dnu_delta, 0, 0, 0, 0, 0);
    }

    if (NO >= 6) {
      if (delta_scl*scl_ksi[1] != 0e0) {
	A[++m][i] = get_a(scl_ksi[1], K_re_delta_scl, 2, 2, 0, 0, 1);
	A[++m][i] = get_a(scl_ksi[1], K_re_delta_scl, 1, 1, 1, 1, 1);
	A[++m][i] = get_a(scl_ksi[1], K_re_delta_scl, 0, 0, 2, 2, 1);
      }

      if (scl_ksi[3] != 0e0) {
	A[++m][i] = get_a(scl_ksi[3], K_re_scl, 1, 1, 0, 0, 3);
	A[++m][i] = get_a(scl_ksi[3], K_re_scl, 0, 0, 1, 1, 3);
      }
    }

    if (NO >= 7) {
      if (scl_dnu[1] != 0e0) {
	A[++m][i] = get_a(scl_dnu[1], K_re_scl, 3, 3, 0, 0, 0);
	A[++m][i] = get_a(scl_dnu[1], K_re_scl, 2, 2, 1, 1, 0);
	A[++m][i] = get_a(scl_dnu[1], K_re_scl, 1, 1, 2, 2, 0);
	A[++m][i] = get_a(scl_dnu[1], K_re_scl, 0, 0, 3, 3, 0);
      }

      if (delta_scl*scl_ksi[2] != 0e0) {
	A[++m][i] = get_a(scl_ksi[2], K_re_delta_scl, 2, 2, 0, 0, 2);
	A[++m][i] = get_a(scl_ksi[2], K_re_delta_scl, 1, 1, 1, 1, 2);
	A[++m][i] = get_a(scl_ksi[2], K_re_delta_scl, 0, 0, 2, 2, 2);
      }

      if (scl_ksi[4] != 0e0) {
	A[++m][i] = get_a(scl_ksi[4], K_re_scl, 1, 1, 0, 0, 4);
	A[++m][i] = get_a(scl_ksi[4], K_re_scl, 0, 0, 1, 1, 4);
      }
    }

    if (NO >= 8) {
      if (delta_scl*scl_ksi[1] != 0e0) {
	A[++m][i] = get_a(scl_ksi[1], K_re_delta_scl, 3, 3, 0, 0, 1);
	A[++m][i] = get_a(scl_ksi[1], K_re_delta_scl, 2, 2, 1, 1, 1);
	A[++m][i] = get_a(scl_ksi[1], K_re_delta_scl, 1, 1, 2, 2, 1);
	A[++m][i] = get_a(scl_ksi[1], K_re_delta_scl, 0, 0, 3, 3, 1);
      }

      if (delta_scl*scl_ksi[3] != 0e0) {
	A[++m][i] = get_a(scl_ksi[3], K_re_delta_scl, 2, 2, 0, 0, 3);
	A[++m][i] = get_a(scl_ksi[3], K_re_delta_scl, 1, 1, 1, 1, 3);
	A[++m][i] = get_a(scl_ksi[3], K_re_delta_scl, 0, 0, 2, 2, 3);
      }
    }

    if (NO >= 9) {
      if (scl_dnu[2] != 0e0) {
	A[++m][i] = get_a(scl_dnu[2], K_re_scl, 4, 4, 0, 0, 0);
	A[++m][i] = get_a(scl_dnu[2], K_re_scl, 3, 3, 1, 1, 0);
	A[++m][i] = get_a(scl_dnu[2], K_re_scl, 2, 2, 2, 2, 0);
	A[++m][i] = get_a(scl_dnu[2], K_re_scl, 1, 1, 3, 3, 0);
	A[++m][i] = get_a(scl_dnu[2], K_re_scl, 0, 0, 4, 4, 0);
      }

      if (delta_scl*scl_ksi[2] != 0e0) {
	A[++m][i] = get_a(scl_ksi[2], K_re_delta_scl, 3, 3, 0, 0, 2);
	A[++m][i] = get_a(scl_ksi[2], K_re_delta_scl, 2, 2, 1, 1, 2);
	A[++m][i] = get_a(scl_ksi[2], K_re_delta_scl, 1, 1, 2, 2, 2);
	A[++m][i] = get_a(scl_ksi[2], K_re_delta_scl, 0, 0, 3, 3, 2);
      }
    }

    if (NO >= 11) {
      if (scl_dnu[3] != 0e0) {
	A[++m][i] = get_a(scl_dnu[3], K_re_scl, 5, 5, 0, 0, 0);
	A[++m][i] = get_a(scl_dnu[3], K_re_scl, 4, 4, 1, 1, 0);
	A[++m][i] = get_a(scl_dnu[3], K_re_scl, 3, 3, 2, 2, 0);
	A[++m][i] = get_a(scl_dnu[3], K_re_scl, 2, 2, 3, 3, 0);
	A[++m][i] = get_a(scl_dnu[3], K_re_scl, 1, 1, 4, 4, 0);
	A[++m][i] = get_a(scl_dnu[3], K_re_scl, 0, 0, 5, 5, 0);
      }
    }

    for (j = 1; j <= m; j++)
      A[j][i] *= bn_prms.bn_scl[i-1];

    bn_prms.clr_prm_dep(i-1);
  }

  printf("\n dnu = %10.3e\n", dnu.cst());
  printf(" dnu = %10.3e\n", dnu_delta.cst());

  m = 0;
  if (scl_h[0] != 0e0) {
    f[++m] = get_b("h_re_10002", scl_h[0], h_re_scl, 1, 0, 0, 0, 2);
    f[++m] = get_b("h_re_20001", scl_h[0], h_re_scl, 2, 0, 0, 0, 1);
    f[++m] = get_b("h_re_00201", scl_h[0], h_re_scl, 0, 0, 2, 0, 1);

    if (!symm) {
      f[++m] = get_b("h_im_10002", scl_h[0], h_im_scl, 1, 0, 0, 0, 2);
      f[++m] = get_b("h_im_20001", scl_h[0], h_im_scl, 2, 0, 0, 0, 1);
      f[++m] = get_b("h_im_00201", scl_h[0], h_im_scl, 0, 0, 2, 0, 1);
    }

    f[++m] = get_b("h_re_10110", scl_h[0], h_re_scl, 1, 0, 1, 1, 0);
    f[++m] = get_b("h_re_21000", scl_h[0], h_re_scl, 2, 1, 0, 0, 0);
    f[++m] = get_b("h_re_30000", scl_h[0], h_re_scl, 3, 0, 0, 0, 0);
    f[++m] = get_b("h_re_10020", scl_h[0], h_re_scl, 1, 0, 0, 2, 0);
    f[++m] = get_b("h_re_10200", scl_h[0], h_re_scl, 1, 0, 2, 0, 0);

    if (!symm) {
      f[++m] = get_b("h_im_10110", scl_h[0], h_im_scl, 1, 0, 1, 1, 0);
      f[++m] = get_b("h_im_21000", scl_h[0], h_im_scl, 2, 1, 0, 0, 0);
      f[++m] = get_b("h_im_30000", scl_h[0], h_im_scl, 3, 0, 0, 0, 0);
      f[++m] = get_b("h_im_10020", scl_h[0], h_im_scl, 1, 0, 0, 2, 0);
      f[++m] = get_b("h_im_10200", scl_h[0], h_im_scl, 1, 0, 2, 0, 0);
    }
  }

  if (scl_h[1] != 0e0) {
    if (NO >= 5) {
      f[++m] = get_b("h_re_20110", scl_h[1], h_re_scl, 2, 0, 1, 1, 0);
      f[++m] = get_b("h_re_31000", scl_h[1], h_re_scl, 3, 1, 0, 0, 0);
      f[++m] = get_b("h_re_40000", scl_h[1], h_re_scl, 4, 0, 0, 0, 0);
      f[++m] = get_b("h_re_20020", scl_h[1], h_re_scl, 2, 0, 0, 2, 0);
      f[++m] = get_b("h_re_20200", scl_h[1], h_re_scl, 2, 0, 2, 0, 0);
      f[++m] = get_b("h_re_00400", scl_h[1], h_re_scl, 0, 0, 4, 0, 0);
      f[++m] = get_b("h_re_00310", scl_h[1], h_re_scl, 0, 0, 3, 1, 0);
      f[++m] = get_b("h_re_11200", scl_h[1], h_re_scl, 1, 1, 2, 0, 0);
  
      if (!symm) {
	f[++m] = get_b("h_im_20110", scl_h[1], h_im_scl, 2, 0, 1, 1, 0);
	f[++m] = get_b("h_im_31000", scl_h[1], h_im_scl, 3, 1, 0, 0, 0);
	f[++m] = get_b("h_im_40000", scl_h[1], h_im_scl, 4, 0, 0, 0, 0);
	f[++m] = get_b("h_im_20020", scl_h[1], h_im_scl, 2, 0, 0, 2, 0);
	f[++m] = get_b("h_im_20200", scl_h[1], h_im_scl, 2, 0, 2, 0, 0);
	f[++m] = get_b("h_im_00400", scl_h[1], h_im_scl, 0, 0, 4, 0, 0);
	f[++m] = get_b("h_im_00310", scl_h[1], h_im_scl, 0, 0, 3, 1, 0);
	f[++m] = get_b("h_im_11200", scl_h[1], h_im_scl, 1, 1, 2, 0, 0);
      }
    }
  }

  if (scl_h[2] != 0e0) {
    if (NO >= 6) {
      f[++m] = get_b("h_re_50000", scl_h[2], h_re_scl, 5, 0, 0, 0, 0);
      f[++m] = get_b("h_re_41000", scl_h[2], h_re_scl, 4, 1, 0, 0, 0);
      f[++m] = get_b("h_re_32000", scl_h[2], h_re_scl, 3, 2, 0, 0, 0);
      f[++m] = get_b("h_re_30200", scl_h[2], h_re_scl, 3, 0, 2, 0, 0);
      f[++m] = get_b("h_re_21200", scl_h[2], h_re_scl, 2, 1, 2, 0, 0);
      f[++m] = get_b("h_re_21020", scl_h[2], h_re_scl, 2, 1, 0, 2, 0);
      f[++m] = get_b("h_re_30020", scl_h[2], h_re_scl, 3, 0, 0, 2, 0);
      f[++m] = get_b("h_re_10400", scl_h[2], h_re_scl, 1, 0, 4, 0, 0);
      f[++m] = get_b("h_re_10040", scl_h[2], h_re_scl, 1, 0, 0, 4, 0);
      f[++m] = get_b("h_re_30110", scl_h[2], h_re_scl, 3, 0, 1, 1, 0);
      f[++m] = get_b("h_re_21110", scl_h[2], h_re_scl, 2, 1, 1, 1, 0);
      f[++m] = get_b("h_re_10310", scl_h[2], h_re_scl, 1, 0, 3, 1, 0);
      f[++m] = get_b("h_re_10130", scl_h[2], h_re_scl, 1, 0, 1, 3, 0);
      f[++m] = get_b("h_re_10220", scl_h[2], h_re_scl, 1, 0, 2, 2, 0);
  
      if (!symm) {
	f[++m] = get_b("h_im_50000", scl_h[2], h_im_scl, 5, 0, 0, 0, 0);
	f[++m] = get_b("h_im_41000", scl_h[2], h_im_scl, 4, 1, 0, 0, 0);
	f[++m] = get_b("h_im_32000", scl_h[2], h_im_scl, 3, 2, 0, 0, 0);
	f[++m] = get_b("h_im_30200", scl_h[2], h_im_scl, 3, 0, 2, 0, 0);
	f[++m] = get_b("h_im_21200", scl_h[2], h_im_scl, 2, 1, 2, 0, 0);
	f[++m] = get_b("h_im_21020", scl_h[2], h_im_scl, 2, 1, 0, 2, 0);
	f[++m] = get_b("h_im_30020", scl_h[2], h_im_scl, 3, 0, 0, 2, 0);
	f[++m] = get_b("h_im_10400", scl_h[2], h_im_scl, 1, 0, 4, 0, 0);
	f[++m] = get_b("h_im_10040", scl_h[2], h_im_scl, 1, 0, 0, 4, 0);
	f[++m] = get_b("h_im_30110", scl_h[2], h_im_scl, 3, 0, 1, 1, 0);
	f[++m] = get_b("h_im_21110", scl_h[2], h_im_scl, 2, 1, 1, 1, 0);
	f[++m] = get_b("h_im_10310", scl_h[2], h_im_scl, 1, 0, 3, 1, 0);
	f[++m] = get_b("h_im_10130", scl_h[2], h_im_scl, 1, 0, 1, 3, 0);
	f[++m] = get_b("h_im_10220", scl_h[2], h_im_scl, 1, 0, 2, 2, 0);
      }
    }
  }

  if (scl_ksi[0] != 0e0) {
    f[++m] = get_b("K_re_11001", scl_ksi[0], K_re_scl, 1, 1, 0, 0, 1);
    f[++m] = get_b("K_re_00111", scl_ksi[0], K_re_scl, 0, 0, 1, 1, 1);
  }

  if (NO >= 5) {
    if (scl_dnu[0] != 0e0) {
      f[++m] = get_b("K_re_22000", scl_dnu[0], K_re_scl, 2, 2, 0, 0, 0);
      f[++m] = get_b("K_re_11110", scl_dnu[0], K_re_scl, 1, 1, 1, 1, 0);
      f[++m] = get_b("K_re_00220", scl_dnu[0], K_re_scl, 0, 0, 2, 2, 0);
    }

    if (scl_ksi[2] != 0e0) {
      f[++m] = get_b("K_re_11002", scl_ksi[2], K_re_scl, 1, 1, 0, 0, 2);
      f[++m] = get_b("K_re_00112", scl_ksi[2], K_re_scl, 0, 0, 1, 1, 2);
    }

    if (scl_dnu_conf != 0e0) {
      f[++m] = get_b("|dnu|", scl_dnu_conf, dnu, 0, 0, 0, 0, 0);
    }
    if (scl_dnu_delta_conf != 0e0) {
      f[++m] =
	get_b("|dnu(delta)|", scl_dnu_delta_conf, dnu_delta, 0, 0, 0, 0, 0);
    }
  }

  if (NO >= 6) {
    if (delta_scl*scl_ksi[1] != 0e0) {
      f[++m] = get_b("K_re_22001", scl_ksi[1], K_re_delta_scl, 2, 2, 0, 0, 1);
      f[++m] = get_b("K_re_11111", scl_ksi[1], K_re_delta_scl, 1, 1, 1, 1, 1);
      f[++m] = get_b("K_re_00221", scl_ksi[1], K_re_delta_scl, 0, 0, 2, 2, 1);
    }

    if (scl_ksi[3] != 0e0) {
      f[++m] = get_b("K_re_11003", scl_ksi[3], K_re_scl, 1, 1, 0, 0, 3);
      f[++m] = get_b("K_re_00113", scl_ksi[3], K_re_scl, 0, 0, 1, 1, 3);
    }
  }

  if (NO >= 7) {
    if (scl_dnu[1] != 0e0) {
      f[++m] = get_b("K_re_33000", scl_dnu[1], K_re_scl, 3, 3, 0, 0, 0);
      f[++m] = get_b("K_re_22110", scl_dnu[1], K_re_scl, 2, 2, 1, 1, 0);
      f[++m] = get_b("K_re_11220", scl_dnu[1], K_re_scl, 1, 1, 2, 2, 0);
      f[++m] = get_b("K_re_00330", scl_dnu[1], K_re_scl, 0, 0, 3, 3, 0);
    }

    if (delta_scl*scl_ksi[2] != 0e0) {
      f[++m] = get_b("K_re_22002", scl_ksi[2], K_re_delta_scl, 2, 2, 0, 0, 2);
      f[++m] = get_b("K_re_11112", scl_ksi[2], K_re_delta_scl, 1, 1, 1, 1, 2);
      f[++m] = get_b("K_re_00222", scl_ksi[2], K_re_delta_scl, 0, 0, 2, 2, 2);
    }

    if (scl_ksi[4] != 0e0) {
      f[++m] = get_b("K_re_11004", scl_ksi[4], K_re_scl, 1, 1, 0, 0, 4);
      f[++m] = get_b("K_re_00114", scl_ksi[4], K_re_scl, 0, 0, 1, 1, 4);
    }
  }

  if (NO >= 8) {
    if (delta_scl*scl_ksi[1] != 0e0) {
      f[++m] = get_b("K_re_33001", scl_ksi[1], K_re_delta_scl, 3, 3, 0, 0, 1);
      f[++m] = get_b("K_re_22111", scl_ksi[1], K_re_delta_scl, 2, 2, 1, 1, 1);
      f[++m] = get_b("K_re_11221", scl_ksi[1], K_re_delta_scl, 1, 1, 2, 2, 1);
      f[++m] = get_b("K_re_00331", scl_ksi[1], K_re_delta_scl, 0, 0, 3, 3, 1);
    }

    if (delta_scl*scl_ksi[3] != 0e0) {
      f[++m] = get_b("K_re_22003", scl_ksi[3], K_re_delta_scl, 2, 2, 0, 0, 3);
      f[++m] = get_b("K_re_11113", scl_ksi[3], K_re_delta_scl, 1, 1, 1, 1, 3);
      f[++m] = get_b("K_re_00223", scl_ksi[3], K_re_delta_scl, 0, 0, 2, 2, 3);
    }
  }

  if (NO >= 9) {
    if (scl_dnu[2] != 0e0) {
      f[++m] = get_b("K_re_44000", scl_dnu[2], K_re_scl, 4, 4, 0, 0, 0);
      f[++m] = get_b("K_re_33110", scl_dnu[2], K_re_scl, 3, 3, 1, 1, 0);
      f[++m] = get_b("K_re_22220", scl_dnu[2], K_re_scl, 2, 2, 2, 2, 0);
      f[++m] = get_b("K_re_11330", scl_dnu[2], K_re_scl, 1, 1, 3, 3, 0);
      f[++m] = get_b("K_re_00440", scl_dnu[2], K_re_scl, 0, 0, 4, 4, 0);
    }

    if (delta_scl*scl_ksi[2] != 0e0) {
      f[++m] = get_b("K_re_33002", scl_ksi[2], K_re_delta_scl, 3, 3, 0, 0, 2);
      f[++m] = get_b("K_re_22112", scl_ksi[2], K_re_delta_scl, 2, 2, 1, 1, 2);
      f[++m] = get_b("K_re_11222", scl_ksi[2], K_re_delta_scl, 1, 1, 2, 2, 2);
      f[++m] = get_b("K_re_00332", scl_ksi[2], K_re_delta_scl, 0, 0, 3, 3, 2);
    }
  }

  if (NO >= 11) {
    if (scl_dnu[3] != 0e0) {
      f[++m] = get_b("K_re_55000", scl_dnu[3], K_re_scl, 5, 5, 0, 0, 0);
      f[++m] = get_b("K_re_44110", scl_dnu[3], K_re_scl, 4, 4, 1, 1, 0);
      f[++m] = get_b("K_re_33220", scl_dnu[3], K_re_scl, 3, 3, 2, 2, 0);
      f[++m] = get_b("K_re_22330", scl_dnu[3], K_re_scl, 2, 2, 3, 3, 0);
      f[++m] = get_b("K_re_11440", scl_dnu[3], K_re_scl, 1, 1, 4, 4, 0);
      f[++m] = get_b("K_re_00550", scl_dnu[3], K_re_scl, 0, 0, 5, 5, 0);
    }
  }

  prt_h_K();

  chi2 = 0e0;
  for (j = 1; j <= m; j++)
    chi2 += sqr(f[j]);
}


void Bubble_Sort2(std::vector<double> &w, std::vector<int> &order)
{
  bool   swapped;
  int    k, ind;
  double val;

  do {
    swapped = false;
    for (k = 0; k < (int)order.size()-1; k++) {
      if (fabs(w[k]) < fabs(w[k+1])) {
	ind = order[k]; order[k] = order[k+1]; order[k+1] = ind;
	val = w[k]; w[k] = w[k+1]; w[k+1] = val;
	swapped = true;
      }
    }
  } while (swapped);
}


void SVD_zero_n(const int n, double *w, const int svd_n)
{
  int                 j, ind;
  std::vector<int>    order;
  std::vector<double> w1;

  for (j = 1; j <= n; j++) {
    order.push_back(j); w1.push_back(w[j]);
  }
  Bubble_Sort2(w1, order);

  printf("\nzeroing singular values:\n");
  for (j = 1; j <= svd_n; j++) {
    ind = order[n-j];
    printf("%3d %3d %11.3e\n", j, ind, w[ind]);
    w[ind] = 0e0;
  }
  std::cout << "\n";
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

#if 0
  SVD_lim(m, n_bn, A, b, bn_prms.bn_lim, bn_prms.svd_n_cut, bn_prms.bn,
  	  bn_prms.dbn);
  // SVD_lim(m, n_bn, A, b, bn_prms.bn_lim, bn_prms.svd_list, bn_prms.bn,
  // 	  bn_prms.dbn);
#else
  double *w, **U, **V;

  w = dvector(1, n_bn);
  U = dmatrix(1, m, 1, n_bn); V = dmatrix(1, n_bn, 1, n_bn);

  dmcopy(A, m, n_bn, U); dsvdcmp(U, m, n_bn, w, V);

  std::cout << "\nsingular values:\n";
  for (i = 1; i <= n_bn; i++) {
    std::cout << std::scientific << std::setprecision(3)
	      << std::setw(11) << w[i];
    if (i % 10 == 0) std::cout << "\n";
  }
  if (n_bn % 10 != 0) std::cout << "\n";

  SVD_zero_n(n_bn, w, bn_prms.svd_n_cut);
  std::cout << "\n";

  dsvbksb(U, w, V, m, n_bn, b, bn_prms.dbn);

  std::cout << "dcorr.:" << std::endl;
  prt_vec(n_bn, bn_prms.dbn);

  free_dvector(w, 1, n_bn);
  free_dmatrix(U, 1, m, 1, n_bn); free_dmatrix(V, 1, n_bn, 1, n_bn);
#endif

  dvcopy(bn_prms.bn, n_bn, bn_ref);
  if (true)
    conj_grad(n_iter, bn_prms.bn, bn_prms.dbn, g_, h_, get_f);
  else
    bn_prms.set_dprm();

  printf("\nbn & dbn:\n");
  for (i = 1; i <= n_bn; i++) {
    printf(" %12.5e", bn_prms.bn_scl[i-1]*bn_prms.bn[i]);
    if (i % n_prt == 0) printf("\n");
  }
  if (n_bn % n_prt != 0) printf("\n");
  printf("\n");
  dbn_max = 0e0;
  for (i = 1; i <= n_bn; i++) {
    dbn_max = max(fabs((bn_prms.bn[i]-bn_ref[i])), dbn_max);
    printf(" %12.5e", bn_prms.bn[i]-bn_ref[i]);
    if (i % n_prt == 0) printf("\n");
  }
  if (n_bn % n_prt != 0) printf("\n");

  prt_bn(bn_prms);
  // bn_prms.set_prm();
  prt_mfile("flat_file.fit");

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

  } while ((dbn_max > bn_prms.bn_tol) && (n_iter < n_iter_max));
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
      get_bn(bn_prms.Fnum[i-1], 1, bn_prms.n[i-1])/bn_prms.bn_scl[i-1];
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

  const int n_data_max = 100;

  n_bn = bn_prms.n_prm;

  ia = ivector(1, n_bn);
  x = dvector(1, n_data_max); y = dvector(1, n_data_max);
  sigma = dvector(1, n_data_max); covar = dmatrix(1, n_bn, 1, n_bn);
  alpha = dmatrix(1, n_bn, 1, n_bn); f_lm = dvector(1, n_data_max);
  A_lm = dmatrix(1, n_data_max, 1, n_bn);

  get_f_grad(n_bn, f_lm, A_lm, chi2, n_data);
  if (n_data > n_data_max) {
    printf("\nn_data exceeded: %d (%d)\n", n_data, n_data_max);
    exit(0);
  }
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

  free_dvector(f_lm, 1, n_data_max); free_dmatrix(A_lm, 1, n_data_max, 1, n_bn);
  free_ivector(ia, 1, n_bn);
  free_dvector(x, 1, n_data_max); free_dvector(y, 1, n_data_max);
  free_dvector(sigma, 1, n_data_max);
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

  int           i, j, m, n_b2;
  double        **A, *b, *b2_lim, *b2, *db2, step;
  double        nu_fract[2], dnu[2];
  ss_vect<tps>  dnus;
  std::ofstream quad_out;

  const bool    debug = true;
  const int     m_max = 2;
  const double  step0 = 0.5, s_cut = 1e-10;

  b2_max = 10e0;

  n_b2 = b2_Fam.size();

  b = dvector(1, m_max); b2_lim = dvector(1, n_b2);
  b2 = dvector(1, n_b2); db2 = dvector(1, n_b2);
  A = dmatrix(1, m_max, 1, n_b2);

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

      m = 0;
      A[++m][i] = h_ijklm_p(nus[3], 0, 0, 0, 0, 0, 7);
      A[++m][i] = h_ijklm_p(nus[4], 0, 0, 0, 0, 0, 7);

      for (j = 1; j <= get_n_Kids(b2_Fam[i-1]); j++)
	clr_bn_par(b2_Fam[i-1], j, Quad);
    }

    m = 0;
    b[++m] = -dnu[X_];
    b[++m] = -dnu[Y_];

    SVD_lim(m, n_b2, A, b, b2_lim, s_cut, b2, db2);

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

  prt_mfile("flat_file.fit");
  prt_b2(b2_Fam);

  free_dvector(b, 1, m_max);
  free_dvector(b2_lim, 1, n_b2); free_dvector(b2, 1, n_b2);
  free_dvector(db2, 1, n_b2);
  free_dmatrix(A, 1, m_max, 1, n_b2);
}


void fit_tune(const double nu_x, const double nu_y)
{
  std::vector<int> b2_Fnum;

  const double eps = 1e-5;

  switch (lat_case) {
  case 6:
    // DIAMOND-II, RB-6-BA.
    // Low dispersion quadrupoles.
    b2_Fnum.push_back(get_Fnum("qd22"));
    b2_Fnum.push_back(get_Fnum("qf11"));
    b2_Fnum.push_back(get_Fnum("qd11"));
    b2_Fnum.push_back(get_Fnum("qd21"));

    b2_Fnum.push_back(get_Fnum("qf6"));
    b2_Fnum.push_back(get_Fnum("qf8"));

    b2_Fnum.push_back(get_Fnum("qd2"));
    b2_Fnum.push_back(get_Fnum("qf1"));
    break;
  }

  fit_tune(nu_x, nu_y, b2_Fnum, eps, true);
}


void lat_select(const int lat_case)
{
  switch (lat_case) {
  case 1:
    // MAX-V.
    n_cell = 1;

    // bn_prms.add_prm("sfh", 3, 5e5, 1.0);
    // bn_prms.add_prm("sd",  3, 5e5, 1.0);

    bn_prms.add_prm("o1",  4, 5e5, 1.0);
    bn_prms.add_prm("o2",  4, 5e5, 1.0);
    bn_prms.add_prm("o3",  4, 5e5, 1.0);
    bn_prms.add_prm("o4",  4, 5e5, 1.0);

    if (false) {
      bn_prms.add_prm("o1", 6, 5e10, 1.0);
      bn_prms.add_prm("o2", 6, 5e10, 1.0);
      bn_prms.add_prm("o3", 6, 5e10, 1.0);
      bn_prms.add_prm("o4", 6, 5e10, 1.0);
    }
    break;
  case 2:
    // SLS-2.
    n_cell = 1;

    // bn_prms.add_prm("sd", 3, 5e5, 1.0);
    // bn_prms.add_prm("sf", 3, 5e5, 1.0);

    if (!fit_ksi) {
      // bn_prms.add_prm("ocxm", 4, 5e5, 1.0);
      // bn_prms.add_prm("ocym", 4, 5e5, 1.0);

      // bn_prms.add_prm("ocxm", 5, 5e5, 1.0);
      // bn_prms.add_prm("ocym", 5, 5e5, 1.0);

      bn_prms.add_prm("oxx",  4, 5e5, 1.0);
      bn_prms.add_prm("oxy",  4, 5e5, 1.0);
      bn_prms.add_prm("oyy",  4, 5e5, 1.0);

      // bn_prms.add_prm("oxx",  6, 5e5, 1.0);
      // bn_prms.add_prm("oxy",  6, 5e5, 1.0);
      // bn_prms.add_prm("oyy",  6, 5e5, 1.0);
    }

    // bn_prms.add_prm("sdh", 3, 5e5, 1.0);
    // bn_prms.add_prm("sfh", 3, 5e5, 1.0);

    // if (!fit_ksi) {
    //   bn_prms.add_prm("sdh",  4, 5e5, 1.0);
    //   bn_prms.add_prm("sfh",  4, 5e5, 1.0);

    //   bn_prms.add_prm("sdh",  5, 5e5, 1.0);
    //   bn_prms.add_prm("sfh",  5, 5e5, 1.0);

    //   bn_prms.add_prm("sxxh", 3, 5e5, 1.0);
    //   bn_prms.add_prm("sxyh", 3, 5e5, 1.0);
    //   bn_prms.add_prm("syyh", 3, 5e5, 1.0);

    //   bn_prms.add_prm("sxxh", 4, 5e5, 1.0);
    //   bn_prms.add_prm("sxyh", 4, 5e5, 1.0);
    //   bn_prms.add_prm("syyh", 4, 5e5, 1.0);

    //   bn_prms.add_prm("sxxh", 6, 5e5, 1.0);
    //   bn_prms.add_prm("sxyh", 6, 5e5, 1.0);
    //   bn_prms.add_prm("syyh", 6, 5e5, 1.0);
    // }
    break;
  case 3 ... 4:
    // DIAMOND.
    n_cell = 1;

    if (true) {
      bn_prms.add_prm("ts1a",  3, 5e5, 1.0);
      bn_prms.add_prm("ts2a",  3, 5e5, 1.0);
      bn_prms.add_prm("ts1b",  3, 5e5, 1.0);
      bn_prms.add_prm("ts2b",  3, 5e5, 1.0);
      bn_prms.add_prm("ts1c",  3, 5e5, 1.0);
      bn_prms.add_prm("ts2c",  3, 5e5, 1.0);
      bn_prms.add_prm("ts1d",  3, 5e5, 1.0);
      bn_prms.add_prm("ts2d",  3, 5e5, 1.0);
    } else {
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
    }

    if (lat_case == 4) {
      // DIAMOND & VMX.
      bn_prms.add_prm("s1", 3, 5e5, 1.0);
      bn_prms.add_prm("s2", 3, 5e5, 1.0);
      bn_prms.add_prm("s3", 3, 5e5, 1.0);
      bn_prms.add_prm("s4", 3, 5e5, 1.0);
      bn_prms.add_prm("s5", 3, 5e5, 1.0);
    }
    break;
  case 5:
    // M-H6BAi.
    n_cell = 2;

    bn_prms.add_prm("sf1", 3, 5e5, 1.0);
    bn_prms.add_prm("sd1", 3, 5e5, 1.0);
    bn_prms.add_prm("sd2", 3, 5e5, 1.0);

    if (!fit_ksi) {
      // bn_prms.add_prm("sf1", 4, 5e5, 1.0);
      // bn_prms.add_prm("sd1", 4, 5e5, 1.0);
      // bn_prms.add_prm("sd2", 4, 5e5, 1.0);

      // bn_prms.add_prm("sh1a", 4, 5e5, 1.0);
      // bn_prms.add_prm("sh2b", 4, 5e5, 1.0);
      // bn_prms.add_prm("of1s", 4, 5e5, 1.0);

      // bn_prms.add_prm("sf1", 5, 5e5, 1.0);
      // bn_prms.add_prm("sd1", 5, 5e5, 1.0);
      // bn_prms.add_prm("sd2", 5, 5e5, 1.0);
    }
    break;
  case 6:
    // M-H6BA.
    n_cell = 1;

    bn_prms.add_prm("sf1", 3, 1e4, 1.0);
    bn_prms.add_prm("sd1", 3, 1e4, 1.0);
    bn_prms.add_prm("sd2", 3, 1e4, 1.0);

    // bn_prms.add_prm("sf1", 4, 1e4, 1e0);
    // bn_prms.add_prm("sd1", 4, 1e4, 1e0);
    // bn_prms.add_prm("sd2", 4, 1e4, 1e0);

    // bn_prms.add_prm("sh2", 4, 1e4, 1e0);
    // bn_prms.add_prm("s",   4, 1e4, 1e0);
    // bn_prms.add_prm("of1", 4, 1e4, 1e0);
    break;
  case 7:
    // ALS-U.
    n_cell = 1;

    bn_prms.add_prm("sf1",  3, 5e5, 1.0);
    bn_prms.add_prm("sf2",  3, 5e5, 1.0);
    bn_prms.add_prm("sf3",  3, 5e5, 1.0);
    bn_prms.add_prm("sd1h", 3, 5e5, 1.0);
    bn_prms.add_prm("sd2h", 3, 5e5, 1.0);
    bn_prms.add_prm("sd3h", 3, 5e5, 1.0);
    bn_prms.add_prm("sd4h", 3, 5e5, 1.0);

    bn_prms.add_prm("sf1",  4, 5e5, 1.0);
    bn_prms.add_prm("sf2",  4, 5e5, 1.0);
    bn_prms.add_prm("sf3",  4, 5e5, 1.0);
    bn_prms.add_prm("sd1h", 4, 5e5, 1.0);
    bn_prms.add_prm("sd2h", 4, 5e5, 1.0);
    bn_prms.add_prm("sd3h", 4, 5e5, 1.0);
    bn_prms.add_prm("sd4h", 4, 5e5, 1.0);

    bn_prms.add_prm("sf1",  5, 5e5, 1.0);
    bn_prms.add_prm("sf2",  5, 5e5, 1.0);
    bn_prms.add_prm("sf3",  5, 5e5, 1.0);
    bn_prms.add_prm("sd1h", 5, 5e5, 1.0);
    bn_prms.add_prm("sd2h", 5, 5e5, 1.0);
    bn_prms.add_prm("sd3h", 5, 5e5, 1.0);
    bn_prms.add_prm("sd4h", 5, 5e5, 1.0);

    bn_prms.add_prm("sh1",  4, 5e5, 1.0);
    bn_prms.add_prm("sh2",  4, 5e5, 1.0);
    bn_prms.add_prm("sh3",  4, 5e5, 1.0);

    // bn_prms.add_prm("sh1",  6, 5e5, 1.0);
    // bn_prms.add_prm("sh2",  6, 5e5, 1.0);
    // bn_prms.add_prm("sh3",  6, 5e5, 1.0);
    break;
  default:
    printf("\nlat_select: unknown Lattice type: %d\n", lat_case);
    break;
  }

  if (false) {
    // fit_tune(57.15/6.0, 22.25/6.0);
    // fit_tune(0.530831725+1e-4, 0.685735574-0*1e-4);
    fit_tune(58.12/6.0, 21.29/6.0);
    get_nu_ksi();
    exit(0);
  }
}

void get_nu_from_lin_map(const double delta, double nu[])
{
  int    k;
  double c;

  danot_(1);
  Map.identity(); Map[delta_] += tps(delta, 5);
  Map.propagate(1, n_elem);
  for (k = 0; k < 2; k++) {
    c = (Map[2*k][2*k]+Map[2*k+1][2*k+1])/2e0;
    nu[k] = acos(c)/(2e0*M_PI); nu[k] = sgn(Map[2*k][2*k+1])*nu[k];
    if (nu[k] < 0e0) nu[k] += 1e0;
  }
}


void get_ksi1(void)
{
  int    k;
  double nu1[2], nu2[2], ksi1[2];

  const double d_delta = 1e-5;

  get_nu_from_lin_map(d_delta, nu2); get_nu_from_lin_map(-d_delta, nu1);
  for (k = 0; k < 2; k++)
    ksi1[k] = (nu2[k]-nu1[k])/(2e0*d_delta);
  printf("\nksi1 = [%7.5f, %7.5f]\n", ksi1[X_], ksi1[Y_]);
}


int main(int argc, char *argv[])
{
  int j;

  rad_on    = false; H_exact        = false; totpath_on   = false;
  cavity_on = false; quad_fringe_on = false; emittance_on = false;
  IBS_on    = false;

  if (argc != 3) {
    printf("\n*** bad command line, no of arguments: %d\n", argc);
    exit(0);
  }

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

  n_cut = atoi(argv[2]);

  lat_select(lat_case);

  printf("\nn_cell:             %1d\n", n_cell);
  printf("scl_h:              %7.1e, %7.1e, %7.1e\n",
	 scl_h[0], scl_h[1], scl_h[2]);
  printf("scl_dnu:            %7.1e, %7.1e, %7.1e\n",
	 scl_dnu[0], scl_dnu[1], scl_dnu[2]);
  printf("scl_ksi:            %7.1e, %7.1e, %7.1e, %7.1e, %7.1e\n",
	 scl_ksi[0], scl_ksi[1], scl_ksi[2], scl_ksi[3], scl_ksi[4]);
  printf("scl_dnu_conf:       %7.1e\n", scl_dnu_conf);
  printf("scl_dnu_delta_conf: %7.1e\n", scl_dnu_delta_conf);
  printf("dnu/dJ:             %d\n", DNU);
  printf("symmetric:          %d\n", symm);
  printf("\n2J:                %10.3e, %10.3e\n", twoJ[X_], twoJ[Y_]);
  printf("delta_max:          %7.1e\n", delta_max);
  if (c_g) {
    printf("Conj. Grad.:        %d\n", n_cut);
    // printf("Conj. Grad. List of Singular Values:\n");
    // for (j = 0; j < (int)bn_prms.svd_list.size(); j++)
    //   printf(" %2d", bn_prms.svd_list[j]);
    // printf("\n");
  } else
    printf("Lev. Marq.\n");

  if (!false) {
    // no_mpoles(Sext);
    get_nu_ksi();
    get_twiss();
    prt_lat("linlat.out", 10);
    prt_lat("linlat1.out");
    // exit(0);
  }

  Id_scl.identity();
  for (j = 0; j < 4; j++)
    Id_scl[j] *= sqrt(twoJ[j/2]);
  Id_scl[delta_] *= delta_max;

  Id_delta_scl.identity();
  for (j = 0; j < 4; j++)
    Id_delta_scl[j] *= sqrt(twoJ_delta[j/2]);
  Id_delta_scl[delta_] *= delta_max;

  if (false) {
    danot_(NO-1);
    cavity_on = true; rad_on = true;
    get_map_n(n_cell);
    // MAX-V:
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
    if (false) {
      // Tweak linear chromaticity.
      bn_prms.step = 1.0;
      bn_prms.ini_prm();
      fit_ksi1(0e0, 0e0);
    }

    danot_(NO-1);
    get_Map();
    danot_(NO);
    K = MapNorm(Map, g, A1, A0, Map_res, 1);
    CtoR(K, K_re, K_im); CtoR(get_h(), h_re, h_im);
    prt_h_K();
    get_ksi1();
    exit(0);
  }

  // Step is 1.0 for conjugated gradient method.
  bn_prms.bn_tol    = 1e-3;
  // bn_prms.svd_n_cut = 0;
  bn_prms.step      = 1.0;

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

  if (c_g) {
    bn_prms.svd_n_cut = n_cut;
    min_conj_grad(true);
  } else
    min_lev_marq();

  prt_h_K();
}
