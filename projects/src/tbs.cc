#include <cfloat>

#define NO 5

#include "thor_lib.h"

int no_tps = NO,

#define DOF_3 0

#if !DOF_3
  ndpt_tps = 5;
#else
  // The cavity must be turned on.
  ndpt_tps = 0;
#endif


// Center of straight.
const double
  beta_inj[]   = {8.1, 2.3},
  twoJ[]       = {sqr(8e-3)/beta_inj[X_], sqr(2e-3)/beta_inj[Y_]},
  delta_max    = 3e-2;


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


tps          h_re, h_im, h_re_scl, h_im_scl, K_re, K_im, K_re_scl;
ss_vect<tps> nus, nus_scl, Id_scl;
param_type   bn_prms;

const int n_prt  = 8;
const bool scale = false;


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


void get_dchi2(double *df)
{
  int    k, loc;
  double eps;

  const bool prt = false;

  for (k = 0; k < bn_prms.n_prm; k++) {

    if (prt) {
      printf("\nget_dchi2: ");
      loc = get_loc(bn_prms.Fnum[k], 1) - 1;
      prt_name(stdout, elem[loc].Name, ":", 6);
      printf(" %2d %10.3e", bn_prms.n[k], eps);
    }

    // constr_dparam(bn_prms.Fnum[k], bn_prms.n[k], eps);
    // df[k+1] = lat_constr.get_chi2();

    // constr_dparam(bn_prms.Fnum[k], bn_prms.n[k], -2e0*eps);
    // df[k+1] -= lat_constr.get_chi2();
    // df[k+1] /= 2e0*eps;

    // constr_dparam(bn_prms.Fnum[k], bn_prms.n[k], eps);
  }

  // Avoid: "warning: deprecated conversion from string constant to ‘char*’".
  // dvdump(stdout,
  // 	 (char *)"\nget_dchi2:", df, bn_prms.n_prm, (char *)" %12.5e");
}


void f_der(double *b2, double *df)
{
  // bn_prms.set_prm(b3);

  // lat_constr.get_dchi2(df);
}


double f_nl(double bn[])
{
  int    i, j;
  double chi2;
  tps    K_re, g2, dnu2;

  const bool prt = true;

  n_iter++;

  for (i = 1; i <= n_prm; i++)
    set_bn(prm[i-1], prm_n[i-1], bn[i]);

  get_dyn(K_re, g2, dnu2);

  chi2 = 0e0;

  chi2 += scl_ksi1*sqr(ksi1[X_]*M_PI+h_ijklm(K_re, 1, 1, 0, 0, 1));
  chi2 += scl_ksi1*sqr(ksi1[Y_]*M_PI+h_ijklm(K_re, 0, 0, 1, 1, 1));
  chi2 += scl_dnu2*dnu2.cst();

  chi2 += scl_res*g2.cst();

  if (chi2 < chi2_min) {
    if (prt) {
      cout << scientific << setprecision(3)
	   << "ksi = "
	   << scl_ksi1*sqr(ksi1[X_]*M_PI+h_ijklm(K_re, 1, 1, 0, 0, 1))
	   << " " << scl_ksi1*sqr(ksi1[Y_]*M_PI+h_ijklm(K_re, 0, 0, 1, 1, 1))
	   << endl;
      cout << scientific << setprecision(3)
	   << "g2 = " << scl_res*g2.cst()
	   << ", dnu2 = " << scl_dnu2*dnu2.cst() << endl;
    }

    chi2_min = min(chi2, chi2_min);

    cout << endl;
    cout << "bnL:";
    for (i = 1; i <= n_prm; i++)
      cout << scientific << setprecision(3)
	   << setw(11) << get_bnL(prm[i-1], 1, prm_n[i-1]);
    cout << endl;

    cout << endl;
    cout << scientific << setprecision(1)
	 << setw(2) << n_iter << ", chi2_min: " << chi2_min << endl;

    sext_out << endl;
    sext_out << "n = " << n_iter << ":" << endl;
    for (i = 1; i <= n_prm; i++)
      for (j = 1; j <= get_n_Kids(prm[i-1]); j++) {
	sext_out << fixed << setprecision(7) 
		 << setw(9) << get_Name(prm[i-1])
		 << "(" << j << ") = "
		 << setw(11) << get_bnL(prm[i-1], 1, prm_n[i-1])
		 << setw(2) << prm_n[i-1] << endl;
      }

    sext_out.flush();
  }

  return chi2;
}


void df_nl(double bn[], double df[])
{
  int    k;
  tps    K_re, g2, dnu2;

  cout << "df_nl" << endl;

  for (k = 1; k <= n_prm; k++)
    set_bn(prm[k-1], prm_n[k-1], bn[k]);

  for (k = 1; k <= n_prm; k++) {
    set_bn_par(prm[k-1], prm_n[k-1], 7);

    get_dyn(K_re, g2, dnu2);

    df[k] = 0e0;

    df[k] +=
      scl_ksi1*2e0*h_ijklm_p(K_re, 1, 1, 0, 0, 1, 7)
      *(ksi1[X_]*M_PI+h_ijklm(K_re, 1, 1, 0, 0, 1));

    df[k] +=
      scl_ksi1*2e0*h_ijklm_p(K_re, 0, 0, 1, 1, 1, 7)
      *(ksi1[Y_]*M_PI+h_ijklm(K_re, 0, 0, 1, 1, 1));

    df[k] += scl_dnu2*h_ijklm_p(dnu2, 0, 0, 0, 0, 0, 7);

    df[k] += scl_res*h_ijklm_p(g2, 0, 0, 0, 0, 0, 7);

    clr_bn_par(prm[k-1], prm_n[k-1]);
  }
}


void fit_conj_grad(param_type &bn_prms, double (*f)(double *))
{
  int          n_b2, iter;
  double       *b2, fret;
  ss_vect<tps> A;

  const double ftol = 1e-8;

  n_b2 = bn_prms.n_prm;

  b2 = dvector(1, n_b2);

  bn_prms.ini_prm();
  f(b2);

  dfrprmn(b2, n_b2, ftol, &iter, &fret, f, f_der);

  printf("\n  iter = %d fret = %12.5e\n", iter, fret);
  printf("b2s:\n");
  // bn_prms.prt_prm(b2);
  // bn_prms.set_prm(b2);
  // eps_x = get_lin_opt(false);
  // f_prt(b2);

  free_dvector(b2, 1, n_b2);
}


void lat_select(const int lat_case)
{
  bn_prms.add_prm("sf1", 3, 1e4, 1.0);
  bn_prms.add_prm("sd1", 3, 1e4, 1.0);
  bn_prms.add_prm("sd2", 3, 1e4, 1.0);

  bn_prms.add_prm("sh1a", 4, 1e4, 1.0);
  bn_prms.add_prm("sh1b", 4, 1e4, 1.0);
  bn_prms.add_prm("sh2",  4, 1e4, 1e2);
  bn_prms.add_prm("s",    4, 1e4, 1e2);

  // bn_prms.add_prm("of1",  4, 1e4, 1e2);
}


void get_nu_k(void)
{
  danot_(NO-1);
  get_Map();
  danot_(NO);
  K = MapNorm(Map, g, A1, A0, Map_res, 1);
  nus = dHdJ(K); nus_scl = nus*Id_scl;
  CtoR(K, K_re, K_im);
  K_re_scl = K_re*Id_scl;
  CtoR(get_h(), h_re, h_im); h_re_scl = h_re*Id_scl; h_im_scl = h_im*Id_scl;
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

}
