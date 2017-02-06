#define NO 5

#include "thor_lib.h"

int no_tps   = NO,
    ndpt_tps = 5;

extern ss_vect<tps> Map;


typedef struct {
  int n;
  int Fnum;
 } bn_type;

int m_cstr, n_prm;

const bool   t566 = true, chrom_trms = false, geom_trms = false;
const int    disp_n = NO-1;
const int    n_prm_max = 20, m_cstr_max = 30, n_delta = 30;
const double delta = -5e-2, T566_ref = -2.34e0, b3L_max = 500e0;
const double scl_T566 = 1e0, scl_ksi = 1e0;
const double scl_h3_chrom = 1e0, scl_h3_geom = 1e0;
const double scl_delta = 1e0;
const double step  = 0.3e0, s_cut = 1e-10;

int     n_iter;
bn_type bn_prm[n_prm_max];


void no_mpoles(void)
{
  int j, k;

  std::cout << std::endl;
  std::cout << "zeroing multipoles" << std::endl;
  std::cout << std::endl;
  for (j = 0; j < n_elem; j++)
    if (elem[j].kind == Mpole)
      for (k = Sext; k < mpole_max; k++) {
	//	std::cout << "zeroing " << elem[j].Name << std::endl;
	set_bn(elem[j].Fnum, elem[j].Knum, k, 0.0);
      }
}


int get_ind(const int k)
{
  int index[] = {x_, px_, y_, py_, ct_, delta_};
  return index[k];
}


void prt_lin_map1(const int n_DOF, const ss_vect<tps> &map)
{
  // Phase-space coordinates in canonical order.
  int  i, j;

  std::cout << std::endl;
  for (i = 0; i < 2*n_DOF; i++) {
    for (j = 0; j < 2*n_DOF; j++) {
      if (true)
	std::cout << std::scientific << std::setprecision(6)
	     << std::setw(14) << map[get_ind(i)][get_ind(j)];
      else
	std::cout << std::scientific << std::setprecision(16)
	     << std::setw(24) << map[get_ind(i)][get_ind(j)];
    }
    std::cout << std::endl;
  }
}


template<typename T>
T atan2(const T y, const T x)
{
  T z;

  if (x == 0e0)
    z = (y > 0e0)? M_PI/2e0 : -M_PI/2e0;
  else
    z = atan(y/x);
  if (x < 0e0) z += (y >= 0e0)? M_PI : -M_PI;
  return z;
}


void get_nu(tps nu[])
{
  int          k;
  ss_vect<tps> A;

  A = elem_tps[0].A1; A.propagate(1, n_elem);

  for (k = 0; k < 2; k++) {
    nu[k] = atan2(Der(A[2*k], 2*k+2), Der(A[2*k], 2*k+1))/(2e0*M_PI);
    if (nu[k] < 0e0) nu[k] += 1e0;
  }
}


void get_ksi(double ksi[])
{
  int k;
  tps nu[2];

  get_nu(nu);
  for (k = 0; k < 2; k++) ksi[k] = h_ijklm(nu[k], 0, 0, 0, 0, 1);

  std::cout << std::endl;
  std::cout << std::fixed << std::setprecision(5)
       << "nu  = [" << std::setw(8) << nu[X_].cst() << ", "
       << std::setw(8) << nu[Y_].cst() << "]" << std::endl;
  std::cout << std::fixed << std::setprecision(5)
       << "ksi = [" << std::setw(8) << ksi[X_] << ", "
       << std::setw(8) << ksi[Y_] << "]" << std::endl;
}


void scan_delta(const int n)
{
  int             k;
  double          d;
  ss_vect<double> ps;
  std::ofstream        outf;

  const std::string file_name = "scan_delta.out";

  file_wr(outf, file_name.c_str());

  for (k = 0; k < n; k++) {
    d = (double)k/double(n-1)*delta;
    ps.zero(); ps[delta_] = d; ps.propagate(1, n_elem);
    outf << std::scientific << std::setprecision(6)
	 << std::setw(14) << d << std::setw(14) << ps << std::endl;
  }
  outf.close();
}


void prt_ct(void)
{
  int          k;
  ss_vect<tps> ps;
  std::ofstream     outf;

  const std::string file_name = "ct.out";

  file_wr(outf, file_name.c_str());

  ps.identity();
  for (k = 1; k <= n_elem; k++) {
    ps.propagate(k, k);

    outf << std::setw(4) << k << " " << std::setw(15) << elem[k-1].Name
	 << std::fixed << std::setprecision(5) << std::setw(9) << elem[k-1].S
	 << std::setprecision(1) << std::setw(5) << get_code(elem[k-1])
	 << std::scientific << std::setprecision(5)
	 << std::setw(13) << h_ijklm(ps[ct_], 0, 0, 0, 0, 1)
	 << std::setw(13) << h_ijklm(ps[ct_], 0, 0, 0, 0, 2) << std::endl;
  }

  outf.close();
}


void prt_eta(ss_vect<tps> &ps)
{
  int          i, k;
  ss_vect<tps> ps_scl, Id_scl;
  std::ofstream     outf;

  const std::string file_name = "eta2.out";

  file_wr(outf, file_name.c_str());

  Id_scl.identity(); Id_scl[delta_] *= fabs(delta);

  for (i = 1; i <= n_elem; i++) {
    ps.propagate(i, i); ps_scl = ps*Id_scl;

    outf << std::setw(4) << i << " " << std::setw(15) << elem[i-1].Name
	 << std::fixed << std::setprecision(5) << std::setw(9) << elem[i-1].S
	 << std::setprecision(1) << std::setw(5) << get_code(elem[i-1]);
    for (k = 1; k <= no_tps; k++)
      if (true)
	outf << std::scientific << std::setprecision(5)
	     << std::setw(13) << h_ijklm(ps_scl[x_], 0, 0, 0, 0, k)
	     << std::setw(13) << h_ijklm(ps_scl[px_], 0, 0, 0, 0, k);
      else
	outf << std::scientific << std::setprecision(5)
	     << std::setw(13) << h_ijklm(ps[x_], 0, 0, 0, 0, k)
	     << std::setw(13) << h_ijklm(ps[px_], 0, 0, 0, 0, k);
    outf << std::endl;
  }

  outf.close();
}


template<typename T>
T trapz(const int n, const double h, const T y[])
{
  int k;
  T   intgrl;

  intgrl = 0e0;
  for (k = 0; k < n-1; k++) {
    intgrl += (y[k+1]+y[k])/2e0*h;
  }

 return intgrl;
}


template<typename T>
void f_opt0(const int n, T &x2d_intgrl, T &px2d_intgrl,
	    double &x_max, double &px_max)
{
  int        k;
  double     h;
  T          x_delta[n], px_delta[n];
  ss_vect<T> ps;

  danot_(1);
  x_max = 0e0; px_max = 0e0; h = delta/(n-1);
  for (k = 0; k < n; k++) {
    ps.zero(); ps[delta_] = k*h; ps.propagate(1, n_elem);
    x_delta[k] = sqr(ps[x_]); px_delta[k] = sqr(ps[px_]);
    x_max = max(fabs(is_double<T>::cst(ps[x_])), x_max);
    px_max = max(fabs(is_double<T>::cst(ps[px_])), px_max);
  }
  x2d_intgrl = trapz(n, h, x_delta); px2d_intgrl = trapz(n, h, px_delta);
}


void h_init(double *bn_max, double *bn)
{
  int    i;
  double L;

  std::cout << std::endl << "b3L0:" << std::endl;
  for (i = 1; i <= n_prm; i++) {
    // Note, Jacobian is a function of the multipole strengths.
    L = get_L(bn_prm[i-1].Fnum, 1);
    if (L == 0e0) L = 1e0;
    bn_max[i] = b3L_max/L;
    bn[i] = get_bn(bn_prm[i-1].Fnum, 1, bn_prm[i-1].n);
    std::cout << std::scientific << std::setprecision(5)
	 << std::setw(13) << get_bnL(bn_prm[i-1].Fnum, 1, bn_prm[i-1].n);
    if (i % 5 == 0) std::cout << std::endl;
  }
  if (n_prm % 5 != 0) std::cout << std::endl;
  std::cout << std::endl;
}


void prt_system(const int m_cstr, const int n_prm, double **A, double *b,
	        const double chi2_old, const double chi2,
		std::ostringstream hs[])
{
  int i, j;

  std::cout << std::scientific << std::setprecision(1)
       << std::endl << n_iter << " Ax = b, chi2: " << chi2_old << " -> " << chi2
       << ":" << std::endl << std::endl;
  for (i = 1; i <= m_cstr; i++) {
    std::cout  << std::setw(3) << i << " " << hs[i-1].str();
    for (j = 1; j <= n_prm; j++)
      std::cout << std::scientific << std::setprecision(2)
		<< std::setw(10) << A[i][j];
    std::cout << std::scientific << std::setprecision(2)
	      << std::setw(10) << -b[i] << std::endl;
  }
}


void prt_sext(std::ofstream &sext_out)
{
  int i;

  sext_out << std::endl;
  sext_out << "n = " << n_iter << ":" << std::endl;
  for (i = 0; i < n_prm; i++)
      sext_out << std::fixed << std::setprecision(5) 
	       << get_Name(bn_prm[i].Fnum)
	       << ": Multipole, L = 0.0, HOM = (" << bn_prm[i].n << ", "
	       << std::setw(10) << get_bn(bn_prm[i].Fnum, 1, bn_prm[i].n)
	       << ", 0.0), N = 1, Method = meth;" << std::endl;
  sext_out.flush();
}


double get_chi2(const int n, const double data[])
{
  int    j;
  double chi2 = 0e0;

  for (j = 1; j <= n; j++) chi2 += sqr(data[j]);

  return chi2;
}


tps get_h(const ss_vect<tps> A0, const ss_vect<tps> A1)
{
  ss_vect<tps> Map1, R;

  // Parallel transport nonlinear kick to start of lattice,
  // assumes left to right evaluation.

  if (true)
    // Dragt-Finn factorization
    return LieFact_DF(Inv(A1)*Map*A0, R)*R;
  else {
    // single Lie exponent
    danot_(1); Map1 = Map; danot_(no_tps);
    return LieFact(Inv(A1)*Map*Inv(Map1)*A0);
  }
}


double f_opt1(double b3s[])
{
  int          i, k;
  double       f, *b;
  tps          nu[2], h_drv;
  ss_vect<tps> Id_scl, map_scl;

  b = dvector(1, m_cstr_max);

  Id_scl.identity(); Id_scl[delta_] *= fabs(delta);

  for (i = 1; i <= n_prm; i++)
    set_bn(bn_prm[i-1].Fnum, bn_prm[i-1].n, b3s[i]);

  danot_(no_tps); get_Map(); danot_(no_tps);
  get_nu(nu);
  h_drv = get_h(elem_tps[0].A1, elem_tps[n_elem-1].A1);

  map_scl = Map*Id_scl;

  m_cstr = 0;
  if (t566)
    b[++m_cstr] = sqrt(scl_T566)*(h_ijklm(Map[ct_], 0, 0, 0, 0, 2)-T566_ref);

  b[++m_cstr] = sqrt(scl_ksi)*h_ijklm(nu[X_], 0, 0, 0, 0, 1);
  b[++m_cstr] = sqrt(scl_ksi)*h_ijklm(nu[Y_], 0, 0, 0, 0, 1);

  if (chrom_trms) {
    b[++m_cstr] = sqrt(scl_h3_chrom)*h_ijklm(h_drv, 1, 0, 0, 0, 2);
    b[++m_cstr] = sqrt(scl_h3_chrom)*h_ijklm(h_drv, 0, 1, 0, 0, 2);
    b[++m_cstr] = sqrt(scl_h3_chrom)*h_ijklm(h_drv, 2, 0, 0, 0, 1);
    b[++m_cstr] = sqrt(scl_h3_chrom)*h_ijklm(h_drv, 0, 2, 0, 0, 1);
    b[++m_cstr] = sqrt(scl_h3_chrom)*h_ijklm(h_drv, 0, 0, 2, 0, 1);
    b[++m_cstr] = sqrt(scl_h3_chrom)*h_ijklm(h_drv, 0, 0, 0, 2, 1);
  }

  if (geom_trms) {
    b[++m_cstr] = sqrt(scl_h3_geom)*h_ijklm(h_drv, 2, 1, 0, 0, 0);
    b[++m_cstr] = sqrt(scl_h3_geom)*h_ijklm(h_drv, 1, 2, 0, 0, 0);
    b[++m_cstr] = sqrt(scl_h3_geom)*h_ijklm(h_drv, 1, 0, 1, 1, 0);
    b[++m_cstr] = sqrt(scl_h3_geom)*h_ijklm(h_drv, 0, 1, 1, 1, 0);
    b[++m_cstr] = sqrt(scl_h3_geom)*h_ijklm(h_drv, 3, 0, 0, 0, 0);
    b[++m_cstr] = sqrt(scl_h3_geom)*h_ijklm(h_drv, 0, 3, 0, 0, 0);
    b[++m_cstr] = sqrt(scl_h3_geom)*h_ijklm(h_drv, 1, 0, 0, 2, 0);
    b[++m_cstr] = sqrt(scl_h3_geom)*h_ijklm(h_drv, 0, 1, 2, 0, 0);
    b[++m_cstr] = sqrt(scl_h3_geom)*h_ijklm(h_drv, 1, 0, 2, 0, 0);
    b[++m_cstr] = sqrt(scl_h3_geom)*h_ijklm(h_drv, 0, 1, 0, 2, 0);
  }

  for (k = 2; k <= disp_n; k++) {
    b[++m_cstr] = sqrt(scl_delta)*h_ijklm(map_scl[x_],  0, 0, 0, 0, k);
    b[++m_cstr] = sqrt(scl_delta)*h_ijklm(map_scl[px_], 0, 0, 0, 0, k);
  }

  f = 0e0;
  for (i = 1; i <= m_cstr; i++)
    f += sqr(b[i]);

  std::cout << std::scientific << std::setprecision(3)
       << std::setw(4) << n_iter << " f = " << std::setw(9) << f << std::endl;

  free_dvector(b, 1, m_cstr_max);

  return f;
}


void h_opt1(double &chi2, const double *bn_max, double *bn,
	    double &dbn_max, double *g, double *h, const bool prt_iter,
	    std::ofstream &sext_out)
{
  // Conjugate gradient method.
  std::ostringstream hs[m_cstr_max];
  int           i, k;
  double        L, chi2_old, fret, g2, gamma, dg2;
  double        **A, *b, *w, **U, **V, *dbn, *bn0;
  tps           nu[2], h_drv;
  ss_vect<tps>  Id_scl, map_scl;

  const bool   cong_grad = true, prt = true;

  b = dvector(1, m_cstr_max); w = dvector(1, n_prm); dbn = dvector(1, n_prm);
  bn0 = dvector(1, n_prm);
  A = dmatrix(1, m_cstr_max, 1, n_prm); U = dmatrix(1, m_cstr_max, 1, n_prm);
  V = dmatrix(1, n_prm, 1, n_prm);

  Id_scl.identity(); Id_scl[delta_] *= fabs(delta);

  std::cout << std::endl;
  for (i = 1; i <= n_prm; i++) {
    set_bn_par(bn_prm[i-1].Fnum, bn_prm[i-1].n, 7);

    danot_(no_tps); get_Map(); danot_(no_tps);
    get_nu(nu);
    h_drv = get_h(elem_tps[0].A1, elem_tps[n_elem-1].A1);
 
    map_scl = Map*Id_scl;

    m_cstr = 0;
    if (t566)
      A[++m_cstr][i] = sqrt(scl_T566)*h_ijklm_p(Map[ct_], 0, 0, 0, 0, 2, 7);

    A[++m_cstr][i] = sqrt(scl_ksi)*h_ijklm_p(nu[X_], 0, 0, 0, 0, 1, 7);
    A[++m_cstr][i] = sqrt(scl_ksi)*h_ijklm_p(nu[Y_], 0, 0, 0, 0, 1, 7);

    if (chrom_trms) {
      A[++m_cstr][i] = sqrt(scl_h3_chrom)*h_ijklm_p(h_drv, 1, 0, 0, 0, 2, 7);
      A[++m_cstr][i] = sqrt(scl_h3_chrom)*h_ijklm_p(h_drv, 0, 1, 0, 0, 2, 7);
      A[++m_cstr][i] = sqrt(scl_h3_chrom)*h_ijklm_p(h_drv, 2, 0, 0, 0, 1, 7);
      A[++m_cstr][i] = sqrt(scl_h3_chrom)*h_ijklm_p(h_drv, 0, 2, 0, 0, 1, 7);
      A[++m_cstr][i] = sqrt(scl_h3_chrom)*h_ijklm_p(h_drv, 0, 0, 2, 0, 1, 7);
      A[++m_cstr][i] = sqrt(scl_h3_chrom)*h_ijklm_p(h_drv, 0, 0, 0, 2, 1, 7);
    }

    if (geom_trms) {
      A[++m_cstr][i] = sqrt(scl_h3_geom)*h_ijklm_p(h_drv, 2, 1, 0, 0, 0, 7);
      A[++m_cstr][i] = sqrt(scl_h3_geom)*h_ijklm_p(h_drv, 1, 2, 0, 0, 0, 7);
      A[++m_cstr][i] = sqrt(scl_h3_geom)*h_ijklm_p(h_drv, 1, 0, 1, 1, 0, 7);
      A[++m_cstr][i] = sqrt(scl_h3_geom)*h_ijklm_p(h_drv, 0, 1, 1, 1, 0, 7);
      A[++m_cstr][i] = sqrt(scl_h3_geom)*h_ijklm_p(h_drv, 3, 0, 0, 0, 0, 7);
      A[++m_cstr][i] = sqrt(scl_h3_geom)*h_ijklm_p(h_drv, 0, 3, 0, 0, 0, 7);
      A[++m_cstr][i] = sqrt(scl_h3_geom)*h_ijklm_p(h_drv, 1, 0, 0, 2, 0, 7);
      A[++m_cstr][i] = sqrt(scl_h3_geom)*h_ijklm_p(h_drv, 0, 1, 2, 0, 0, 7);
      A[++m_cstr][i] = sqrt(scl_h3_geom)*h_ijklm_p(h_drv, 1, 0, 2, 0, 0, 7);
      A[++m_cstr][i] = sqrt(scl_h3_geom)*h_ijklm_p(h_drv, 0, 1, 0, 2, 0, 7);
    }

    for (k = 2; k <= disp_n; k++) {
      A[++m_cstr][i] = sqrt(scl_delta)*h_ijklm_p(map_scl[x_],  0, 0, 0, 0, k, 7);
      A[++m_cstr][i] = sqrt(scl_delta)*h_ijklm_p(map_scl[px_], 0, 0, 0, 0, k, 7);
    }

    if (i == n_prm) {
      m_cstr = 0;
      if (t566) {
	hs[m_cstr++] << "T566   ";
	b[m_cstr] = -sqrt(scl_T566)*(h_ijklm(Map[ct_], 0, 0, 0, 0, 2)-T566_ref);
      }

      hs[m_cstr++] << "ksi_x  ";
      b[m_cstr] = -sqrt(scl_ksi)*h_ijklm(nu[X_], 0, 0, 0, 0, 1);
      hs[m_cstr++] << "ksi_y  ";
      b[m_cstr] = -sqrt(scl_ksi)*h_ijklm(nu[Y_], 0, 0, 0, 0, 1);

      if (chrom_trms) {
	hs[m_cstr++] << "h_10002";
	b[m_cstr] = -sqrt(scl_h3_chrom)*h_ijklm(h_drv, 1, 0, 0, 0, 2);
	hs[m_cstr++] << "h_01002";
	b[m_cstr] = -sqrt(scl_h3_chrom)*h_ijklm(h_drv, 0, 1, 0, 0, 2);
	hs[m_cstr++] << "h_20001";
	b[m_cstr] = -sqrt(scl_h3_chrom)*h_ijklm(h_drv, 2, 0, 0, 0, 1);
	hs[m_cstr++] << "h_02001";
	b[m_cstr] = -sqrt(scl_h3_chrom)*h_ijklm(h_drv, 0, 2, 0, 0, 1);
	hs[m_cstr++] << "h_00201";
	b[m_cstr] = -sqrt(scl_h3_chrom)*h_ijklm(h_drv, 0, 0, 2, 0, 1);
	hs[m_cstr++] << "h_00021";
	b[m_cstr] = -sqrt(scl_h3_chrom)*h_ijklm(h_drv, 0, 0, 0, 2, 1);
      }

      if (geom_trms) {
	hs[m_cstr++] << "h_21000";
	b[m_cstr] = -sqrt(scl_h3_geom)*h_ijklm(h_drv, 2, 1, 0, 0, 0);
	hs[m_cstr++] << "h_12000";
	b[m_cstr] = -sqrt(scl_h3_geom)*h_ijklm(h_drv, 1, 2, 0, 0, 0);
	hs[m_cstr++] << "h_10110";
	b[m_cstr] = -sqrt(scl_h3_geom)*h_ijklm(h_drv, 1, 0, 1, 1, 0);
	hs[m_cstr++] << "h_01110";
	b[m_cstr] = -sqrt(scl_h3_geom)*h_ijklm(h_drv, 0, 1, 1, 1, 0);
	hs[m_cstr++] << "h_30000";
	b[m_cstr] = -sqrt(scl_h3_geom)*h_ijklm(h_drv, 3, 0, 0, 0, 0);
	hs[m_cstr++] << "h_03000";
	b[m_cstr] = -sqrt(scl_h3_geom)*h_ijklm(h_drv, 0, 3, 0, 0, 0);
	hs[m_cstr++] << "h_10020";
	b[m_cstr] = -sqrt(scl_h3_geom)*h_ijklm(h_drv, 1, 0, 0, 2, 0);
	hs[m_cstr++] << "h_10020";
	b[m_cstr] = -sqrt(scl_h3_geom)*h_ijklm(h_drv, 0, 1, 2, 0, 0);
	hs[m_cstr++] << "h_01200";
	b[m_cstr] = -sqrt(scl_h3_geom)*h_ijklm(h_drv, 1, 0, 2, 0, 0);
	hs[m_cstr++] << "h_10200";
	b[m_cstr] = -sqrt(scl_h3_geom)*h_ijklm(h_drv, 0, 1, 0, 2, 0);
      }

      for (k = 2; k <= disp_n; k++) {
	hs[m_cstr++] << "eta" << k << "   ";
	b[m_cstr] = -sqrt(scl_delta)*h_ijklm(map_scl[x_],  0, 0, 0, 0, k);
	hs[m_cstr++] << "eta'" << k << "  ";
	b[m_cstr] = -sqrt(scl_delta)*h_ijklm(map_scl[px_], 0, 0, 0, 0, k);
      }
    }

    clr_bn_par(bn_prm[i-1].Fnum, bn_prm[i-1].n);
  }

  chi2_old = chi2; chi2 = get_chi2(m_cstr, b);

  if (prt) prt_system(m_cstr, n_prm, A, b, chi2_old, chi2, hs);
    
  SVD_lim(m_cstr, n_prm, A, b, bn_max, s_cut, bn, dbn);

  if (cong_grad) {
    dvcopy(bn, n_prm, bn0);

    if (n_iter == 1) {
      for (i = 1; i <= n_prm; i++) {
	g[i] = dbn[i]; h[i] = g[i];
      }
    } else {
      dg2 = g2 = 0e0;
      for (i = 1; i <= n_prm; i++) {
	g2 += sqr(g[i]); dg2 += (dbn[i]-g[i])*dbn[i];
      }
      if (g2 != 0e0) {
	gamma = dg2/g2;
	for (i = 1; i <= n_prm; i++) {
	  g[i] = dbn[i]; dbn[i] = h[i] = g[i] + gamma*h[i];
	}
      } else {
	std::cout << "g.g = 0" << std::endl; exit(0);
      }
    }

    std::cout << std::endl;
    d_linmin(bn, dbn, n_prm, &fret, f_opt1);
  }

  std::cout << std::endl << "db3L (dcorr*L):" << std::endl;
  dbn_max = 0e0;
  for (i = 1; i <= n_prm; i++) {
    L = get_L(bn_prm[i-1].Fnum, 1);
    if (L == 0e0) L = 1e0;
    if (cong_grad) {
      dbn_max = max(fabs((bn[i]-bn0[i])*L), dbn_max);
      std::cout << std::scientific << std::setprecision(3)
		<< std::setw(11) << (bn[i]-bn0[i])*L;
    } else {
      dbn_max = max(fabs(step*dbn[i]*L), dbn_max);
      std::cout << std::scientific << std::setprecision(3)
		<< std::setw(11) << dbn[i]*L;
    }
  }
  std::cout << std::endl;

  std::cout << "b3L:" << std::endl;
  for (i = 1; i <= n_prm; i++) {
    if (cong_grad)
      set_bn(bn_prm[i-1].Fnum, bn_prm[i-1].n, bn[i]);
    else
      set_dbn(bn_prm[i-1].Fnum, bn_prm[i-1].n, step*dbn[i]);
    bn[i] = get_bn(bn_prm[i-1].Fnum, 1, bn_prm[i-1].n);
    std::cout << std::scientific << std::setprecision(6)
	 << std::setw(14) << get_bnL(bn_prm[i-1].Fnum, 1, bn_prm[i-1].n);
    if (i % 5 == 0) std::cout << std::endl;
  }
  if (n_prm % 5 != 0) std::cout << std::endl;

  if (prt_iter) prt_sext(sext_out);

  free_dvector(b, 1, m_cstr_max); free_dvector(w, 1, n_prm);
  free_dvector(dbn, 1, n_prm); free_dvector(bn0, 1, n_prm);
  free_dmatrix(A, 1, m_cstr_max, 1, n_prm);
  free_dmatrix(U, 1, m_cstr_max, 1, n_prm); free_dmatrix(V, 1, n_prm, 1, n_prm);
}


void h_zero(void)
{
  std::string       str;
  double       *g, *h, dbn_max, chi2, *bn, *bn_max;
  ss_vect<tps> ps;
  std::ofstream     sext_out;

  const bool   prt_iter   = true;
  const int    n_iter_max = 1000;
  const double bn_tol     = 1e-2;

  bn = dvector(1, n_prm); bn_max = dvector(1, n_prm);
  g = dvector(1, n_prm); h = dvector(1, n_prm);

  file_wr(sext_out, "sext.dat");

  h_init(bn_max, bn);

  n_iter = 0; chi2 = 0e0;
  do {
    n_iter++;
    h_opt1(chi2, bn_max, bn, dbn_max, g, h, prt_iter, sext_out);

    prt_ct();
    ps.identity();
    prt_eta(ps);
    prt_mfile("flat_file.dat");
  } while ((dbn_max > bn_tol) && (n_iter < n_iter_max));

  if (!prt_iter) prt_sext(sext_out);

  free_dvector(bn, 1, n_prm); free_dvector(bn_max, 1, n_prm);
  free_dvector(g, 1, n_prm); free_dvector(h, 1, n_prm);
}


void get_prms(void)
{
  int k;

  // Arc3.
  n_prm = 0;
  // bn_prm[n_prm].Fnum = get_Fnum("s0"); bn_prm[n_prm++].n = Sext;
  bn_prm[n_prm].Fnum = get_Fnum("s1"); bn_prm[n_prm++].n = Sext;
  bn_prm[n_prm].Fnum = get_Fnum("s2"); bn_prm[n_prm++].n = Sext;
  bn_prm[n_prm].Fnum = get_Fnum("s3"); bn_prm[n_prm++].n = Sext;
  bn_prm[n_prm].Fnum = get_Fnum("s4"); bn_prm[n_prm++].n = Sext;
  bn_prm[n_prm].Fnum = get_Fnum("s5"); bn_prm[n_prm++].n = Sext;
  bn_prm[n_prm].Fnum = get_Fnum("s6"); bn_prm[n_prm++].n = Sext;
  // bn_prm[n_prm].Fnum = get_Fnum("s7"); bn_prm[n_prm++].n = Sext;
  // bn_prm[n_prm].Fnum = get_Fnum("s8"); bn_prm[n_prm++].n = Sext;
  // bn_prm[n_prm].Fnum = get_Fnum("s9"); bn_prm[n_prm++].n = Sext;
  // bn_prm[n_prm].Fnum = get_Fnum("s10"); bn_prm[n_prm++].n = Sext;
  // bn_prm[n_prm].Fnum = get_Fnum("s11"); bn_prm[n_prm++].n = Sext;
  // bn_prm[n_prm].Fnum = get_Fnum("s12"); bn_prm[n_prm++].n = Sext;

  // bn_prm[n_prm].Fnum = get_Fnum("oct1"); bn_prm[n_prm++].n = Oct;
  // bn_prm[n_prm].Fnum = get_Fnum("oct2"); bn_prm[n_prm++].n = Oct;
  // bn_prm[n_prm].Fnum = get_Fnum("oct3"); bn_prm[n_prm++].n = Oct;
  // bn_prm[n_prm].Fnum = get_Fnum("oct4"); bn_prm[n_prm++].n = Oct;
  // bn_prm[n_prm].Fnum = get_Fnum("oct5"); bn_prm[n_prm++].n = Oct;

  // bn_prm[n_prm].Fnum = get_Fnum("dec1"); bn_prm[n_prm++].n = Dec;
  // bn_prm[n_prm].Fnum = get_Fnum("dec2"); bn_prm[n_prm++].n = Dec;
  // bn_prm[n_prm].Fnum = get_Fnum("dec3"); bn_prm[n_prm++].n = Dec;
  // bn_prm[n_prm].Fnum = get_Fnum("dec4"); bn_prm[n_prm++].n = Dec;
  // bn_prm[n_prm].Fnum = get_Fnum("dec5"); bn_prm[n_prm++].n = Dec;

  if (true) {
    std::cout << std::endl;
    for (k = 0; k < n_prm; k++) {
      std::cout << std::setw(3) << k+1
	   << ", " << std::setw(6) << get_Name(bn_prm[k].Fnum)
	   << ", " << std::setw(2) << bn_prm[k].Fnum
	   << ", n = " << bn_prm[k].n << std::endl;
    }
  }
}


template<typename T>
ss_vect<T> x_px2x_xp(const ss_vect<T> map)
{
  // Transform from [x, px, y, py, ct, delta] to [x, x', y, y', ct, delta]
  int        k;
  ss_vect<T> Id, Id_scl, map1;

  Id.identity(); Id_scl.identity(); map1 = map;
  for (k = 0; k < 2; k++)
    Id_scl[2*k+1] /= 1e0 + Id[delta_];
  map1 = map*Id_scl;
  for (k = 0; k < 2; k++)
    map1[2*k+1] *= 1e0 + Id[delta_];

  return map1;
}


void get_twiss(const double alpha[], const double beta[],
	       const double eta[], const double etap[])
{
  int      j, k;
  double   alpha1[2], beta1[2], eta1[2], etap1[2], dnu1[2], dnu2[2];
  std::ofstream outf;

  const std::string file_name = "linlat.out";

  // Crucial; to only store linear part of A.
  danot_(1);

  file_wr(outf, file_name.c_str());

  for (k = 0; k < 2; k++) dnu1[k] = 0e0;
  A1 = get_A(alpha, beta, eta, etap);
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
    for (k = 0; k < 2; k++) dnu1[k] = dnu2[k];

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


void fit_ksi(const double ksi_x, const double ksi_y,
	     const int n_b3, const int b3s[])
{
  int    i, j;
  double L, ksi[2], *b, *b3, *db3, *b3_max, **A;
  tps    nu[2];

  b = dvector(1, 2); b3 = dvector(1, n_b3); db3 = dvector(1, n_b3);
  b3_max = dvector(1, n_b3); A = dmatrix(1, 2, 1, n_b3);

  danot_(3);

  std::cout << std::endl;
  for (i = 1; i <= n_b3; i++) {
    std::cout << "b3: " << b3s[i-1] << std::endl;
    L = get_L(b3s[i-1], 1);
    if (L == 0e0) L = 1e0;
    b3_max[i] = b3L_max/L; b3[i] = get_bn(b3s[i-1], 1, Sext);
  }

  get_ksi(ksi);

  for (i = 1; i <= n_b3; i++) {
    set_bn_par(b3s[i-1], Sext, 7);

    get_nu(nu);
    A[1][i] = h_ijklm_p(nu[X_], 0, 0, 0, 0, 1, 7);
    A[2][i] = h_ijklm_p(nu[Y_], 0, 0, 0, 0, 1, 7);
    if (i == n_b3) {
      b[1] = -h_ijklm(nu[X_], 0, 0, 0, 0, 1);
      b[2] = -h_ijklm(nu[Y_], 0, 0, 0, 0, 1);
    }

    clr_bn_par(b3s[i-1], Sext);
  }

  std::cout << std::endl;
  for (j = 1; j <= 2; j++) {
    for (i = 1; i <= n_b3; i++)
      std::cout << std::scientific << std::setprecision(3)
		<< std::setw(11) << A[j][i];
    std::cout << std::scientific << std::setprecision(3)
	      << std::setw(11) << b[j] << std::endl;
  }

  SVD_lim(2, n_b3, A, b, b3_max, s_cut, b3, db3);

  std::cout << "b3L:" << std::endl;
  for (i = 1; i <= n_b3; i++) {
    set_dbn(b3s[i-1], Sext, db3[i]);
    std::cout << std::scientific << std::setprecision(6)
	 << std::setw(14) << get_bn(b3s[i-1], 1, Sext);
  }

  get_ksi(ksi);

  free_dvector(b, 1, 2); free_dvector(b3, 1, n_b3); free_dvector(db3, 1, n_b3);
  free_dvector(b3_max, 1, n_b3); free_dmatrix(A, 1, 2, 1, n_b3);
}


void get_ct(void)
{
  ss_vect<tps> ps;

  ps.identity(); ps.propagate(1, n_elem);
  std::cout << std::fixed << std::setprecision(5)
	    << std::endl << "R_56 = " << h_ijklm(ps[ct_], 0, 0, 0, 0, 1) << ", "
	    << "T_566 = " << h_ijklm(ps[ct_], 0, 0, 0, 0, 2) << std::endl;
}


void fit_ct(const double ct, const int n_bn, const int fam_bn[],
	    const int order_bn[])
{
  const int m = 5;

  int          i, j, m1;
  double       L, *b, *bn, *dbn, *bn_max, **A;
  ss_vect<tps> nus;

  const double s_cut = 1e-8;

  b = dvector(1, m); bn = dvector(1, n_bn); dbn = dvector(1, n_bn);
  bn_max = dvector(1, n_bn); A = dmatrix(1, m, 1, n_bn);

  danot_(3);

  std::cout << std::endl;
  for (i = 1; i <= n_bn; i++) {
    std::cout << "bn: " << fam_bn[i-1] << std::endl;
    L = get_L(fam_bn[i-1], 1);
    if (L == 0e0) L = 1e0;
    bn_max[i] = b3L_max/L; bn[i] = get_bn(fam_bn[i-1], 1, order_bn[i-1]);
  }

  for (i = 1; i <= n_bn; i++) {
    set_bn_par(fam_bn[i-1], order_bn[i-1], 7);

    Map.identity(); Map.propagate(1, n_elem);

    K = MapNorm(Map, g, A1, A0, Map_res, 1); nus = dHdJ(K);

    m1 = 0;
    A[++m1][i] = h_ijklm_p(nus[3],       0, 0, 0, 0, 1, 7);
    A[++m1][i] = h_ijklm_p(nus[4],       0, 0, 0, 0, 1, 7);
    A[++m1][i] = h_ijklm_p(Map[ct_],     0, 0, 0, 0, 2, 7);
    if (i == n_bn) {
      m1 = 0;
      b[++m1] = -h_ijklm(nus[3],       0, 0, 0, 0, 1);
      b[++m1] = -h_ijklm(nus[4],       0, 0, 0, 0, 1);
      b[++m1] = ct - h_ijklm(Map[ct_], 0, 0, 0, 0, 2);
    }

    clr_bn_par(fam_bn[i-1], order_bn[i-1]);
  }

  std::cout << std::endl;
  for (j = 1; j <= m1; j++) {
    for (i = 1; i <= n_bn; i++)
      std::cout << std::scientific << std::setprecision(3)
		<< std::setw(11) << A[j][i];
    std::cout << std::scientific << std::setprecision(3)
	      << std::setw(11) << b[j] <<std:: endl;
  }

  SVD_lim(m1, n_bn, A, b, bn_max, s_cut, bn, dbn);

  std::cout << "bnL:" << std::endl;
  for (i = 1; i <= n_bn; i++) {
    set_dbn(fam_bn[i-1], order_bn[i-1], dbn[i]);
    std::cout << std::scientific << std::setprecision(6)
	 << std::setw(14) << get_bn(fam_bn[i-1], 1, order_bn[i-1]);
  }
  std::cout << std::endl;

  free_dvector(b, 1, m); free_dvector(bn, 1, n_bn); free_dvector(dbn, 1, n_bn);
  free_dvector(bn_max, 1, n_bn); free_dmatrix(A, 1, m, 1, n_bn);
}


void fit_eta(const int n_bn, const int fam_bn[], const int order_bn[])
{
  const int m = 7;

  int          i, j, k, m1;
  double       L, *b, *bn, *dbn, *bn_max, **A;
  tps          K_scl;
  ss_vect<tps> Id_scl, map_scl;

  const double s_cut = 1e-15, bnL_max = 1e7;

  b = dvector(1, m); bn = dvector(1, n_bn); dbn = dvector(1, n_bn);
  bn_max = dvector(1, n_bn); A = dmatrix(1, m, 1, n_bn);

  Id_scl.identity(); Id_scl[delta_] *= fabs(delta);

  std::cout << std::endl;
  for (i = 1; i <= n_bn; i++) {
    std::cout << "bn: " << fam_bn[i-1] << std::endl;
    L = get_L(fam_bn[i-1], 1);
    if (L == 0e0) L = 1e0;
    bn_max[i] = bnL_max/L; bn[i] = get_bn(fam_bn[i-1], 1, order_bn[i-1]);
  }

  for (k = 1; k <= 3; k++) {
    for (i = 1; i <= n_bn; i++) {
      set_bn_par(fam_bn[i-1], order_bn[i-1], 7);

      danot_(no_tps-1);
      Map.identity(); Map.propagate(1, n_elem);
      danot_(no_tps);

      K = MapNorm(Map, g, A1, A0, Map_res, 1);

      Map.identity();
      Map[x_] += h_ijklm(A0[x_], 0, 0, 0, 0, 1)*tps(0e0, delta_+1);
      Map.propagate(1, n_elem);

      map_scl = Map*Id_scl; K_scl = K*Id_scl;

      m1 = 0;
      A[++m1][i] = h_ijklm_p(map_scl[x_], 0, 0, 0, 0, 3, 7);
      A[++m1][i] = h_ijklm_p(K_scl,       2, 2, 0, 0, 0, 7);
      A[++m1][i] = h_ijklm_p(K_scl,       0, 0, 2, 2, 0, 7);
      A[++m1][i] = h_ijklm_p(K_scl,       1, 1, 1, 1, 0, 7);
      if (i == n_bn) {
	m1 = 0;
	b[++m1] = -h_ijklm(map_scl[x_], 0, 0, 0, 0, 3);
	b[++m1] = -h_ijklm(K_scl,       2, 2, 0, 0, 0);
	b[++m1] = -h_ijklm(K_scl,       0, 0, 2, 2, 0);
	b[++m1] = -h_ijklm(K_scl,       1, 1, 1, 1, 0);
      }

      clr_bn_par(fam_bn[i-1], order_bn[i-1]);
    }

    std::cout << std::endl;
    for (j = 1; j <= m1; j++) {
      for (i = 1; i <= n_bn; i++)
	std::cout << std::scientific << std::setprecision(3)
		  << std::setw(11) << A[j][i];
      std::cout << std::scientific << std::setprecision(3)
		<< std::setw(11) << b[j] << std::endl;
    }

    SVD_lim(m1, n_bn, A, b, bn_max, s_cut, bn, dbn);

    std::cout << "bnL:" << std::endl;
    for (i = 1; i <= n_bn; i++) {
      set_dbn(fam_bn[i-1], order_bn[i-1], dbn[i]);
      std::cout << std::scientific << std::setprecision(6)
	   << std::setw(14) << get_bn(fam_bn[i-1], 1, order_bn[i-1]);
    }
    std::cout << std::endl;
  }

  free_dvector(b, 1, m); free_dvector(bn, 1, n_bn); free_dvector(dbn, 1, n_bn);
  free_dvector(bn_max, 1, n_bn); free_dmatrix(A, 1, m, 1, n_bn);
}


int main(int argc, char *argv[])
{
  const int n_b3_max = 5, n_bn_max = 8;

  int             k, n_b3, b3s[n_b3_max];
  int             n_bn, fam_bn[n_bn_max], order_bn[n_bn_max];
  double          alpha[2], beta[2], eta[2], etap[2], ksi[2];
  ss_vect<double> ps;
  ss_vect<tps>    ps1, nus;

  rad_on    = false; H_exact        = false; totpath_on   = false;
  cavity_on = false; quad_fringe_on = false; emittance_on = false;
  IBS_on    = false;

  std::string home_dir = "/home/johan/projects/src/";

  // E0 contains p_0 [GeV].
  E0 = 750e-3; // Energy in ARC3.
  gamma0 = sqrt(sqr(m_e)+sqr(1e9*E0))/m_e;
  beta0  = sqrt(1e0-1e0/sqr(gamma0));

  printf("\np0 = %12.5e, beta0 = %12.5e, gamma0 = %12.5e\n",
	 1e9*E0, beta0, gamma0);

  rd_mfile((home_dir+"flat_file.dat").c_str(), elem);
  rd_mfile((home_dir+"flat_file.dat").c_str(), elem_tps);
  
  // Initialize the symplectic integrator after energy has been defined
  ini_si();

  // disable from TPSALib and LieLib log messages
  idprset(-1);

  for (k = 0; k < 2; k++) {
    eta[k] = 0e0; etap[k] = 0e0; alpha[k] = 0e0; beta[k] = 0e0;
  }
  switch (1) {
  case 1:
    // Arc3.
    // beta[X_] = 5.9942; beta[Y_] = 1.8373;
    // beta[X_] = 0.35751; beta[Y_] = 5.64898; eta[X_] = 0.0698;
    if (false) {
      // beta[X_] = 4.73529; beta[Y_] = 1.55552; eta[X_] = 0.31861;
      beta[X_] = 1.92895; beta[Y_] = 1.01853; eta[X_] = 0.6;
    } else {
      // JB_1.lat.
      // Fine tuned optics.
      // beta[X_] = 6.01442; beta[Y_] =1.78576 ;
      if (false) {
      // Dispersion wave.
	beta[X_] = 8.70141; beta[Y_] = 2.00042;
      } else {
	// Matching section.
	// beta[X_]  = 1.52376; beta[Y_]  = 12.34581;
	// alpha[X_] = 3.35393; alpha[Y_] =  5.27033;
	beta[X_]  = 2.70756; beta[Y_]  = 10.31726;
	alpha[X_] = 1.09794; alpha[Y_] =  3.45210;
      }
      // 8-unit cell structure.
      // alpha[X_] = 0.92093; alpha[Y_] = -0.61452;
      // beta[X_]  = 1.00644; beta[Y_]  =  2.70480;
      // 4-unit cell structure.
      // alpha[X_] = 1.20744; alpha[Y_] = -0.18567;
      // beta[X_]  = 2.84388; beta[Y_]  =  3.90905;
    }
    break;
   case 2:
     // Chicane2arc3combiner2.
    alpha[X_] = -0.25; alpha[Y_] = 0.00; beta[X_] =  2.50;  beta[Y_] = 5.00;
    break;
   case 3:
     // Arc1.
    beta[X_] = 1.2413; beta[Y_] = 1.4246;
    break;
  }

  danot_(1);
  get_Map();
  prt_lin_map(3, Map);
  get_twiss(alpha, beta, eta, etap);

  prt_lat("linlat.out", 10);
  prt_lat("linlat1.out");

  danot_(2);
  get_ksi(ksi); get_ct();
  prt_ct();

  if (true) {
    danot_(3);
    get_Map();
    K = MapNorm(Map, g, A1, A0, Map_res, 1); nus = dHdJ(K);
    ksi[X_] = h_ijklm(nus[3], 0, 0, 0, 0, 1);
    ksi[Y_] = h_ijklm(nus[4], 0, 0, 0, 0, 1);

    std::cout << std::fixed << std::setprecision(5)
	 << "\nksi = [" << ksi[X_] << ", " << ksi[Y_] << "]" << "\n";
  }

  if (false) {
    n_bn = 0;
    fam_bn[n_bn] = get_Fnum("s1"); order_bn[n_bn++] = 3;
    fam_bn[n_bn] = get_Fnum("s2"); order_bn[n_bn++] = 3;
    fam_bn[n_bn] = get_Fnum("s3"); order_bn[n_bn++] = 3;
    fam_bn[n_bn] = get_Fnum("s4"); order_bn[n_bn++] = 3;

    fit_ct(-2.5, n_bn, fam_bn, order_bn);

    // exit(0);
  }

  if (false) {
    n_bn = 0;
    fam_bn[n_bn] = get_Fnum("s1"); order_bn[n_bn++] = 4;
    fam_bn[n_bn] = get_Fnum("s2"); order_bn[n_bn++] = 4;
    fam_bn[n_bn] = get_Fnum("s3"); order_bn[n_bn++] = 4;
    fam_bn[n_bn] = get_Fnum("s4"); order_bn[n_bn++] = 4;

    fit_eta(n_bn, fam_bn, order_bn);

    // exit(0);
  }

  danot_(no_tps);
  ps1.identity();
  if (false) {
    // Get periodic solution for one cell.
    danot_(no_tps);
    Map.identity();
    if (false)
      Map.propagate(1, 19);
    else
      Map.propagate(1, n_elem);
    GoFix(Map, A0, A0_inv, no_tps);
    // A0[x_]  -= h_ijklm(A0[x_],  0, 0, 0, 0, 1)*tps(0e0, 5);
    // A0[px_] -= h_ijklm(A0[px_], 0, 0, 0, 0, 1)*tps(0e0, 5);
    printf("\neta1  = %13.5e, etap1 = %13.5e\n",
	   h_ijklm(A0[x_],  0, 0, 0, 0, 1), h_ijklm(A0[px_],  0, 0, 0, 0, 1));
    printf("eta2  = %13.5e, etap2 = %13.5e\n",
	   h_ijklm(A0[x_],  0, 0, 0, 0, 2), h_ijklm(A0[px_],  0, 0, 0, 0, 2));
    for (k = 0; k < 2; k++)
      ps1[k] = A0[k];
  }
  prt_eta(ps1);

  if (false) {
    n_b3 = 0;
    b3s[n_b3++] = get_Fnum("s2"); b3s[n_b3++] = get_Fnum("s3");

    if (n_b3 > n_b3_max) {
      std::cout << std::endl << "n_b3_max exceeded: "
		<< n_b3 << " (" << n_b3_max << ")" << std::endl;
      exit(1);
    }

    if (false)
      fit_ksi(0e0, 0e0, n_b3, b3s);
    else
      // Periodic solution.
      fit_chrom(0e0, 0e0, n_b3, b3s, true);
  
    prt_eta(ps1);
    get_ct();
    prt_ct();
  }

  if (false) {
    get_prms();

    if (false) no_mpoles();

    get_twiss(alpha, beta, eta, etap);

    h_zero();

    get_ct(); get_ksi(ksi);
    prt_ct();
    danot_(no_tps);
    ps1.identity();
    prt_eta(ps1);

    prt_mfile("flat_file.dat");
  }
}
