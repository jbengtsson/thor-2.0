#define NO 4

#include "thor_lib.h"

int no_tps   = NO,
    ndpt_tps = 5;

extern ss_vect<tps> Map;


int    n_prt, Fnums[25];
double f_ref;


typedef struct {
  int n;
  int Fnum;
 } bn_type;

const int    n_prm = 4;
const double f_scl = 1e13;
int     n_iter;
bn_type bn_prm[n_prm];


int get_ind(const int k)
{
  int index[] = {x_, px_, y_, py_, ct_, delta_};
  return index[k];
}


void prt_lin_map(const int n_DOF, const ss_vect<tps> &map)
{
  int  i, j;

  cout << endl;
  for (i = 0; i < 2*n_DOF; i++) {
    for (j = 0; j < 2*n_DOF; j++) {
      if (true)
	cout << scientific << setprecision(6)
	     << setw(14) << map[get_ind(i)][get_ind(j)];
      else
	cout << scientific << setprecision(16)
	     << setw(24) << map[get_ind(i)][get_ind(j)];
    }
    cout << endl;
  }
}


void scan_delta(const int n, const double delta)
{
  int             k;
  double          d;
  ss_vect<double> ps;
  ofstream        outf;

  const string file_name = "scan_delta.out";

  file_wr(outf, file_name.c_str());

  for (k = 0; k < n; k++) {
    d = (double)k/double(n-1)*delta;
    ps.zero(); ps[delta_] = d; ps.propagate(1, n_elem);
    outf << scientific << setprecision(6)
	 << setw(14) << d << setw(14) << ps << endl;
  }
  outf.close();
}


void h_init(const int n_bn,  const bn_type bns[], const double bnL_max[],
	    double *bn_max, double *bn)
{
  int    i;
  double L;

  cout << endl << "b3L0:" << endl;
  for (i = 1; i <= n_bn; i++) {
    // Note, Jacobian is a function of the multipole strengths.
    L = get_L(bns[i-1].Fnum, 1);
    if (L == 0e0) L = 1e0;
    bn_max[i] = bnL_max[bns[i-1].n]/L;
    if (true)
      bn[i] = get_bn(bns[i-1].Fnum, 1, bns[i-1].n);
    else {
      // Zero sextupoles.
      // bn[i] = 0e0;
      bn[i] = 81.6515/2.0;
      set_bn(bns[i-1].Fnum, bns[i-1].n, bn[i]);
    }
    cout << scientific << setprecision(5)
	 << setw(13) << get_bnL(bns[i-1].Fnum, 1, bns[i-1].n);
  }
  cout << endl;
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


void prt_system(const int m, const int n_bn, double **A, double *b,
		const double T566, const int n_iter, const double chi2_old,
		const double chi2, ostringstream hs[])
{
  int i, j;

  cout << scientific << setprecision(1)
       << endl << n_iter << " Ax = b, chi2: " << chi2_old << " -> " << chi2
       << ":" << endl << endl;
  for (i = 1; i <= m; i++) {
    cout  << setw(3) << i << " " << hs[i-1].str();
    for (j = 1; j <= n_bn; j++)
      cout << scientific << setprecision(2) << setw(10) << A[i][j];
    if (i != m)
      cout << scientific << setprecision(2) << setw(10) << -b[i] << endl;
    else
      cout << scientific << setprecision(2) << setw(10) << -(b[i]-T566) << endl;
  }
}


void prt_sext(ofstream &sext_out, const int n_iter, const int n_bn,
	      const bn_type bns[])
{
  int i, j;

  sext_out << endl;
  sext_out << "n = " << n_iter << ":" << endl;
  for (i = 1; i <= n_bn; i++)
    for (j = 1; j <= get_n_Kids(bns[i-1].Fnum); j++)
      sext_out << fixed << setprecision(7) 
	       << setw(9) << elem[bns[i-1].Fnum].Name << "(" << j << ") = "
	       << setw(11) << get_bnL(bns[i-1].Fnum, j, bns[i-1].n)
	       << setw(2) << bns[i-1].n << endl;
}


double get_chi2(const int n, const double data[])
{
  int    j;
  double chi2 = 0e0;

  for (j = 1; j <= n; j++) chi2 += sqr(data[j]);

  return chi2;
}


template<typename T>
void f_opt1(const int n, const double delta, T &x2d_intgrl, T &px2d_intgrl)
{
  int        k;
  double     h, x_max, px_max;
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

  cout << scientific << setprecision(3)
       << "x_max = " << x_max << ", px_max = " << px_max<< endl;

  x2d_intgrl = trapz(n, h, x_delta); px2d_intgrl = trapz(n, h, px_delta);

  danot_(no_tps);
}


void h_opt1(const int m, const int n_iter, double &chi2, const int n_bn,
	    const bn_type bns[], const double *bn_max, double *bn, double **A,
	    double *b, double *dbn, double &dbn_max, const double T566,
	    const double delta, const bool prt_iter, ofstream &sext_out)
{
  ostringstream hs[m];
  int           i, k;
  double        L, chi2_old;
  tps           x2d_intgrl, px2d_intgrl;

  const bool   prt = true;
  const int    n_delta = 30;
  const double step = 0.3, scl_T566 = 1e-6, svd_tol = 1e-10;

  cout << endl;
  for (i = 1; i <= n_bn; i++) {
    set_bn_par(bns[i-1].Fnum, bns[i-1].n, 7);

    f_opt1(n_delta, -delta, x2d_intgrl, px2d_intgrl);

    danot_(no_tps-1); get_Map(); danot_(no_tps);

    k = 1;
    A[k++][i] = h_ijklm_p(x2d_intgrl, 0, 0, 0, 0, 0, 7);
    A[k++][i] = h_ijklm_p(px2d_intgrl, 0, 0, 0, 0, 0, 7);
    A[k++][i] = scl_T566*h_ijklm_p(Map[ct_], 0, 0, 0, 0, 2, 7);

    if (i == n_bn) {
      k = 1;
      hs[k-1] << " x-delta area"; b[k++] = -h_ijklm(x2d_intgrl, 0, 0, 0, 0, 0);
      hs[k-1] << "px-delta area"; b[k++] = -h_ijklm(px2d_intgrl, 0, 0, 0, 0, 0);
      hs[k-1] << "T566         ";
      b[k++] = -scl_T566*(h_ijklm(Map[ct_], 0, 0, 0, 0, 2)-T566);
    }

    clr_bn_par(bns[i-1].Fnum, bns[i-1].n);
  }

  chi2_old = chi2; chi2 = get_chi2(m, b);

  if (prt) prt_system(m, n_bn, A, b, T566, n_iter, chi2_old, chi2, hs);
    
  SVD_lim(m, n_bn, A, b, bn_max, svd_tol, bn, dbn);

  cout << "db3L (step*dcorr*L):" << endl;
  dbn_max = 0e0;
  for (i = 1; i <= n_bn; i++) {
    L = get_L(bns[i-1].Fnum, 1);
    if (L == 0e0) L = 1e0;
    dbn_max = max(fabs(step*dbn[i]*L), dbn_max);
    cout << scientific << setprecision(3) << setw(11) << step*dbn[i]*L;
  }
  cout << endl;

  cout << "b3L:" << endl;
  for (i = 1; i <= n_bn; i++) {
    set_dbn(bns[i-1].Fnum, bns[i-1].n, step*dbn[i]);
    bn[i] = get_bn(bns[i-1].Fnum, 1, bns[i-1].n);
    cout << scientific << setprecision(3)
	 << setw(11) << get_bnL(bns[i-1].Fnum, 1, bns[i-1].n);
  }
  cout << endl;

  if (prt_iter) prt_sext(sext_out, n_iter, n_bn, bns);
}


void h_zero1(const int m, const int n_bn, const bn_type bns[],
	     const double bn_tol, const int n_iter_max, const double T566,
	     const double delta, const bool prt_iter)
{
  string   str;
  int      i, n_iter;
  double   dbn_max, bnL_max[mpole_max], chi2;
  double   **A, *b, *w, **U, **V, *dbn, *bn, *bn_max, **A_inv;
  ofstream sext_out;

  b = dvector(1, m); w = dvector(1, n_bn); dbn = dvector(1, n_bn);
  bn = dvector(1, n_bn); bn_max = dvector(1, n_bn);
  A = dmatrix(1, m, 1, n_bn); U = dmatrix(1, m, 1, n_bn);
  V = dmatrix(1, n_bn, 1, n_bn);
  A_inv = dmatrix(1, n_bn, 1, m);

  // Max sextupole strength.
  bnL_max[Sext] = 50e0;

  file_wr(sext_out, "sext.dat");

  h_init(n_bn, bns, bnL_max, bn_max, bn);

  n_iter = 0; chi2 = 0e0;
  do {
    n_iter++;
    h_opt1(m, n_iter, chi2, n_bn, bns, bn_max, bn, A, b, dbn, dbn_max,
	   T566, delta, prt_iter, sext_out);

   scan_delta(20, -delta);
  } while ((dbn_max > bn_tol) && (n_iter < n_iter_max));

  if (!prt_iter) {
    for (i = 1; i <= n_bn; i++)

    for (i = 1; i <= n_bn; i++)
      sext_out << fixed << setprecision(7) 
	       << setw(9) << get_Name(bns[i-1].Fnum)
	       << "(" << elem[bns[i-1].Fnum].Knum << ") = "
	       << setw(11) << get_bnL(bns[i-1].Fnum, 1, bns[i-1].n)
	       << setw(2) << bns[i-1].n << endl;
    sext_out.flush();
  }

  free_dvector(b, 1, m); free_dvector(w, 1, n_bn);
  free_dvector(dbn, 1, n_bn); free_dvector(bn, 1, n_bn);
  free_dvector(bn_max, 1, n_bn);
  free_dmatrix(A, 1, m, 1, n_bn);
  free_dmatrix(U, 1, m, 1, n_bn); free_dmatrix(V, 1, n_bn, 1, n_bn);
  free_dmatrix(A_inv, 1, n_bn, 1, m);
}


double f_opt2(double b3s[])
{
  int          i;
  double       x2d_intgrl, px2d_intgrl, f, dT566;

  const int    n_delta = 25;
  const double delta = -5e-2, T566 = -1.125e0, scl_T566 = 1e-8;

  for (i = 1; i <= n_prm; i++)
    set_bn(bn_prm[i-1].Fnum, bn_prm[i-1].n, b3s[i]);

  f_opt1(n_delta, delta, x2d_intgrl, px2d_intgrl);

  danot_(no_tps-1); get_Map(); danot_(no_tps);
  dT566 = h_ijklm(Map[ct_], 0, 0, 0, 0, 2) - T566;
  cout << scientific << setprecision(3) << "dT566 = " << dT566 << endl;

  f = scl_T566*sqr(dT566) + sqr(x2d_intgrl) + sqr(px2d_intgrl);

  f *= f_scl;

  cout << endl << scientific << setprecision(3) << setw(4) << n_iter
       << " f = " << f << setw(11) << dT566 << setw(11) << x2d_intgrl
       << setw(11) << px2d_intgrl << endl;
  cout << "b3s:" << endl;
  for (i = 1; i <= n_prm; i++)
    cout << scientific << setprecision(10) << setw(18) << b3s[i];
  cout << endl;

  return f;
}


void df_opt2(double b3s[], double df[])
{
  int i;
  tps x2d_intgrl, px2d_intgrl, dT566, f;

  const int    n_delta = 25;
  const double delta = -5e-2, T566 = -1.125e0, scl_T566 = 1e-12;

  for (i = 1; i <= n_prm; i++)
    set_bn(bn_prm[i-1].Fnum, bn_prm[i-1].n, b3s[i]);

  for (i = 1; i <= n_prm; i++) {
    set_bn_par(bn_prm[i-1].Fnum, bn_prm[i-1].n, 7);

    f_opt1(n_delta, delta, x2d_intgrl, px2d_intgrl);
    danot_(no_tps-1); get_Map(); danot_(no_tps);
    dT566 =
      h_ijklm(Map[ct_], 0, 0, 0, 0, 2) - T566
      + h_ijklm_p(Map[ct_], 0, 0, 0, 0, 2, 7)*tps(0e0, 7);

    f = scl_T566*sqr(dT566) + sqr(x2d_intgrl) + sqr(px2d_intgrl);

    f *= f_scl;

    df[i] = h_ijklm_p(f, 0, 0, 0, 0, 0, 7);

    clr_bn_par(bn_prm[i-1].Fnum, bn_prm[i-1].n);
  }

  cout << scientific << setprecision(3) << "dT566 = " << dT566 << endl;
  cout << endl << setw(4) << n_iter << " df: ";
  for (i = 1; i <= n_prm; i++)
      cout << scientific << setprecision(3) << setw(11) << df[i];
  cout << endl;
  cout << "b3s: ";
  for (i = 1; i <= n_prm; i++)
    cout << scientific << setprecision(10) << setw(18) << b3s[i];
  cout << endl;

}


void h_zero2(const int n_bn, const bn_type bns[])
{
  int    i;
  double *b3s, fret;

  const double ftol = 1e-20;

  b3s = dvector(1, n_bn);

  for (i = 1; i <= n_bn; i++) {
    bn_prm[i-1] = bns[i-1];
    b3s[i] = get_bn(bns[i-1].Fnum, 1, bns[i-1].n);
  }

  n_iter = 0;
  cout << endl;
  dfrprmn(b3s, n_bn, ftol, &n_iter, &fret, f_opt2, df_opt2);
  // d_dfpmin(b3s, n_bn, ftol, &n_iter, &fret, f_opt2, df_opt2);

  free_dvector(b3s, 1, n_bn);
}


void h_zero3(const int n_bn, const bn_type bns[])
{
  int    i, j;
  double *b3s, **xi, fret;

  const double ftol = 1e-5;

  b3s = dvector(1, n_bn); xi = dmatrix(1, n_prm, 1, n_prm);

  for (i = 1; i <= n_bn; i++) {
    bn_prm[i-1] = bns[i-1];
    b3s[i] = get_bn(bns[i-1].Fnum, 1, bns[i-1].n);
  }

  for (i = 1; i <= n_prm; i++) {
    for (j = 1; j <= n_prm; j++)
      if (i == j)
        xi[i][j] = 0.1e0;
      else
        xi[i][j] = 0e0;
  }

  cout << endl;
  n_iter = 0; f_ref = 1e30;
  dpowell(b3s, xi, n_prm, ftol, &n_iter, &fret, f_opt2);

  free_dvector(b3s, 1, n_bn); free_dmatrix(xi, 1, n_prm, 1, n_prm);
}


void get_prms(const int n_bn, bn_type bns[])
{
  int k;

  // Arc3.
  k = 0;
  bns[k].Fnum = get_Fnum("a3s1"); bns[k++].n = Sext;
  bns[k].Fnum = get_Fnum("a3s2"); bns[k++].n = Sext;
  bns[k].Fnum = get_Fnum("a3s3"); bns[k++].n = Sext;
  bns[k].Fnum = get_Fnum("a3s4"); bns[k++].n = Sext;

  if (true) {
    cout << endl;
    for (k = 0; k < n_bn; k++) {
      cout << setw(3) << k+1
	   << ", " << setw(4) << get_Name(bns[k].Fnum)
	   << ", " << setw(2) << bns[k].Fnum
	   << ", n = " << bns[k].n << endl;
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


void analyze_nl_dyn(const double delta)
{
  int             k;
  tps             h;
  ss_vect<double> ps;
  ss_vect<tps>    Id, Id_scl, lin_map;

  Id.identity();
  Id_scl.identity(); Id_scl[delta_] *= delta;

  danot_(no_tps-1); get_Map(); danot_(no_tps);
  danot_(1); lin_map = Map; danot_(no_tps);
  h = LieFact(Map*Inv(lin_map));

  prt_lin_map(3, Map);
  cout << endl << scientific << setprecision(6)
       << "R_56:    " << setw(14) << h_ijklm(Map[ct_], 0, 0, 0, 0, 1) << endl;
  cout << scientific << setprecision(6)
       << "T_566:   " << setw(14) << h_ijklm(Map[ct_], 0, 0, 0, 0, 2) << endl;

  ps.zero(); ps[delta_] = -delta; ps.propagate(1, n_elem);
  cout << endl << scientific << setprecision(3) << setw(11) << ps << endl;

  cout << endl << fixed << setprecision(1)
       << "x0(" << -1e2*delta << "%) [m]:" << endl;
  for (k = 2; k <= no_tps-1; k++) {
    cout << scientific << setprecision(6) << "  x0(delta^" << k << ") = "
	 << setw(14) << h_ijklm(Map[x_]*Id_scl, 0, 0, 0, 0, k) << endl;
  }

  cout << endl << fixed << setprecision(1)
       << "px0(" << -1e2*delta << "%) [rad]:" << endl;
  for (k = 2; k <= no_tps-1; k++) {
    cout << scientific << setprecision(6) << "px0(delta^" << k << ") = "
	 << setw(14) << h_ijklm(Map[px_]*Id_scl, 0, 0, 0, 0, k) << endl;
  }
  // cout << scientific << setprecision(6) << setw(13) << h*Id_scl << endl;
  // danot_(no_tps-1);
  // cout << scientific << setprecision(6) << setw(13) << x_px2x_xp(Map)
  //      << endl;
  // cout << scientific << setprecision(6) << setw(13) << Map << endl;
}


int main(int argc, char *argv[])
{
  const int n_bn = 4, m = 3;

  bn_type bns[n_bn];
  ss_vect<double> ps;

  const bool   prt_iter = true;
  const int    n_iter_max = 1000;
  const double bn_tol = 1e-3, T566 = -1.125e0, delta = 5e-2;

  rad_on    = false; H_exact        = false;  totpath_on   = false;
  cavity_on = false; quad_fringe_on = false; emittance_on = false;
  IBS_on    = false;

  string home_dir = "/home/johan/projects/src/";

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

  if (false) {
    analyze_nl_dyn(delta);
    scan_delta(20, -delta);
    exit(0);
  }

  // opt_nl_disp();

  if (true) {
    get_prms(n_bn, bns);

    if (false)
      h_zero1(m, n_bn, bns, bn_tol, n_iter_max, T566, delta, prt_iter);
    else
      h_zero2(n_bn, bns);

    analyze_nl_dyn(delta);
    scan_delta(20, -delta);
  }
}
