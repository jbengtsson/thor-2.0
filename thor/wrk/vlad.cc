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

const int    n_prm = 4, m = 3, n_delta = 30;
const double delta = -5e-2, T566_ref = -1.125e0, scl_T566 = 1e-19;
const double f_scl = 1e13;

int     n_iter;
bn_type bn_prm[n_prm];


void no_mpoles(void)
{
  int j, k;

  cout << endl;
  cout << "zeroing multipoles" << endl;
  cout << endl;
  for (j = 0; j < n_elem; j++)
    if (elem[j].kind == Mpole)
      for (k = Sext; k < mpole_max; k++) {
	//	cout << "zeroing " << elem[j].Name << endl;
	set_bn(elem[j].Fnum, elem[j].Knum, k, 0.0);
      }
}


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


void scan_delta(const int n)
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


void h_init(const double bnL_max[], double *bn_max, double *bn)
{
  int    i;
  double L;

  cout << endl << "b3L0:" << endl;
  for (i = 1; i <= n_prm; i++) {
    // Note, Jacobian is a function of the multipole strengths.
    L = get_L(bn_prm[i-1].Fnum, 1);
    if (L == 0e0) L = 1e0;
    bn_max[i] = bnL_max[bn_prm[i-1].n]/L;
    if (true)
      bn[i] = get_bn(bn_prm[i-1].Fnum, 1, bn_prm[i-1].n);
    else {
      // Zero sextupoles.
      // bn[i] = 0e0;
      bn[i] = 81.6515/2.0;
      set_bn(bn_prm[i-1].Fnum, bn_prm[i-1].n, bn[i]);
    }
    cout << scientific << setprecision(5)
	 << setw(13) << get_bnL(bn_prm[i-1].Fnum, 1, bn_prm[i-1].n);
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


void prt_system(const int m, const int n_prm, double **A, double *b,
	        const double chi2_old, const double chi2, ostringstream hs[])
{
  int i, j;

  cout << scientific << setprecision(1)
       << endl << n_iter << " Ax = b, chi2: " << chi2_old << " -> " << chi2
       << ":" << endl << endl;
  for (i = 1; i <= m; i++) {
    cout  << setw(3) << i << " " << hs[i-1].str();
    for (j = 1; j <= n_prm; j++)
      cout << scientific << setprecision(2) << setw(10) << A[i][j];
    if (i != m)
      cout << scientific << setprecision(2) << setw(10) << -b[i] << endl;
    else
      cout << scientific << setprecision(2) << setw(10)
	   << -(b[i]-T566_ref) << endl;
  }
}


void prt_sext(ofstream &sext_out)
{
  int i, j;

  sext_out << endl;
  sext_out << "n = " << n_iter << ":" << endl;
  for (i = 1; i <= n_prm; i++)
    for (j = 1; j <= get_n_Kids(bn_prm[i-1].Fnum); j++)
      sext_out << fixed << setprecision(7) 
	       << setw(9) << elem[bn_prm[i-1].Fnum].Name << "(" << j << ") = "
	       << setw(11) << get_bnL(bn_prm[i-1].Fnum, j, bn_prm[i-1].n)
	       << setw(2) << bn_prm[i-1].n << endl;
}


double get_chi2(const int n, const double data[])
{
  int    j;
  double chi2 = 0e0;

  for (j = 1; j <= n; j++) chi2 += sqr(data[j]);

  return chi2;
}


template<typename T>
void f_opt1(const int n, T &x2d_intgrl, T &px2d_intgrl,
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


double f_opt2(double b3s[])
{
  int    i;
  double x2d_intgrl, px2d_intgrl, x_max, px_max, f, T566, dT566;

  for (i = 1; i <= n_prm; i++)
    set_bn(bn_prm[i-1].Fnum, bn_prm[i-1].n, b3s[i]);

  f_opt1(n_delta, x2d_intgrl, px2d_intgrl, x_max, px_max);
  danot_(no_tps-1); get_Map(); danot_(no_tps);

  T566 = h_ijklm(Map[ct_], 0, 0, 0, 0, 2); dT566 = T566 - T566_ref;

  f = scl_T566*sqr(dT566) + sqr(x2d_intgrl) + sqr(px2d_intgrl);

  cout << scientific << setprecision(3)
       << setw(4) << n_iter << " f = " << setw(9) << f
       << ", T566 = " << setw(10) << T566
       << ", x_max = " << setw(9) << x_max
       << ", px_max = " << setw(9) << px_max << endl;

  return f;
}


void h_opt1(const int m, double &chi2, const double *bn_max, double *bn,
	    double &dbn_max, double *g, double *h, const bool prt_iter,
	    ofstream &sext_out)
{
  // Conjugate gradient method.
  ostringstream hs[m];
  int           i, k;
  double        L, chi2_old, fret, x_max, px_max, g2, gamma, dg2;
  double        **A, *b, *w, **U, **V, *dbn;
  tps           x2d_intgrl, px2d_intgrl;

  const bool   prt = true;
  const double svd_tol = 1e-12;

  b = dvector(1, m); w = dvector(1, n_prm); dbn = dvector(1, n_prm);
  A = dmatrix(1, m, 1, n_prm); U = dmatrix(1, m, 1, n_prm);
  V = dmatrix(1, n_prm, 1, n_prm);

  cout << endl;
  for (i = 1; i <= n_prm; i++) {
    set_bn_par(bn_prm[i-1].Fnum, bn_prm[i-1].n, 7);

    f_opt1(n_delta, x2d_intgrl, px2d_intgrl, x_max, px_max);
    danot_(no_tps-1); get_Map(); danot_(no_tps);

    k = 1;
    A[k++][i] = h_ijklm_p(x2d_intgrl, 0, 0, 0, 0, 0, 7);
    A[k++][i] = h_ijklm_p(px2d_intgrl, 0, 0, 0, 0, 0, 7);
    A[k++][i] = sqrt(scl_T566)*h_ijklm_p(Map[ct_], 0, 0, 0, 0, 2, 7);

    if (i == n_prm) {
      k = 1;
      hs[k-1] << " x-delta area"; b[k++] = -h_ijklm(x2d_intgrl, 0, 0, 0, 0, 0);
      hs[k-1] << "px-delta area"; b[k++] = -h_ijklm(px2d_intgrl, 0, 0, 0, 0, 0);
      hs[k-1] << "T566         ";
      b[k++] = -sqrt(scl_T566)*(h_ijklm(Map[ct_], 0, 0, 0, 0, 2)-T566_ref);
    }

    clr_bn_par(bn_prm[i-1].Fnum, bn_prm[i-1].n);
  }

  chi2_old = chi2; chi2 = get_chi2(m, b);

  if (prt) prt_system(m, n_prm, A, b, chi2_old, chi2, hs);
    
  SVD_lim(m, n_prm, A, b, bn_max, svd_tol, bn, dbn);

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
      cout << "g.g = 0" << endl; exit(0);
    }
  }

  cout << endl;
  d_linmin(bn, dbn, n_prm, &fret, f_opt2);

  cout << endl << "db3L (dcorr*L):" << endl;
  dbn_max = 0e0;
  for (i = 1; i <= n_prm; i++) {
    L = get_L(bn_prm[i-1].Fnum, 1);
    if (L == 0e0) L = 1e0;
    dbn_max = max(fabs(dbn[i]*L), dbn_max);
    cout << scientific << setprecision(3) << setw(11) << dbn[i]*L;
  }
  cout << endl;

  cout << "b3L:" << endl;
  for (i = 1; i <= n_prm; i++) {
    set_bn(bn_prm[i-1].Fnum, bn_prm[i-1].n, bn[i]);
    bn[i] = get_bn(bn_prm[i-1].Fnum, 1, bn_prm[i-1].n);
    cout << scientific << setprecision(6)
	 << setw(14) << get_bnL(bn_prm[i-1].Fnum, 1, bn_prm[i-1].n);
  }
  cout << endl;

  if (prt_iter) prt_sext(sext_out);

  free_dvector(b, 1, m); free_dvector(w, 1, n_prm); free_dvector(dbn, 1, n_prm);
  free_dmatrix(A, 1, m, 1, n_prm); free_dmatrix(U, 1, m, 1, n_prm);
  free_dmatrix(V, 1, n_prm, 1, n_prm);
}


void h_zero1(void)
{
  string   str;
  int      i;
  double   *g, *h, dbn_max, bnL_max[mpole_max], chi2, *bn, *bn_max;
  ofstream sext_out;

  const bool   prt_iter   = true;
  const int    n_iter_max = 1000;
  const double bn_tol     = 1e-3;

  bn = dvector(1, n_prm); bn_max = dvector(1, n_prm);
  g = dvector(1, n_prm); h = dvector(1, n_prm);

  // Max sextupole strength.
  bnL_max[Sext] = 50e0;

  file_wr(sext_out, "sext.dat");

  h_init(bnL_max, bn_max, bn);

  n_iter = 0; chi2 = 0e0;
  do {
    n_iter++;
    h_opt1(m, chi2, bn_max, bn, dbn_max, g, h, prt_iter, sext_out);

    scan_delta(20);
  } while ((dbn_max > bn_tol) && (n_iter < n_iter_max));

  if (!prt_iter) {
    for (i = 1; i <= n_prm; i++)

    for (i = 1; i <= n_prm; i++)
      sext_out << fixed << setprecision(7) 
	       << setw(9) << get_Name(bn_prm[i-1].Fnum)
	       << "(" << elem[bn_prm[i-1].Fnum].Knum << ") = "
	       << setw(11) << get_bnL(bn_prm[i-1].Fnum, 1, bn_prm[i-1].n)
	       << setw(2) << bn_prm[i-1].n << endl;
    sext_out.flush();
  }

  free_dvector(bn, 1, n_prm); free_dvector(bn_max, 1, n_prm);
  free_dvector(g, 1, n_prm); free_dvector(h, 1, n_prm);
}


void get_prms(void)
{
  int k;

  // Arc3.
  k = 0;
  bn_prm[k].Fnum = get_Fnum("a3s1"); bn_prm[k++].n = Sext;
  bn_prm[k].Fnum = get_Fnum("a3s2"); bn_prm[k++].n = Sext;
  bn_prm[k].Fnum = get_Fnum("a3s3"); bn_prm[k++].n = Sext;
  bn_prm[k].Fnum = get_Fnum("a3s4"); bn_prm[k++].n = Sext;

  if (true) {
    cout << endl;
    for (k = 0; k < n_prm; k++) {
      cout << setw(3) << k+1
	   << ", " << setw(4) << get_Name(bn_prm[k].Fnum)
	   << ", " << setw(2) << bn_prm[k].Fnum
	   << ", n = " << bn_prm[k].n << endl;
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


void analyze_nl_dyn(void)
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


double get_code(elem_type<double> &elem)
{
  double  code;

  switch (elem.kind) {
  case Drift:
    code = 0e0;
    break;
  case Mpole:
    if (elem.mpole->h_bend != 0e0)
      code = 0.5e0;
    else if (elem.mpole->bn[Quad-1] != 0)
      code = sgn(elem.mpole->bn[Quad-1]);
    else if (elem.mpole->bn[Sext-1] != 0)
      code = 1.5*sgn(elem.mpole->bn[Sext-1]);
    else
      code = 0e0;
    break;
  default:
    code = 0e0;
    break;
  }

  return code;
}

#if 0

void prt_lat(const char *fname)
{
  long int i;
  FILE     *outf;

  outf = file_write(fname);
  fprintf(outf, "#        name           s   code"
	        "  alphax  betax   nux   etax   etapx");
  fprintf(outf, "  alphay  betay   nuy   etay   etapy    I5\n");
  fprintf(outf, "#                      [m]"
	        "                 [m]           [m]");
  fprintf(outf, "                   [m]           [m]\n");
  fprintf(outf, "#\n");

  for (i = 0; i <= n_elem; i++) {
    fprintf(outf, "%4ld %15s %9.5f %4.1f"
	    " %9.5f %8.5f %8.5f %8.5f %8.5f"
	    " %9.5f %8.5f %8.5f %8.5f %8.5f\n",
	    i, elem[i].Name, elem[i].S, get_code(elem[i]),
	    elem[i].Alpha[X_], elem[i].Beta[X_], elem[i].Nu[X_],
	    elem[i].Eta[X_], elem[i].Etap[X_],
	    elem[i].Alpha[Y_], elem[i].Beta[Y_], elem[i].Nu[Y_],
	    elem[i].Eta[Y_], elem[i].Etap[Y_]);
  }

  fclose(outf);
}

#endif

void get_twiss(const double alpha[], const double beta[],
	       const double eta[], const double etap[])
{
  int          j;
  ss_vect<tps> A1_A1tp;

  danot_(1);

  get_A1(alpha[X_], beta[X_], alpha[Y_], beta[Y_]);
  for (j = 1; j <= n_elem; j++) {
    A1.propagate(j, j); elem_tps[j].A1 = A1;
    A1_A1tp = A1*tp_S(2, A1);

    elem[j-1].Alpha[X_] = -h_ijklm(A1_A1tp[x_], 0, 1, 0, 0, 0);
    elem[j-1].Alpha[Y_] = -h_ijklm(A1_A1tp[y_], 0, 0, 0, 1, 0);
    elem[j-1].Beta[X_]  =  h_ijklm(A1_A1tp[x_], 1, 0, 0, 0, 0);
    elem[j-1].Beta[Y_]  =  h_ijklm(A1_A1tp[y_], 0, 0, 1, 0, 0);
  }
}


void ctrl_lin_opt()
{
  int    k;
  double alpha[2], beta[2], eta[2], etap[2];

  for (k = 0; k < 2; k++) {
    alpha[k] = 0e0; beta[k] = 0e0; eta[k] = 0e0; etap[k] = 0e0;
  }

  // ARC3.
  beta[X_] = 5.9942; beta[Y_] = 1.8373;

  get_twiss(alpha, beta, eta, etap);

  prt_lat("linlat.out");
}


int main(int argc, char *argv[])
{
  ss_vect<double> ps;

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
    analyze_nl_dyn();
    scan_delta(20);
    exit(0);
  }

  ctrl_lin_opt();

  // opt_nl_disp();

  if (false) {
    get_prms();

    if (false) no_mpoles();

    h_zero1();

    analyze_nl_dyn();
    scan_delta(20);
  }
}
