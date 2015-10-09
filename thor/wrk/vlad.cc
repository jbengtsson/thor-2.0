#define NO 6

#include "thor_lib.h"

int no_tps   = NO,
    ndpt_tps = 5;

extern ss_vect<tps> Map;


typedef struct {
  int n;
  int Fnum;
 } bn_type;


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


void h_init(bool &first, const int n_bn,  const bn_type bns[],
	    const double bnL_max[], double *bn_max, double b3[], double *bn)
{
  int    i;
  double L;

  cout << endl;
  for (i = 1; i <= n_bn; i++) {
    // Note, Jacobian is a function of the multipole strengths.
    L = get_L(bns[i-1].Fnum, 1);
    if (L == 0e0) L = 1e0;
    bn_max[i] = bnL_max[bns[i-1].n]/L;
    if (true)
      // bn[i] = get_bn(bns[i-1].Fnum, 1, bns[i-1].n);
      bn[i] = 81.6515/2.0;
    else
      // Zero sextupoles.
      bn[i] = 0e0;
    cout << scientific << setprecision(5)
	 << setw(4) << get_Name(bns[0].Fnum) << setw(13) << bn[i]*L << endl;
  }

  if (first) {
    // Store initial sextupole strengths.
    first = false;
    cout << endl << "initial b3s:" << endl;
    for (i = 0; i < n_bn; i++) {
      b3[i] = get_bnL(elem[bns[i].Fnum].Fnum, 1, bns[i].n);
      cout << scientific << setprecision(3)
	   << setw(11) << b3[i] << setw(2) << bns[i].n << endl;
    }
  } else {
    // Initialize sextupoles.
    for (i = 0; i < n_bn; i++)
      set_bnL(elem[bns[i].Fnum].Fnum, bns[i].n, bn[i+1]);
  }
}


void h_opt(const int m, const int n_iter, double &chi2, const int n_bn,
	   const bn_type bns[], const double *bn_max, double *bn, double **A,
	   double *b, double *dbn, double &dbn_max, const double T566,
	   const double delta, const bool prt_iter, ofstream &sext_out)
{
  ostringstream hs[m];
  int           i, j, k;
  double        L, chi2_old;
  tps           h;
  ss_vect<tps>  Id_scl, lin_map, map_scl;

  const bool   prt = true;
  const int    n_delta = 4, n_prt = 9;
  const double step = 0.7, scl_T566 = 1e-2;

  Id_scl.identity(); Id_scl[delta_] *= delta;

  cout << endl;
  for (i = 1; i <= n_bn; i++) {
    set_bn_par(bns[i-1].Fnum, bns[i-1].n, 7);

    danot_(no_tps-1); get_Map(); danot_(no_tps);
    danot_(1); lin_map = Map; danot_(no_tps);
    h = LieFact(Map*Inv(lin_map));

    map_scl = Map*Id_scl;

    k = 1;
    for (j = 2; j <= n_delta; j++)
      A[k++][i] = h_ijklm_p(map_scl[x_], 0, 0, 0, 0, j, 7);
    for (j = 2; j <= n_delta; j++)
      A[k++][i] = h_ijklm_p(map_scl[px_], 0, 0, 0, 0, j, 7);
    A[k++][i] = scl_T566*h_ijklm_p(Map[ct_], 0, 0, 0, 0, 2, 7);

    if (i == n_bn) {
      k = 1;
      for (j = 2; j <= n_delta; j++) {
	hs[k-1] << "  x(delta^" << j << ")";
	b[k++] = -h_ijklm(map_scl[x_], 0, 0, 0, 0, j);
      }
      for (j = 2; j <= n_delta; j++) {
	hs[k-1] << "p_x(delta^" << j << ")";
	b[k++] = -h_ijklm(map_scl[px_], 0, 0, 0, 0, j);
      }
      hs[k-1] << "T566        ";
      b[k++] = -scl_T566*(h_ijklm(Map[ct_], 0, 0, 0, 0, 2)-T566);
    }

    clr_bn_par(bns[i-1].Fnum, bns[i-1].n);
  }

  chi2_old = chi2; chi2 = 0e0;
  for (j = 1; j <= m; j++)
    chi2 += sqr(b[j]);

  if (prt) {
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
	  cout << scientific << setprecision(2) << setw(10) << -(b[i]-T566)
	       << endl;
    }
  }
    
  SVD_lim(m, n_bn, A, b, bn_max, 1e-10, bn, dbn);

  cout << "db_3L:" << endl;
  dbn_max = 0e0;
  for (i = 1; i <= n_bn; i++) {
    L = get_L(bns[i-1].Fnum, 1);
    if (L == 0e0) L = 1e0;
    dbn_max = max(fabs(step*dbn[i]*L), dbn_max);
    cout << scientific << setprecision(3) << setw(11) << step*dbn[i]*L;
    if (i % n_prt == 0) cout << endl;
  }
  if (n_bn % n_prt != 0) cout << endl;

  cout << "b_3L:" << endl;
  for (i = 1; i <= n_bn; i++) {
    set_dbn(bns[i-1].Fnum, bns[i-1].n, step*dbn[i]);
    bn[i] = get_bn(bns[i-1].Fnum, 1, bns[i-1].n);
    cout << scientific << setprecision(3)
	 << setw(11) << get_bnL(bns[i-1].Fnum, 1, bns[i-1].n);
    if (i % n_prt == 0) cout << endl;
  }
  if (n_bn % n_prt != 0) cout << endl;

  if (prt_iter) {
    sext_out << endl;
    sext_out << "n = " << n_iter << ":" << endl;
    for (i = 1; i <= n_bn; i++)
      for (j = 1; j <= get_n_Kids(bns[i-1].Fnum); j++)
	sext_out << fixed << setprecision(7) 
		 << setw(9) << elem[bns[i-1].Fnum].Name << "(" << j << ") = "
		 << setw(11) << get_bnL(bns[i-1].Fnum, j, bns[i-1].n)
		 << setw(2) << bns[i-1].n << endl;
  }
}


void h_zero(const int m, const int n_bn, const bn_type bns[],
	    const double bn_tol, const int n_iter_max, const double T566,
	    const double delta, const bool prt_iter)
{
  string   str;
  bool     first;
  int      i, n_iter;
  double   dbn_max, bnL_max[mpole_max], b3[n_bn], chi2;
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

  first = false; h_init(first, n_bn, bns, bnL_max, bn_max, b3, bn);

  n_iter = 0; chi2 = 0e0;
  do {
    n_iter++;
    h_opt(m, n_iter, chi2, n_bn, bns, bn_max, bn, A, b, dbn, dbn_max,
	  T566, delta, prt_iter, sext_out);
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
    ps.zero(); ps[delta_] = d;
    ps.propagate(1, n_elem);
    outf << scientific << setprecision(6)
	 << setw(14) << d << setw(14) << ps << endl;
  }
  outf.close();
}


void analyze_nl_dyn(const double delta)
{
  int             k;
  tps             h;
  ss_vect<double> ps;
  ss_vect<tps>    Id_scl, lin_map;

  Id_scl.identity(); Id_scl[delta_] *= delta;

  danot_(no_tps-1); get_Map(); danot_(no_tps);
  danot_(1); lin_map = Map; danot_(no_tps);
  h = LieFact(Map*Inv(lin_map));

  prt_lin_map(3, Map);
  cout << endl << scientific << setprecision(6)
       << "R_56:    " << setw(14) << h_ijklm(Map[ct_], 0, 0, 0, 0, 1) << endl;
  cout << scientific << setprecision(6)
       << "T_566:   " << setw(14) << h_ijklm(Map[ct_], 0, 0, 0, 0, 2) << endl;

  ps.zero(); ps[delta_] = -delta;
  ps.propagate(1, n_elem);
  cout << endl << scientific << setprecision(3)
       << setw(11) << ps << endl;

  ps.zero(); ps[delta_] = -delta;
  cout << endl << scientific << setprecision(3)
       << setw(11) << (Map*ps).cst() << endl;
  cout << endl
       << fixed << setprecision(1) << "x0(" << -1e2*delta << "%) = "
       << scientific << setprecision(3) << (Map*ps).cst()[x_]
       << fixed << setprecision(1) << ", p_x0(" << -1e2*delta << "%) = "
       << scientific << setprecision(3) << (Map*ps).cst()[px_]  << endl;
 
  cout << endl << fixed << setprecision(1)
       << "x0(" << -1e2*delta << "%) [m]:" << endl;
  for (k = 1; k <= no_tps-1; k++) {
    cout << scientific << setprecision(6)
	 << "h_1000" << k+1 << ": "
	 << setw(14) << h_ijklm(Map[x_]*Id_scl, 0, 0, 0, 0, k) << endl;
  }

  cout << endl << fixed << setprecision(1)
       << "p_x0(" << -1e2*delta << "%) [rad]:" << endl;
  for (k = 1; k <= no_tps-1; k++) {
    cout << scientific << setprecision(6)
	 << "h_0100" << k+1 << ": "
	 << setw(14) << h_ijklm(Map[px_]*Id_scl, 0, 0, 0, 0, k) << endl;
  }
  // cout << scientific << setprecision(6) << setw(13) << h*Id_scl << endl;
  // cout << scientific << setprecision(6) << setw(13) << Map[x_] << endl;
}


int main(int argc, char *argv[])
{
  const int n_bn = 4, m = 7;

  bn_type bns[n_bn];

  const bool   prt_iter = true;
  const int    n_iter_max = 100;
  const double bn_tol = 1e-2, T566 = -1.125e0, delta = 5e-2;

  rad_on    = false; H_exact        = false; totpath_on   = false;
  cavity_on = false; quad_fringe_on = false; emittance_on = false;
  IBS_on    = false;

  string home_dir = "/home/johan/projects/src/";

  // E0 contains kinetic energy [eV].
  // E0 = 10e3;

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

  if (true) {
    get_prms(n_bn, bns);
    h_zero(m, n_bn, bns, bn_tol, n_iter_max, T566, delta, prt_iter);
    analyze_nl_dyn(delta);
    scan_delta(20, -delta);
  }
}
