#define NO 3

#include "thor_lib.h"

int no_tps   = NO,
    ndpt_tps = 5;

extern ss_vect<tps> Map;

// Program to design matching cell an 8-cell higher order achromat.
// Use A script to propagate initial conditions for dispersion.

typedef struct {
  int      n;
  int      Fnum, Fnum2;
  elemkind type;
 } bn_type;

int m_cstr, n_prm;

const int    n_prm_max = 8, m_cstr_max = 11, n_delta = 30;
const double b2L_max = 10e0, L_max = 1.0, tr_cut = 1.95;


// Match ARC3.
const double alpha_[]  = {0e0, 0e0}, beta_[]  = {4.83023, 1.55324};
// const double eta_[] = {0.31861, 0e0}, etap_[]  = {0e0, 0e0};
const double eta_[] = {0.55140, 0e0}, etap_[]  = {0e0, 0e0};

int     n_iter;
bn_type bn_prm[n_prm_max];


void set_s_par1(const int Fnum, const int Knum, const int j)
{
  // Set s-dependence.
  long int  k;
  double    L;

  k = get_loc(Fnum, Knum) - 1;
  L = elem_tps[k].L.cst(); elem_tps[k].L = tps(L, j);
}


void set_s_par1(const int Fnum, const int j)
{
  // Set s-dependence.
  int  k;

  for (k = 1; k <= get_n_Kids(Fnum); k++)
    set_s_par1(Fnum, k, j);
}


void clr_s_par1(const int Fnum, const int Knum)
{
  // Clear s-dependence.
  int     k;
  double  L;

  k = get_loc(Fnum, Knum) - 1;
  L = elem_tps[k].L.cst(); elem_tps[k].L = L;
}


void clr_s_par1(const int Fnum)
{
  // Clear s-dependence.
  int  k;

  for (k = 1; k <= get_n_Kids(Fnum); k++)
    clr_s_par1(Fnum, k);
}


double get_bn_s1(const int Fnum, const int Knum, const int n)
{
  return (n != 0)? get_bn(Fnum, Knum, n) : elem[get_loc(Fnum, Knum)-1].L;
}


void set_bn_s1(const int Fnum, const int Knum, const int n, const double bn)
{
  long int  k = 0;

  if (n != 0)
    set_bn(Fnum, Knum, n, bn);
  else {
    k = get_loc(abs(Fnum), Knum) - 1;
    set_L(elem[k].Fnum, elem[k].Knum, bn);
  }
}


void set_bn_s1(const int Fnum, const int n, const double bn)
{
  int  k;

  for (k = 1; k <= get_n_Kids(abs(Fnum)); k++)
    set_bn_s1(Fnum, k, n, bn);
}


void prt_ct()
{
  int          k;
  ss_vect<tps> ps;
  std::ofstream     outf;

  const std::string file_name = "ct.out";

  file_wr(outf, file_name.c_str());

  std::cout << std::endl;
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


void get_ct(void)
{
  ss_vect<tps> ps;

  ps.identity(); ps.propagate(1, n_elem);
  std::cout << std::fixed << std::setprecision(5)
       << std::endl << "R_56 = " << h_ijklm(ps[ct_], 0, 0, 0, 0, 1) << ", "
       << "T_566 = " << h_ijklm(ps[ct_], 0, 0, 0, 0, 2) << std::endl;
}


void get_twiss(const double alpha[], const double beta[],
	       const double eta[], const double etap[])
{
  int          j, k;
  double       alpha1[2], beta1[2], eta1[2], etap1[2], dnu1[2], dnu2[2];

  // Crucial; to only store linear part of A.
  danot_(1);

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


double get_chi2(const int n, const double data[])
{
  int    j;
  double chi2 = 0e0;

  for (j = 1; j <= n; j++) chi2 += sqr(data[j]);

  return chi2;
}


void h_init(double *bn_max, double *bn)
{
  int    i;
  double L;

  std::cout << std::endl << "Initial b2:" << std::endl;
  for (i = 1; i <= n_prm; i++) {
    // Note, Jacobian is a function of the multipole strengths.
    if (bn_prm[i-1].n > 0) {
      L = get_L(bn_prm[i-1].Fnum, 1);
      if (L == 0e0) L = 1e0;
      bn_max[i] = b2L_max/L;
      bn[i] = get_bn(bn_prm[i-1].Fnum, 1, bn_prm[i-1].n);
    } else {
      bn_max[i] = L_max;
      bn[i] = get_bn_s1(bn_prm[i-1].Fnum, 1, bn_prm[i-1].n);
    }
    std::cout << std::scientific << std::setprecision(5)
	 << std::setw(13) << get_bn_s1(bn_prm[i-1].Fnum, 1, bn_prm[i-1].n);
  }
  std::cout << std::endl;
}


double f_opt(double b2s[])
{
  int          i;
  double       f, *b;
  tps          Tr[2];
  ss_vect<tps> A1, A1_A1tp;

  b = dvector(1, m_cstr_max);

  for (i = 1; i <= n_prm; i++)
    set_bn_s1(bn_prm[i-1].Fnum, bn_prm[i-1].n, b2s[i]);

  danot_(no_tps);

  A1 = get_A(alpha_, beta_, eta_, etap_); A1.propagate(1, n_elem);
  A1_A1tp = A1*tp_S(3, A1);

  danot_(no_tps-1); get_Map(); danot_(no_tps);
 
  for (i = 0; i < 2; i++)
    Tr[i] = Der(Map[2*i], 2*i+1) + Der(Map[2*i+1], 2*i+2);

  m_cstr = 0;
  b[++m_cstr] = h_ijklm(A1[x_],  0, 0, 0, 0, 1);
  b[++m_cstr] = h_ijklm(A1[px_], 0, 0, 0, 0, 1);
  // if (fabs(Tr[X_].cst()) > tr_cut)
  //   b[++m_cstr] = h_ijklm(Tr[X_], 0, 0, 0, 0, 0) - sgn(Tr[X_].cst())*tr_cut;
  // if (fabs(Tr[Y_].cst()) > tr_cut)
  //   b[++m_cstr] = h_ijklm(Tr[Y_], 0, 0, 0, 0, 0) - sgn(Tr[Y_].cst())*tr_cut;

  f = 0e0;
  for (i = 1; i <= m_cstr; i++)
    f += sqr(b[i]);

  std::cout << std::scientific << std::setprecision(3)
       << std::setw(4) << n_iter << " f = " << std::setw(9) << f << std::endl;

  free_dvector(b, 1, m_cstr_max);

  return f;
}


void h_opt(double &chi2, const double *bn_max, double *bn,
	   double &dbn_max, double *g, double *h, const bool prt_iter,
	   std::ofstream &sext_out)
{
  // Conjugate gradient method.
  std::ostringstream hs[m_cstr_max];
  int           i, k;
  double        L, chi2_old, fret, g2, gamma, dg2;
  double        **A, *b, *w, **U, **V, *dbn, *bn0;
  tps           Tr[2];
  ss_vect<tps>  A1, A1_A1tp;

  const bool   prt = true;
  const double svd_tol = 1e-12;

  b = dvector(1, m_cstr_max); w = dvector(1, n_prm); dbn = dvector(1, n_prm);
  bn0 = dvector(1, n_prm);
  A = dmatrix(1, m_cstr_max, 1, n_prm); U = dmatrix(1, m_cstr_max, 1, n_prm);
  V = dmatrix(1, n_prm, 1, n_prm);

  std::cout << std::endl;
  for (i = 1; i <= n_prm; i++) {
    if (bn_prm[i-1].n > 0)
      set_bn_par(bn_prm[i-1].Fnum, bn_prm[i-1].n, 7);
    else
      set_s_par1(bn_prm[i-1].Fnum, 7);

    std::cout << "Family: " << std::setw(6) << get_Name(bn_prm[i-1].Fnum)
	      << std::setw(4) << bn_prm[i-1].Fnum
	      << std::setw(3) << bn_prm[i-1].n << std::endl;

    danot_(no_tps);

    A1 = get_A(alpha_, beta_, eta_, etap_); A1.propagate(1, n_elem);
    A1_A1tp = A1*tp_S(3, A1);

    danot_(no_tps-1); get_Map(); danot_(no_tps);

    for (k = 0; k < 2; k++)
      Tr[k] = Der(Map[2*k], 2*k+1) + Der(Map[2*k+1], 2*k+2);

    m_cstr = 0;
    A[++m_cstr][i] = h_ijklm_p(A1[x_],  0, 0, 0, 0, 1, 7);
    A[++m_cstr][i] = h_ijklm_p(A1[px_], 0, 0, 0, 0, 1, 7);
    // if (fabs(Tr[X_].cst()) > tr_cut)
    //   A[++m_cstr][i] = h_ijklm_p(Tr[X_], 0, 0, 0, 0, 0, 7);
    // if (fabs(Tr[Y_].cst()) > tr_cut)
    //   A[++m_cstr][i] = h_ijklm_p(Tr[Y_], 0, 0, 0, 0, 0, 7);

    if (i == n_prm) {
      m_cstr = 0;
      hs[m_cstr++] << "eta_x  ";
      b[m_cstr] = -h_ijklm(A1[x_],  0, 0, 0, 0, 1);
      hs[m_cstr++] << "eta'_x ";
      b[m_cstr] = -h_ijklm(A1[px_], 0, 0, 0, 0, 1);
      // if (fabs(Tr[X_].cst()) > tr_cut) {
      // 	hs[m_cstr++] << "Tr{M_x}";
      // 	b[m_cstr] = -(h_ijklm(Tr[X_], 0, 0, 0, 0, 0)-sgn(Tr[X_].cst())*tr_cut);
      // }
      // if (fabs(Tr[Y_].cst()) > tr_cut) {
      // 	hs[m_cstr++] << "Tr{M_y}";
      // 	b[m_cstr] = -(h_ijklm(Tr[Y_], 0, 0, 0, 0, 0)-sgn(Tr[Y_].cst())*tr_cut);
      // }
    }

    if (bn_prm[i-1].n > 0)
      clr_bn_par(bn_prm[i-1].Fnum, bn_prm[i-1].n);
    else
      clr_s_par1(bn_prm[i-1].Fnum);
  }

  chi2_old = chi2; chi2 = get_chi2(m_cstr, b);

  if (prt) prt_system(m_cstr, n_prm, A, b, chi2_old, chi2, hs);

  SVD_lim(m_cstr, n_prm, A, b, bn_max, svd_tol, bn, dbn);

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
  d_linmin(bn, dbn, n_prm, &fret, f_opt);

  std::cout << std::endl << "db2L (dcorr*L):" << std::endl;
  dbn_max = 0e0;
  for (i = 1; i <= n_prm; i++) {
    if (bn_prm[i-1].n > 0) {
      L = get_L(bn_prm[i-1].Fnum, 1);
      if (L == 0e0) L = 1e0;
      dbn_max = max(fabs((bn[i]-bn0[i])*L), dbn_max);
    } else
      dbn_max = max(fabs((bn[i]-bn0[i])), dbn_max);
    if (bn_prm[i-1].n > 0) 
      std::cout << std::scientific << std::setprecision(3)
		<< std::setw(11) << (bn[i]-bn0[i])*L;
    else
      std::cout << std::scientific << std::setprecision(3)
		<< std::setw(11) << bn[i]-bn0[i];
  }
  std::cout << std::endl;

  std::cout << "b2:" << std::endl;
  for (i = 1; i <= n_prm; i++) {
    set_bn_s1(bn_prm[i-1].Fnum, bn_prm[i-1].n, bn[i]);
    bn[i] = get_bn_s1(bn_prm[i-1].Fnum, 1, bn_prm[i-1].n);
    std::cout << std::scientific << std::setprecision(6)
	 << std::setw(14) << get_bn_s1(bn_prm[i-1].Fnum, 1, bn_prm[i-1].n);
  }
  std::cout << std::endl;

  free_dvector(b, 1, m_cstr_max); free_dvector(w, 1, n_prm);
  free_dvector(dbn, 1, n_prm); free_dvector(bn0, 1, n_prm);
  free_dmatrix(A, 1, m_cstr_max, 1, n_prm);
  free_dmatrix(U, 1, m_cstr_max, 1, n_prm); free_dmatrix(V, 1, n_prm, 1, n_prm);
}


void h_zero(void)
{
  std::string   str;
  int      i;
  double   *g, *h, dbn_max, chi2, *bn, *bn_max;
  std::ofstream sext_out;

  const bool   prt_iter   = true;
  const int    n_iter_max = 1000;
  const double bn_tol     = 1e-5;

  bn = dvector(1, n_prm); bn_max = dvector(1, n_prm);
  g = dvector(1, n_prm); h = dvector(1, n_prm);

  file_wr(sext_out, "sext.dat");

  h_init(bn_max, bn);

  n_iter = 0; chi2 = 0e0;
  do {
    n_iter++;
    h_opt(chi2, bn_max, bn, dbn_max, g, h, prt_iter, sext_out);

    prt_mfile("flat_file.dat");
  } while ((dbn_max > bn_tol) && (n_iter < n_iter_max));

  if (!prt_iter) {
    for (i = 1; i <= n_prm; i++)

    for (i = 1; i <= n_prm; i++)
      sext_out << std::fixed << std::setprecision(7) 
	       << std::setw(9) << get_Name(bn_prm[i-1].Fnum)
	       << "(" << elem[bn_prm[i-1].Fnum].Knum << ") = "
	       << std::setw(11) << get_bnL_s(bn_prm[i-1].Fnum, 1, bn_prm[i-1].n)
	       << std::setw(2) << bn_prm[i-1].n << std::endl;
    sext_out.flush();
  }

  free_dvector(bn, 1, n_prm); free_dvector(bn_max, 1, n_prm);
  free_dvector(g, 1, n_prm); free_dvector(h, 1, n_prm);
}


void get_prms(void)
{
  int k;

  // Arc3.
  n_prm = 0;
  bn_prm[n_prm].Fnum = get_Fnum("qm1");  bn_prm[n_prm++].n = Quad;
  bn_prm[n_prm].Fnum = get_Fnum("qm2");  bn_prm[n_prm++].n = Quad;
  // bn_prm[n_prm].Fnum = get_Fnum("dm"); bn_prm[n_prm++].n = 0;
  // bn_prm[n_prm].Fnum = get_Fnum("bm"); bn_prm[n_prm++].n = 0;

  std::cout << std::endl;
  for (k = 0; k < n_prm; k++) {
    std::cout << std::setw(3) << k+1
	 << ", " << std::setw(6) << get_Name(bn_prm[k].Fnum)
	 << ", " << std::setw(2) << bn_prm[k].Fnum
	 << ", n = " << bn_prm[k].n << std::endl;
  }
}


int main(int argc, char *argv[])
{

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

  danot_(1);
  get_twiss(alpha_, beta_, eta_, etap_);

  prt_lat("linlat.out", 10);
  prt_lat("linlat1.out");

  get_ct();
  prt_ct();

  if (true) {
    get_prms();

    h_zero();

    prt_lin_map(3, Map);

    prt_ct();
    prt_mfile("flat_file.dat");
  }
}
