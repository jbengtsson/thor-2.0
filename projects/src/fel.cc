#define NO 5

#include "thor_lib.h"

int  no_tps   = NO,
     ndpt_tps = 5;


long int      rseed0, rseed;
double        normcut_, chi2 = 0e0;
ss_vect<tps>  Id_scl;


double           alpha0[2], beta0[2], gamma0[2];
ss_vect<double>  eta0;

// const int  n_cell = 20;
const int  n_cell = 209;


#define k 19
#define c 656329L
#define m 100000001

void iniranf(const long i)
{
  rseed0 = i; rseed = i;
}

void newseed(void)
{
  rseed0 = (k*rseed0+c) % m; rseed = (rseed0+54321) % m;
}

double ranf(void)
{
  // Random number [0, 1] with rectangular distribution.
  rseed = (k*rseed+c) % m; return (rseed/1e8);
}

#undef k
#undef c
#undef m  


void setrancut(const double cut)
{

  printf("\n");
  printf("setrancut: cut set to %3.1f\n", cut);

  normcut_ = cut;
}


double normranf(void)
{
  int     i, j;
  double  f, w;

  const int  maxiter = 100, n = 12;

  j = 0;
  do {
    j++;
    w = 0.0;
    for (i = 1; i <= n; i++)
      w += ranf();
    f = w - 6.0;
  } while (fabs(f) > fabs(normcut_) && j <= maxiter);

  if (j > maxiter)
    fprintf(stdout,"*** fatal error in normranf\n");
  return f;
}


double gaussian(void)
{
  // Guassian distribution by the Box-Muller method.
  static bool    even = false;
  static double  v, fac;
  double         s, z, u;

  if (even)
    z = v*fac;
  else {
    do {
      u = 2e0*(double)rand()/(double)RAND_MAX - 1e0;
      v = 2e0*(double)rand()/(double)RAND_MAX - 1e0;
      s = sqr(u) + sqr(v);
    } while(s >= 1e0 || s == 0e0);
    fac = sqrt(-2e0*log(s)/s); z = u*fac;
  }

  even = !even;

  return z;
}


tps dacfu1(const tps &a, double (*func)(const int []))
{
  char    name[11];
  int     j, n, jj[ss_dim], ibuf1[bufsize], ibuf2[bufsize];
  double  rbuf[bufsize];
  tps     b;

  a.exprt(rbuf, ibuf1, ibuf2, name); n = (int)rbuf[0];

  for (j = 0; j < n; j++) {
    dehash_(no_tps, ss_dim, ibuf1[j], ibuf2[j], jj);

     rbuf[j+1] *= (*func)(jj);

    hash_(no_tps, ss_dim, jj, ibuf1[j], ibuf2[j]);
  }

  b.imprt(n, rbuf, ibuf1, ibuf2);

  return b;
}


void prt_lin_map(const int n_DOF, const ss_vect<tps> &map)
{
  int  i, j;

  cout << endl;
  for (i = 0; i < 2*n_DOF; i++) {
    for (j = 0; j < 2*n_DOF; j++)
      if (true)
        cout << scientific << setprecision(5) << setw(13) << map[i][j];
      else
        cout << scientific << setprecision(16) << setw(24) << map[i][j];
    cout << endl;
  }
}


ss_vect<tps> pow(const ss_vect<tps> &a, const int n)
{
  ss_vect<tps>  b;

  if (n < 0) {
    cout << "pow: undefined exponent " << n << endl;
    exit(1);
  }

  if (n == 0)
    b.identity();
  else
    b = a*pow(a, n-1);

  return b;
}


void no_mpoles(const int n)
{
  int j, k;

  cout << endl;
  cout << "zeroing multipoles" << endl;
  cout << endl;
  for (j = 0; j < n_elem; j++)
    if (elem[j].kind == Mpole)
      for (k = n; k < mpole_max; k++) {
	set_bn(elem[j].Fnum, elem[j].Knum, k, 0.0);
      }
}


void fit_isoch(const double m_11[], const int n_b2, const int b2s[],
	       const double eps)
{
  int       i, j, n;
  double    **A, *b, *b2_lim, *b2, *db2;
  ofstream  quad_out;

  const int     m     = 5;
  const double  s_cut = 1e-9, step = 0.7, scl_ct = 1e0, R_56 = 0e0;

  b = dvector(1, m); b2_lim = dvector(1, n_b2);
  b2 = dvector(1, n_b2); db2 = dvector(1, n_b2);
  A = dmatrix(1, m, 1, n_b2);

  danot_(2); 

  quad_out.open("fit_isoch.dat");

  for (i = 1; i <= n_b2; i++) {
    b2_lim[i] = b2_max; b2[i] = get_bn(b2s[i-1], 1, Quad);
  }

  n = 0; b[1] = 1e30;
  while ((fabs(b[1]) > eps) || (fabs(b[2]) > eps) || (fabs(b[3]) > eps)) {
    n++;
    for (i = 1; i <= n_b2; i++) {
      for (j = 1; j <= get_n_Kids(b2s[i-1]); j++)
	set_bn_par(b2s[i-1], j, Quad, 7);

      get_Map(); Map = pow(Map, n_cell);

      A[1][i] = h_ijklm_p(Map[x_], 1, 0, 0, 0, 0, 7);
      A[2][i] = h_ijklm_p(Map[x_], 0, 1, 0, 0, 0, 7);
      A[3][i] = h_ijklm_p(Map[y_], 0, 0, 1, 0, 0, 7);
      A[4][i] = h_ijklm_p(Map[y_], 0, 0, 0, 1, 0, 7);
      A[5][i] = scl_ct*h_ijklm_p(Map[ct_], 0, 0, 0, 0, 1, 7);

      for (j = 1; j <= get_n_Kids(b2s[i-1]); j++)
	clr_bn_par(b2s[i-1], j, Quad);
    }

    b[1] = -(h_ijklm(Map[x_], 1, 0, 0, 0, 0)-m_11[X_]);
    b[2] = -h_ijklm(Map[x_], 0, 1, 0, 0, 0);
    b[3] = -(h_ijklm(Map[y_], 0, 0, 1, 0, 0)-m_11[Y_]);
    b[4] = -h_ijklm(Map[y_], 0, 0, 0, 1, 0);
    b[5] = -scl_ct*(h_ijklm(Map[ct_], 0, 0, 0, 0, 1)-R_56);

    cout << endl;
    cout << " Ax = b:" << endl;
    cout << endl;
    for (i = 1; i <= m; i++) {
      for (j = 1; j <= n_b2; j++)
	cout << scientific << setprecision(3) << setw(11) << A[i][j];
      cout << scientific << setprecision(3) << setw(11) << b[i] << endl;
    }
	
    SVD_lim(m, n_b2, A, b, b2_lim, s_cut, b2, db2);

    for (i = 1; i <= n_b2; i++) {
      set_dbn(b2s[i-1], Quad, step*db2[i]);
      b2[i] = get_bn(b2s[i-1], 1, Quad);
    }

    get_Map(); Map = pow(Map, n_cell);

    cout << endl;
    cout << scientific << setprecision(3)
	 << "m_11 - 1:" << setw(11) << h_ijklm(Map[x_], 1, 0, 0, 0, 0)-1e0
	 << "," << setw(11) << h_ijklm(Map[y_], 0, 0, 1, 0, 0)-1e0 << endl;
    cout << scientific << setprecision(3)
	 << "m_12:    " << setw(11) << h_ijklm(Map[x_], 0, 1, 0, 0, 0)
	 << "," << setw(11) << h_ijklm(Map[y_], 0, 0, 0, 1, 0) << endl;
    cout << scientific << setprecision(3)
	 << "R_56 =   " << setw(11) << h_ijklm(Map[ct_], 0, 0, 0, 0, 1)
	 << endl;

    quad_out << endl;
    quad_out << "n = " << n << endl;
    for (i = 1; i <= n_b2; i++)
      for (j = 1; j <= get_n_Kids(b2s[i-1]); j++)
	quad_out << fixed << setprecision(16) 
		 << setw(6) << get_Name(b2s[i-1]) << "(" << j << ") = "
		 << setw(24) << get_bnL(b2s[i-1], j, Quad)
		 << setw(2) << Quad << endl;
  }

  quad_out.close();

  free_dvector(b, 1, m);
  free_dvector(b2_lim, 1, n_b2); free_dvector(b2, 1, n_b2);
  free_dvector(db2, 1, n_b2);
  free_dmatrix(A, 1, m, 1, n_b2);
}


void fit_isoch1(const double m_11[], const int n_b2, const int b2s[],
		const double eps)
{
  int       i, j, n;
  float     **A, *b, *w, **U, **V, *db2;
  double    s_max;
  ofstream  quad_out;

  const int     m     = 5;
  const double  s_cut = 1e-9, step = 0.7, scl_ct = 1e0, R_56 = 0e0;

  b = vector(1, m); w = vector(1, n_b2); db2 = vector(1, n_b2);
  A = matrix(1, m, 1, n_b2); U = matrix(1, m, 1, n_b2);
  V = matrix(1, n_b2, 1, n_b2);

  danot_(2); 

  quad_out.open("fit_isoch.dat");

  n = 0; b[1] = 1e30;
  while ((fabs(b[1]) > eps) || (fabs(b[2]) > eps) || (fabs(b[3]) > eps)) {
    n++;
    for (i = 1; i <= n_b2; i++) {
      for (j = 1; j <= get_n_Kids(b2s[i-1]); j++)
	set_bn_par(b2s[i-1], j, Quad, 7);

      get_Map(); Map = pow(Map, n_cell);

      A[1][i] = h_ijklm_p(Map[x_], 1, 0, 0, 0, 0, 7);
      A[2][i] = h_ijklm_p(Map[x_], 0, 1, 0, 0, 0, 7);
      A[3][i] = h_ijklm_p(Map[y_], 0, 0, 1, 0, 0, 7);
      A[4][i] = h_ijklm_p(Map[y_], 0, 0, 0, 1, 0, 7);
      A[5][i] = scl_ct*h_ijklm_p(Map[ct_], 0, 0, 0, 0, 1, 7);

      for (j = 1; j <= get_n_Kids(b2s[i-1]); j++)
	clr_bn_par(b2s[i-1], j, Quad);
    }

    b[1] = -(h_ijklm(Map[x_], 1, 0, 0, 0, 0)-m_11[X_]);
    b[2] = -h_ijklm(Map[x_], 0, 1, 0, 0, 0);
    b[3] = -(h_ijklm(Map[y_], 0, 0, 1, 0, 0)-m_11[Y_]);
    b[4] = -h_ijklm(Map[y_], 0, 0, 0, 1, 0);
    b[5] = -scl_ct*(h_ijklm(Map[ct_], 0, 0, 0, 0, 1)-R_56);

    cout << endl;
    cout << " Ax = b:" << endl;
    cout << endl;
    for (i = 1; i <= m; i++) {
      for (j = 1; j <= n_b2; j++)
	cout << scientific << setprecision(3) << setw(11) << A[i][j];
      cout << scientific << setprecision(3) << setw(11) << b[i] << endl;
    }
	
    for (i = 1; i <= m; i++)
      for (j = 1; j <= n_b2; j++)
	U[i][j] = A[i][j];

    svdcmp(U, m, n_b2, w, V);

    s_max = -1e30;
    for (i = 1; i <= n_b2; i++)
      s_max = max(w[i], s_max);
  
    cout << endl;
    cout << "singular values:" << endl;
    for (i = 1; i <= n_b2; i++) {
      cout << scientific << setprecision(3) << setw(10) << w[i];
      if (w[i]/s_max < s_cut) {
	w[i] = 0.0;
	cout << " (zeroed)";
      }
      cout << endl;
    }

    svbksb(U, w, V, m, n_b2, b, db2);

    for (i = 1; i <= n_b2; i++)
      set_dbn(b2s[i-1], Quad, step*db2[i]);

    get_Map(); Map = pow(Map, n_cell);

    cout << endl;
    cout << scientific << setprecision(3)
	 << "m_11 - 1:" << setw(11) << h_ijklm(Map[x_], 1, 0, 0, 0, 0)-1e0
	 << "," << setw(11) << h_ijklm(Map[y_], 0, 0, 1, 0, 0)-1e0 << endl;
    cout << scientific << setprecision(3)
	 << "m_12:    " << setw(11) << h_ijklm(Map[x_], 0, 1, 0, 0, 0)
	 << "," << setw(11) << h_ijklm(Map[y_], 0, 0, 0, 1, 0) << endl;
    cout << scientific << setprecision(3)
	 << "R_56 =   " << setw(11) << h_ijklm(Map[ct_], 0, 0, 0, 0, 1)
	 << endl;

    quad_out << endl;
    quad_out << "n = " << n << endl;
    for (i = 1; i <= n_b2; i++)
      for (j = 1; j <= get_n_Kids(b2s[i-1]); j++)
	quad_out << fixed << setprecision(16) 
		 << setw(6) << get_Name(b2s[i-1]) << "(" << j << ") = "
		 << setw(24) << get_bnL(b2s[i-1], j, Quad)
		 << setw(2) << Quad << endl;
  }

  quad_out.close();

  free_vector(b, 1, m); free_vector(w, 1, n_b2);
  free_vector(db2, 1, n_b2);
  free_matrix(A, 1, m, 1, n_b2); free_matrix(U, 1, m, 1, n_b2);
  free_matrix(V, 1, n_b2, 1, n_b2);
}


double f_prm(const int jj[])
{
  int     k, n;
  double  f;

  n = 0;
  for (k = 0; k < 2*nd_tps; k++)
    n += jj[k];

  f = ((n == no_tps) && (jj[ss_dim-1] == 0))? 0e0 : 1e0;

 return f;
}


void fit_alpha(const int n_bn, const int bn_Fam[], const int bn_type[],
	       const double eps, const int n_max)
{
  const int     m_max = 200;  // max no of constraints

  int            i, j, k, l, m, m1 = 0, i1, n;
  float          **A, *b, *w, **U, **V, *dbn;
  double         s_max;
  ss_vect<tps>   nus;
  ostringstream  hs[m_max];
  ofstream       bn_out;

  const double  s_cut = 1e-15, step = 0.3;
  const double  R_566_scl   = 1e0, R_5666_scl  = 1e0;
  const double  T_10002_scl = 1e0, T_20001_scl = 1e0;

  b = vector(1, m_max); w = vector(1, n_bn); dbn = vector(1, n_bn);
  A = matrix(1, m_max, 1, n_bn); U = matrix(1, m_max, 1, n_bn);
  V = matrix(1, n_bn, 1, n_bn);

  bn_out.open("fit_alpha.dat", ios::out);

  danot_(no_tps);

  n = 0;
  do {
    n++;
    for (i1 = 1; i1 <= n_bn; i1++) {
      for (j = 1; j <= get_n_Kids(bn_Fam[i1-1]); j++)
	set_bn_par(bn_Fam[i1-1], j, bn_type[i1-1], 7);

      get_Map(); Map = pow(Map, n_cell);
      K = MapNorm(Map, g, A1, A0, Map_res, 1); nus = dHdJ(K);

      Map[ct_] = dacfu1(Map[ct_], f_prm); Map[ct_] = Map[ct_]*Id_scl;

      m1 = 0;

      for (i = 0; i <= no_tps-1; i++)
	for (j = 0; j <= no_tps-1; j++)
	  for (k = 0; k <= no_tps-1; k++)
	    for (l = 0; l <= no_tps-1; l++)
	      for (m = 0; m <= no_tps-1; m++) {
		if ((0 < i+j+k+l+m) && (i+j+k+l+m <= no_tps-1)) {
		  if (m1 == m_max-1) {
		    cout << "m_max reached " << m1 << "(" << m_max << ")"
			 << endl;
		    exit(1);
		  }

		  if (!((i == 1) && (j == 0) && (k == 0) && (l == 0)
			&& (m == 0)) &&
		      !((i == 0) && (j == 1) && (k == 0) && (l == 0)
			&& (m == 0)) &&
		      !((i == 0) && (j == 0) && (k == 0) && (l == 0)
			&& (m == 1)) &&
		      (fabs(h_ijklm(Map[ct_], i, j, k, l, m)) > 0e0)) {
		    A[++m1][i1] = h_ijklm_p(Map[ct_], i, j, k, l, m, 7);
		    if (i1 == n_bn) {
		      hs[m1-1].str("");
		      hs[m1-1] << "ct_" << i << j << k << l << m;
		      b[m1] = -h_ijklm(Map[ct_], i, j, k, l, m);
		    }

		    if ((i == 0) && (j == 0) && (k == 0) && (l == 0)
			&& (m == 2)) {
		      A[m1][i1] *= R_566_scl; b[m1] *= R_566_scl;
		    } else if ((i == 0) && (j == 0) && (k == 0) && (l == 0)
			       && (m == 3)) {
		      A[m1][i1] *= R_5666_scl; b[m1] *= R_5666_scl;
		    } else if ((i == 1) && (j == 0) && (k == 0) && (l == 0)
			       && (m == 2)) {
		      A[m1][i1] *= T_10002_scl; b[m1] *= T_10002_scl;
		    } else if ((i == 2) && (j == 0) && (k == 0) && (l == 0)
			       && (m == 1)) {
		      A[m1][i1] *= T_20001_scl; b[m1] *= T_20001_scl;
		    }
		  }
		}
	      }

      for (j = 1; j <= get_n_Kids(bn_Fam[i1-1]); j++)
	clr_bn_par(bn_Fam[i1-1], j, bn_type[i1-1]);
    }

    cout << endl;
    cout << n << " Ax = b:" << endl;
    cout << endl;
    for (i = 1; i <= m1; i++) {
      cout  << setw(3) << i << " " << hs[i-1].str();
      for (j = 1; j <= n_bn; j++)
	cout << scientific << setprecision(2) << setw(10) << A[i][j];
      cout << scientific << setprecision(2) << setw(10) << b[i] << endl;
    }
    
    for (i = 1; i <= m1; i++)
      for (j = 1; j <= n_bn; j++)
	U[i][j] = A[i][j];

    svdcmp(U, m1, n_bn, w, V);

    s_max = -1e30;
    for (i = 1; i <= n_bn; i++)
      s_max = max(w[i], s_max);
  
    cout << endl;
    cout << "singular values:" << endl;
    for (i = 1; i <= n_bn; i++) {
      cout << scientific << setprecision(3) << setw(10) << w[i];
      if (w[i]/s_max < s_cut) {
	w[i] = 0.0;
	cout << " (zeroed)";
      }
      cout << endl;
    }

    svbksb(U, w, V, m1, n_bn, b, dbn);

    for (i = 1; i <= n_bn; i++)
      for (j = 1; j <= get_n_Kids(bn_Fam[i-1]); j++)
	set_dbn(bn_Fam[i-1], j, bn_type[i-1], step*dbn[i]);

    get_Map(); Map = pow(Map, n_cell);
    K = MapNorm(Map, g, A1, A0, Map_res, 1); nus = dHdJ(K);

    Map[ct_] = Map[ct_]*Id_scl;

    cout << endl;
    cout << scientific << setprecision(3)
	 << "ksi:      "
	 << setw(11) << h_ijklm(nus[3], 0, 0, 0, 0, 1) << ","
	 << setw(11) << h_ijklm(nus[4], 0, 0, 0, 0, 1) << endl;
    cout << scientific << setprecision(5)
	 << "R_566 =   "
	 << setw(13) << h_ijklm(Map[ct_], 0, 0, 0, 0, 2) << endl;
    cout << scientific << setprecision(5)
	 << "R_5666 =  "
	 << setw(13) << h_ijklm(Map[ct_], 0, 0, 0, 0, 3) << endl;
    cout << scientific << setprecision(5)
	 << "R_56666 = "
	 << setw(13) << h_ijklm(Map[ct_], 0, 0, 0, 0, 4) << endl;

    bn_out << endl;
    bn_out << "n = " << n << endl;
    for (i = 1; i <= n_bn; i++)
      for (j = 1; j <= get_n_Kids(bn_Fam[i-1]); j++)
	bn_out << fixed << setprecision(16) 
	       << setw(6) << get_Name(bn_Fam[i-1]) << "("
	       << j << ") = " << setw(24)
	       << get_bnL(bn_Fam[i-1], j, bn_type[i-1])
	       << setw(2) << bn_type[i-1] << endl;

    cout << endl;
    cout << scientific << setprecision(1)
	 << setw(2) << n << " chi2: " << chi2;

    chi2 = 0e0;
    for (i = 1; i <= m1; i++)
      chi2 += sqr(b[i]);

    cout << scientific << setprecision(1)
	 << " -> " << chi2 << endl;

  } while ((sqrt(chi2) > eps) && (n < n_max));

  bn_out.close();

  free_vector(b, 1, m_max); free_vector(w, 1, n_bn);
  free_vector(dbn, 1, n_bn);
  free_matrix(A, 1, m_max, 1, n_bn); free_matrix(U, 1, m_max, 1, n_bn);
  free_matrix(V, 1, n_bn, 1, n_bn);
}


void prt_ps_long(const int n, const double eps[], const double delta)
{
  int              i, k;
  double           twoJ[2], phi[2], ct_sum, ct_sum2, ct_mean, ct_sigma;
  double           delta_sum, delta_sum2, delta_mean, delta_sigma;
  ss_vect<double>  ps, ps_Fl;
  ofstream         outf;

  outf.open("ps_long.out");

  delta_sum = 0e0; delta_sum2 = 0e0;
  ct_sum = 0e0; ct_sum2 = 0e0;

  ps_Fl.zero();
  for (i = 1; i <= n; i++) {
    for (k = 0; k < 2; k++) {
      phi[k] = 2e0*M_PI*ranf();
      twoJ[k] = eps[k]*normranf();
      while (twoJ[k] < 0e0)
	twoJ[k] = eps[k]*normranf();

      ps_Fl[2*k] = sqrt(twoJ[k])*cos(phi[k]);
      ps_Fl[2*k+1] = -sqrt(twoJ[k])*sin(phi[k]);
    }

    ps_Fl[delta_] = delta*normranf();

    delta_sum += ps_Fl[delta_]; delta_sum2 += sqr(ps_Fl[delta_]);

    ps = (A1*ps_Fl).cst();

    if (false)
      for (k = 1; k <= n_cell; k++)
	ps.propagate(1, n_elem);
    else
      ps = (Map*ps).cst();

    ct_sum += ps[ct_]; ct_sum2 += sqr(ps[ct_]);

    outf << scientific << setprecision(5)
	 << setw(4) << i << setw(13) << ps << endl;
  }

  delta_mean = delta_sum/n;
  delta_sigma = sqrt((n*delta_sum2-sqr(delta_sum))/(n*(n-1e0)));  
  ct_mean = ct_sum/n; ct_sigma = sqrt((n*ct_sum2-sqr(ct_sum))/(n*(n-1e0)));  

  cout << endl;
  cout << scientific << setprecision(3)
       << "delta = " << setw(11) << delta_mean
       << " +/-" << setw(10) << delta_sigma << endl;
  cout << scientific << setprecision(3)
       << "ct    = " << setw(11) << ct_mean
       << " +/-" << setw(10) << ct_sigma << endl;

  outf.close();
}


void get_orbit(void)
{
  int              j;
  ss_vect<double>  ps;
  ofstream         outf;

  outf.open("cod.out");

  no_mpoles(Quad);

  ps.zero();

  for (j = 1; j < n_elem; j++) {
    ps.propagate(j, j); 
    if ((elem[j-1].kind == Mpole) && (elem[j-1].mpole->n_design >= Quad)) {
      outf << scientific << setprecision(3)
	   << setw(3) << j
	   << " " << left << setw(6) << elem[j-1].Name << right
	   << setw(3) << elem[j-1].Fnum << setw(3) << elem[j-1].Knum
	   << setw(11) << ps[x_] << setw(11) << ps[y_] << endl;
    }
  }

  outf.close();
}


void mpole_align(void)
{
  string    line;
  int       j;
  double    dx[2];
  ifstream  inf;

  inf.open("cod.out");

  while (getline(inf, line) != NULL) {
    sscanf(line.c_str(), "%d %*6s %*d %*d %lf %lf", &j, &dx[X_], &dx[Y_]);
    elem[j-1].dx[X_] = -dx[X_]; elem[j-1].dx[Y_] = -dx[Y_];
    
    cout << scientific << setprecision(3)
	 << setw(3) << j
	 << setw(11) << dx[X_] << setw(11) << dx[Y_]
	 << endl;
  }
}


void get_moments(const int n, const double eps[], const double delta,
		 double sigma[][ps_dim])
{
  // Note, eps = <J>.
  int              i, j, k;
  double           phi[2], twoJ[2], delta_rnf = 0e0, eps1[2], delta1;
  ss_vect<double>  ps, ps_Fl, ps_mean;

  const int  seed = 1121;

  // Compute averages.

  srand(seed);

  for (k = 0; k < 2; k++)
    eps1[k] = 0e0;
  delta1 = 0e0;

  ps_mean.zero();
  for (i = 1; i <= n; i++) {
    ps_Fl.zero();
    for (k = 0; k < 2; k++) {
      phi[k] = 2e0*M_PI*ranf();
      twoJ[k] = 2e0*eps[k]*gaussian();
      while (twoJ[k] < 0e0)
	twoJ[k] = 2e0*eps[k]*gaussian();

      ps_Fl[2*k] = sqrt(twoJ[k])*cos(phi[k]);
      ps_Fl[2*k+1] = -sqrt(twoJ[k])*sin(phi[k]);

      eps1[k] += sqr(twoJ[k]/2e0);
    }

    delta_rnf = delta*gaussian();

    ps_Fl[delta_] = delta_rnf; ps = (A1*ps_Fl).cst();

    for (j = 0; j < ps_dim; j++)
      ps_mean[j] += ps[j];

    delta1 += sqr(delta_rnf);
  }

  for (j = 0; j < ps_dim; j++)
    ps_mean[j] /= n;

  delta1 /= n;
  for (k = 0; k < 2; k++)
    eps1[k] /= n;

  cout << endl;
  cout << scientific << setprecision(3)
       << "eps rms   =  [" << setw(9) << sqrt(eps1[X_]) << ", "
       << setw(9) << sqrt(eps1[Y_]) << "]" << endl;
  cout << scientific << setprecision(3)
       << "delta rms = " << setw(11) << sqrt(delta1) << endl;

  cout << endl;
  cout << scientific << setprecision(3)
       << "x_max:"
       << setw(11) << sqrt(beta0[X_]*2e0*eps[X_])
       << setw(11) << sqrt(gamma0[X_]*2e0*eps[X_])
       << setw(11) << sqrt(beta0[Y_]*2e0*eps[Y_])
       << setw(11) << sqrt(gamma0[Y_]*2e0*eps[Y_])
       << setw(11) << delta << setw(11) << 0e0 << endl;
  cout << scientific << setprecision(3)
       << "<x>:  " << setw(11) << ps_mean << endl;

  // Compute covariances.

  srand(seed);

  for (j = 0; j < ps_dim; j++)
    for (k = 0; k < ps_dim; k++)
      sigma[j][k] = 0e0;

  for (i = 1; i <= n; i++) {
    ps_Fl.zero();
    for (k = 0; k < 2; k++) {
      phi[k] = 2e0*M_PI*ranf();
      twoJ[k] = 2e0*eps[k]*gaussian();
      while (twoJ[k] < 0e0)
	twoJ[k] = 2e0*eps[k]*gaussian();

      ps_Fl[2*k] = sqrt(twoJ[k])*cos(phi[k]);
      ps_Fl[2*k+1] = -sqrt(twoJ[k])*sin(phi[k]);
    }

    delta_rnf = delta*gaussian();

    ps_Fl[delta_] = delta_rnf; ps = (A1*ps_Fl).cst();

    for (j = 0; j < ps_dim; j++) {
      for (k = 0; k < ps_dim; k++)
	sigma[j][k] += (ps[j]-ps_mean[j])*(ps[k]-ps_mean[k]);
    }
  }

  cout << endl;
  for (j = 0; j < ps_dim; j++) {
    for (k = 0; k < ps_dim; k++) {
      sigma[j][k] /= n - 1;

      cout << scientific << setprecision(3)
	   << setw(11) << sigma[j][k];
    }
    cout << endl;
  }

  cout << endl;
  cout << scientific << setprecision(3)
       << "sigma_x,x =      " << setw(11) << sigma[x_][x_]
       << setw(11) << beta0[X_]*eps[X_] << endl;
  cout << scientific << setprecision(3)
       << "sigma_x,p_x =    " << setw(11) << sigma[x_][px_]
       << setw(11) << -alpha0[X_]*eps[X_]
       << endl;
  cout << scientific << setprecision(3)
       << "sigma_x,delta =  " << setw(11) << sigma[x_][delta_]
       << setw(11) << eta0[x_]*sqr(delta) << endl;

  cout << endl;
  cout << scientific << setprecision(3)
       << "sigma_p_x,p_x =  " << setw(11) << sigma[px_][px_]
       << setw(11) << gamma0[X_]*eps[X_] << endl;
  cout << scientific << setprecision(3)
       << "sigma_p_x,delta =" << setw(11) << sigma[px_][delta_]
       << setw(11) << eta0[px_]*sqr(delta) << endl;

  cout << endl;
  cout << scientific << setprecision(3)
       << "sigma_y,y =      " << setw(11) << sigma[y_][y_]
       << setw(11) << beta0[Y_]*eps[Y_] << endl;
  cout << scientific << setprecision(3)
       << "sigma_y,p_y =    " << setw(11) << sigma[y_][py_]
       << setw(11) << -alpha0[Y_]*eps[Y_] << endl;
  cout << scientific << setprecision(3)
       << "sigma_y,delta =  " << setw(11) << sigma[y_][delta_]
       << setw(11) << eta0[y_]*sqr(delta) << endl;

  cout << endl;
  cout << scientific << setprecision(3)
       << "sigma_p_y,p_y =  " << setw(11) << sigma[py_][py_]
       << setw(11) << gamma0[Y_]*eps[Y_] << endl;
  cout << scientific << setprecision(3)
       << "sigma_p_y,delta =" << setw(11) << sigma[py_][delta_]
       << setw(11) << eta0[py_]*sqr(delta) << endl;
}


int main(int argc, char *argv[])
{
  const int  b2s_max = 10, bns_max = 20;

  int              k;
  int              n_b2, b2s[b2s_max], n_bn, bn_Fam[bns_max], bn_n[bns_max];
  double           s[2][3];
  double           sigma0[ps_dim][ps_dim];
  tps              qf, ct2;
  ss_vect<double>  ps;
  ss_vect<tps>     Id, nus, sigma;
  ofstream         outf;

  const long int  seed   = 1121;

  const double    nu[]   = {0.95, 0.45};
  const double    m_11[] = {1.0,  1.0};
//   const double    m_11[] = {cos(nu[X_]*2e0*M_PI), cos(nu[Y_]*2e0*M_PI)};

  const double    eps[]  = {5*5e-12, 5*5e-12};
  const double    delta  = 1e-4;

  rad_on    = false; H_exact        = false; totpath_on   = false;
  cavity_on = false; quad_fringe_on = false; emittance_on = false;
  IBS_on    = false;

  iniranf(seed); setrancut(3e0);

  rd_mfile(argv[1], elem); rd_mfile(argv[1], elem_tps);

  // initialize the symplectic integrator after energy has been defined
  ini_si();

  // disable from TPSALib- and LieLib log messages
  idprset(-1);

//   rad_on = true; cavity_on = true;

  if (false) {
    get_orbit();
    exit(0);
  }

//   mpole_align();

  get_Map(); Map = pow(Map, n_cell);

  cout << endl;
  cout << "M:" << endl;
  prt_lin_map(3, Map);

  if (false) {
    outf.open("map.dat");
    outf << Map;
    outf.close();
    prt_lin_map(3, Map);

    ps.zero();
    ps[x_]     =  0.1e0;   ps[px_] =  0.15e0;
    ps[y_]     = -0.05e0;  ps[py_] = -0.02e0;
    ps[delta_] =  0.05e0;  ps[ct_] =  0.1e0;

    cout << scientific << setprecision(16)
	 << setw(24) << ps << endl;
    cout << scientific << setprecision(16)
	 << setw(24) << (Map*ps).cst() << endl;

    exit(0);
  }

  K = MapNorm(Map, g, A1, A0, Map_res, 1); nus = dHdJ(K);

  if (false) {
    cout << endl;
    cout << "A0:" << endl;
    prt_lin_map(3, A0);

    cout << endl;
    cout << "A1:" << endl;
    prt_lin_map(3, A1);
  }

  get_ab(alpha0, beta0, 0);

  // Single pass system.
  for (k = 0; k < 4; k++)
    eta0[k] = h_ijklm(Map[k], 0, 0, 0, 0, 1);

  cout << endl;
  cout << scientific << setprecision(3)
       << "eps: [" << setw(9) << eps[X_] << ", "
       << setw(11) << eps[Y_] << "]" << endl;

  cout << endl;
  cout << scientific << setprecision(3)
       << "alpha: [" << setw(11) << alpha0[X_] << ", "
       << setw(11) << alpha0[Y_] << "]" << endl;
  cout << scientific << setprecision(3)
       << "beta:  [" << setw(11) << beta0[X_] << ", "
       << setw(11) << beta0[Y_] << "]" << endl;
  cout << scientific << setprecision(3)
       << "eta:   [" << setw(11) << eta0[x_] << ", "
       << setw(11) << eta0[px_] << ", " << eta0[y_] << ", "
       << setw(11) << eta0[py_] << "]" << endl;

  for (k = 0; k < 2; k++) {
    gamma0[k] = (1e0+sqr(alpha0[k]))/beta0[k];

    s[k][0] = beta0[k]*eps[k]; s[k][1] = -alpha0[k]*eps[k];
    s[k][2] = gamma0[k]*eps[k];
  }

  if (false) {
    get_moments(10000, eps, delta, sigma0);
    exit(0);
  }

  cout << endl;
  cout << scientific << setprecision(16)
       << "sigma_trans: " << setw(24) << s[X_][0] << setw(24) << s[X_][1]
       << setw(24) << s[X_][2] << endl;
  cout << scientific << setprecision(16)
       << "sigma_trans: " << setw(24) << s[Y_][0] << setw(24) << s[Y_][1]
       << setw(24) << s[Y_][2] << endl;
  cout << endl;
  cout << scientific << setprecision(3)
       << "sigma_x = " << sqrt(s[X_][0])
       << ", sigma_px = " << sqrt(s[X_][2]) << endl;
  cout << scientific << setprecision(3)
       << "sigma_y = " << sqrt(s[Y_][0])
       << ", sigma_py = " << sqrt(s[Y_][2]) << endl;

  // 3 sigma
  Id_scl.identity();
  Id_scl[x_] *= 3e0*sqrt(s[X_][0]); Id_scl[px_] *= 3e0*sqrt(s[X_][2]);
  Id_scl[y_] *= 3e0*sqrt(s[Y_][0]); Id_scl[py_] *= 3e0*sqrt(s[Y_][2]);
  Id_scl[delta_] *= 3e0*delta;

  if (false) Map[ct_] = Map[ct_]*Id_scl;

  cout << endl;
  cout << fixed << setprecision(10)
       << "nu:      " << setw(13) << nus[0].cst() << "," 
       << setw(13) << nus[1].cst() << endl;
  cout << scientific << setprecision(3)
       << "ksi:     " << setw(11) << h_ijklm(nus[0], 0, 0, 0, 0, 1) << ","
       << setw(11) << h_ijklm(nus[1], 0, 0, 0, 0, 1) << endl;
  cout << scientific << setprecision(3)
       << "R_56:    " << setw(11) << h_ijklm(Map[ct_], 0, 0, 0, 0, 1) << endl;
  cout << scientific << setprecision(3)
       << "R_566:   " << setw(11) << h_ijklm(Map[ct_], 0, 0, 0, 0, 2) << endl;
  cout << scientific << setprecision(3)
       << "R_5666:  " << setw(11) << h_ijklm(Map[ct_], 0, 0, 0, 0, 3) << endl;
  cout << scientific << setprecision(3)
       << "R_56666: " << setw(11) << h_ijklm(Map[ct_], 0, 0, 0, 0, 4) << endl;

  if (true) {
    n_b2 = 0;
    b2s[n_b2++] = get_Fnum("qf2"); b2s[n_b2++] = get_Fnum("qd2");
    b2s[n_b2++] = get_Fnum("qf3"); b2s[n_b2++] = get_Fnum("qd3");

    b2_max = 30e0;

    switch (0) {
    case 1:
      fit_tune(nu[X_], nu[Y_], n_b2, b2s, 1e-10, true);
      break;
    case 2:
      fit_isoch(m_11, n_b2, b2s, 1e-12);
      break;
    case 3:
      fit_isoch1(m_11, n_b2, b2s, 1e-12);
      break;
    }

//     no_mpoles(Sext);

    n_bn = 0;
    bn_Fam[n_bn] = get_Fnum("s1"); bn_n[n_bn++] = Sext;
    bn_Fam[n_bn] = get_Fnum("s2"); bn_n[n_bn++] = Sext;
    bn_Fam[n_bn] = get_Fnum("s3"); bn_n[n_bn++] = Sext;
    bn_Fam[n_bn] = get_Fnum("s4"); bn_n[n_bn++] = Sext;
    bn_Fam[n_bn] = get_Fnum("s5"); bn_n[n_bn++] = Sext;
    bn_Fam[n_bn] = get_Fnum("s1"); bn_n[n_bn++] = Oct;
    bn_Fam[n_bn] = get_Fnum("s2"); bn_n[n_bn++] = Oct;
    bn_Fam[n_bn] = get_Fnum("s3"); bn_n[n_bn++] = Oct;
    bn_Fam[n_bn] = get_Fnum("s4"); bn_n[n_bn++] = Oct;
    bn_Fam[n_bn] = get_Fnum("s5"); bn_n[n_bn++] = Oct;
    bn_Fam[n_bn] = get_Fnum("s1"); bn_n[n_bn++] = Dec;
    bn_Fam[n_bn] = get_Fnum("s2"); bn_n[n_bn++] = Dec;
    bn_Fam[n_bn] = get_Fnum("s3"); bn_n[n_bn++] = Dec;
    bn_Fam[n_bn] = get_Fnum("s4"); bn_n[n_bn++] = Dec;
    bn_Fam[n_bn] = get_Fnum("s5"); bn_n[n_bn++] = Dec;

    fit_alpha(n_bn, bn_Fam, bn_n, 1e-12, 10);
  }

  get_Map(); Map = pow(Map, n_cell);

  K = MapNorm(Map, g, A1, A0, Map_res, 1); nus = dHdJ(K);

  prt_ps_long(1000, eps, delta);

//   header = true;

  outf.open("map.dat");
  outf << Map;
  outf.close();

  Map[ct_] = Map[ct_]*Id_scl;

  outf.open("ct_Taylor.dat");
  outf << Map[ct_];
  outf.close();

  prt_lin_map(3, Map);

  cout << endl;
  cout << fixed << setprecision(10)
       << "nu:       " << setw(13) << nus[0].cst() << "," 
       << setw(13) << nus[1].cst() << endl;
  cout << scientific << setprecision(3)
       << "m_11 - 1:" << setw(11) << h_ijklm(Map[x_], 1, 0, 0, 0, 0)-1e0
       << "," << setw(11) << h_ijklm(Map[y_], 0, 0, 1, 0, 0)-1e0 << endl;
  cout << scientific << setprecision(3)
       << "m_12:    " << setw(11) << h_ijklm(Map[x_], 0, 1, 0, 0, 0)
       << "," << setw(11) << h_ijklm(Map[y_], 0, 0, 0, 1, 0) << endl;
  cout << scientific << setprecision(3)
       << "ksi:     " << setw(11) << h_ijklm(nus[0], 0, 0, 0, 0, 1) << ","
       << setw(11) << h_ijklm(nus[1], 0, 0, 0, 0, 1) << endl;
  cout << scientific << setprecision(3)
       << "R_56:    " << setw(11) << h_ijklm(Map[ct_], 0, 0, 0, 0, 1) << endl;
  cout << scientific << setprecision(3)
       << "R_566:   " << setw(11) << h_ijklm(Map[ct_], 0, 0, 0, 0, 2) << endl;
  cout << scientific << setprecision(3)
       << "R_5666:  " << setw(11) << h_ijklm(Map[ct_], 0, 0, 0, 0, 3) << endl;
  cout << scientific << setprecision(3)
       << "R_56666: " << setw(11) << h_ijklm(Map[ct_], 0, 0, 0, 0, 4) << endl;

  get_Map(); Map = pow(Map, n_cell);

  Id.identity();

  sigma.zero();
  for (k = 0; k < 2; k++) {
    sigma[2*k]   = s[k][0]*Id[2*k] + s[k][1]*Id[2*k+1];
    sigma[2*k+1] = s[k][1]*Id[2*k] + s[k][2]*Id[2*k+1];
  }
  sigma[delta_] = sqr(delta)*Id[delta_];

  prt_lin_map(3, sigma);

  sigma = Map*sigma*tp_S(ps_dim/2, Map);

  prt_lin_map(3, sigma);

  qf =
    s[0][0]*sqr(Id[x_]) + 2e0*s[0][1]*Id[x_]*Id[px_] + s[0][2]*sqr(Id[px_])
    + s[1][0]*sqr(Id[y_]) + 2e0*s[1][1]*Id[y_]*Id[py_] + s[1][2]*sqr(Id[py_])
    + sqr(delta)*sqr(Id[delta_]);

  cout << qf << endl;

  ct2 = sqr(Map[ct_]);

  cout << ct2;
}
