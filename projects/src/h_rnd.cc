#define NO 4

#include "thor_lib.h"

int  no_tps   = NO,
     ndpt_tps = 5;


extern double        b2_max;
extern tps           K, g;
extern ss_vect<tps>  Map, A0, A1, Map_res;

const int   max_n_femto = 7;

const int   n_track     = 512;
const int   n_aper      = 25;

const bool  tune_scn = false;

const int   max_ind = 10;

bool          h_sel[max_ind][max_ind][max_ind][max_ind][max_ind], fit_chrm;
int           n_prm, prms[n_prm_max], bns[n_prm_max], n_cell, n_iter;
int           check_range, adj_tune, adj_chrom, n_steps, n_b3, b3s[n_prm_max];
long int      beta_loc1, beta_loc2, beta_loc3, rseed;
float         **A, *b, *w, **U, **V, *db3;
double        ksi[2], dksi[2], Jx, Jy, delta, ksi1[2], nu0_x, nu0_y, eps_nu;
double        beta1[2], beta2[2], beta3[2], beta4[2];
double        bnL_max[mpole_max], eps_ksi, L;
double        scl_dnu, scl_ksi_nl, scl_ksi3, step;
double        alpha_femto[2], beta_femto[2], nu_femto[2];
double        nu_x_min, nu_x_max, nu_y_min, nu_y_max;
ss_vect<tps>  Id_scl;


const double  max_Ax = 20e-3, max_Ay = 10e-3, max_delta = 3e-2;

const int     m_ksi  = 2; // (ksi_x, ksi_y)

void ini_ranf(const int i) { rseed = i; }


double ranf(void)
{

  const int       k = 19;
  const long int  c = 656329, m = 100000001;

  rseed = (k*rseed+c) % m;

  return rseed/1e8;
}


void no_mpoles(void)
{
  int j, k;

  for (j = 0; j < n_elem; j++)
    if (elem[j].kind == Mpole)
      for (k = Sext; k < mpole_max; k++) {
//	cout << "zeroing " << elem[j].Name << endl;
	set_bn(elem[j].Fnum, elem[j].Knum, k, 0.0);
      }
}


void fit_chrom_D(const double ksi_x, const double ksi_y, const double DeltaL,
		 const int n_b3, const int b3s[], const double eps,
		 const bool first)
{
  int           i, j;
  double        s_max, x, y, Delta;
  ss_vect<tps>  nus;
  ofstream      sext_out;

  const bool    debug = false;
  const double  s_cut = 1e-7;

  if (first) {
    danot_(3); get_Map();
    danot_(4); K = MapNorm(Map, g, A1, A0, Map_res, 1); nus = dHdJ(K);
    ksi[X_] = h_ijklm(nus[3], 0, 0, 0, 0, 1);
    ksi[Y_] = h_ijklm(nus[4], 0, 0, 0, 0, 1);
    dksi[X_] = ksi[X_] - ksi_x; dksi[Y_] = ksi[Y_] - ksi_y;

    cout << endl;
    cout << fixed << setprecision(5)
	 << "ksi_x = " << setw(8) << ksi[X_]
	 << ", ksi_y = " << setw(8) << ksi[Y_] << endl;

    for (i = 1; i <= n_b3; i++) {
      for (j = 1; j <= get_n_Kids(b3s[i-1]); j++)
	set_bn_par(b3s[i-1], j, Sext, 7);

      danot_(3); get_Map();
      danot_(4); K = MapNorm(Map, g, A1, A0, Map_res, 1); nus = dHdJ(K);

      A[1][i] = -h_ijklm_p(nus[3], 0, 0, 0, 0, 1, 7);
      A[2][i] = -h_ijklm_p(nus[4], 0, 0, 0, 0, 1, 7);

      for (j = 1; j <= get_n_Kids(b3s[i-1]); j++)
	clr_bn_par(b3s[i-1], j, Sext);
    }

    L = get_L(b3s[0], 1);
    for (i = 1; i <= 2; i++) {
      x = A[i][1] + A[i][2]; y = A[i][1] - A[i][2];
      A[i][1] = x; A[i][2] = A[i][3]; A[i][3] = y;
    }

    if (debug) {
      cout << endl;
      cout << " Ax = b:" << endl;
      cout << endl;
      for (i = 1; i <= m_ksi; i++) {
	for (j = 1; j <= n_b3; j++)
	  cout << scientific << setprecision(3) << setw(11) << A[i][j] << endl;
      }
    }
    
    for (i = 1; i <= m_ksi; i++)
      for (j = 1; j <= 2; j++)
	U[i][j] = A[i][j];

    svdcmp(U, m_ksi, 2, w, V);

    s_max = -1e30;
    for (i = 1; i <= 2; i++)
      s_max = max(w[i], s_max);
  
    cout << endl;
    cout << "singular values:" << endl;
    for (i = 1; i <= 2; i++) {
      cout << scientific << setprecision(3) << setw(10) << w[i];
      if (w[i]/s_max < s_cut) {
	w[i] = 0.0;
	cout << " (zeroed)";
      }
      cout << endl;
    }
  } else {
    Delta = DeltaL/L;
    b[1] = dksi[X_] - A[1][3]*Delta/2.0; b[2] = dksi[Y_] - A[2][3]*Delta/2.0;

    svbksb(U, w, V, m_ksi, 2, b, db3);

    if (debug) {
      cout << endl;
      cout << "b3s:" << endl;
    }
    for (i = 1; i <= n_b3; i++)
      for (j = 1; j <= get_n_Kids(b3s[i-1]); j++) {
	if (i == 1)
	  set_dbn(b3s[0], j, Sext, db3[1]+Delta/2.0);
	else if (i == 2)
	  set_dbn(b3s[1], j, Sext, db3[1]-Delta/2.0);
	else
	  set_dbn(b3s[2], j, Sext, db3[2]);
	if (debug) {
	  cout << fixed << setprecision(5)
	       << setw(9) << get_bnL(b3s[i-1], j, Sext) << endl;
	}
      }

    danot_(3); get_Map();
    danot_(4); K = MapNorm(Map, g, A1, A0, Map_res, 1); nus = dHdJ(K);
    ksi[X_] = h_ijklm(nus[3], 0, 0, 0, 0, 1);
    ksi[Y_] = h_ijklm(nus[4], 0, 0, 0, 0, 1);
    cout << scientific << setprecision(1)
	 << "ksi_x = " << setw(8) << ksi[X_]
	 << ", ksi_y = " << setw(8) << ksi[Y_];
  }
}


void select_h(void)
{
  int  i, j, k, l, m;

  for (i = 0; i <= no_tps-1; i++)
    for (j = 0; j <= no_tps-1; j++)
      for (k = 0; k <= no_tps-1; k++)
	for (l = 0; l <= no_tps-1; l++)
	  for (m = 0; m <= no_tps-1; m++)
	    if (false)
	      h_sel[i][j][k][l][m] = true;
	    else
	      if (i+j+k+l <= 4)
		// enable all non-chromatic driving terms
		h_sel[i][j][k][l][m] = (m == 0)? true : false;

  // linear chromaticity
//  h_sel[1][1][0][0][1] = true; h_sel[0][0][1][1][1] = true;

  // 2nd order amplitude dependent tune shift
  h_sel[2][2][0][0][0] = true; h_sel[1][1][1][1][0] = true;
  h_sel[0][0][2][2][0] = true;

  // 4th order amplitude dependent tune shift
  h_sel[3][3][0][0][0] = true; h_sel[0][0][3][3][0] = true;
  h_sel[2][2][1][1][0] = true; h_sel[1][1][2][2][0] = true;

  // 6th order amplitude dependent tune shift
  h_sel[4][4][0][0][0] = true; h_sel[0][0][4][4][0] = true;
  h_sel[3][3][1][1][0] = true; h_sel[1][1][3][3][0] = true;
  h_sel[2][2][2][2][0] = true;

  if (true) {
    // nonlinear chromaticity
    for (m = 2; m <= no_tps-3; m++) {
      h_sel[1][1][0][0][m] = true; h_sel[0][0][1][1][m] = true;
    }

    // amplitude dependent chromaticity

    for (m = 1; m <= no_tps-5; m++) {
      h_sel[2][2][0][0][m] = true; h_sel[0][0][2][2][m] = true;
      h_sel[1][1][1][1][m] = true;
    }

    for (m = 1; m <= no_tps-7; m++) {
      h_sel[3][3][0][0][m] = true; h_sel[0][0][3][3][m] = true;
      h_sel[2][2][1][1][m] = true; h_sel[1][1][2][2][m] = true;
    }
  }

  // exclude tune
  h_sel[1][1][0][0][0] = false; h_sel[0][0][1][1][0] = false;

  // exclude delta dependence of pathlength
  for (m = 2; m <= no_tps-1; m++)
    h_sel[0][0][0][0][m] = false;

  if (true) {
    // higher order dispersion
    for (m = 2; m <= no_tps-2; m++) {
      h_sel[1][0][0][0][m] = true; h_sel[0][1][0][0][m] = true;
    }

    // delta dependance of the beta functions
    for (m = 1; m <= no_tps-3; m++) {
      h_sel[2][0][0][0][m] = true; h_sel[0][2][0][0][m] = true;
      h_sel[0][0][2][0][m] = true; h_sel[0][0][0][2][m] = true;
    }
  }
}


void H_rnd(const int n_seed, const double DeltaL_max)
{
  const int     m_max = 500;  // max no of constraints

  char          hs[m_max][max_str];
  int           i, j, k, l, m, i1, m1 = 0;
  double        ksi2 = 0.0, *h, DeltaL;
  tps           r, K_re, K_im, g_re, g_im;
  ss_vect<tps>  nus;
  ofstream      sext_out;

  const bool    prt = false, mirror_sym = false;
  const double  scl_ksi1[]    = { 1e4, 1e4 };

  b = vector(1, m_ksi); w = vector(1, n_b3); db3 = vector(1, n_b3);
  A = matrix(1, m_ksi, 1, n_b3); U = matrix(1, m_ksi, 1, n_b3);
  V = matrix(1, n_b3, 1, n_b3);

  h = dvector(1, m_max);

  select_h();

  file_wr(sext_out, "sext_rnd.dat");
 
  no_mpoles();
  fit_chrom_D(ksi1[X_], ksi1[Y_], 0.0, n_b3, b3s, eps_ksi, true);

  cout << endl;
  for (i1 = 1; i1 <= n_seed; i1++) {
    DeltaL = DeltaL_max*2.0*(ranf()-0.5);

    no_mpoles();

    cout << setw(3) << i1 << ", ";
    fit_chrom_D(ksi1[X_], ksi1[Y_], DeltaL, n_b3, b3s, eps_ksi, false);

    for (i = 1; i <= n_prm; i++)
      set_bnL(prms[i-1], bns[i-1], bnL_max[Sext]*2.0*(ranf()-0.5));

    danot_(no_tps-1); get_Map();
    danot_(no_tps); K = MapNorm(Map, g, A1, A0, Map_res, 1);
    CtoR(K, K_re, K_im); K_re = K_re*Id_scl; nus = dHdJ(K);

    g = g*Id_scl;

    CtoR(g, g_re, g_im);

    // mirror symmetric cell => g_re = 0
    m1 = 0;
    for (i = 0; i <= no_tps-1; i++)
      for (j = 0; j <= no_tps-1; j++)
	for (k = 0; k <= no_tps-1; k++)
	  for (l = 0; l <= no_tps-1; l++)
	    for (m = 0; m <= no_tps-1; m++) {
	      if ((0 < i+j+k+l+m) && (i+j+k+l+m <= no_tps-1) &&
		  h_sel[i][j][k][l][m] &&
		  ((fabs(h_ijklm(g_im, i, j, k, l, m)) > 0.0) ||
		   (fabs(h_ijklm(K_re, i, j, k, l, m)) > 0.0))) {
		m1++;
		sprintf(hs[m1-1], "h_%d%d%d%d%d", i, j, k, l, m);
		h[m1] = h_ijklm(g_im, i, j, k, l, m);

		if (is_h_ijklm(1, 1, 0, 0, 1, i, j, k, l, m)) {
		  // horizontal linear chromaticity
		  h[m1] = scl_ksi1[X_]*(n_cell*ksi1[X_]*M_PI*2.0*Jx*delta
					+h_ijklm(K_re, i, j, k, l, m));
		} else if (is_h_ijklm(0, 0, 1, 1, 1, i, j, k, l, m)) {
		  // vertical linear chromaticity
		  h[m1] = scl_ksi1[Y_]*(n_cell*ksi1[Y_]*M_PI*2.0*Jy*delta
					+h_ijklm(K_re, i, j, k, l, m));
		} else if (is_h_ijklm(2, 2, 0, 0, 0, i, j, k, l, m) ||
			   is_h_ijklm(0, 0, 2, 2, 0, i, j, k, l, m) ||
			   is_h_ijklm(1, 1, 1, 1, 0, i, j, k, l, m) ||
			   is_h_ijklm(3, 3, 0, 0, 0, i, j, k, l, m) ||
			   is_h_ijklm(2, 2, 1, 1, 0, i, j, k, l, m) ||
			   is_h_ijklm(1, 1, 2, 2, 0, i, j, k, l, m) ||
			   is_h_ijklm(0, 0, 3, 3, 0, i, j, k, l, m) ||
			   is_h_ijklm(4, 4, 0, 0, 0, i, j, k, l, m) ||
			   is_h_ijklm(0, 0, 4, 4, 0, i, j, k, l, m) ||
			   is_h_ijklm(3, 3, 1, 1, 0, i, j, k, l, m) ||
			   is_h_ijklm(1, 1, 3, 3, 0, i, j, k, l, m) ||
			   is_h_ijklm(2, 2, 2, 2, 0, i, j, k, l, m)) {
		  // amplitude dependent tune shift
		  h[m1] = scl_dnu*h_ijklm(K_re, i, j, k, l, m);
		} else if (is_h_ijklm(1, 1, 0, 0, 2, i, j, k, l, m) ||
			   is_h_ijklm(0, 0, 1, 1, 2, i, j, k, l, m) ||
			   is_h_ijklm(1, 1, 0, 0, 3, i, j, k, l, m) ||
			   is_h_ijklm(0, 0, 1, 1, 3, i, j, k, l, m) ||
			   is_h_ijklm(1, 1, 0, 0, 4, i, j, k, l, m) ||
			   is_h_ijklm(0, 0, 1, 1, 4, i, j, k, l, m) ||
			   is_h_ijklm(1, 1, 0, 0, 5, i, j, k, l, m) ||
			   is_h_ijklm(0, 0, 1, 1, 5, i, j, k, l, m) ||
			   is_h_ijklm(1, 1, 0, 0, 6, i, j, k, l, m) ||
			   is_h_ijklm(0, 0, 1, 1, 6, i, j, k, l, m)) {
		  // nonlinear chromaticity
		  h[m1] = scl_ksi_nl*h_ijklm(K_re, i, j, k, l, m);
		} else if (is_h_ijklm(2, 2, 0, 0, 1, i, j, k, l, m) ||
			   is_h_ijklm(0, 0, 2, 2, 1, i, j, k, l, m) ||
			   is_h_ijklm(1, 1, 1, 1, 1, i, j, k, l, m) ||
			   is_h_ijklm(2, 2, 0, 0, 2, i, j, k, l, m) ||
			   is_h_ijklm(0, 0, 2, 2, 2, i, j, k, l, m) ||
			   is_h_ijklm(1, 1, 1, 1, 2, i, j, k, l, m) ||
			   is_h_ijklm(2, 2, 0, 0, 3, i, j, k, l, m) ||
			   is_h_ijklm(0, 0, 2, 2, 3, i, j, k, l, m) ||
			   is_h_ijklm(1, 1, 1, 1, 3, i, j, k, l, m) ||
			   is_h_ijklm(2, 2, 0, 0, 4, i, j, k, l, m) ||
			   is_h_ijklm(0, 0, 2, 2, 4, i, j, k, l, m) ||
			   is_h_ijklm(1, 1, 1, 1, 4, i, j, k, l, m) ||
			   is_h_ijklm(3, 3, 0, 0, 1, i, j, k, l, m) ||
			   is_h_ijklm(2, 2, 1, 1, 1, i, j, k, l, m) ||
			   is_h_ijklm(1, 1, 2, 2, 1, i, j, k, l, m) ||
			   is_h_ijklm(0, 0, 3, 3, 1, i, j, k, l, m) ||
			   is_h_ijklm(3, 3, 0, 0, 2, i, j, k, l, m) ||
			   is_h_ijklm(2, 2, 1, 1, 2, i, j, k, l, m) ||
			   is_h_ijklm(1, 1, 2, 2, 2, i, j, k, l, m) ||
			   is_h_ijklm(0, 0, 3, 3, 2, i, j, k, l, m)) {
		  // amplitude dependent chromaticity
		  h[m1] = scl_ksi_nl*h_ijklm(K_re, i, j, k, l, m);
		} else if (!mirror_sym) {
		  m1++;
		  sprintf(hs[m1-2], "i_%d%d%d%d%d", i, j, k, l, m);
		  h[m1-1] = h_ijklm(g_im, i, j, k, l, m);
		  sprintf(hs[m1-1], "r_%d%d%d%d%d", i, j, k, l, m);
		  h[m1] = h_ijklm(g_re, i, j, k, l, m);
		}
	      }
	    }

    if (prt) {
      cout << endl;
      cout << "h:" << endl;
      for (i = 1; i <= m1; i++) {
	cout  << setw(3) << i << " " << hs[i-1]
	      << scientific << setprecision(2) << setw(10) << h[i] << endl;
      }
    }
    
    ksi2 = 0.0;
    for (i = 1; i <= m1; i++)
      ksi2 += sqr(h[i]);

    cout << scientific << setprecision(1)
	 << ", ksi2 = " << ksi2 << endl;

    sext_out << endl;
    sext_out << scientific << setprecision(1)
	     << setw(3) << i1 << " ksi2 = " << ksi2;
    for (i = 1; i <= n_prm; i++)
      sext_out << fixed << setprecision(5) 
	       << setw(5) << get_Name(abs(prms[i-1])) << " = "
	       << setw(8) << get_bnL(prms[i-1], 1, bns[i-1]);
    sext_out << fixed << setprecision(5) 
	     << setw(5) << "sm1a" << " = "
	     << setw(8) << get_bnL(get_Fnum("sm1a"), 1, Sext)
	     << setw(5) << "sm1b" << " = "
	     << setw(8) << get_bnL(get_Fnum("sm1b"), 1, Sext)
	     << setw(5) << "sm2" << " = "
	     << setw(8) << get_bnL(get_Fnum("sm2"), 1, Sext);
    sext_out << fixed << setprecision(5)
	     << setw(5) << "DL" << " = " << setw(8) << DeltaL << endl;
    sext_out.flush();
  }

  sext_out.close();

  free_dvector(h, 1, m_max);

  free_vector(b, 1, m_ksi); free_vector(w, 1, n_b3);
  free_vector(db3, 1, n_b3);
  free_matrix(A, 1, m_ksi, 1, n_b3); free_matrix(U, 1, m_ksi, 1, n_b3);
  free_matrix(V, 1, n_b3, 1, n_b3);
}


void get_prm(char *file_name)
{
  char      line[max_str];      
  ifstream  prm_in;

  file_rd(prm_in, file_name);

  do
    prm_in.getline(line, max_str);
  while (strstr(line, "#") != NULL);

  sscanf(line, "%*s %d %lf %lf %lf", &adj_tune, &nu0_x, &nu0_y, &eps_nu);

  prm_in.getline(line, max_str);
  sscanf(line, "%*s %d %lf %lf %lf",
	 &adj_chrom, &ksi1[X_], &ksi1[Y_], &eps_ksi);
  fit_chrm = eps_ksi < 0.0; eps_ksi = fabs(eps_ksi);

  prm_in.getline(line, max_str);
  sscanf(line, "%*s %d", &check_range);

  prm_in.getline(line, max_str);
  sscanf(line, "%*s %lf", &nu_x_min);

  prm_in.getline(line, max_str);
  sscanf(line, "%*s %lf", &nu_x_max);

  prm_in.getline(line, max_str);
  sscanf(line, "%*s %lf", &nu_y_min);

  prm_in.getline(line, max_str);
  sscanf(line, "%*s %lf", &nu_y_max);

  prm_in.getline(line, max_str);
  sscanf(line, "%*s %d", &n_steps);

  prm_in.getline(line, max_str);
  sscanf(line, "%*s %d", &n_cell);

  prm_in.getline(line, max_str);
  sscanf(line, "%*s %lf", &ds_max);

  prm_in.getline(line, max_str);
  sscanf(line, "%*s %lf", &b2_max);

  prm_in.getline(line, max_str);
  sscanf(line, "%*s %lf", &bnL_max[Sext]);

  prm_in.getline(line, max_str);
  sscanf(line, "%*s %lf", &bnL_max[Oct]);

  prm_in.getline(line, max_str);
  sscanf(line, "%*s %lf", &bnL_max[Dec]);

  prm_in.getline(line, max_str);
  sscanf(line, "%*s %lf", &scl_dnu);

  prm_in.getline(line, max_str);
  sscanf(line, "%*s %lf", &scl_ksi_nl);

  prm_in.getline(line, max_str);
  sscanf(line, "%*s %lf", &scl_ksi3);

  prm_in.getline(line, max_str);
  sscanf(line, "%*s %lf", &step);

  cout << endl;
  cout << fixed << setprecision(6)
       << "fit_tune    = " << adj_tune
       << ", nu0_x  = " << nu0_x << ", nu0_y  = " << nu0_y
       << scientific << setprecision(1) << ", eps_nu = " << eps_nu << endl;
  cout << fixed << setprecision(6)
       << "fit_chrom   = " << adj_chrom
       << ", ksi0_x = " << ksi1[X_] << ", ksi0_y = " << ksi1[Y_]
       << scientific << setprecision(1) << ", eps_ksi = " << eps_ksi << endl;
  cout << "check_range = " << check_range << endl;
  cout << endl;
  cout << fixed << setprecision(6)
       << "n_steps = " << n_steps
       << ", nu_x_min = " << nu_x_min << ", nu_x_max = " << nu_x_max
       << ", nu_y_min = " << nu_y_min << ", nu_y_max = " << nu_y_max << endl;
  cout << endl;
  cout << "n_cell      = " << n_cell << endl;
  cout << fixed << setprecision(2)
       << "ds_max      = " << ds_max << endl;
  cout << fixed << setprecision(1)
       << "b2_max      = " << b2_max << endl;
  cout << fixed << setprecision(1)
       << "b3L_max     = " << bnL_max[Sext] << endl;
  cout << fixed << setprecision(1)
       << "b4L_max     = " << bnL_max[Oct] << endl;
  cout << fixed << setprecision(1)
       << "b5L_max     = " << bnL_max[Dec] << endl;
  cout << fixed << setprecision(1)
       << "scl_dnu     = " << scl_dnu << endl;
  cout << fixed << setprecision(1)
       << "scl_ksi_nl  = " << scl_ksi_nl << endl;
  cout << scientific << setprecision(1)
       << "scl_ksi3    = " << scl_ksi3   << endl;
  cout << fixed << setprecision(2)
       << "step        = " << step << endl;
}


void get_b2s(int &n_b2, int b2_Fams[])
{

  n_b2 = 0;
  b2_Fams[n_b2++] = get_Fnum("ql1");
  b2_Fams[n_b2++] = get_Fnum("ql2");
  b2_Fams[n_b2++] = get_Fnum("ql3");

  b2_Fams[n_b2++] = get_Fnum("qh1");
  b2_Fams[n_b2++] = get_Fnum("qh2");
  b2_Fams[n_b2++] = get_Fnum("qh3");
}


void get_prms()
{

  n_prm = 0;

  prms[n_prm] = get_Fnum("sl1"); bns[n_prm++] = Sext;
  prms[n_prm] = get_Fnum("sl2"); bns[n_prm++] = Sext;
  prms[n_prm] = get_Fnum("sl3"); bns[n_prm++] = Sext;

//  prms[n_prm] = get_Fnum("sm1a"); bns[n_prm++] = Sext;
//  prms[n_prm] = get_Fnum("sm1b"); bns[n_prm++] = Sext;
//  prms[n_prm] = get_Fnum("sm2"); bns[n_prm++] = Sext;

  prms[n_prm] = get_Fnum("sh1"); bns[n_prm++] = Sext;
  prms[n_prm] = get_Fnum("sh3"); bns[n_prm++] = Sext;
  prms[n_prm] = get_Fnum("sh4"); bns[n_prm++] = Sext;

  if (n_prm > n_prm_max) {
    cout << "get_prms: n_prm_max exceeded " << n_prm << "(" << n_prm_max
	 << ")" << endl;
    exit(0);
  }

  cout << endl;
  cout << "get_prms: no of multipole families " << n_prm << endl;
}


void chk_lat(double nu[], double ksi[])
{
  double        alpha1[2];
  ss_vect<tps>  nus;

//  get_Map();
  get_COD(10, 1e-10, 0.0, true);
  K = MapNorm(Map, g, A1, A0, Map_res, 1);
  nus = dHdJ(K); get_nu_ksi(nus, nu, ksi); get_ab(alpha1, beta1, 0);
  cout << endl;
  cout << fixed << setprecision(3)
       << "alpha_x  = " << alpha1[X_] << ", alpha_y = " << alpha1[Y_]
       << ", beta_x = " << beta1[X_] << ", beta_y  = " << beta1[Y_] << endl;
  prt_nu(nus);
}


void get_locs()
{
  double  alpha1[2], alpha2[2], alpha3[2];

  beta_loc1 = get_loc(get_Fnum("mp"), 1); get_ab(alpha1, beta1, beta_loc1);
  beta_loc2 = get_loc(get_Fnum("ss"), 1); get_ab(alpha2, beta2, beta_loc2);
//  beta2[X_] = 1.0; beta2[Y_] = 1.0;
  beta_loc3 = get_loc(get_Fnum("ls"), 1); get_ab(alpha3, beta3, beta_loc3);
//  beta3[X_] = 15.0; beta3[Y_] = 3.0;
}


void get_sxt(const char file_name[], const int n)
{
  const int  line_max = 200;

  char      line[line_max];
  int       j;
  double    b3L[n_prm_max], ksi2 = 0.0;
  tps       K_re, K_im;
  ifstream  inf;

  file_rd(inf, file_name);

  // skip blank line
  inf.getline(line, line_max);
  while (!inf.getline(line, line_max).eof()) {
    sscanf(line,
	   "%d %*s = %lf %*s = %lf %*s = %lf %*s = %lf %*s = %lf %*s = %lf"
	   " %*s = %lf %*s = %lf %*s = %lf %*s = %lf",
	   &j, &ksi2, &b3L[0], &b3L[1], &b3L[2], &b3L[3], &b3L[4], &b3L[5],
	   &b3L[6], &b3L[7], &b3L[8]);

    if (j == n) {
      printf("\n");
      printf("%3d %7.1e %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f"
	     " %8.5f\n",
	     j, ksi2, b3L[0], b3L[1], b3L[2], b3L[3], b3L[4], b3L[5], b3L[6],
	     b3L[7], b3L[8]);

      set_bnL(get_Fnum("sl1"), Sext, b3L[0]);
      set_bnL(get_Fnum("sl2"), Sext, b3L[1]);
      set_bnL(get_Fnum("sl3"), Sext, b3L[2]);
      set_bnL(get_Fnum("sh1"), Sext, b3L[3]);
      set_bnL(get_Fnum("sh3"), Sext, b3L[4]);
      set_bnL(get_Fnum("sh4"), Sext, b3L[5]);
      set_bnL(get_Fnum("sm1a"), Sext, b3L[6]);
      set_bnL(get_Fnum("sm1b"), Sext, b3L[7]);
      set_bnL(get_Fnum("sm2"), Sext, b3L[8]);

      break;
    }

    inf.getline(line, line_max);
  }

  danot_(no_tps-1); get_Map();
  danot_(no_tps); K = MapNorm(Map, g, A1, A0, Map_res, 1);
  CtoR(K, K_re, K_im); K_re = K_re*Id_scl; g = g*Id_scl;

  if (false) {
    cout << endl;
    cout << "g:" << endl;
    cout << g;
    cout << endl;
    cout << "K:" << endl;
    cout << K_re;
  }
}


int main()
{
  int              b2_Fams[n_prm_max], n_b2 = 0;
  double           nu[2], ksi[2];
  double           alpha[2], beta[2];
  tps              Hnl, H2, gn, h, h_re, h_im, H_num, dH, H, H_re, H_im;
  tps              K_re, K_im;
  tps              g_re, g_im;
  ss_vect<double>  x;
  ss_vect<tps>     A_inv, ps, ps_re, ps_im, nus, R, R_inv, Map2;
  ss_vect<tps>     Map_Fl, Mn, Mk, J, eta, M1, M2, M3, M4, Mp;
  ofstream         outf, K_out, nus_out, A1_out, J_out;
  ifstream         inf;

  const int     seed = 1005;

  danot_(no_tps-1);

  rad_on    = false; H_exact        = false; totpath_on   = false;
  cavity_on = false; quad_fringe_on = false; emittance_on = false;
  IBS_on    = false;

  ini_ranf(seed);

  rd_mfile("flat_file.dat", elem); rd_mfile("flat_file.dat", elem_tps);

  // initialize the symplectic integrator after energy has been defined
  ini_si();

  // disable from TPSALib- and LieLib log messages
  idprset(-1);

  danot_(3);

  if (true) chk_lat(nu, ksi);

  get_ab(alpha, beta, 0);
  cout << endl;
  cout << fixed << setprecision(3)
       << "alpha_x  = " << setw(6) << alpha[X_]
       << ", alpha_y = " << setw(6) << alpha[Y_]
       << ", beta_x = " << setw(6) << beta[X_]
       << ", beta_y  = " << setw(6) << beta[Y_] << endl;


  if (true) prt_alphac();

  Jx = sqr(max_Ax)/(2.0*beta1[X_]); Jy = sqr(max_Ay)/(2.0*beta1[Y_]);
  delta = max_delta;

  Id_scl.identity();
  Id_scl[x_] *= sqrt(2.0*Jx); Id_scl[px_] *= sqrt(2.0*Jx);
  Id_scl[y_] *= sqrt(2.0*Jy); Id_scl[py_] *= sqrt(2.0*Jy);
  Id_scl[delta_] *= delta;


  get_prm("tune_scan.prm");

  if (adj_tune) {
    get_locs();

    get_b2s(n_b2, b2_Fams);

    fit_tune(nu0_x, nu0_y,
	     beta1[X_], beta1[Y_], beta_loc1,
	     beta2[X_], beta2[Y_], beta_loc2,
	     beta3[X_], beta3[Y_], beta_loc3,
	     n_b2, b2_Fams, eps_nu, true);
  }

  n_b3 = 0;
  b3s[n_b3++] = get_Fnum("sm1a"); b3s[n_b3++] = get_Fnum("sm1b");
  b3s[n_b3++] = get_Fnum("sm2");

  get_prms();

  H_rnd(10000, 1.0);

  if (false) {
    get_sxt("sext_rnd.dat", 5805);
    get_sxt("sext_rnd.dat", 1360);
    get_sxt("sext_rnd.dat", 3416);
    get_sxt("sext_rnd.dat", 9879);
    get_sxt("sext_rnd.dat", 5707);
  }
}
