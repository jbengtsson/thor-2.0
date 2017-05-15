#define NO 5

#include "thor_lib.h"

int  no_tps   = NO,
     ndpt_tps = 5;


void get_b3s(int &n, int b3s[])
{

  n = 0;

  b3s[n++] = get_Fnum("ss1a_1"); b3s[n++] = get_Fnum("ss2a_1");
  b3s[n++] = get_Fnum("ss1b_1"); b3s[n++] = get_Fnum("ss2b_1");
  b3s[n++] = get_Fnum("ss1c_1"); b3s[n++] = get_Fnum("ss2c_1");
  b3s[n++] = get_Fnum("ss1d_1"); b3s[n++] = get_Fnum("ss2d_1");

  b3s[n++] = get_Fnum("ss1a_2"); b3s[n++] = get_Fnum("ss2a_2");
  b3s[n++] = get_Fnum("ss1b_2"); b3s[n++] = get_Fnum("ss2b_2");
  b3s[n++] = get_Fnum("ss1c_2"); b3s[n++] = get_Fnum("ss2c_2");
  b3s[n++] = get_Fnum("ss1d_2"); b3s[n++] = get_Fnum("ss2d_2");

  b3s[n++] = get_Fnum("ss1a_3"); b3s[n++] = get_Fnum("ss2a_3");
  b3s[n++] = get_Fnum("ss1b_3"); b3s[n++] = get_Fnum("ss2b_3");
  b3s[n++] = get_Fnum("ss1c_3"); b3s[n++] = get_Fnum("ss2c_3");
  b3s[n++] = get_Fnum("ss1d_3"); b3s[n++] = get_Fnum("ss2d_3");

  cout << endl;
  cout << "get_b3s: " << n << " parameters" << endl;
}


void get_SRM(const int n, const int b3s[], double **A)
{
  int     i, j, k, m;
  double  *w, **U, **V, **A_inv, **Id;
  tps     g_re, g_im, K_re, K_im;

  const double  eps = 1e-5;

  for (k = 1; k <= n; k++) {
//    if (k % 8 != 1) set_bn_par(b3s[k-1], Sext, 7);
    set_bn_par(b3s[k-1], Sext, 7);

    danot_(no_tps-1); get_Map();
    danot_(no_tps); K = MapNorm(Map, g, A1, A0, Map_res, 1);
    CtoR(K, K_re, K_im); CtoR(g, g_re, g_im);

    m = 0;

    // Linear chromaticity
    A[++m][k] = h_ijklm_p(K_re, 1, 1, 0, 0, 1, 7);
    A[++m][k] = h_ijklm_p(K_re, 0, 0, 1, 1, 1, 7);

    // dbeta/ddelta
    A[++m][k] = h_ijklm_p(g_im, 2, 0, 0, 0, 1, 7);
    A[++m][k] = h_ijklm_p(g_re, 2, 0, 0, 0, 1, 7);

    A[++m][k] = h_ijklm_p(g_im, 0, 0, 2, 0, 1, 7);
    A[++m][k] = h_ijklm_p(g_re, 0, 0, 2, 0, 1, 7);

    // 1st order geometric terms
    A[++m][k] = h_ijklm_p(g_im, 2, 1, 0, 0, 0, 7);
    A[++m][k] = h_ijklm_p(g_re, 2, 1, 0, 0, 0, 7);

    A[++m][k] = h_ijklm_p(g_im, 1, 0, 1, 1, 0, 7);
    A[++m][k] = h_ijklm_p(g_re, 1, 0, 1, 1, 0, 7);

    A[++m][k] = h_ijklm_p(g_im, 3, 0, 0, 0, 0, 7);
    A[++m][k] = h_ijklm_p(g_re, 3, 0, 0, 0, 0, 7);

    A[++m][k] = h_ijklm_p(g_im, 1, 0, 2, 0, 0, 7);
    A[++m][k] = h_ijklm_p(g_re, 1, 0, 2, 0, 0, 7);

    A[++m][k] = h_ijklm_p(g_im, 1, 0, 0, 2, 0, 7);
    A[++m][k] = h_ijklm_p(g_re, 1, 0, 0, 2, 0, 7);

    // 2nd order momentum compaction
//    A[++m][k] = h_ijklm_p(Map[ct_], 0, 0, 0, 0, 2, 7);

    // 3rd order momentum compaction
//    A[++m][k] = h_ijklm_p(Map[ct_], 0, 0, 0, 0, 3, 7);

    if (false) {
      // 2nd order geometric terms
      A[++m][k] = h_ijklm_p(g_im, 2, 0, 1, 1, 0, 7);
      A[++m][k] = h_ijklm_p(g_re, 2, 0, 1, 1, 0, 7);

      A[++m][k] = h_ijklm_p(g_im, 1, 1, 2, 0, 0, 7);
      A[++m][k] = h_ijklm_p(g_re, 1, 1, 2, 0, 0, 7);

      A[++m][k] = h_ijklm_p(g_im, 3, 1, 0, 0, 0, 7);
      A[++m][k] = h_ijklm_p(g_re, 3, 1, 0, 0, 0, 7);

      A[++m][k] = h_ijklm_p(g_im, 0, 0, 3, 1, 0, 7);
      A[++m][k] = h_ijklm_p(g_re, 0, 0, 3, 1, 0, 7);

      A[++m][k] = h_ijklm_p(g_im, 4, 0, 0, 0, 0, 7);
      A[++m][k] = h_ijklm_p(g_re, 4, 0, 0, 0, 0, 7);

      A[++m][k] = h_ijklm_p(g_im, 0, 0, 4, 0, 0, 7);
      A[++m][k] = h_ijklm_p(g_re, 0, 0, 4, 0, 0, 7);

      A[++m][k] = h_ijklm_p(g_im, 2, 0, 2, 0, 0, 7);
      A[++m][k] = h_ijklm_p(g_re, 2, 0, 2, 0, 0, 7);

      A[++m][k] = h_ijklm_p(g_im, 2, 0, 0, 2, 0, 7);
      A[++m][k] = h_ijklm_p(g_re, 2, 0, 0, 2, 0, 7);

      // Amplitude dependent tune shift
      A[++m][k] = h_ijklm_p(K_re, 2, 2, 0, 0, 0, 7);
      A[++m][k] = h_ijklm_p(K_re, 1, 1, 1, 1, 0, 7);
      A[++m][k] = h_ijklm_p(K_re, 0, 0, 2, 2, 0, 7);
    }

//    if (k % 8 != 1) clr_bn_par(b3s[k-1], Sext);
    clr_bn_par(b3s[k-1], Sext);
  }

  cout << endl;
  cout << scientific << setprecision(3)
       << "alphac = " << h_ijklm(Map[ct_], 0, 0, 0, 0, 1)/elem[n_elem-1].S
       << endl;

  cout << endl;
  dWriteMatrix(stdout, "A:", A, m, n, "%11.3e");

  w = dvector(1, n);
  A_inv = dmatrix(1, n, 1, m);
  U = dmatrix(1, m, 1, n); V = dmatrix(1, n, 1, n);
  Id = dmatrix(1, m, 1, m);

  dmcopy(A, m, n, U); dsvdcmp(U, m, n, w, V);

  cout << endl;
  cout << "singular values:" << endl;
  for (k = 1; k <= n; k++) {
    cout << scientific << setprecision(3) << setw(11) << w[k];
    if (k % 10 == 0) cout << endl;
    if (fabs(w[k]) < eps) {
      w[k] = 0.0; cout << " (zeroed)  ";
    }
  }
  if (n % 10 != 0) cout << endl;

  for (i = 1; i <= n; i++)
    for (j = 1; j <= m; j++) {
      A_inv[i][j] = 0.0;
      for (k = 1; k <= n; k++)
	if (w[k] != 0.0) A_inv[i][j] += V[i][k]/w[k]*U[j][k];
  }

  cout << endl;
  for (i = 0; i < n; i++)
    cout << setw(11) << get_Name(b3s[i]);
  cout << endl;

  dWriteMatrix(stdout, "A^-1:", A_inv, n, m, "%11.3e");

  dmmult(A, m, n, A_inv, n, m, Id);

  cout << endl;
  dWriteMatrix(stdout, "AxA^-1:", Id, m, m, "%11.3e");

  free_dvector(w, 1, n);
  free_dmatrix(A_inv, 1, n, 1, m);
  free_dmatrix(U, 1, m, 1, n); free_dmatrix(V, 1, n, 1, n);
  free_dmatrix(Id, 1, m, 1, m);
}


int main()
{
  const int  m_max = 18, n_max = 24;

  int     n_b3, b3s[n_max];
  double  **A;


  danot_(no_tps-1);

  rad_on    = false; H_exact        = false; totpath_on   = false;
  cavity_on = false; quad_fringe_on = false; emittance_on = false;
  IBS_on    = false;


  rd_mfile("flat_file.dat", elem); rd_mfile("flat_file.dat", elem_tps);

  // initialize the symplectic integrator after energy has been defined
  ini_si();

  // disable from TPSALib- and LieLib log messages
  idprset(-1);

  A = dmatrix(1, m_max, 1, n_max);

  get_b3s(n_b3, b3s); get_SRM(n_b3, b3s, A);

  free_dmatrix(A, 1, m_max, 1, n_max);
}
