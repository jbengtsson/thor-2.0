#define NO 8
#include "mpi.h"
#include "thor_lib.h"

#define MASTER_PROC 0
#define WORKTAG     1
#define DIETAG      0

#define DBG_LOG(msg, args...) \
  printf("[%d, %s]: "msg, mpi_rank, __func__, ##args)

// Example:
//
//   DBG_LOG("%f %f\n", min_realvar[0], max_realvar[1]);

int no_tps   = NO,
    ndpt_tps = 5;

extern double       b2_max;
extern tps          K, g;
extern ss_vect<tps> Map, A0, A1, Map_res;

const int max_ind = 10, n_max = 50;

bool         h[max_ind][max_ind][max_ind][max_ind][max_ind], fit_chrm;
int          n_prm, prms[n_prm_max], bns[n_prm_max];
int          adj_tune, adj_chrom, n_iter, sext_scheme;
long int     beta_loc1, beta_loc2, beta_loc3, beta_loc4, rseed;
double       Jx, Jy, delta, beta1[2], beta2[2], beta3[2], beta4[2];
double       b3s[n_prm_max], chi2 = 0e0;
double       nu0[2], eps_nu, ksi1[2], eps_ksi;
double       nu_x_min, nu_x_max, nu_y_min, nu_y_max;
double       bnL_max[mpole_max];
double       scl_res, scl_dnu, scl_ksi1, scl_ksi_nl, scl_dnu2, scl_L = 0.0;
double       step;
ss_vect<tps> Id_scl, map1, map2;

// const double A_max[] = {20e-3, 7.5e-3}, delta_max = 3e-2;
const double A_max[] = {30e-3, 10e-3}, delta_max = 3.5e-2;

const int n_index = 5;

typedef int unit_of_work_t[n_index];

int  mpi_rank = 0, mpi_size = 1;


void select_h(void)
{
  int  i, j, k, l, m;

  for (i = 0; i <= no_tps-1; i++)
    for (j = 0; j <= no_tps-1; j++)
      for (k = 0; k <= no_tps-1; k++)
	for (l = 0; l <= no_tps-1; l++)
	  for (m = 0; m <= no_tps-1; m++)
	    if ((3 <= i+j+k+l)&& (i+j+k+l <= 4))
		// enable all non-chromatic driving terms
		h[i][j][k][l][m] = (m == 0)? true : false;

  // linear chromaticity
  h[1][1][0][0][1] = true; h[0][0][1][1][1] = true;

  // 2nd order amplitude dependent tune shift
  h[2][2][0][0][0] = true; h[1][1][1][1][0] = true; h[0][0][2][2][0] = true;

  // 4th order amplitude dependent tune shift
  h[3][3][0][0][0] = true; h[2][2][1][1][0] = true;
  h[0][0][3][3][0] = true; h[1][1][2][2][0] = true;

  // 6th order amplitude dependent tune shift
  h[4][4][0][0][0] = true; h[0][0][4][4][0] = true;
  h[3][3][1][1][0] = true; h[1][1][3][3][0] = true;
  h[2][2][2][2][0] = true;

  if (true) {
    // nonlinear chromaticity
    for (m = 2; m <= no_tps-3; m++) {
      h[1][1][0][0][m] = true; h[0][0][1][1][m] = true;
    }
  }

  if (true) {
    // amplitude dependent chromaticity

    for (m = 1; m <= no_tps-5; m++) {
      h[2][2][0][0][m] = true; h[0][0][2][2][m] = true;
      h[1][1][1][1][m] = true;
    }

    for (m = 1; m <= no_tps-7; m++) {
      h[3][3][0][0][m] = true; h[0][0][3][3][m] = true;
      h[2][2][1][1][m] = true; h[1][1][2][2][m] = true;
    }
  }

  // exclude tune
  h[1][1][0][0][0] = false; h[0][0][1][1][0] = false;

  // exclude delta dependence of pathlength
  for (m = 2; m <= no_tps-1; m++)
    h[0][0][0][0][m] = false;

  if (true) {
    // higher order dispersion
    for (m = 2; m <= no_tps-2; m++) {
      h[1][0][0][0][m] = true; h[0][1][0][0][m] = true;
    }

    // delta dependance of the beta functions
    for (m = 1; m <= no_tps-3; m++) {
      h[2][0][0][0][m] = true; h[0][2][0][0][m] = true;
      h[0][0][2][0][m] = true; h[0][0][0][2][m] = true;
    }
  }
}

tps Int(const tps &a, const int k)
{
  char    name[11];
  int     j, n, jj[ss_dim], ibuf1[bufsize], ibuf2[bufsize];
  double  rbuf[bufsize];
  tps     b;

  a.exprt(rbuf, ibuf1, ibuf2, name); n = (int)rbuf[0];

  for (j = 0; j < n; j++) {
    dehash_(no_tps, ss_dim, ibuf1[j], ibuf2[j], jj);

    if (jj[k-1] != 0) rbuf[j+1] /= jj[k-1] + 1e0;

    hash_(no_tps, ss_dim, jj, ibuf1[j], ibuf2[j]);
  }

  b.imprt(n, rbuf, ibuf1, ibuf2);

  return b*tps(0e0, k);
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


double f_prm(const int jj[])
{
  int     k, n;
  double  f;

  n = 0;
  for (k = 0; k < 2*nd_tps; k++)
    n += jj[k];

  // O{K} = no, O{nu ~ dK/dJ} = no-2.
  f =
    ((jj[ss_dim-1] > 1) || (jj[delta_] > 3) ||
     ((n == no_tps-2) && (jj[ss_dim-1] == 0)))?
    0e0 : 1e0;

 return f;
}


void get_dnu2(const tps &K, tps dnu2[])
{
  // no must be 8 (for delta^3).
  int            i, j, k;
  double         twoJ[2], twoJ1[2], delta1;
  ss_vect<tps>   Id_tr, dnus, ps;
  ostringstream  str;
  ofstream       outf;

  const bool  prt    = true;
  const int   n_step = 25;

  twoJ[X_] = 2e0*Jx; twoJ[Y_] = 2e0*Jy;

  dnus = dHdJ(K);
  for (k = 0; k < 2; k++)
    dnus[3+k] -= dnus[3+k].cst();

  Id_tr.identity(); Id_tr[px_] = 1e0; Id_tr[py_] = 1e0;
  for (k = 0; k < 2; k++)
    dnu2[k] = dacfu1(dnus[3+k], f_prm)*Id_tr;

  for (i = 0; i < 2; i++) {
    dnu2[i] = sqr(dnu2[i]);

    dnu2[i] = Int(dnu2[i], 1);
    ps.identity(); ps[0] = twoJ[X_]; dnu2[i] = dnu2[i]*ps;

    dnu2[i] = Int(dnu2[i], 3);
    ps.identity(); ps[2] = twoJ[Y_]; dnu2[i] = dnu2[i]*ps;

    dnu2[i] = Int(dnu2[i], 5);
    ps.identity(); ps[delta_] = delta_max; dnu2[i] = dnu2[i]*ps;

    dnu2[i] /= twoJ[X_]*twoJ[Y_]*delta_max;
  }

  if (prt) {
    str.str("get_dnu2_1.out"); outf.open(str.str().c_str());

    ps.zero(); 
    for (j = 0; j <= n_step; j++) {
      twoJ1[X_] = j*twoJ[X_]/n_step;
      ps[x_] = sqrt(twoJ1[X_]); ps[px_] = sqrt(twoJ1[X_]);
      for (k = 0; k <= n_step; k++) {
	twoJ1[Y_] = k*twoJ[Y_]/n_step;
	ps[y_] = sqrt(twoJ1[Y_]); ps[py_] = sqrt(twoJ1[Y_]);

	outf << scientific << setprecision(5)
	     << setw(13) << twoJ1[X_] << setw(13) << twoJ1[Y_]
	     << setw(13) << (dnus[3]*ps).cst()
	     << setw(13) << (dnus[4]*ps).cst() << endl;
      }
      outf << endl;
    }

    outf.close();

    str.str("get_dnu2_2.out"); outf.open(str.str().c_str());

    ps.zero(); 
    for (j = -n_step; j <= n_step; j++) {
      delta1 = j*delta_max/n_step;
      ps[delta_] = delta1;
      for (k = 0; k <= n_step; k++) {
	twoJ1[X_] = k*twoJ[X_]/n_step;
	ps[x_] = sqrt(twoJ1[X_]); ps[px_] = sqrt(twoJ1[X_]);

	outf << scientific << setprecision(5)
	     << setw(13) << delta1 << setw(13) << twoJ1[X_]
	     << setw(13) << (dnus[3]*ps).cst()
	     << setw(13) << (dnus[4]*ps).cst() << endl;
      }
      outf << endl;
    }

    outf.close();
  }
}


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


 void get_prm(const char *file_name)
{
  char      line[max_str];      
  ifstream  prm_in;

  file_rd(prm_in, file_name);

  do
    prm_in.getline(line, max_str);
  while (strstr(line, "#") != NULL);

  sscanf(line, "%*s %d %lf %lf %lf", &adj_tune, &nu0[X_], &nu0[Y_], &eps_nu);

  prm_in.getline(line, max_str);
  sscanf(line, "%*s %d %lf %lf %lf",
	 &adj_chrom, &ksi1[X_], &ksi1[Y_], &eps_ksi);
  fit_chrm = eps_ksi < 0.0; eps_ksi = fabs(eps_ksi);

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
  sscanf(line, "%*s %lf", &scl_res);

  prm_in.getline(line, max_str);
  sscanf(line, "%*s %lf", &scl_dnu);

  prm_in.getline(line, max_str);
  sscanf(line, "%*s %lf", &scl_ksi1);

  prm_in.getline(line, max_str);
  sscanf(line, "%*s %lf", &scl_ksi_nl);

  prm_in.getline(line, max_str);
  sscanf(line, "%*s %lf", &scl_dnu2);

  prm_in.getline(line, max_str);
  sscanf(line, "%*s %lf", &step);

  cout << endl;
  cout << fixed << setprecision(6)
       << "fit_tune      = " << adj_tune
       << ", nu0_x  = " << nu0[X_] << ", nu0_y  = " << nu0[Y_]
       << scientific << setprecision(1) << ", eps_nu = " << eps_nu << endl;
  cout << fixed << setprecision(6)
       << "fit_chrom     = " << adj_chrom
       << ", ksi0_x = " << ksi1[X_] << ", ksi0_y = " << ksi1[Y_]
       << scientific << setprecision(1) << ", eps_ksi = " << eps_ksi
       << ", fit_chrm = " << fit_chrm << endl;
  cout << endl;
  cout << fixed << setprecision(2)
       << "ds_max        = " << ds_max << endl;
  cout << fixed << setprecision(1)
       << "b2_max        = " << b2_max << endl;
  cout << fixed << setprecision(1)
       << "b3L_max       = " << bnL_max[Sext] << endl;
  cout << fixed << setprecision(1)
       << "b4L_max       = " << bnL_max[Oct] << endl;
  cout << fixed << setprecision(1)
       << "b5L_max       = " << bnL_max[Dec] << endl;
  cout << fixed << setprecision(1)
       << "scl_res       = " << scl_res << endl;
  cout << fixed << setprecision(1)
       << "scl_dnu       = " << scl_dnu << endl;
  cout << scientific << setprecision(1)
       << "scl_ksi1      = " << scl_ksi1 << endl;
  cout << fixed << setprecision(1)
       << "scl_ksi_nl    = " << scl_ksi_nl << endl;
  cout << fixed << setprecision(1)
       << "scl_dnu2      = " << scl_dnu2 << endl;
  cout << fixed << setprecision(2)
       << "step          = " << step << endl;
}


void get_prms()
{

  n_prm = 0;

  prms[n_prm] = get_Fnum("sl1"); bns[n_prm++] = Sext;
  prms[n_prm] = get_Fnum("sl2"); bns[n_prm++] = Sext;
  prms[n_prm] = get_Fnum("sl3"); bns[n_prm++] = Sext;

  prms[n_prm] = get_Fnum("sh1"); bns[n_prm++] = Sext;
  prms[n_prm] = get_Fnum("sh3"); bns[n_prm++] = Sext;
  prms[n_prm] = get_Fnum("sh4"); bns[n_prm++] = Sext;

  switch (sext_scheme) {
  case 0:
    prms[n_prm] = get_Fnum("sm1"); bns[n_prm++] = Sext;
    prms[n_prm] = get_Fnum("sm2"); bns[n_prm++] = Sext;
    prms[n_prm] = get_Fnum("sm3"); bns[n_prm++] = Sext;
    break;
  case 1:
    prms[n_prm] = get_Fnum("sm1a"); bns[n_prm++] = Sext;
    prms[n_prm] = get_Fnum("sm1b"); bns[n_prm++] = Sext;
    prms[n_prm] = get_Fnum("sm2"); bns[n_prm++] = Sext;
    break;
  case 2:
    prms[n_prm] = get_Fnum("sm1a"); bns[n_prm++] = Sext;
    prms[n_prm] = get_Fnum("sm1b"); bns[n_prm++] = Sext;
    prms[n_prm] = get_Fnum("sm2");  bns[n_prm++] = Sext;

    prms[n_prm] = get_Fnum("sm3");  bns[n_prm++] = Sext;
    break;
  case 3:
    prms[n_prm] = get_Fnum("sm1a"); bns[n_prm++] = Sext;
    prms[n_prm] = get_Fnum("sm1b"); bns[n_prm++] = Sext;
    prms[n_prm] = get_Fnum("sm2");  bns[n_prm++] = Sext;

    prms[n_prm] = get_Fnum("sm3");  bns[n_prm++] = Oct;
    break;
  case 4:
    prms[n_prm] = get_Fnum("sm1a"); bns[n_prm++] = Sext;
    prms[n_prm] = get_Fnum("sm1b"); bns[n_prm++] = Sext;
    prms[n_prm] = get_Fnum("sm2");  bns[n_prm++] = Sext;

    prms[n_prm] = get_Fnum("sm3a");  bns[n_prm++] = Oct;
    prms[n_prm] = get_Fnum("sm3a");  bns[n_prm++] = Dec;
    prms[n_prm] = get_Fnum("sm3b");  bns[n_prm++] = Oct;
    prms[n_prm] = get_Fnum("sm3b");  bns[n_prm++] = Dec;
    break;
  case 5:
    prms[n_prm] = get_Fnum("sm1a"); bns[n_prm++] = Sext;
    prms[n_prm] = get_Fnum("sm1b"); bns[n_prm++] = Sext;
    prms[n_prm] = get_Fnum("sm2");  bns[n_prm++] = Sext;

    prms[n_prm] = get_Fnum("sm3a");  bns[n_prm++] = Sext;
    prms[n_prm] = get_Fnum("sm3a");  bns[n_prm++] = Oct;
    prms[n_prm] = get_Fnum("sm3a");  bns[n_prm++] = Dec;
    prms[n_prm] = get_Fnum("sm3b");  bns[n_prm++] = Sext;
    prms[n_prm] = get_Fnum("sm3b");  bns[n_prm++] = Oct;
    prms[n_prm] = get_Fnum("sm3b");  bns[n_prm++] = Dec;
    break;
  }

  if (n_prm > n_prm_max) {
    cout << "get_prms: n_prm_max exceeded " << n_prm << "(" << n_prm_max
	 << ")" << endl;
    exit(0);
  }

  cout << endl;
  cout << "get_prms: no of multipole families " << n_prm << endl;
}


void chk_lat(void)
{
  double        nu[2], ksi[2], alpha1[2];
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


 void fit_chrom(void)
{
  int k, n_b3, b3s[n_prm_max];

  n_b3 = 0;
  for (k = 0; k < n_prm; k++) {
    if (bns[k] == Sext) {
      n_b3++; b3s[n_b3-1] = prms[k];
    }
  }

  // fit chromaticity
  cavity_on = false;
  fit_chrom(ksi1[X_], ksi1[Y_], n_b3, b3s, true);
  //  fit_chrom1(0.0, 0.0, n_prm, prms, eps_ksi, true);
}


void model_init(void)
{
  double  alpha[2], beta[2];

  danot_(3);

  if (true) chk_lat();

  get_ab(alpha, beta, 0);
  cout << endl;
  cout << fixed << setprecision(3)
       << "alpha_x  = " << setw(6) << alpha[X_]
       << ", alpha_y = " << setw(6) << alpha[Y_]
       << ", beta_x = " << setw(6) << beta[X_]
       << ", beta_y  = " << setw(6) << beta[Y_] << endl;

  if (true) prt_alphac();

  if (false) {
    // no scaling
    Jx = 0.5, Jy = 0.5, delta = 1.0;
  } else {
    Jx = sqr(A_max[X_])/(2.0*beta1[X_]); Jy = sqr(A_max[Y_])/(2.0*beta1[Y_]);
    delta = delta_max;
  }

  Id_scl.identity();
  Id_scl[x_] *= sqrt(2e0*Jx); Id_scl[px_] *= sqrt(2e0*Jx);
  Id_scl[y_] *= sqrt(2e0*Jy); Id_scl[py_] *= sqrt(2e0*Jy);
  Id_scl[delta_] *= delta;

  get_prm("nl_opt.dat");

  n_prm = 0;
  switch (sext_scheme) {
  case 0:
    prms[n_prm] = get_Fnum("sm1"); bns[n_prm++] = Sext;
    prms[n_prm] = get_Fnum("sm2"); bns[n_prm++] = Sext;
    break;
  case 1:
  case 2:
  case 3:
  case 4:
  case 5:
    prms[n_prm] = get_Fnum("sm1a"); bns[n_prm++] = Sext;
    prms[n_prm] = get_Fnum("sm1b"); bns[n_prm++] = Sext;
    prms[n_prm] = get_Fnum("sm2");  bns[n_prm++] = Sext;
    break;
  }

  if (fit_chrm) {
    danot_(3);
    no_mpoles(); fit_chrom();
  }

  get_prms();
}


void prep_in_buf(unit_of_work_t in_buf[], int &n_work)
{
  int i, j, k, l, m;

  select_h();

  n_work = 0;
  for (i = 0; i <= no_tps-1; i++)
    for (j = 0; j <= no_tps-1-i; j++)
      for (k = 0; k <= no_tps-1-i-j; k++)
	for (l = 0; l <= no_tps-1-i-j-k; l++)
	  for (m = 0; m <= no_tps-1-i-j-k-l; m++)
	    if (h[i][j][k][l][m]) {
	      n_work++;
	      in_buf[n_work-1][0] = i; in_buf[n_work-1][1] = j;
	      in_buf[n_work-1][2] = k; in_buf[n_work-1][3] = l;
	      in_buf[n_work-1][4] = m;
	    }
}


void get_next_work_item(int &n_sent)
{
  int rank;

  n_sent++;

  if (n_sent > n_max) {
    printf("*** error get_next_work_item: n_sebt (%d) > n_max (%d)\n",
	   n_sent, n_max);
    for (rank = 1; rank < mpi_size; rank++)
      MPI_Send(0, 0, MPI_INT, rank, DIETAG, MPI_COMM_WORLD);
    exit(1);
  }
}


void process_results(const double result[], double out_buf[])
{
  int k;

  for (k = 0; k < n_prm; k++)
    out_buf[k] = result[k];
}


void do_work(const int i1, int &m1,
	    const int i, const int j, const int k, const int l, const int m,
	    char hs[][max_str], double *A[], double b[])
{
  tps K_re, K_im, K_re_scl, g_re, g_im, g_re_scl, g_im_scl, dnu2[2];

  const bool mirror_sym = false;

  if (prms[i1-1] > 0)
    set_bn_par(prms[i1-1], bns[i1-1], 7);
  else
    set_s_par(abs(prms[i1-1]), 7);

  danot_(no_tps-1);

  get_Map();

  danot_(no_tps);

  K = MapNorm(Map, g, A1, A0, Map_res, 1);

  CtoR(K, K_re, K_im); CtoR(g, g_re, g_im);

  K_re_scl = K_re*Id_scl;

  g_re_scl = scl_res*(g_re*Id_scl);
  g_im_scl = scl_res*(g_im*Id_scl);

  get_dnu2(K, dnu2);

  m1 = 0;

  A[++m1][i1] = scl_L*h_ijklm_p(scl_dnu2*dnu2[0], 0, 0, 0, 0, 0, 7);
  if (i1 == n_prm) {
    sprintf(hs[m1-1], "Ax     ");
    b[m1] = -(h_ijklm(scl_dnu2*dnu2[0], 0, 0, 0, 0, 0));
  }

  A[++m1][i1] = scl_L*h_ijklm_p(scl_dnu2*dnu2[1], 0, 0, 0, 0, 0, 7);
  if (i1 == n_prm) {
    sprintf(hs[m1-1], "Ay     ");
    b[m1] = -(h_ijklm(scl_dnu2*dnu2[1], 0, 0, 0, 0, 0));
  }

  if (h[i][j][k][l][m] &&
      (fabs(h_ijklm(scl_ksi1*K_re_scl, i, j, k, l, m))
       > 0e0)) {
    if (is_h_ijklm(1, 1, 0, 0, 1, i, j, k, l, m)) {
      // horizontal linear chromaticity
      A[++m1][i1] =
	scl_L*scl_ksi1*h_ijklm_p(K_re_scl, i, j, k, l, m, 7);
      if (i1 == n_prm) {
	sprintf(hs[m1-1], "K_%d%d%d%d%d", i, j, k, l, m);
	b[m1] =
	  -scl_ksi1*(ksi1[X_]*M_PI*2.0*Jx*delta
		     +h_ijklm(K_re_scl, i, j, k, l, m));
      }
    } else if (is_h_ijklm(0, 0, 1, 1, 1, i, j, k, l, m)) {
      // vertical linear chromaticity
      A[++m1][i1] =
	scl_L*scl_ksi1*h_ijklm_p(K_re_scl, i, j, k, l, m, 7);
      if (i1 == n_prm) {
	sprintf(hs[m1-1], "K_%d%d%d%d%d", i, j, k, l, m);
	b[m1] =
	  -scl_ksi1*(ksi1[Y_]*M_PI*2.0*Jy*delta
		     +h_ijklm(K_re_scl, i, j, k, l, m));
      }
    }
  }

  if (h[i][j][k][l][m] &&
      (fabs(h_ijklm(scl_ksi_nl*K_re_scl, i, j, k, l, m))
       > 0e0)) {
    if (is_h_ijklm(1, 1, 0, 0, 2, i, j, k, l, m) ||
	is_h_ijklm(0, 0, 1, 1, 2, i, j, k, l, m) ||
	is_h_ijklm(1, 1, 0, 0, 3, i, j, k, l, m) ||
	is_h_ijklm(0, 0, 1, 1, 3, i, j, k, l, m)) {
      // nonlinear chromaticity
      A[++m1][i1] =
	scl_L*scl_ksi_nl
	*h_ijklm_p(K_re_scl, i, j, k, l, m, 7);
      if (i1 == n_prm) {
	sprintf(hs[m1-1], "K_%d%d%d%d%d", i, j, k, l, m);
	b[m1] = -scl_ksi_nl*h_ijklm(K_re_scl, i, j, k, l, m);
      }
    }
  }

  if (h[i][j][k][l][m] &&
      (fabs(h_ijklm(scl_dnu*K_re_scl, i, j, k, l, m))
       > 0e0)) {
    if (is_h_ijklm(2, 2, 0, 0, 0, i, j, k, l, m) ||
	is_h_ijklm(0, 0, 2, 2, 0, i, j, k, l, m) ||
	is_h_ijklm(1, 1, 1, 1, 0, i, j, k, l, m)) {
      // amplitude dependent tune shift
      A[++m1][i1] =
	scl_L*scl_dnu*h_ijklm_p(K_re_scl, i, j, k, l, m, 7);
      if (i1 == n_prm) {
	sprintf(hs[m1-1], "K_%d%d%d%d%d", i, j, k, l, m);
	b[m1] = -scl_dnu*h_ijklm(K_re_scl, i, j, k, l, m);
      }
    }
  }

  if (h[i][j][k][l][m] &&
      (fabs(h_ijklm(g_im_scl, i, j, k, l, m)) > 0e0)) {
    A[++m1][i1] = scl_L*h_ijklm_p(g_im_scl, i, j, k, l, m, 7);
    if (i1 == n_prm) {
      sprintf(hs[m1-1], "i_%d%d%d%d%d", i, j, k, l, m);
      b[m1] = -h_ijklm(g_im_scl, i, j, k, l, m);
    }

    if (!mirror_sym) {
      A[++m1][i1] =
	scl_L*h_ijklm_p(g_re_scl, i, j, k, l, m, 7);
      if (i1 == n_prm) {
	sprintf(hs[m1-1], "r_%d%d%d%d%d", i, j, k, l, m);
	b[m1] = -h_ijklm(g_re_scl, i, j, k, l, m);
      }
    }
  }
}


void master(int &n_work)
{
  typedef double unit_result_t[n_prm];

  unit_of_work_t *in_buf;
  unit_result_t  *out_buf;

  int            rank, n_sent;
  unit_result_t  result;
  MPI_Status     status;

  // nok(no+nv, nv)
  const int n_mn = nok(no_tps-1+n_index, n_index);

  in_buf  = (unit_of_work_t*) malloc(n_mn*sizeof(unit_of_work_t));
  out_buf = (unit_result_t*) malloc(n_mn*sizeof(unit_result_t));

  prep_in_buf(in_buf, n_work);

  if (n_work <= n_mn)
    printf("n_work = %d (%d)\n", n_work, n_mn);
  else {
    printf("*** error prep_in_buf: n_work (%d) > n_mn (%d)\n", n_work, n_mn);
    for (rank = 1; rank < mpi_size; rank++)
      MPI_Send(0, 0, MPI_INT, rank, DIETAG, MPI_COMM_WORLD);
    exit(1);
  }

  n_sent = 0;
  for (rank = 1; rank < min(mpi_size, n_work); rank++) {
    get_next_work_item(n_sent);
    MPI_Send(in_buf[n_sent-1], n_index, MPI_INT, rank, n_sent, MPI_COMM_WORLD);
  }

  while (n_sent < n_work) {
    MPI_Recv(result, n_prm, MPI_DOUBLE, MPI_ANY_SOURCE, MPI_ANY_TAG,
	     MPI_COMM_WORLD, &status);
    process_results(result, out_buf[status.MPI_TAG]);

    get_next_work_item(n_sent);
    MPI_Send(in_buf[n_sent-1], n_index, MPI_INT, status.MPI_SOURCE, n_sent,
	     MPI_COMM_WORLD);
  }

  for (rank = 1; rank < mpi_size; rank++) {
    MPI_Recv(result, n_prm, MPI_DOUBLE, MPI_ANY_SOURCE, MPI_ANY_TAG,
	     MPI_COMM_WORLD, &status);
    process_results(result, out_buf[status.MPI_TAG]);
  }


  free(in_buf); free(out_buf);
}


void slave(const int i1, int &m1,
	   char hs[][max_str], double **A, double b[])
{
  unit_of_work_t work;
  double         result[n_prm];
  MPI_Status     status;

  while (true) {
    MPI_Recv(work, n_index, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

    if (status.MPI_TAG == DIETAG) return;

    DBG_LOG("%d: %d %d %d %d %d\n",
	    status.MPI_TAG, work[0], work[1], work[2], work[3], work[4]);

    do_work(i1, m1,
	    work[0], work[1], work[2], work[3], work[4],
	    hs, A, b);

    MPI_Send(result, n_prm, MPI_DOUBLE, 0, status.MPI_TAG, MPI_COMM_WORLD);
  }
}


void H_zero(const double prm_tol, const int n_max, const bool prt_iter)
{
  // Minimize tune footprint.

  const int m_max = pow(no_tps-1+5, 5); // nok(no+nv, nv)

  bool     first = true;
  char     hs[m_max][max_str], str[max_str];
  int      i, j, n, i1, m1 = 0, rank, n_work;
  double   L, dprm_max;
  double   **A, *b, *w, **U, **V, *dbn, *bn, *bn_max, **A_inv;
  tps      K_re, K_im, K_re_scl, g_re, g_im, g_re_scl, g_im_scl, dnu2[2];
  ofstream sext_out;

  const bool prt   = true;
  const int  n_prt = 9;

  if (mpi_rank == MASTER_PROC) {
    b = dvector(1, m_max); w = dvector(1, n_prm); dbn = dvector(1, n_prm);
    bn = dvector(1, n_prm); bn_max = dvector(1, n_prm);
    A = dmatrix(1, m_max, 1, n_prm); U = dmatrix(1, m_max, 1, n_prm);
    V = dmatrix(1, n_prm, 1, n_prm);
    A_inv = dmatrix(1, n_prm, 1, m_max);
  }

  // Initiate all processes.
  for (i = 1; i <= n_prm; i++) {
    if (prms[i-1] > 0) {
      // Note, Jacobian is a function of multipole strengths
      L = get_L(prms[i-1], 1);
      if (L == 0.0) L = 1.0;
      if (false && (strcmp(get_Name(prms[i-1]), "om1") == 0)) {
	cout << endl;
	cout << get_Name(prms[i-1]) << ": bnL_max set to 1.0" << endl;
	bn_max[i] = 1.0/L; 
      } else
	bn_max[i] = bnL_max[bns[i-1]]/L; 
    } else
      bn_max[i] = ds_max;

    bn[i] = get_bn_s(prms[i-1], 1, bns[i-1]);
  }

  if (first) {
    // store initial sextupole strengths
    first = false;
    cout << endl;
    cout << "initial b3s:" << endl;
    for (i = 1; i <= n_prm; i++) {
      b3s[i-1] = get_bnL_s(prms[i-1], 1, bns[i-1]);
      cout << scientific << setprecision(3)
	   << setw(11) << b3s[i-1] << setw(2) << bns[i-1] << endl;
    }
  } else if (true) {
    // initialize sextupoles
    for (i = 1; i <= n_prm; i++)
      set_bnL_s(prms[i-1], bns[i-1], b3s[i-1]);
  }

  n = 0;
  do {
    n++; n_iter++;
    cout << endl;
    for (i1 = 1; i1 <= n_prm; i1++) {
      scl_L = (prms[i1-1] > 0)? 1e0 : scl_ds;

//      Jx = 0.5; Jy = 0.5; delta = 1.0;

      // mirror symmetric cell => g_re = 0


      if ((i == 1) && (mpi_rank == MASTER_PROC))
	master(n_work);
      else
	slave(m1, i1, hs, A, b);


      if (prms[i1-1] > 0)
	clr_bn_par(prms[i1-1], bns[i1-1]);
      else
	clr_s_par(abs(prms[i1-1]));
    }

    if (prt) {
      cout << endl;
      cout << n << " Ax = b:" << endl;
      cout << endl;
      for (i = 1; i <= m1; i++) {
	cout  << setw(3) << i << " " << hs[i-1];
	for (j = 1; j <= n_prm; j++)
	  cout << scientific << setprecision(2) << setw(10) << A[i][j];
	cout << scientific << setprecision(2) << setw(10) << b[i] << endl;
      }
    }
    
    SVD_lim(m1, n_prm, A, b, bn_max, 1e-11, bn, dbn);

    cout << endl;
    cout << "dcorr. (int):" << endl;
    dprm_max = 0.0;
    for (i = 1; i <= n_prm; i++) {
      if (prms[i-1] > 0) {
	L = get_L(prms[i-1], 1);
	if (L == 0.0) L = 1.0;
	dprm_max = max(fabs(step*dbn[i]*L), dprm_max);
	cout << scientific << setprecision(3) << setw(11) << step*dbn[i]*L;
      } else {
	dprm_max = max(fabs(step*scl_ds*dbn[i]), dprm_max);
	cout << scientific << setprecision(3)
	     << setw(11) << step*scl_ds*dbn[i];
      }
      if (i % n_prt == 0) cout << endl;
    }
    if (n_prm % n_prt != 0) cout << endl;

    cout << endl;
    cout << "corr.:" << endl;
    for (i = 1; i <= n_prm; i++) {
      set_dbn_s(prms[i-1], bns[i-1], step*dbn[i]);
      bn[i] = get_bn_s(prms[i-1], 1, bns[i-1]);
      cout << scientific << setprecision(3)
	   << setw(11) << get_bnL_s(prms[i-1], 1, bns[i-1]);
      if (i % n_prt == 0) cout << endl;
    }
    if (n_prm % n_prt != 0) cout << endl;

    if (prt_iter) {
      sext_out << endl;
      sext_out << "n = " << n_iter << ":" << endl;
      for (i = 1; i <= n_prm; i++)
	for (j = 1; j <= get_n_Kids(abs(prms[i-1])); j++) {
	  if (prms[i-1] > 0) 
	    sext_out << fixed << setprecision(7) 
		     << setw(9) << get_Name(abs(prms[i-1]))
		     << "(" << j << ") = "
		     << setw(11) << get_bnL_s(prms[i-1], 1, bns[i-1])
		     << setw(2) << bns[i-1] << endl;
	  else {
	    sprintf(str, "du_%s", get_Name(abs(prms[i-1])));
	    sext_out << fixed << setprecision(7) 
		     << setw(9) << str << "(" << j << ") = "
		     << setw(11) << get_bnL_s(prms[i-1], 1, bns[i-1])
		     << setw(2) << 0 << endl;
	    sprintf(str, "dd_%s", get_Name(abs(prms[i-1])));
	    sext_out << fixed << setprecision(7) 
		      << setw(9) << str << "(" << j << ") = "
		     << setw(11) << -get_bnL_s(prms[i-1], 1, bns[i-1])
		     << setw(2) << 0 << endl;
	  }
	}
      sext_out.flush();
    }

    cout << endl;
    cout << scientific << setprecision(1)
	 << setw(2) << n_iter << " chi2: " << chi2;

    chi2 = 0.0;
    for (i = 1; i <= m1; i++)
      chi2 += sqr(b[i]);

    cout << scientific << setprecision(1)
	 << " -> " << chi2 << endl;
  } while ((dprm_max > prm_tol) && (n < n_max));

  if (!prt_iter) {
    for (i = 1; i <= n_prm; i++)
      for (j = 1; j <= get_n_Kids(abs(prms[i-1])); j++) {
	sext_out << fixed << setprecision(7) 
		 << setw(6) << get_Name(abs(prms[i-1])) << "(" << j << ") = "
		 << setw(11) << get_bnL_s(prms[i-1], 1, bns[i-1])
		 << setw(2) << bns[i-1] << endl;
      }
    sext_out.flush();
  }

  for (rank = 1; rank < min(mpi_size, n_work); rank++)
    MPI_Send(0, 0, MPI_INT, rank, DIETAG, MPI_COMM_WORLD);

  if (mpi_rank == MASTER_PROC) {
    free_dvector(b, 1, m_max); free_dvector(w, 1, n_prm);
    free_dvector(dbn, 1, n_prm); free_dvector(bn, 1, n_prm);
    free_dvector(bn_max, 1, n_prm);
    free_dmatrix(A, 1, m_max, 1, n_prm);
    free_dmatrix(U, 1, m_max, 1, n_prm); free_dmatrix(V, 1, n_prm, 1, n_prm);
    free_dmatrix(A_inv, 1, n_prm, 1, m_max);
  }
}


int main(int argc, char *argv[])
{
  tps           h_re, h_im, H, H_re, H_im, K_re, K_im, g_re, g_im;
  ss_vect<tps>  nus;
  ofstream      outf, K_out, nus_out;

  MPI_Init(&argc, &argv); MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);

  danot_(no_tps-1);

  rad_on    = false; H_exact        = false; totpath_on   = false;
  cavity_on = false; quad_fringe_on = false; emittance_on = false;
  IBS_on    = false;

  sscanf(argv[1], "%d", &sext_scheme);

  rd_mfile("flat_file.dat", elem); rd_mfile("flat_file.dat", elem_tps);

  // initialize the symplectic integrator after energy has been defined
  ini_si();

  // disable from TPSALib- and LieLib log messages
  idprset(-1);

  model_init();

  H_zero(eps_ksi, 50, true);

  if (mpi_rank == MASTER_PROC) {
    danot_(no_tps-1);

    get_Map();

    file_wr(outf, "map.dat"); outf << Map; outf.close();

    danot_(no_tps);

    CtoR(get_h()*Id_scl, h_re, h_im);
    outf.open("h.dat"); outf << h_re; outf.close();

    K = MapNorm(Map, g, A1, A0, Map_res, no_tps); CtoR(K, K_re, K_im);

    outf.open("K.dat");
    outf << K_re << K_re*Id_scl;
    outf.close();

    nus = dHdJ(K);
    outf.open("nus.dat");
    outf << nus[3] << nus[3]*Id_scl << nus[4] << nus[4]*Id_scl;
    outf.close();

    CtoR(get_H()*Id_scl, H_re, H_im);
    outf.open("H.dat");
    outf << H_re << H_re*Id_scl;
    outf.close();

    CtoR(g*Id_scl, g_re, g_im);
    outf.open("g.dat");
    outf << g_im << g_im*Id_scl;
    outf.close();
  }

  MPI_Finalize();
}
