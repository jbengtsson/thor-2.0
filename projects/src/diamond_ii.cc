
#include <cfloat>

#define NO 7

#include "thor_lib.h"

int no_tps = NO,

#define DOF_3 0

#if !DOF_3
  ndpt_tps = 5;
#else
// Requires that cavity is turned on.
ndpt_tps = 0;
#endif


extern tps          K, g;
extern ss_vect<tps> Map, A0, A1, Map_res;

double chi2 = 0e0, *f_lm, **A_lm;

const bool symm = true, tune_conf = false;
const int  n_cell = 2;

// MAX-IV                      1,
// SLS-2                       2,
// DIAMOND                     3,
// DIAMOND with VMX            4,
// DIAMOND-II 4-BA             5,
// DIAMOND-II 6-BA             6,
// DIAMOND-II 6-BA_jb          7,
// DIAMOND-II 8-BA             8,
// DIAMOND-II 8-BA by Hossein  9.
const int lat_case = 9, n_prt = 8;

// Center of straight.
const double
  beta_inj[][2] = {{3.0, 3.0}, {3.4, 1.9}, {9.9, 5.4}, {9.9, 5.4}, {10.6, 8.6},
		   {10.9, 2.9}, {5.3, 2.0}, {3.4, 1.9}, {10.5, 5.2}},
  A_max[][2] =
    {{1.2e-3, 1.2e-3}, {6e-3, 4e-3}, {15e-3, 8e-3}, {15e-3, 8e-3}, {5e-3, 3e-3},
     {6e-3, 4e-3}, {7e-3, 4e-3}, {2e-3, 1e-3}, {2e-3, 1e-3}},
  delta_max[] = {3e-2, 5e-2, 3e-2, 3e-2, 3e-2, 3e-2, 3e-2, 3e-2, 3e-2};


// const double scl_h[] = {1e0, 1e0}, scl_dnu[] = {1e5, 1e0, 1e-1, 1e-9};
const double scl_h[]   = {1e0, 1e0, 1e0},
             scl_dnu[] = {1e-1, 1e-1},
             scl_ksi[] = {1e5, 1e-1, 1e-1};


struct param_type {
private:

public:
  int                 m_constr, n_prm;
  double              bn_tol, svd_cut, step;
  double              *bn_lim, *bn, *dbn;
  std::vector<double> bn_max, bn_scl;
  std::vector<int>    Fnum, n;

  void add_prm(const std::string Fname, const int n,
	       const double bn_max, const double bn_scl);
  void ini_prm(void);
  void set_prm_dep(const int k) const;
  void clr_prm_dep(const int k) const;
  double set_dprm(void) const;
  void set_prm(void) const;
};


param_type   bn_prms;
int          n_iter, n_powell;
double       twoJ[2];
ss_vect<tps> Id_scl;

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

  n_prm = Fnum.size();

  bn_prms.bn_lim = dvector(1, n_prm); bn_prms.bn = dvector(1, n_prm);
  bn_prms.dbn = dvector(1, n_prm);

  printf("\nInitial bn; (incl. scaling) (%d):\n", n_prm);
  for (i = 1; i <= n_prm; i++) {
    bn_lim[i] = bn_max[i-1];
    if (n[i-1] > 0)
      // Multipole.
      bn[i] = get_bn(Fnum[i-1], 1, n[i-1])/bn_scl[i-1];
    else if (n[i-1] == -1)
      // Drift.
      bn[i] = get_L(Fnum[i-1], 1)/bn_scl[i-1];
    else if (n[i-1] == -2)
      // Location.
      bn[i] = get_bn_s(-Fnum[i-1], 1, n[i-1])/bn_scl[i-1];
    printf(" %12.5e", bn[i]);
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


void no_mpoles(const int n)
{
  int j;

  printf("\nzeroing multipoles: %d\n", n);
  for (j = 0; j < n_elem; j++)
    if (elem[j].kind == Mpole)
      set_bn(elem[j].Fnum, elem[j].Knum, n, 0e0);
}


void get_map_n(const int n)
{
  int          k;
  ss_vect<tps> map_n;

  get_Map(); map_n.identity();
  for (k = 0; k < n; k++)
    map_n = map_n*Map;
  Map = map_n;
}


void get_S(void)
{
  int    j;
  double S;

  S = 0e0;
  for (j = 0; j < n_elem; j++) {
    S += elem[j].L;
    elem[j].S = S; elem_tps[j].S = S;
  }
}


void get_nu_ksi(void)
{
  ss_vect<tps> nus;

  danot_(2);
  get_map_n(n_cell);
  danot_(3);
  K = MapNorm(Map, g, A1, A0, Map_res, 1); nus = dHdJ(K);

  printf("\nnu  = [%8.5f, %8.5f]\n",
	 nus[0].cst(), nus[1].cst());
  printf("ksi = [%8.5f, %8.5f]\n",
	 h_ijklm(nus[0], 0, 0, 0, 0, 1),
	 h_ijklm(nus[1], 0, 0, 0, 0, 1));
}


void prt_map(void)
{
  tps           h, h_re, h_im, K_re, K_im;
  ss_vect<tps>  nus;
  std::ofstream outf;

  danot_(NO-1);
  get_map_n(n_cell);

  file_wr(outf, "map.out");
  outf << Map[ct_];
  outf.close();
}


void prt_h_K(void)
{
  int           k;
  double        mu[2];
  tps           h_re, h_im, K_re, K_im;
  ss_vect<tps>  nus, R_half;
  std::ofstream outf;

  danot_(NO-1);
  get_map_n(n_cell);
  danot_(NO);
  K = MapNorm(Map, g, A1, A0, Map_res, no_tps); nus = dHdJ(K);
  CtoR(K, K_re, K_im);

#if 0
  // Compute h at symmetry point of lattice; for a lattice with two super
  // periods.
  printf("\nprt_h_K: nu/2 = [%7.5f, %7.5f]\n",
	 nus[0].cst()/2e0, nus[1].cst()/2e0);
  R_half.identity();
  mu[0] = 2e0*M_PI*nus[0].cst()/2e0; mu[1] = 2e0*M_PI*nus[1].cst()/2e0;
  for (k = 0; k < 2; k++) {
    R_half[2*k]   =  cos(mu[k])*tps(0e0, 2*k+1) + sin(mu[k])*tps(0e0, 2*k+2);
    R_half[2*k+1] = -sin(mu[k])*tps(0e0, 2*k+1) + cos(mu[k])*tps(0e0, 2*k+2);
  }
  CtoR(get_h()*Inv(R_half), h_re, h_im);
#else
  CtoR(get_h(), h_re, h_im);
#endif

  file_wr(outf, "h.out");
  outf << h_re*Id_scl << h_im*Id_scl;
  outf.close();

  file_wr(outf, "K.out");
  outf << K_re*Id_scl << K_im*Id_scl;
  outf.close();

  nus[3] = nus[3]*Id_scl; nus[4] = nus[4]*Id_scl;
  daeps_(1e-5);
  file_wr(outf, "nus.out");
  // Remove numeric noise.
  outf << 1e0*nus[3] << 1e0*nus[4];
  outf.close();

  daeps_(1e-30);
}


void prt_system(const int m, const int n_b2, double **A, double *b)
{
  int i, j, d, n_h;

  printf("\n Ax = b:\n");
  for (j = 1; j <= n_b2; j++)
    if (j == 1)
      printf("%11d", j);
    else
      printf("%11d", j);
  printf("\n");
  d = (symm)? 0 : 16;
  for (i = 1; i <= m; i++) {
    if (i-1 == 0)
      printf("1st order chromatic\n");
    else if (i-1 == 3)
      printf("1st order geometric\n");
    else if ((NO >= 5) && (i-1 == 3+5))
      printf("2nd order geometric\n");
    else if ((NO >= 6) && (i-1 == 3+5+8))
      printf("3rd order geometric\n");

    n_h = 0;
    if (NO >= 3+1) n_h += 3 + 5;   // 8.
    if (NO >= 4+1) n_h += 8;       // 16.
    if (NO >= 5+1) n_h += 14;      // 30.
    if (!symm)     n_h += 16 + 14;

    if (i-1 == n_h)
      printf("linear chromaticity\n");
    else if (i-1 == n_h+2)
      printf("ampl. dependant tune shift\n");
    else if (i-1 == n_h+2+3)
      printf("2nd order chromaticity\n");
    else if (i-1 == n_h+2+3+2)
      printf("cross terms\n");
    else if (i-1 == n_h+2+3+2+3)
      printf("3rd order chromaticity\n");
    else if (i-1 == n_h+2+3+2+3+2) {
      if (!tune_conf)
	printf("ampl. dependant tune shift\n");
      else
	printf("tune confinement\n");
    }

    printf("%4d", i);
    for (j = 1; j <= n_b2; j++)
      printf("%11.3e", A[i][j]);
    printf("%11.3e\n", b[i]);
  }
}


void prt_dnu(tps &K)
{
  tps          K_re, K_im, nu_scl[2];
  ss_vect<tps> nus;

  CtoR(K, K_re, K_im); nus = dHdJ(K);
  nu_scl[0] = nus[3]*Id_scl; nu_scl[1] = nus[4]*Id_scl;

  printf("\ndnu:\n");
  printf(" %8.5f",   h_ijklm(nu_scl[0], 1, 1, 0, 0, 0));
  printf(" %8.5f,",  h_ijklm(nu_scl[0], 0, 0, 1, 1, 0));

  printf(" %8.5f",   h_ijklm(nu_scl[0], 2, 2, 0, 0, 0));
  printf(" %8.5f",   h_ijklm(nu_scl[0], 1, 1, 1, 1, 0));
  printf(" %8.5f,",  h_ijklm(nu_scl[0], 0, 0, 2, 2, 0));

  printf(" %8.5f",   h_ijklm(nu_scl[0], 3, 3, 0, 0, 0));
  printf(" %8.5f",   h_ijklm(nu_scl[0], 2, 2, 1, 1, 0));
  printf(" %8.5f",   h_ijklm(nu_scl[0], 1, 1, 2, 2, 0));
  printf(" %8.5f\n", h_ijklm(nu_scl[0], 0, 0, 3, 3, 0));

  printf(" %8.5f",   h_ijklm(nu_scl[1], 1, 1, 0, 0, 0));
  printf(" %8.5f,",  h_ijklm(nu_scl[1], 0, 0, 1, 1, 0));

  printf(" %8.5f",   h_ijklm(nu_scl[1], 2, 2, 0, 0, 0));
  printf(" %8.5f",   h_ijklm(nu_scl[1], 1, 1, 1, 1, 0));
  printf(" %8.5f,",  h_ijklm(nu_scl[1], 0, 0, 2, 2, 0));

  printf(" %8.5f",   h_ijklm(nu_scl[1], 3, 3, 0, 0, 0));
  printf(" %8.5f",   h_ijklm(nu_scl[1], 2, 2, 1, 1, 0));
  printf(" %8.5f",   h_ijklm(nu_scl[1], 1, 1, 2, 2, 0));
  printf(" %8.5f\n", h_ijklm(nu_scl[1], 0, 0, 3, 3, 0));

  printf("\n %8.5f", h_ijklm(nu_scl[0], 0, 0, 0, 0, 2));
  printf(" %8.5f",   h_ijklm(nu_scl[0], 1, 1, 0, 0, 2));
  printf(" %8.5f\n", h_ijklm(nu_scl[0], 0, 0, 1, 1, 2));

  printf(" %8.5f",   h_ijklm(nu_scl[1], 0, 0, 0, 0, 2));
  printf(" %8.5f",   h_ijklm(nu_scl[1], 1, 1, 0, 0, 2));
  printf(" %8.5f\n", h_ijklm(nu_scl[1], 0, 0, 1, 1, 2));

  printf("\nTune confinement:\n");
  printf(" %11.3e %11.3e\n",
	 h_ijklm(K_re, 3, 3, 0, 0, 0),
	 h_ijklm(K_re/(3e0*twoJ[X_]), 2, 2, 0, 0, 0));
  printf(" %11.3e %11.3e\n",
	 h_ijklm(K_re, 0, 0, 3, 3, 0),
	 h_ijklm(K_re/(3e0*twoJ[Y_]), 0, 0, 2, 2, 0));
  printf(" %11.3e %11.3e %11.3e\n",
	 h_ijklm(K_re*Id_scl, 1, 1, 1, 1, 0),
	 h_ijklm(K_re*Id_scl, 2, 2, 1, 1, 0),
	 h_ijklm(K_re*Id_scl, 1, 1, 2, 2, 0));
}


void prt_bn_1(const param_type &bn_prms)
{
  int  k;
  FILE *outf;

  const std::string file_name = "dnu.out";

  outf = file_write(file_name.c_str());

  k = 0;
  fprintf(outf, "\ns1b: sextupole, l = 0.25, k = %12.5e, n = nsext"
	  ", Method = Meth;\n", bn_prms.bn_scl[k]*bn_prms.bn[k+1]);
  k++;
  fprintf(outf, "s1d: sextupole, l = 0.25, k = %12.5e, n = nsext"
	  ", Method = Meth;\n", bn_prms.bn_scl[k]*bn_prms.bn[k+1]);
  k++;
  fprintf(outf, "s2b: sextupole, l = 0.25, k = %12.5e, n = nsext"
	  ", Method = Meth;\n", bn_prms.bn_scl[k]*bn_prms.bn[k+1]);
  k++;
  fprintf(outf, "s2d: sextupole, l = 0.25, k = %12.5e, n = nsext"
	  ", Method = Meth;\n", bn_prms.bn_scl[k]*bn_prms.bn[k+1]);
  k++;
  fprintf(outf, "sx1: sextupole, l = 0.25, k = %12.5e, n = nsext"
	  ", Method = Meth;\n", bn_prms.bn_scl[k]*bn_prms.bn[k+1]);
  k++;
  fprintf(outf, "sy1: sextupole, l = 0.2,  k = %12.5e, n = nsext"
	  ", Method = Meth;\n", bn_prms.bn_scl[k]*bn_prms.bn[k+1]);
  k++;
  fprintf(outf, "s3:  sextupole, l = 0.0,  k = %12.5e, n = nsext"
	  ", Method = Meth;\n", bn_prms.bn_scl[k]*bn_prms.bn[k+1]);

  fclose(outf);
}


void prt_bn_3(const param_type &bn_prms)
{
  int  k;
  FILE *outf;

  const std::string file_name = "dnu.out";

  outf = file_write(file_name.c_str());

  k = 0;
  fprintf(outf, "\nts1a:  sextupole, l = 0.29,  k = %12.5e, n = nsext"
	  ", Method = Meth;\n", bn_prms.bn_scl[k]*bn_prms.bn[k+1]);
  k++;
  fprintf(outf, "ts1ab: sextupole, l = 0.29,  k = %12.5e, n = nsext"
	  ", Method = Meth;\n", bn_prms.bn_scl[k]*bn_prms.bn[k+1]);
  k++;
  fprintf(outf, "ts2a:  sextupole, l = 0.19,  k = %12.5e, n = nsext"
	  ", Method = Meth;\n", bn_prms.bn_scl[k]*bn_prms.bn[k+1]);
  k++;
  fprintf(outf, "ts2ab: sextupole, l = 0.19,  k = %12.5e, n = nsext"
	  ", Method = Meth;\n", bn_prms.bn_scl[k]*bn_prms.bn[k+1]);
  k++;
  fprintf(outf, "ts1b:  sextupole, l = 0.29,  k = %12.5e, n = nsext"
	  ", Method = Meth;\n", bn_prms.bn_scl[k]*bn_prms.bn[k+1]);
  k++;
  fprintf(outf, "ts2b:  sextupole, l = 0.29,  k = %12.5e, n = nsext"
	  ", Method = Meth;\n", bn_prms.bn_scl[k]*bn_prms.bn[k+1]);
  k++;
  fprintf(outf, "ts1c:  sextupole, l = 0.29,  k = %12.5e, n = nsext"
	  ", Method = Meth;\n", bn_prms.bn_scl[k]*bn_prms.bn[k+1]);
  k++;
  fprintf(outf, "ts2c:  sextupole, l = 0.29,  k = %12.5e, n = nsext"
	  ", Method = Meth;\n", bn_prms.bn_scl[k]*bn_prms.bn[k+1]);
  // k++;
  // fprintf(outf, "ts1d:  sextupole, l = 0.29,  k = %12.5e, n = nsext"
  // 	  ", Method = Meth;\n", bn_prms.bn_scl[k]*bn_prms.bn[k+1]);
  k++;
  fprintf(outf, "ts2d:  sextupole, l = 0.29,  k = %12.5e, n = nsext"
	  ", Method = Meth;\n", bn_prms.bn_scl[k]*bn_prms.bn[k+1]);
  k++;
  fprintf(outf, "ts1e:  sextupole, l = 0.29,  k = %12.5e, n = nsext"
	  ", Method = Meth;\n", bn_prms.bn_scl[k]*bn_prms.bn[k+1]);
  k++;
  fprintf(outf, "ts2e:  sextupole, l = 0.29,  k = %12.5e, n = nsext"
	  ", Method = Meth;\n", bn_prms.bn_scl[k]*bn_prms.bn[k+1]);

  if (lat_case == 4) {
    k++;
    fprintf(outf, "s1:    sextupole, l = 0.175, k = %12.5e, n = nsext"
	    ", Method = Meth;\n", bn_prms.bn_scl[k]*bn_prms.bn[k+1]);
    k++;
    fprintf(outf, "s2:    sextupole, l = 0.175, k = %12.5e, n = nsext"
	    ", Method = Meth;\n", bn_prms.bn_scl[k]*bn_prms.bn[k+1]);
    k++;
    fprintf(outf, "s3:    sextupole, l = 0.175, k = %12.5e, n = nsext"
	    ", Method = Meth;\n", bn_prms.bn_scl[k]*bn_prms.bn[k+1]);
    k++;
    fprintf(outf, "s4:    sextupole, l = 0.175, k = %12.5e, n = nsext"
	    ", Method = Meth;\n", bn_prms.bn_scl[k]*bn_prms.bn[k+1]);
    k++;
    fprintf(outf, "s5:    sextupole, l = 0.175, k = %12.5e, n = nsext"
	    ", Method = Meth;\n", bn_prms.bn_scl[k]*bn_prms.bn[k+1]);
  }

  fclose(outf);
}


void prt_bn_6(const param_type &bn_prms)
{
  int  k;
  FILE *outf;

  const std::string file_name = "dnu.out";

  outf = file_write(file_name.c_str());

  k = 0;
  fprintf(outf, "\nsd1:  sextupole, l = 0.14, k = %12.5e, n = nsext"
	  ", Method = Meth;\n", bn_prms.bn_scl[k]*bn_prms.bn[k+1]);
  k++;
  fprintf(outf, "sd2:  sextupole, l = 0.14, k = %12.5e, n = nsext"
	  ", Method = Meth;\n", bn_prms.bn_scl[k]*bn_prms.bn[k+1]);
  k++;
  fprintf(outf, "sd3:  sextupole, l = 0.14, k = %12.5e, n = nsext"
	  ", Method = Meth;\n", bn_prms.bn_scl[k]*bn_prms.bn[k+1]);
  k++;
  fprintf(outf, "sf21: sextupole, l = 0.07, k = %12.5e, n = nsext"
	  ", Method = Meth;\n", bn_prms.bn_scl[k]*bn_prms.bn[k+1]);
  k++;
  fprintf(outf, "sd31: sextupole, l = 0.14, k = %12.5e, n = nsext"
	  ", Method = Meth;\n", bn_prms.bn_scl[k]*bn_prms.bn[k+1]);
  k++;
  fprintf(outf, "sf1:  sextupole, l = 0.07, k = %12.5e, n = nsext"
	  ", Method = Meth;\n", bn_prms.bn_scl[k]*bn_prms.bn[k+1]);
  k++;
  fprintf(outf, "sh1a: sextupole, l = 0.14, k = %12.5e, n = nsext"
	  ", Method = Meth;\n", bn_prms.bn_scl[k]*bn_prms.bn[k+1]);
  k++;
  fprintf(outf, "sh1e: sextupole, l = 0.14, k = %12.5e, n = nsext"
	  ", Method = Meth;\n", bn_prms.bn_scl[k]*bn_prms.bn[k+1]);

  fclose(outf);
}


void prt_bn_7(const param_type &bn_prms)
{
  int  k;
  FILE *outf;

  const std::string file_name = "dnu.out";

  outf = file_write(file_name.c_str());

  k = 0;
  fprintf(outf, "\nsf:  sextupole, l = 0.14, k = %12.5e, n = nsext"
	  ", Method = Meth;\n", bn_prms.bn_scl[k]*bn_prms.bn[k+1]);
  k++;
  fprintf(outf, "sd:  sextupole, l = 0.14, k = %12.5e, n = nsext"
	  ", Method = Meth;\n", bn_prms.bn_scl[k]*bn_prms.bn[k+1]);

  k++;
  fprintf(outf, "\ns1:  sextupole, l = 0.0, k = %12.5e, n = nsext"
	  ", Method = Meth;\n", bn_prms.bn_scl[k]*bn_prms.bn[k+1]);
  k++;
  fprintf(outf, "s2: sextupole, l = 0.0, k = %12.5e, n = nsext"
	  ", Method = Meth;\n", bn_prms.bn_scl[k]*bn_prms.bn[k+1]);
  k++;
  fprintf(outf, "s3:  sextupole, l = 0.0, k = %12.5e, n = nsext"
	  ", Method = Meth;\n", bn_prms.bn_scl[k]*bn_prms.bn[k+1]);
  k++;
  fprintf(outf, "s4:  sextupole, l = 0.0, k = %12.5e, n = nsext"
	  ", Method = Meth;\n", bn_prms.bn_scl[k]*bn_prms.bn[k+1]);
  k++;
  fprintf(outf, "s5:  sextupole, l = 0.0, k = %12.5e, n = nsext"
	  ", Method = Meth;\n", bn_prms.bn_scl[k]*bn_prms.bn[k+1]);
  k++;
  fprintf(outf, "s6:  sextupole, l = 0.0, k = %12.5e, n = nsext"
	  ", Method = Meth;\n", bn_prms.bn_scl[k]*bn_prms.bn[k+1]);

  fclose(outf);
}


void prt_bn_8(const param_type &bn_prms)
{
  int  k;
  FILE *outf;

  const std::string file_name = "dnu.out";

  outf = file_write(file_name.c_str());

  k = 0;
  fprintf(outf, "\nsfh:  sextupole, l = 0.05, k = %12.5e, n = nsext"
	  ", Method = Meth;\n", bn_prms.bn_scl[k]*bn_prms.bn[k+1]);
  k++;
  fprintf(outf, "sdh:  sextupole, l = 0.05, k = %12.5e, n = nsext"
	  ", Method = Meth;\n", bn_prms.bn_scl[k]*bn_prms.bn[k+1]);
  k++;
  fprintf(outf, "sfmh: sextupole, l = 0.05, k = %12.5e, n = nsext"
	  ", Method = Meth;\n", bn_prms.bn_scl[k]*bn_prms.bn[k+1]);
  k++;
  fprintf(outf, "sdmh: sextupole, l = 0.05, k = %12.5e, n = nsext"
	  ", Method = Meth;\n", bn_prms.bn_scl[k]*bn_prms.bn[k+1]);

  k++;
  fprintf(outf, "\nsxxh: sextupole, l = 0.05, k = %12.5e, n = nsext"
	  ", Method = Meth;\n", bn_prms.bn_scl[k]*bn_prms.bn[k+1]);
  k++;
  fprintf(outf, "sxyh: sextupole, l = 0.05, k = %12.5e, n = nsext"
	  ", Method = Meth;\n", bn_prms.bn_scl[k]*bn_prms.bn[k+1]);
  k++;
  fprintf(outf, "syyh: sextupole, l = 0.05, k = %12.5e, n = nsext"
	  ", Method = Meth;\n", bn_prms.bn_scl[k]*bn_prms.bn[k+1]);

  fclose(outf);
}


void prt_bn_9(const param_type &bn_prms)
{
  int  k;
  FILE *outf;

  const std::string file_name = "dnu.out";

  outf = file_write(file_name.c_str());

  k = 0;
  fprintf(outf, "\ns1:  sextupole, l = 0.14, k = %12.5e, n = nsext"
	  ", Method = Meth;\n", bn_prms.bn_scl[k]*bn_prms.bn[k+1]);
  k++;
  fprintf(outf, "s2:  sextupole, l = 0.14, k = %12.5e, n = nsext"
	  ", Method = Meth;\n", bn_prms.bn_scl[k]*bn_prms.bn[k+1]);
  k++;
  fprintf(outf, "s3:  sextupole, l = 0.07, k = %12.5e, n = nsext"
	  ", Method = Meth;\n", bn_prms.bn_scl[k]*bn_prms.bn[k+1]);
  k++;
  fprintf(outf, "s4:  sextupole, l = 0.14, k = %12.5e, n = nsext"
	  ", Method = Meth;\n", bn_prms.bn_scl[k]*bn_prms.bn[k+1]);
  k++;
  fprintf(outf, "s5:  sextupole, l = 0.05, k = %12.5e, n = nsext"
	  ", Method = Meth;\n", bn_prms.bn_scl[k]*bn_prms.bn[k+1]);
  k++;
  fprintf(outf, "s6:  sextupole, l = 0.05, k = %12.5e, n = nsext"
	  ", Method = Meth;\n", bn_prms.bn_scl[k]*bn_prms.bn[k+1]);
  k++;
  fprintf(outf, "s7:  sextupole, l = 0.14, k = %12.5e, n = nsext"
	  ", Method = Meth;\n", bn_prms.bn_scl[k]*bn_prms.bn[k+1]);
  k++;
  fprintf(outf, "s8:  sextupole, l = 0.07, k = %12.5e, n = nsext"
	  ", Method = Meth;\n", bn_prms.bn_scl[k]*bn_prms.bn[k+1]);
  k++;
  fprintf(outf, "s9:  sextupole, l = 0.14, k = %12.5e, n = nsext"
	  ", Method = Meth;\n", bn_prms.bn_scl[k]*bn_prms.bn[k+1]);
  k++;
  fprintf(outf, "s10: sextupole, l = 0.14, k = %12.5e, n = nsext"
	  ", Method = Meth;\n", bn_prms.bn_scl[k]*bn_prms.bn[k+1]);

  fclose(outf);
}


void prt_bn(const param_type &bn_prms)
{

  switch (lat_case) {
  case 3 ... 4:
    prt_bn_3(bn_prms);
    break;
  case 6:
    prt_bn_6(bn_prms);
    break;
  case 7:
    prt_bn_7(bn_prms);
    break;
  case 8:
    prt_bn_8(bn_prms);
    break;
  case 9:
    prt_bn_9(bn_prms);
    break;
  }
}


void conj_grad(const int n_iter, double bn[], double dbn[],
	       double *g, double *h, double (*f)(double *))
{
  // Conjugate gradient method.
  int    n_bn, i;
  double fret, g2, gamma, dg2;

  n_bn = bn_prms.n_prm;

  if (n_iter == 1) {
    for (i = 1; i <= n_bn; i++) {
      g[i] = dbn[i]; h[i] = g[i];
    }
  } else {
    dg2 = g2 = 0e0;
    for (i = 1; i <= n_bn; i++) {
      g2 += sqr(g[i]); dg2 += (dbn[i]-g[i])*dbn[i];
    }
    if (g2 != 0e0) {
      gamma = dg2/g2;
      for (i = 1; i <= n_bn; i++) {
	g[i] = dbn[i]; dbn[i] = h[i] = g[i] + gamma*h[i];
      }
    } else {
      printf("g.g = 0\n");
      exit(0);
    }
  }

  printf("\n");
  d_linmin(bn, dbn, n_bn, &fret, f);
}


double get_a(const double scale, const tps &t,
	     const int i, const int j, const int k, const int l, const int m)
{
  return scale*(h_ijklm_p(t, i, j, k, l, m, 7));
}


double get_b(const double scale, const tps &t,
	     const int i, const int j, const int k, const int l, const int m)
{
  return scale*(h_ijklm(t, i, j, k, l, m));
}


void fit_ksi1(const double ksi_x, const double ksi_y)
{
  int          n_bn, i, m;
  double       **A, *b;
  ss_vect<tps> nus;

  const int m_max = 2;

  n_bn = bn_prms.n_prm;

  b = dvector(1, m_max); A = dmatrix(1, m_max, 1, n_bn);

  no_mpoles(Sext); no_mpoles(Oct);

  printf("\n");
  for (i = 1; i <= n_bn; i++) {
    bn_prms.set_prm_dep(i-1);

    danot_(3);
    get_map_n(n_cell);
    danot_(4);
    K = MapNorm(Map, g, A1, A0, Map_res, 1); nus = dHdJ(K);

    m = 0;
    A[++m][i] = get_a(1e0, nus[3], 0, 0, 0, 0, 1);
    A[++m][i] = get_a(1e0, nus[4], 0, 0, 0, 0, 1);

    bn_prms.clr_prm_dep(i-1);
  }

  m = 0;
  b[++m] = -(get_b(1e0, nus[3], 0, 0, 0, 0, 1)-ksi_x);
  b[++m] = -(get_b(1e0, nus[4], 0, 0, 0, 0, 1)-ksi_y);

  prt_system(m, n_bn, A, b);

  SVD_lim(m, n_bn, A, b, bn_prms.bn_lim, bn_prms.svd_cut, bn_prms.bn,
	  bn_prms.dbn);

  bn_prms.set_dprm();

  for (i = 1; i <= n_bn; i++)
    bn_prms.bn[i] = get_bn(bn_prms.Fnum[i-1], 1, bn_prms.n[i-1]);

  printf("\nfit ksi:\n");
  for (i = 1; i <= n_bn; i++) {
    printf(" %12.5e", bn_prms.bn[i]);
    if (i % n_prt == 0) printf("\n");
  }
  if (n_bn % n_prt != 0) printf("\n");

  prt_mfile("flat_file.fit");
  prt_bn(bn_prms);

  free_dvector(b, 1, m_max); free_dmatrix(A, 1, m_max, 1, n_bn);
}


double get_f(double *bns)
{
  int                 i;
  static double       chi2_ref = 1e30;
  double              chi2;
  tps                 h_re, h_im, K_re, K_im, K_re_scl;
  std::vector<double> b;

  const bool prt = false;

  n_powell++;

  // Do not change parameters.
  if (prt) printf("get_f (incl. scaling):\n");
  for (i = 1; i <= bn_prms.n_prm; i++) {
    set_bn(bn_prms.Fnum[i-1], bn_prms.n[i-1], bn_prms.bn_scl[i-1]*bns[i]);
    if (prt) printf(" %12.5e", bn_prms.bn_scl[i-1]*bns[i]);
  }
  if (prt) printf("\n");

  danot_(NO-1);
  get_map_n(n_cell);
  danot_(NO);
  K = MapNorm(Map, g, A1, A0, Map_res, 1);
  CtoR(K, K_re, K_im); CtoR(get_h(), h_re, h_im);
  K_re_scl = K_re*Id_scl; h_re = h_re*Id_scl; h_im = h_im*Id_scl;

  b.push_back(get_b(scl_h[0], h_re, 1, 0, 0, 0, 2));
  b.push_back(get_b(scl_h[0], h_re, 2, 0, 0, 0, 1));
  b.push_back(get_b(scl_h[0], h_re, 0, 0, 2, 0, 1));

  b.push_back(get_b(scl_h[0], h_re, 1, 0, 1, 1, 0));
  b.push_back(get_b(scl_h[0], h_re, 2, 1, 0, 0, 0));
  b.push_back(get_b(scl_h[0], h_re, 3, 0, 0, 0, 0));
  b.push_back(get_b(scl_h[0], h_re, 1, 0, 0, 2, 0));
  b.push_back(get_b(scl_h[0], h_re, 1, 0, 2, 0, 0));

  if (NO >= 5) {
    b.push_back(get_b(scl_h[1], h_re, 2, 0, 1, 1, 0));
    b.push_back(get_b(scl_h[1], h_re, 3, 1, 0, 0, 0));
    b.push_back(get_b(scl_h[1], h_re, 4, 0, 0, 0, 0));
    b.push_back(get_b(scl_h[1], h_re, 2, 0, 0, 2, 0));
    b.push_back(get_b(scl_h[1], h_re, 2, 0, 2, 0, 0));
    b.push_back(get_b(scl_h[1], h_re, 0, 0, 4, 0, 0));
    b.push_back(get_b(scl_h[1], h_re, 0, 0, 3, 1, 0));
    b.push_back(get_b(scl_h[1], h_re, 1, 1, 2, 0, 0));
  }

  if (NO >= 6) {
    b.push_back(get_b(scl_h[2], h_re, 5, 0, 0, 0, 0));
    b.push_back(get_b(scl_h[2], h_re, 4, 1, 0, 0, 0));
    b.push_back(get_b(scl_h[2], h_re, 3, 2, 0, 0, 0));
    b.push_back(get_b(scl_h[2], h_re, 3, 0, 2, 0, 0));
    b.push_back(get_b(scl_h[2], h_re, 2, 1, 2, 0, 0));
    b.push_back(get_b(scl_h[2], h_re, 2, 1, 0, 2, 0));
    b.push_back(get_b(scl_h[2], h_re, 3, 0, 0, 2, 0));
    b.push_back(get_b(scl_h[2], h_re, 1, 0, 4, 0, 0));
    b.push_back(get_b(scl_h[2], h_re, 1, 0, 0, 4, 0));
    b.push_back(get_b(scl_h[2], h_re, 3, 0, 1, 1, 0));
    b.push_back(get_b(scl_h[2], h_re, 2, 1, 1, 1, 0));
    b.push_back(get_b(scl_h[2], h_re, 1, 0, 3, 1, 0));
    b.push_back(get_b(scl_h[2], h_re, 1, 0, 1, 3, 0));
    b.push_back(get_b(scl_h[2], h_re, 1, 0, 2, 2, 0));
  }

  if (!symm) {
    b.push_back(get_b(scl_h[0], h_im, 1, 0, 0, 0, 2));
    b.push_back(get_b(scl_h[0], h_im, 2, 0, 0, 0, 1));
    b.push_back(get_b(scl_h[0], h_im, 0, 0, 2, 0, 1));

    b.push_back(get_b(scl_h[0], h_im, 1, 0, 1, 1, 0));
    b.push_back(get_b(scl_h[0], h_im, 2, 1, 0, 0, 0));
    b.push_back(get_b(scl_h[0], h_im, 3, 0, 0, 0, 0));
    b.push_back(get_b(scl_h[0], h_im, 1, 0, 0, 2, 0));
    b.push_back(get_b(scl_h[0], h_im, 1, 0, 2, 0, 0));

    b.push_back(get_b(scl_h[1], h_im, 2, 0, 1, 1, 0));
    b.push_back(get_b(scl_h[1], h_im, 3, 1, 0, 0, 0));
    b.push_back(get_b(scl_h[1], h_im, 4, 0, 0, 0, 0));
    b.push_back(get_b(scl_h[1], h_im, 2, 0, 0, 2, 0));
    b.push_back(get_b(scl_h[1], h_im, 2, 0, 2, 0, 0));
    b.push_back(get_b(scl_h[1], h_im, 0, 0, 4, 0, 0));
    b.push_back(get_b(scl_h[1], h_im, 0, 0, 3, 1, 0));
    b.push_back(get_b(scl_h[1], h_im, 1, 1, 2, 0, 0));

    if (NO >= 6) {
      b.push_back(get_b(scl_h[2], h_im, 5, 0, 0, 0, 0));
      b.push_back(get_b(scl_h[2], h_im, 4, 1, 0, 0, 0));
      b.push_back(get_b(scl_h[2], h_im, 3, 2, 0, 0, 0));
      b.push_back(get_b(scl_h[2], h_im, 3, 0, 2, 0, 0));
      b.push_back(get_b(scl_h[2], h_im, 2, 1, 2, 0, 0));
      b.push_back(get_b(scl_h[2], h_im, 2, 1, 0, 2, 0));
      b.push_back(get_b(scl_h[2], h_im, 3, 0, 0, 2, 0));
      b.push_back(get_b(scl_h[2], h_im, 1, 0, 4, 0, 0));
      b.push_back(get_b(scl_h[2], h_im, 1, 0, 0, 4, 0));
      b.push_back(get_b(scl_h[2], h_im, 3, 0, 1, 1, 0));
      b.push_back(get_b(scl_h[2], h_im, 2, 1, 1, 1, 0));
      b.push_back(get_b(scl_h[2], h_im, 1, 0, 3, 1, 0));
      b.push_back(get_b(scl_h[2], h_im, 1, 0, 1, 3, 0));
      b.push_back(get_b(scl_h[2], h_im, 1, 0, 2, 2, 0));
    }
  }

  b.push_back(get_b(scl_ksi[0], K_re_scl, 1, 1, 0, 0, 1));
  b.push_back(get_b(scl_ksi[0], K_re_scl, 0, 0, 1, 1, 1));

  if (NO >= 5) {
    b.push_back(get_b(scl_dnu[0], K_re_scl, 2, 2, 0, 0, 0));
    b.push_back(get_b(scl_dnu[0], K_re_scl, 0, 0, 2, 2, 0));
    b.push_back(get_b(scl_dnu[0], K_re_scl, 1, 1, 1, 1, 0));

    b.push_back(get_b(scl_ksi[1], K_re_scl, 1, 1, 0, 0, 2));
    b.push_back(get_b(scl_ksi[1], K_re_scl, 0, 0, 1, 1, 2));
  }

  if (NO >= 6) {
    b.push_back(get_b(scl_dnu[0], K_re_scl, 2, 2, 0, 0, 1));
    b.push_back(get_b(scl_dnu[0], K_re_scl, 0, 0, 2, 2, 1));
    b.push_back(get_b(scl_dnu[0], K_re_scl, 1, 1, 1, 1, 1));

    b.push_back(get_b(scl_ksi[2], K_re_scl, 1, 1, 0, 0, 3));
    b.push_back(get_b(scl_ksi[2], K_re_scl, 0, 0, 1, 1, 3));
  }

  if (NO >= 7) {
    if (!tune_conf) {
      b.push_back(get_b(scl_dnu[0], K_re_scl, 3, 3, 0, 0, 0));
      b.push_back(get_b(scl_dnu[0], K_re_scl, 2, 2, 1, 1, 0));
      b.push_back(get_b(scl_dnu[0], K_re_scl, 1, 1, 2, 2, 0));
      b.push_back(get_b(scl_dnu[0], K_re_scl, 0, 0, 3, 3, 0));
    } else {
      b.push_back(get_b(scl_dnu[0], K_re_scl, 2, 2, 1, 1, 0));
      b.push_back(get_b(scl_dnu[0], K_re_scl, 1, 1, 2, 2, 0));
      b.push_back(scl_dnu[1]
		  *(get_b(1e0, K_re, 3, 3, 0, 0, 0)
		    +get_b(1e0/(3e0*twoJ[X_]), K_re, 2, 2, 0, 0, 0)));
      b.push_back(scl_dnu[1]
		  *(get_b(1e0, K_re, 0, 0, 3, 3, 0)
		    +get_b(1e0/(3e0*twoJ[Y_]), K_re, 0, 0, 2, 2, 0)));
    }
  }

  if (NO >= 9) {
    b.push_back(get_b(scl_dnu[0], K_re_scl, 4, 4, 0, 0, 0));
    b.push_back(get_b(scl_dnu[0], K_re_scl, 3, 3, 1, 1, 0));
    b.push_back(get_b(scl_dnu[0], K_re_scl, 2, 2, 2, 2, 0));
    b.push_back(get_b(scl_dnu[0], K_re_scl, 1, 1, 3, 3, 0));
    b.push_back(get_b(scl_dnu[0], K_re_scl, 0, 0, 4, 4, 0));
  }

  chi2 = 0e0;
  for (i = 0; i < (int)b.size(); i++)
    chi2 += sqr(b[i]);

  if (chi2 < chi2_ref) {
    prt_bn(bn_prms);

    printf("\n%3d %12.5e -> %12.5e\n", n_powell, chi2_ref, chi2);
    printf("b & bn:\n");
    for (i = 0; i < (int)b.size(); i++)
      printf("%11.3e", b[i]);
    printf("\n");
    for (i = 1; i <= bn_prms.n_prm; i++) 
      printf("%11.3e", bns[i]);
    printf("\n");
  }

  chi2_ref = min(chi2, chi2_ref);

  return chi2;
}


void get_f_grad(const int n_bn, double *f, double **A, double &chi2, int &m)
{
  int i, j;
  tps h_re, h_im, K_re, K_im, K_re_scl;

  // printf("\n");
  for (i = 1; i <= n_bn; i++) {
    bn_prms.set_prm_dep(i-1);

    danot_(NO-1);
    get_map_n(n_cell);
    danot_(NO);
    K = MapNorm(Map, g, A1, A0, Map_res, 1);
    CtoR(K, K_re, K_im); CtoR(get_h(), h_re, h_im);
    K_re_scl = K_re*Id_scl; h_re = h_re*Id_scl; h_im = h_im*Id_scl;

    m = 0;
    A[++m][i] = get_a(scl_h[0], h_re, 1, 0, 0, 0, 2);
    A[++m][i] = get_a(scl_h[0], h_re, 2, 0, 0, 0, 1);
    A[++m][i] = get_a(scl_h[0], h_re, 0, 0, 2, 0, 1);

    A[++m][i] = get_a(scl_h[0], h_re, 1, 0, 1, 1, 0);
    A[++m][i] = get_a(scl_h[0], h_re, 2, 1, 0, 0, 0);
    A[++m][i] = get_a(scl_h[0], h_re, 3, 0, 0, 0, 0);
    A[++m][i] = get_a(scl_h[0], h_re, 1, 0, 0, 2, 0);
    A[++m][i] = get_a(scl_h[0], h_re, 1, 0, 2, 0, 0);

    if (NO >= 5) {
      A[++m][i] = get_a(scl_h[1], h_re, 2, 0, 1, 1, 0);
      A[++m][i] = get_a(scl_h[1], h_re, 3, 1, 0, 0, 0);
      A[++m][i] = get_a(scl_h[1], h_re, 4, 0, 0, 0, 0);
      A[++m][i] = get_a(scl_h[1], h_re, 2, 0, 0, 2, 0);
      A[++m][i] = get_a(scl_h[1], h_re, 2, 0, 2, 0, 0);
      A[++m][i] = get_a(scl_h[1], h_re, 0, 0, 4, 0, 0);
      A[++m][i] = get_a(scl_h[1], h_re, 0, 0, 3, 1, 0);
      A[++m][i] = get_a(scl_h[1], h_re, 1, 1, 2, 0, 0);
    }

    if (NO >= 6) {
      A[++m][i] = get_a(scl_h[2], h_re, 5, 0, 0, 0, 0);
      A[++m][i] = get_a(scl_h[2], h_re, 4, 1, 0, 0, 0);
      A[++m][i] = get_a(scl_h[2], h_re, 3, 2, 0, 0, 0);
      A[++m][i] = get_a(scl_h[2], h_re, 3, 0, 2, 0, 0);
      A[++m][i] = get_a(scl_h[2], h_re, 2, 1, 2, 0, 0);
      A[++m][i] = get_a(scl_h[2], h_re, 2, 1, 0, 2, 0);
      A[++m][i] = get_a(scl_h[2], h_re, 3, 0, 0, 2, 0);
      A[++m][i] = get_a(scl_h[2], h_re, 1, 0, 4, 0, 0);
      A[++m][i] = get_a(scl_h[2], h_re, 1, 0, 0, 4, 0);
      A[++m][i] = get_a(scl_h[2], h_re, 3, 0, 1, 1, 0);
      A[++m][i] = get_a(scl_h[2], h_re, 2, 1, 1, 1, 0);
      A[++m][i] = get_a(scl_h[2], h_re, 1, 0, 3, 1, 0);
      A[++m][i] = get_a(scl_h[2], h_re, 1, 0, 1, 3, 0);
      A[++m][i] = get_a(scl_h[2], h_re, 1, 0, 2, 2, 0);
    }

    if (!symm) {
      A[++m][i] = get_a(scl_h[0], h_im, 1, 0, 0, 0, 2);
      A[++m][i] = get_a(scl_h[0], h_im, 2, 0, 0, 0, 1);
      A[++m][i] = get_a(scl_h[0], h_im, 0, 0, 2, 0, 1);

      A[++m][i] = get_a(scl_h[0], h_im, 1, 0, 1, 1, 0);
      A[++m][i] = get_a(scl_h[0], h_im, 2, 1, 0, 0, 0);
      A[++m][i] = get_a(scl_h[0], h_im, 3, 0, 0, 0, 0);
      A[++m][i] = get_a(scl_h[0], h_im, 1, 0, 0, 2, 0);
      A[++m][i] = get_a(scl_h[0], h_im, 1, 0, 2, 0, 0);

      A[++m][i] = get_a(scl_h[1], h_im, 2, 0, 1, 1, 0);
      A[++m][i] = get_a(scl_h[1], h_im, 3, 1, 0, 0, 0);
      A[++m][i] = get_a(scl_h[1], h_im, 4, 0, 0, 0, 0);
      A[++m][i] = get_a(scl_h[1], h_im, 2, 0, 0, 2, 0);
      A[++m][i] = get_a(scl_h[1], h_im, 2, 0, 2, 0, 0);
      A[++m][i] = get_a(scl_h[1], h_im, 0, 0, 4, 0, 0);
      A[++m][i] = get_a(scl_h[1], h_im, 0, 0, 3, 1, 0);
      A[++m][i] = get_a(scl_h[1], h_im, 1, 1, 2, 0, 0);

      if (NO >= 6) {
	A[++m][i] = get_a(scl_h[2], h_im, 5, 0, 0, 0, 0);
	A[++m][i] = get_a(scl_h[2], h_im, 4, 1, 0, 0, 0);
	A[++m][i] = get_a(scl_h[2], h_im, 3, 2, 0, 0, 0);
	A[++m][i] = get_a(scl_h[2], h_im, 3, 0, 2, 0, 0);
	A[++m][i] = get_a(scl_h[2], h_im, 2, 1, 2, 0, 0);
	A[++m][i] = get_a(scl_h[2], h_im, 2, 1, 0, 2, 0);
	A[++m][i] = get_a(scl_h[2], h_im, 3, 0, 0, 2, 0);
	A[++m][i] = get_a(scl_h[2], h_im, 1, 0, 4, 0, 0);
	A[++m][i] = get_a(scl_h[2], h_im, 1, 0, 0, 4, 0);
	A[++m][i] = get_a(scl_h[2], h_im, 3, 0, 1, 1, 0);
	A[++m][i] = get_a(scl_h[2], h_im, 2, 1, 1, 1, 0);
	A[++m][i] = get_a(scl_h[2], h_im, 1, 0, 3, 1, 0);
	A[++m][i] = get_a(scl_h[2], h_im, 1, 0, 1, 3, 0);
	A[++m][i] = get_a(scl_h[2], h_im, 1, 0, 2, 2, 0);
      }
    }

    A[++m][i] = get_a(scl_ksi[0], K_re_scl, 1, 1, 0, 0, 1);
    A[++m][i] = get_a(scl_ksi[0], K_re_scl, 0, 0, 1, 1, 1);

    if (NO >= 5) {
      A[++m][i] = get_a(scl_dnu[0], K_re_scl, 2, 2, 0, 0, 0);
      A[++m][i] = get_a(scl_dnu[0], K_re_scl, 0, 0, 2, 2, 0);
      A[++m][i] = get_a(scl_dnu[0], K_re_scl, 1, 1, 1, 1, 0);

      A[++m][i] = get_a(scl_ksi[1], K_re_scl, 1, 1, 0, 0, 2);
      A[++m][i] = get_a(scl_ksi[1], K_re_scl, 0, 0, 1, 1, 2);
    }

    if (NO >= 6) {
      A[++m][i] = get_a(scl_dnu[0], K_re_scl, 2, 2, 0, 0, 1);
      A[++m][i] = get_a(scl_dnu[0], K_re_scl, 0, 0, 2, 2, 1);
      A[++m][i] = get_a(scl_dnu[0], K_re_scl, 1, 1, 1, 1, 1);

      A[++m][i] = get_a(scl_ksi[2], K_re_scl, 1, 1, 0, 0, 3);
      A[++m][i] = get_a(scl_ksi[2], K_re_scl, 0, 0, 1, 1, 3);
    }

    if (NO >= 7) {
      if (!tune_conf) {
	A[++m][i] = get_a(scl_dnu[0], K_re_scl, 3, 3, 0, 0, 0);
	A[++m][i] = get_a(scl_dnu[0], K_re_scl, 2, 2, 1, 1, 0);
	A[++m][i] = get_a(scl_dnu[0], K_re_scl, 1, 1, 2, 2, 0);
	A[++m][i] = get_a(scl_dnu[0], K_re_scl, 0, 0, 3, 3, 0);
      } else {
	A[++m][i] = get_a(scl_dnu[0], K_re_scl, 2, 2, 1, 1, 0);
	A[++m][i] = get_a(scl_dnu[0], K_re_scl, 1, 1, 2, 2, 0);
	A[++m][i] =
	  scl_dnu[1]
	  *(get_a(1e0, K_re, 3, 3, 0, 0, 0)
	    +get_a(1e0/(3e0*twoJ[X_]), K_re, 2, 2, 0, 0, 0));
	A[++m][i] =
	  scl_dnu[1]
	  *(get_a(1e0, K_re, 0, 0, 3, 3, 0)
	    +get_a(1e0/(3e0*twoJ[Y_]), K_re, 0, 0, 2, 2, 0));
      }
    }

    if (NO >= 9) {
      A[++m][i] = get_a(scl_dnu[0], K_re_scl, 4, 4, 0, 0, 0);
      A[++m][i] = get_a(scl_dnu[0], K_re_scl, 3, 3, 1, 1, 0);
      A[++m][i] = get_a(scl_dnu[0], K_re_scl, 2, 2, 2, 2, 0);
      A[++m][i] = get_a(scl_dnu[0], K_re_scl, 1, 1, 3, 3, 0);
      A[++m][i] = get_a(scl_dnu[0], K_re_scl, 0, 0, 4, 4, 0);
    }

    for (j = 1; j <= m; j++)
      A[j][i] *= bn_prms.bn_scl[i-1];

    bn_prms.clr_prm_dep(i-1);
  }

  m = 0;
  f[++m] = get_b(scl_h[0], h_re, 1, 0, 0, 0, 2);
  f[++m] = get_b(scl_h[0], h_re, 2, 0, 0, 0, 1);
  f[++m] = get_b(scl_h[0], h_re, 0, 0, 2, 0, 1);

  f[++m] = get_b(scl_h[0], h_re, 1, 0, 1, 1, 0);
  f[++m] = get_b(scl_h[0], h_re, 2, 1, 0, 0, 0);
  f[++m] = get_b(scl_h[0], h_re, 3, 0, 0, 0, 0);
  f[++m] = get_b(scl_h[0], h_re, 1, 0, 0, 2, 0);
  f[++m] = get_b(scl_h[0], h_re, 1, 0, 2, 0, 0);

  if (NO >= 5) {
    f[++m] = get_b(scl_h[1], h_re, 2, 0, 1, 1, 0);
    f[++m] = get_b(scl_h[1], h_re, 3, 1, 0, 0, 0);
    f[++m] = get_b(scl_h[1], h_re, 4, 0, 0, 0, 0);
    f[++m] = get_b(scl_h[1], h_re, 2, 0, 0, 2, 0);
    f[++m] = get_b(scl_h[1], h_re, 2, 0, 2, 0, 0);
    f[++m] = get_b(scl_h[1], h_re, 0, 0, 4, 0, 0);
    f[++m] = get_b(scl_h[1], h_re, 0, 0, 3, 1, 0);
    f[++m] = get_b(scl_h[1], h_re, 1, 1, 2, 0, 0);
  }

  if (NO >= 6) {
    f[++m] = get_b(scl_h[2], h_re, 5, 0, 0, 0, 0);
    f[++m] = get_b(scl_h[2], h_re, 4, 1, 0, 0, 0);
    f[++m] = get_b(scl_h[2], h_re, 3, 2, 0, 0, 0);
    f[++m] = get_b(scl_h[2], h_re, 3, 0, 2, 0, 0);
    f[++m] = get_b(scl_h[2], h_re, 2, 1, 2, 0, 0);
    f[++m] = get_b(scl_h[2], h_re, 2, 1, 0, 2, 0);
    f[++m] = get_b(scl_h[2], h_re, 3, 0, 0, 2, 0);
    f[++m] = get_b(scl_h[2], h_re, 1, 0, 4, 0, 0);
    f[++m] = get_b(scl_h[2], h_re, 1, 0, 0, 4, 0);
    f[++m] = get_b(scl_h[2], h_re, 3, 0, 1, 1, 0);
    f[++m] = get_b(scl_h[2], h_re, 2, 1, 1, 1, 0);
    f[++m] = get_b(scl_h[2], h_re, 1, 0, 3, 1, 0);
    f[++m] = get_b(scl_h[2], h_re, 1, 0, 1, 3, 0);
    f[++m] = get_b(scl_h[2], h_re, 1, 0, 2, 2, 0);
  }

  if (!symm) {
    f[++m] = get_b(scl_h[0], h_im, 1, 0, 0, 0, 2);
    f[++m] = get_b(scl_h[0], h_im, 2, 0, 0, 0, 1);
    f[++m] = get_b(scl_h[0], h_im, 0, 0, 2, 0, 1);

    f[++m] = get_b(scl_h[0], h_im, 1, 0, 1, 1, 0);
    f[++m] = get_b(scl_h[0], h_im, 2, 1, 0, 0, 0);
    f[++m] = get_b(scl_h[0], h_im, 3, 0, 0, 0, 0);
    f[++m] = get_b(scl_h[0], h_im, 1, 0, 0, 2, 0);
    f[++m] = get_b(scl_h[0], h_im, 1, 0, 2, 0, 0);

    f[++m] = get_b(scl_h[1], h_im, 2, 0, 1, 1, 0);
    f[++m] = get_b(scl_h[1], h_im, 3, 1, 0, 0, 0);
    f[++m] = get_b(scl_h[1], h_im, 4, 0, 0, 0, 0);
    f[++m] = get_b(scl_h[1], h_im, 2, 0, 0, 2, 0);
    f[++m] = get_b(scl_h[1], h_im, 2, 0, 2, 0, 0);
    f[++m] = get_b(scl_h[1], h_im, 0, 0, 4, 0, 0);
    f[++m] = get_b(scl_h[1], h_im, 0, 0, 3, 1, 0);
    f[++m] = get_b(scl_h[1], h_im, 1, 1, 2, 0, 0);

    if (NO >= 6) {
      f[++m] = get_b(scl_h[2], h_im, 5, 0, 0, 0, 0);
      f[++m] = get_b(scl_h[2], h_im, 4, 1, 0, 0, 0);
      f[++m] = get_b(scl_h[2], h_im, 3, 2, 0, 0, 0);
      f[++m] = get_b(scl_h[2], h_im, 3, 0, 2, 0, 0);
      f[++m] = get_b(scl_h[2], h_im, 2, 1, 2, 0, 0);
      f[++m] = get_b(scl_h[2], h_im, 2, 1, 0, 2, 0);
      f[++m] = get_b(scl_h[2], h_im, 3, 0, 0, 2, 0);
      f[++m] = get_b(scl_h[2], h_im, 1, 0, 4, 0, 0);
      f[++m] = get_b(scl_h[2], h_im, 1, 0, 0, 4, 0);
      f[++m] = get_b(scl_h[2], h_im, 3, 0, 1, 1, 0);
      f[++m] = get_b(scl_h[2], h_im, 2, 1, 1, 1, 0);
      f[++m] = get_b(scl_h[2], h_im, 1, 0, 3, 1, 0);
      f[++m] = get_b(scl_h[2], h_im, 1, 0, 1, 3, 0);
      f[++m] = get_b(scl_h[2], h_im, 1, 0, 2, 2, 0);
    }
  }

  f[++m] = get_b(scl_ksi[0], K_re_scl, 1, 1, 0, 0, 1);
  f[++m] = get_b(scl_ksi[0], K_re_scl, 0, 0, 1, 1, 1);

  if (NO >= 5) {
    f[++m] = get_b(scl_dnu[0], K_re_scl, 2, 2, 0, 0, 0);
    f[++m] = get_b(scl_dnu[0], K_re_scl, 0, 0, 2, 2, 0);
    f[++m] = get_b(scl_dnu[0], K_re_scl, 1, 1, 1, 1, 0);

    f[++m] = get_b(scl_ksi[1], K_re_scl, 1, 1, 0, 0, 2);
    f[++m] = get_b(scl_ksi[1], K_re_scl, 0, 0, 1, 1, 2);
  }

  if (NO >= 6) {
    f[++m] = get_b(scl_dnu[0], K_re_scl, 2, 2, 0, 0, 1);
    f[++m] = get_b(scl_dnu[0], K_re_scl, 0, 0, 2, 2, 1);
    f[++m] = get_b(scl_dnu[0], K_re_scl, 1, 1, 1, 1, 1);

    f[++m] = get_b(scl_ksi[2], K_re_scl, 1, 1, 0, 0, 3);
    f[++m] = get_b(scl_ksi[2], K_re_scl, 0, 0, 1, 1, 3);
  }

  if (NO >= 7) {
    if (!tune_conf) {
      f[++m] = get_b(scl_dnu[0], K_re_scl, 3, 3, 0, 0, 0);
      f[++m] = get_b(scl_dnu[0], K_re_scl, 2, 2, 1, 1, 0);
      f[++m] = get_b(scl_dnu[0], K_re_scl, 1, 1, 2, 2, 0);
      f[++m] = get_b(scl_dnu[0], K_re_scl, 0, 0, 3, 3, 0);
    } else {
      f[++m] = get_b(scl_dnu[0], K_re_scl, 2, 2, 1, 1, 0);
      f[++m] = get_b(scl_dnu[0], K_re_scl, 1, 1, 2, 2, 0);
      f[++m] =
	scl_dnu[1]
	*(get_b(1e0, K_re, 3, 3, 0, 0, 0)
	  +get_b(1e0/(3e0*twoJ[X_]), K_re, 2, 2, 0, 0, 0));
      f[++m] =
	scl_dnu[1]
	*(get_b(1e0, K_re, 0, 0, 3, 3, 0)
	  +get_b(1e0/(3e0*twoJ[Y_]), K_re, 0, 0, 2, 2, 0));
    }
  }

  if (NO >= 9) {
    f[++m] = get_b(scl_dnu[0], K_re_scl, 4, 4, 0, 0, 0);
    f[++m] = get_b(scl_dnu[0], K_re_scl, 3, 3, 1, 1, 0);
    f[++m] = get_b(scl_dnu[0], K_re_scl, 2, 2, 2, 2, 0);
    f[++m] = get_b(scl_dnu[0], K_re_scl, 1, 1, 3, 3, 0);
    f[++m] = get_b(scl_dnu[0], K_re_scl, 0, 0, 4, 4, 0);
  }

  chi2 = 0e0;
  for (j = 1; j <= m; j++)
    chi2 += sqr(f[j]);
}


void min_conj_grad(double &chi2, double &dbn_max, double *g_, double *h_,
		   const bool cg_meth)
{
  int    n_bn, i, m;
  double chi2_ref, **A, *f, *b, *bn_ref;

  const int m_max = 100;

  n_bn = bn_prms.n_prm;

  bn_ref = dvector(1, n_bn); f = dvector(1, m_max); b = dvector(1, m_max);
  A = dmatrix(1, m_max, 1, n_bn);

  chi2_ref = chi2;

  get_f_grad(n_bn, f, A, chi2, m);
  prt_system(m, n_bn, A, f);

  printf("\n%4d chi2: %12.5e -> %12.5e\n", n_iter, chi2_ref, chi2);

  for (i = 1; i <= m; i++)
    b[i] = -f[i];

  SVD_lim(m, n_bn, A, b, bn_prms.bn_lim, bn_prms.svd_cut, bn_prms.bn,
	  bn_prms.dbn);

  dvcopy(bn_prms.bn, n_bn, bn_ref);
  if (cg_meth)
    conj_grad(n_iter, bn_prms.bn, bn_prms.dbn, g_, h_, get_f);
  else
    bn_prms.set_dprm();

  for (i = 1; i <= n_bn; i++)
    bn_prms.bn[i] =
      get_bn(bn_prms.Fnum[i-1], 1, bn_prms.n[i-1])/bn_prms.bn_scl[i-1];

  printf("\nbn & dbn (incl. scaling):\n");
  for (i = 1; i <= n_bn; i++)
    printf(" %12.5e", bn_prms.bn[i]);
  printf("\n");
  dbn_max = 0e0;
  for (i = 1; i <= n_bn; i++) {
    dbn_max = max(fabs((bn_prms.bn[i]-bn_ref[i])), dbn_max);
    printf(" %12.5e", bn_prms.bn[i]-bn_ref[i]);
  }
  printf("\n");

  free_dvector(bn_ref, 1, n_bn); free_dvector(f, 1, m_max);
  free_dvector(b, 1, m_max); free_dmatrix(A, 1, m_max, 1, n_bn);
}


void min_conj_grad(const bool cg_meth)
{
  // Control tune foot print; conjugate gradient method.
  std::string str;
  int         n_bn;
  double      *g, *h, dbn_max, chi2;

  const int n_iter_max = 1000;

  n_bn = bn_prms.n_prm;

  g = dvector(1, n_bn); h = dvector(1, n_bn);

  n_iter = 0; chi2 = 0e0;
  do {
    n_iter++;

    min_conj_grad(chi2, dbn_max, g, h, cg_meth);

    prt_mfile("flat_file.fit");
    prt_bn(bn_prms);
  } while ((dbn_max >  bn_prms.bn_tol) && (n_iter < n_iter_max));
}


void min_powell(void)
{
  int    n_bn, i, j, iter;
  double **xi, fret;

  n_bn = bn_prms.n_prm;

  xi = dmatrix(1, n_bn, 1, n_bn);

  // Set initial directions (unit vectors).
  for (i = 1; i <= n_bn; i++)
    for (j = 1; j <= n_bn; j++)
      xi[i][j] = (i == j)? 1e0 : 0e0;

  dpowell(bn_prms.bn, xi, n_bn, bn_prms.bn_tol, &iter, &fret, get_f);

  free_dmatrix(xi, 1, n_bn, 1, n_bn);
}


void prt_lev_marq(const int m, const int n)
{
  int i;

  prt_system(m, n, A_lm, f_lm);
  prt_dnu(K);

  prt_bn(bn_prms);

  printf("\n%d bn:\n", n_powell);
  for (i = 1; i <= bn_prms.n_prm; i++) {
    bn_prms.bn[i] = get_bn(bn_prms.Fnum[i-1], 1, bn_prms.n[i-1]);
    printf("%11.3e", bn_prms.bn[i]);
    if (i % n_prt == 0) printf("\n");
  }
  if (bn_prms.n_prm % n_prt != 0) printf("\n");
  prt_mfile("flat_file.fit");
}


void get_f_der(double x, double *bn, double *yfit, double *dyda, int n)
{
  int        i, m1;
  static int m;
  double     chi2;

  const bool prt = false;

  m1 = (int)(x+0.5);
  if (prt) printf(" %d", m1);

  if (m1 == 1) {
    n_powell++;
    for (i = 1; i <= n; i++)
      set_bn(bn_prms.Fnum[i-1], bn_prms.n[i-1], bn_prms.bn_scl[i-1]*bn[i]);
    get_f_grad(n, f_lm, A_lm, chi2, m);
  }

  if (prt && (m1 == m)) {
    printf("\n");
    for (i = 1; i <= n; i++) {
      printf(" %12.5e", bn_prms.bn_scl[i-1]*bn[i]);
      if (i % n_prt == 0) printf("\n");
    }
    if (n % n_prt != 0) printf("\n");
  }

  *yfit = f_lm[m1];

  for (i = 1; i <= n; i++)
    dyda[i] = A_lm[m1][i];
}


void min_lev_marq(void)
{
  int    n_data, n_bn, i, n, *ia;
  double *x, *y, *sigma, **covar, **alpha, chisq, alambda, alambda0;

  n_data = 0;
  if (NO >= 3+1) n_data += 3 + 5 + 2;  // 10.
  if (NO >= 4+1) n_data += 8 + 3 + 2;  // 23.
  if (NO >= 5+1) n_data += 14 + 3 + 2; // 42.
  if (NO >= 6+1) n_data += 4;          // 46.
  if (NO >= 8+1) n_data += 5;          // 51.
  if (!symm)     n_data += 16 + 14;

  n_bn = bn_prms.n_prm;

  ia = ivector(1, n_bn);
  x = dvector(1, n_data); y = dvector(1, n_data); sigma = dvector(1, n_data);
  covar = dmatrix(1, n_bn, 1, n_bn); alpha = dmatrix(1, n_bn, 1, n_bn);
  f_lm = dvector(1, n_data); A_lm = dmatrix(1, n_data, 1, n_bn);

  get_f_grad(n_bn, f_lm, A_lm, chi2, n_data);
  prt_system(n_data, n_bn, A_lm, f_lm);

  for (i = 1; i <= n_bn; i++)
    ia[i] = 1;

  for (i = 1; i <= n_data; i++) {
    sigma[i] = 1e0; x[i] = i; y[i] = 0e0;
  }

  alambda = -1e0; alambda0 = 1e-3;
  dmrqmin(x, y, sigma, n_data, bn_prms.bn, ia, n_bn, covar, alpha, &chisq,
	  get_f_der, &alambda);
  printf("\nalambda = %7.1e, chi2 = %9.3e\n", alambda, chisq);
  if (alambda < alambda0) prt_lev_marq(n_data, n_bn);
  alambda0 = alambda;

  n = 0;
  do {
    n++;
    dmrqmin(x, y, sigma, n_data, bn_prms.bn, ia, n_bn,  covar, alpha, &chisq,
	    get_f_der, &alambda);
    printf("\nalambda = %7.1e, chi2 = %9.3e\n", alambda, chisq);
    if (alambda < alambda0) prt_lev_marq(n_data, n_bn);
    alambda0 = alambda;
  } while (n < 25);

  alambda = 0e0;
  dmrqmin(x, y, sigma, n_data, bn_prms.bn, ia, n_bn,  covar, alpha, &chisq,
	  get_f_der, &alambda);

  free_dvector(f_lm, 1, n_data); free_dmatrix(A_lm, 1, n_data, 1, n_bn);
  free_ivector(ia, 1, n_bn);
  free_dvector(x, 1, n_data); free_dvector(y, 1, n_data);
  free_dvector(sigma, 1, n_data);
  free_dmatrix(covar, 1, n_bn, 1, n_bn); free_dmatrix(alpha, 1, n_bn, 1, n_bn);
}


void prt_ct(const int n, const double delta)
{
  int             k;
  double          delta1;
  ss_vect<double> ps;
  FILE            *outf;

  outf = file_write("ct.out");

  // cavity_on = true;

  for (k = -n; k <= n; k++) {
    delta1 = (double)k/(double)n*delta;
    ps.zero();
    ps[x_] = 0*2.6e-3; ps[delta_] = delta1;
    ps.propagate(1, n_elem);

    fprintf(outf, "%3d %12.5e %12.5e\n", k, delta1, ps[ct_]);
  }

  fclose(outf);
}


void fit_tune(const double nu_x, const double nu_y,
	      const std::vector<int> &b2_Fam,
	      const double eps, const bool prt)
{
  // Periodic solution: [nu_x, nu_y, beta_x, beta_y]

  int           i, j, n_b2;
  double        **A, *b, *b2_lim, *b2, *db2, step;
  double        nu_fract[2], dnu[2];
  ss_vect<tps>  nus, dnus;
  std::ofstream quad_out;

  const bool    debug = true;
  const int     m     = 2;
  const double  s_cut = 1e-7, step0 = 1.0;

  n_b2 = b2_Fam.size();

  b = dvector(1, m); b2_lim = dvector(1, n_b2);
  b2 = dvector(1, n_b2); db2 = dvector(1, n_b2);
  A = dmatrix(1, m, 1, n_b2);

  nu_fract[X_] = fract(nu_x); nu_fract[Y_] = fract(nu_y);
  std::cout << std::endl;
  std::cout << std::fixed << std::setprecision(5)
       << "fit_tune: nu_x = " << nu_fract[X_] << ", nu_y = " << nu_fract[Y_]
       << std::endl;

  for (i = 1; i <= n_b2; i++) {
    if (b2_Fam[i-1] > 0) {
      b2_lim[i] = b2_max; b2[i] = get_bn(b2_Fam[i-1], 1, Quad);
    } else {
      b2_lim[i] = ds_max; b2[i] = get_L(abs(b2_Fam[i-1]), 1);
    }
  }

  danot_(3);

  get_map_n(n_cell); K = MapNorm(Map, g, A1, A0, Map_res, 1); nus = dHdJ(K);

  dnu[X_] = nus[3].cst() - nu_fract[X_];
  dnu[Y_] = nus[4].cst() - nu_fract[Y_];
  
  while ((fabs(dnu[X_]) > eps) || (fabs(dnu[Y_]) > eps)) {
    step = step0;
    for (i = 1; i <= n_b2; i++) {
      for (j = 1; j <= get_n_Kids(abs(b2_Fam[i-1])); j++)
	if (b2_Fam[i-1] > 0)
	  set_bn_par(b2_Fam[i-1], j, Quad, 7);
	else
	  set_s_par(abs(b2_Fam[i-1]), j, 7);

      get_map_n(n_cell); K = MapNorm(Map, g, A1, A0, Map_res, 1); nus = dHdJ(K);

      A[1][i] = -h_ijklm_p(nus[3], 0, 0, 0, 0, 0, 7);
      A[2][i] = -h_ijklm_p(nus[4], 0, 0, 0, 0, 0, 7);

      if (b2_Fam[i-1] < 0)
	for (j = 1; j <= m; j++)
	  A[j][i] *= scl_ds;

      for (j = 1; j <= get_n_Kids(abs(b2_Fam[i-1])); j++)
	if (b2_Fam[i-1] > 0)
	  clr_bn_par(b2_Fam[i-1], j, Quad);
	else
	  clr_s_par(abs(b2_Fam[i-1]), j);
    }

    b[1] = dnu[X_]; b[2] = dnu[Y_];

    SVD_lim(m, n_b2, A, b, b2_lim, s_cut, b2, db2);

    for (i = 1; i <= n_b2; i++) {
      set_dbn_s(b2_Fam[i-1], Quad, step*db2[i]);
      b2[i] = get_bn_s(b2_Fam[i-1], 1, Quad);
    }

    get_map_n(n_cell); K = MapNorm(Map, g, A1, A0, Map_res, 1);

    while (!stable) {
      // roll back
      for (i = 1; i <= n_b2; i++) {
	set_dbn_s(b2_Fam[i-1], Quad, -step*db2[i]);
	b2[i] = get_bn_s(b2_Fam[i-1], 1, Quad);
      }

      step /= 2.0;
      std::cout << std::endl;
      std::cout << std::scientific << std::setprecision(3)
		<< "step = " << step << std::endl;

      for (i = 1; i <= n_b2; i++) {
	set_dbn_s(b2_Fam[i-1], Quad, step*db2[i]);
	b2[i] = get_bn_s(b2_Fam[i-1], 1, Quad);
      }
	
      get_map_n(n_cell); K = MapNorm(Map, g, A1, A0, Map_res, 1);
    }
      
    nus = dHdJ(K);

    dnu[X_] = nus[3].cst() - nu_fract[X_];
    dnu[Y_] = nus[4].cst() - nu_fract[Y_];
    
    if (debug) {
      std::cout << std::endl;
      std::cout << " Ax = b:" << std::endl;
      std::cout << std::endl;
      for (i = 1; i <= m; i++) {
	for (j = 1; j <= n_b2; j++)
	  std::cout << std::scientific << std::setprecision(3)
		    << std::setw(11) << A[i][j];
	std::cout << std::scientific << std::setprecision(3)
		  << std::setw(11) << b[i] << std::endl;
      }
	
      std::cout << std::endl;
      std::cout << std::fixed << std::setprecision(5)
	   << "dnu_x = " << dnu[X_] << ", dnu_y = " << dnu[Y_]
	   << std::endl;
    }
   
    std::cout << std::endl;
    std::cout << std::fixed << std::setprecision(5)
	      << "nu_x = " << nus[3].cst() << ", nu_y = " << nus[4].cst()
	      << std::endl;
  }

  if (prt) {
    quad_out.open("fit_tune.dat", std::ios::out);
    quad_out << std::endl;
    quad_out << "n = 1:" << std::endl;
    for (i = 1; i <= n_b2; i++)
      for (j = 1; j <= get_n_Kids(abs(b2_Fam[i-1])); j++)
	if (b2_Fam[i-1] > 0)
	  quad_out << std::fixed << std::setprecision(7) 
		   << std::setw(6) << get_Name(b2_Fam[i-1]) << "(" << j
		   << ") = " << std::setw(11) << get_bn(b2_Fam[i-1], j, Quad)
		   << std::setw(2) << Quad << std::endl;
	else {
	  quad_out << std::fixed << std::setprecision(7) 
		   << std::setw(6) << get_Name(abs(b2_Fam[i-1])) << "("
		   << j << ") = "
		   << std::setw(11) << get_L(abs(b2_Fam[i-1]), j)
		   << std::setw(3) << -Quad << std::endl;
	}
    quad_out.close();
  }

  free_dvector(b, 1, m);
  free_dvector(b2_lim, 1, n_b2); free_dvector(b2, 1, n_b2);
  free_dvector(db2, 1, n_b2);
  free_dmatrix(A, 1, m, 1, n_b2);
}


void fit_tune(const double nu_x, const double nu_y)
{
  std::vector<int> b2_Fam;

  const double eps = 1e-5;

  // DIAMOND.
  b2_Fam.push_back(get_Fnum("q1d"));
  b2_Fam.push_back(get_Fnum("q2d"));
  b2_Fam.push_back(get_Fnum("q3d"));
  b2_Fam.push_back(get_Fnum("q1b"));
  b2_Fam.push_back(get_Fnum("q2b"));
  b2_Fam.push_back(get_Fnum("q3b"));
 
  fit_tune(nu_x, nu_y, b2_Fam, eps, true);
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

    // printf("\nDBL_EPSILON = %9.2e\n", DBL_EPSILON);
    // daeps_(DBL_EPSILON);
    daeps_(1e-30);

    danot_(1);

    printf("\nscl_h:     %7.1e, %7.1e, %7.1e\n", scl_h[0], scl_h[1], scl_h[2]);
    printf("scl_dnu:   %7.1e, %7.1e\n", scl_dnu[0], scl_dnu[1]);
    printf("scl_ksi:   %7.1e, %7.1e, %7.1e\n",
	   scl_ksi[0], scl_ksi[1], scl_ksi[2]);
    printf("symmetric: %d\n", symm);
    printf("\nA_max:     %7.1e, %7.1e\n",
	   A_max[lat_case-1][X_], A_max[lat_case-1][Y_]);
    printf("delta_max: %7.1e\n", delta_max[lat_case-1]);
    printf("beta_inj:  %7.1e, %7.1e\n",
	   beta_inj[lat_case-1][X_], beta_inj[lat_case-1][Y_]);

    get_nu_ksi();

    if (false) {
      fit_tune(28.15, 13.22);
      get_nu_ksi();
      exit(0);
    }

    for (j = 0; j < 2; j++)
      twoJ[j] =	sqr(A_max[lat_case-1][j])/beta_inj[lat_case-1][j];

    Id_scl.identity();
    Id_scl[x_] *= sqrt(twoJ[X_]); Id_scl[px_] *= sqrt(twoJ[X_]);
    Id_scl[y_] *= sqrt(twoJ[Y_]); Id_scl[py_] *= sqrt(twoJ[Y_]);
    Id_scl[delta_] *= delta_max[lat_case-1];

    if (false) {
      danot_(NO-1);
      cavity_on = true; rad_on = true;
      get_map_n(n_cell);
      // MAX-VI:
      // prt_H_long(10, M_PI, 10e-2, -405.6e3, false);
      // SLS-2:
      prt_H_long(10, M_PI, 10e-2, -544.7e3, true);
      prt_alphac();
      exit(0);
    }

    if (false) {
      if (false) {
	no_mpoles(Sext);
	no_mpoles(Oct);
      }

      danot_(NO-1);
      get_map_n(n_cell);
      prt_alphac();

      prt_ct(100, 8e-2);
      exit(0);
    }

    if (false) {
      prt_map();
      exit(0);
    }

    if (false) {
      prt_h_K();
      exit(0);
    }

    switch (lat_case) {
    case 1:
      // MAX VI:
      bn_prms.add_prm("o1", 4, 5e5, 1.0);
      bn_prms.add_prm("o2", 4, 5e5, 1.0);
      bn_prms.add_prm("o3", 4, 5e5, 1.0);
      bn_prms.add_prm("o4", 4, 5e5, 1.0);

      bn_prms.add_prm("o1", 6, 5e10, 1.0);
      bn_prms.add_prm("o2", 6, 5e10, 1.0);
      bn_prms.add_prm("o3", 6, 5e10, 1.0);
      bn_prms.add_prm("o4", 6, 5e10, 1.0);
      break;
    case 2:
      // SLS-2:
      if (true) {
	bn_prms.add_prm("sfh",  3, 5e5, 1.0);
	bn_prms.add_prm("sdh",  3, 5e5, 1.0);
	bn_prms.add_prm("sfmh", 3, 5e5, 1.0);
	bn_prms.add_prm("sdmh", 3, 5e5, 1.0);

	bn_prms.add_prm("sxxh", 3, 5e5, 1.0);
	bn_prms.add_prm("sxyh", 3, 5e5, 1.0);
	bn_prms.add_prm("syyh", 3, 5e5, 1.0);
      } else {
	bn_prms.add_prm("ocx",  4, 5e10, 1.0);
	bn_prms.add_prm("ocxm", 4, 5e10, 1.0);
	bn_prms.add_prm("ocy",  4, 5e10, 1.0);
	bn_prms.add_prm("ocym", 4, 5e10, 1.0);

	bn_prms.add_prm("oxx",  4, 5e10, 1.0);
	bn_prms.add_prm("oxy",  4, 5e10, 1.0);
	bn_prms.add_prm("oyy",  4, 5e10, 1.0);
      }
      break;
    case 3 ... 4:
      // DIAMOND:
      if (true) {
	bn_prms.add_prm("ts1a",  3, 5e5, 1.0);
	bn_prms.add_prm("ts1ab", 3, 5e5, 1.0);
	bn_prms.add_prm("ts2a",  3, 5e5, 1.0);
	bn_prms.add_prm("ts2ab", 3, 5e5, 1.0);
	bn_prms.add_prm("ts1b",  3, 5e5, 1.0);
	bn_prms.add_prm("ts2b",  3, 5e5, 1.0);
	bn_prms.add_prm("ts1c",  3, 5e5, 1.0);
	bn_prms.add_prm("ts2c",  3, 5e5, 1.0);
	// bn_prms.add_prm("ts1d",  3, 5e5, 1.0);
	bn_prms.add_prm("ts2d",  3, 5e5, 1.0);
	bn_prms.add_prm("ts1e",  3, 5e5, 1.0);
	bn_prms.add_prm("ts2e",  3, 5e5, 1.0);

	if (lat_case == 4) {
	  // VMX.
	  bn_prms.add_prm("s1", 3, 5e5, 1.0);
	  bn_prms.add_prm("s2", 3, 5e5, 1.0);
	  bn_prms.add_prm("s3", 3, 5e5, 1.0);
	  bn_prms.add_prm("s4", 3, 5e5, 1.0);
	  bn_prms.add_prm("s5", 3, 5e5, 1.0);
	}
      } else {
	bn_prms.add_prm("ts1a",  3, 5e5, 1.0);
	bn_prms.add_prm("ts2a",  3, 5e5, 1.0);
	bn_prms.add_prm("ts1b",  3, 5e5, 1.0);
	bn_prms.add_prm("ts2b",  3, 5e5, 1.0);
	bn_prms.add_prm("ts1c",  3, 5e5, 1.0);
	bn_prms.add_prm("ts2c",  3, 5e5, 1.0);

	// bn_prms.add_prm("s1",    3, 5e5, 1.0);
	// bn_prms.add_prm("s2",    3, 5e5, 1.0);
	// bn_prms.add_prm("s3",    3, 5e5, 1.0);
	// bn_prms.add_prm("s4",    3, 5e5, 1.0);
	// bn_prms.add_prm("s5",    3, 5e5, 1.0);
      }
      break;
    case 5:
      // DIAMOND-II, 4-BA:
      bn_prms.add_prm("s1b", 3, 5e5, 1.0);
      bn_prms.add_prm("s1d", 3, 5e5, 1.0);
      bn_prms.add_prm("s2b", 3, 5e5, 1.0);
      bn_prms.add_prm("s2d", 3, 5e5, 1.0);
      bn_prms.add_prm("sx1", 3, 5e5, 1.0);
      bn_prms.add_prm("sy1", 3, 5e5, 1.0);

      bn_prms.add_prm("s3",  3, 5e5, 1.0);
      break;
    case 6:
      // DIAMOND-II, 6-BA:
      bn_prms.add_prm("sd1",  3, 5e5, 1.0);
      bn_prms.add_prm("sd2",  3, 5e5, 1.0);
      bn_prms.add_prm("sd3",  3, 5e5, 1.0);
      bn_prms.add_prm("sf21", 3, 5e5, 1.0);
      bn_prms.add_prm("sd31", 3, 5e5, 1.0);
      bn_prms.add_prm("sf1",  3, 5e5, 1.0);
      bn_prms.add_prm("sh1a", 3, 5e5, 1.0);
      bn_prms.add_prm("sh1e", 3, 5e5, 1.0);
      break;
    case 7:
      if (true) {
	// DIAMOND-II, 6-BA_Beni:
	bn_prms.add_prm("sfa", 3, 5e5, 1.0);
	bn_prms.add_prm("sfb", 3, 5e5, 1.0);
	bn_prms.add_prm("sda", 3, 5e5, 1.0);
	bn_prms.add_prm("sdb", 3, 5e5, 1.0);
	bn_prms.add_prm("s1",  3, 5e5, 1.0);
	bn_prms.add_prm("s2a", 3, 5e5, 1.0);
	bn_prms.add_prm("s2b", 3, 5e5, 1.0);
	bn_prms.add_prm("s3",  3, 5e5, 1.0);
	bn_prms.add_prm("s4",  3, 5e5, 1.0);
	bn_prms.add_prm("s5",  3, 5e5, 1.0);
	bn_prms.add_prm("s6",  3, 5e5, 1.0);
      } else {
	bn_prms.add_prm("sf", 3, 5e5, 1.0);
	bn_prms.add_prm("sd", 3, 5e5, 1.0);
	bn_prms.add_prm("s1", 3, 5e5, 1.0);
	bn_prms.add_prm("s2", 3, 5e5, 1.0);
	bn_prms.add_prm("s3", 3, 5e5, 1.0);
	bn_prms.add_prm("s4", 3, 5e5, 1.0);
	bn_prms.add_prm("s5", 3, 5e5, 1.0);
	bn_prms.add_prm("s6", 3, 5e5, 1.0);
      }
     break;
    case 8:
      // DIAMOND-II, 8-BA:
      bn_prms.add_prm("sfh",  3, 5e5, 1.0);
      bn_prms.add_prm("sdh",  3, 5e5, 1.0);
      bn_prms.add_prm("sfmh", 3, 5e5, 1.0);
      bn_prms.add_prm("sdmh", 3, 5e5, 1.0);
      if (true) {
	bn_prms.add_prm("sxxh", 3, 5e5, 1.0);
	bn_prms.add_prm("sxyh", 3, 5e5, 1.0);
	bn_prms.add_prm("syyh", 3, 5e5, 1.0);
      }
      break;
    case 9:
      // DIAMOND-II, 8-BA by Hossein:
      bn_prms.add_prm("s1",  3, 5e5, 1.0);
      bn_prms.add_prm("s2",  3, 5e5, 1.0);
      bn_prms.add_prm("s3",  3, 5e5, 1.0);
      bn_prms.add_prm("s4",  3, 5e5, 1.0);
      bn_prms.add_prm("s5",  3, 5e5, 1.0);
      bn_prms.add_prm("s6",  3, 5e5, 1.0);
      bn_prms.add_prm("s7",  3, 5e5, 1.0);
      bn_prms.add_prm("s8",  3, 5e5, 1.0);
      bn_prms.add_prm("s9",  3, 5e5, 1.0);
      bn_prms.add_prm("s10", 3, 5e5, 1.0);
      break;
    }

    // Step is 1.0 for conjugated gradient method.
    bn_prms.bn_tol = 1e-1; bn_prms.svd_cut = 1e-15; bn_prms.step = 1.0;

    // no_mpoles(Sext); no_mpoles(Oct); no_mpoles(Dodec);

    bn_prms.ini_prm();

    prt_bn(bn_prms);

    if (true) {
      fit_ksi1(0e0, 0e0);
      // exit(0);
    }

    if (false)
      min_conj_grad(true);
    else
      min_lev_marq();

    prt_h_K();
}
