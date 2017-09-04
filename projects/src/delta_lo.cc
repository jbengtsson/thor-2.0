#define NO 2

#include "thor_lib.h"

int no_tps = NO,

#define DOF_3 0

#if !DOF_3
  ndpt_tps = 5;
#else
  // Requires that cavity is turned on.
  ndpt_tps = 0;
#endif

// Initial conditions: alpha, beta, eta, etap.
// Provided, roughly periodic:
// const double ic[][2] =
//   {{1.05266, -0.25384}, {0.62733, 5.60502}, {0.06552, 0.0}, {-0.10478, 0.0}};
const double ic[][2] =
  {{1.15199, -0.22236}, {0.65878, 5.53043 }, {0.03741, 0.0}, {-0.04304, 0.0}};

int loc[10], n;

double get_bn_s1(const int Fnum, const int Knum, const int n);
void set_bn_s1(const int Fnum, const int n, const double dbn);
void get_S(void);


struct param_type {
private:

public:
  int                 n_prm;
  double              bn_tol, step;
  std::vector<double> bn_min, bn_max, bn_scl;
  std::vector<int>    Fnum, n;

  void add_prm(const std::string Fname, const int n,
	       const double bn_min, const double bn_max,
	       const double bn_scl);
  void ini_prm(double *bn, double *bn_lim);
  void set_prm(double *bn) const;
  void prt_prm(double *bn) const;
};


param_type b2_prms;


void param_type::add_prm(const std::string Fname, const int n,
			 const double bn_min, const double bn_max,
			 const double bn_scl)
{
  Fnum.push_back(get_Fnum(Fname.c_str()));
  this->n.push_back(n);
  this->bn_min.push_back(bn_min);
  this->bn_max.push_back(bn_max);
  this->bn_scl.push_back(bn_scl);
  n_prm = Fnum.size();
}


void param_type::ini_prm(double *bn, double *bn_lim)
{
  int    i;

  n_prm = Fnum.size();
  for (i = 1; i <= n_prm; i++) {
    bn_lim[i] = bn_max[i-1];
    if (n[i-1] > 0)
      // Multipole.
      bn[i] = get_bn(Fnum[i-1], 1, n[i-1]);
    else if (n[i-1] == -1)
      // Drift.
      bn[i] = get_L(Fnum[i-1], 1);
    else if (n[i-1] == -2)
      // Placement.
      bn[i] = get_bn_s1(-Fnum[i-1], 1, n[i-1]);
  }
}


void param_type::set_prm(double *bn) const
{
  int i;

  for (i = 1; i <= n_prm; i++) {
    if (n[i-1] > 0)
	set_bn(Fnum[i-1], n[i-1], bn[i]);
    else if (n[i-1] == -1) {
      set_L(Fnum[i-1], bn[i]); get_S();
    } else if (n[i-1] == -2)
      set_bn_s1(-Fnum[i-1], n[i-1], bn[i]);
  }
}


void param_type::prt_prm(double *bn) const
{
  int i;

  const int n_prt = 8;

  for (i = 1; i <= n_prm; i++) {
    printf(" %9.5f", bn[i]);
    if (i % n_prt == 0) printf("\n");
  }
  if (n_prm % n_prt != 0) printf("\n");
}


double get_bn_s1(const int Fnum, const int Knum, const int n)
{
  int    k;
  double bn;

  if (Fnum > 0)
    bn = get_bn(Fnum, Knum, n);
  else {
    k = get_loc(abs(Fnum), Knum);

    switch (elem[k-1].Name[1]) {
    case 'u':
      bn = elem[k-1].L;
      break;
    case 'd':
      bn = elem[k+1].L;
      break;
    default:
      printf("get_bn_s1: configuration error %s (%d)\n",
	     elem[k-1].Name, k);
      exit(1);
      break;
    }
  }

  return bn;
}


void set_bn_s1(const int Fnum, const int Knum, const int n, const double bn)
{
  char name[name_length];
  int  loc, loc_d;

  if (Fnum > 0)
    set_bn(Fnum, Knum, n, bn);
  else {
    // Point to multipole.
    loc = get_loc(abs(Fnum), Knum);

    if (elem[loc-1].Name[1] == 'u') {
      strcpy(name, elem[loc-1].Name); name[1] = 'd';
      loc_d = get_Fnum(name);
      set_L(elem[loc-1].Fnum, Knum, bn);
      set_L(elem[loc_d].Fnum, Knum, -bn);
    } else if (elem[loc-1].Name[1] == 'd') {
      strcpy(name, elem[loc-1].Name); name[1] = 'u';
      loc_d = get_Fnum(name);
      set_L(elem[loc-1].Fnum, Knum, -bn);
      set_L(elem[loc_d].Fnum, Knum, bn);
    } else if (elem[loc+1].Name[1] == 'd') {
      strcpy(name, elem[loc+1].Name); name[1] = 'u';
      loc_d = get_Fnum(name);
      set_L(elem[loc+1].Fnum, Knum, -bn);
      set_L(elem[loc_d].Fnum, Knum, bn);
    } else if (elem[loc+1].Name[1] == 'u') {
      strcpy(name, elem[loc+1].Name); name[1] = 'd';
      loc_d = get_Fnum(name);
      set_L(elem[loc+1].Fnum, Knum, bn);
      set_L(elem[loc_d].Fnum, Knum, -bn);
    }
  }
}


void set_bn_s1(const int Fnum, const int n, const double bn)
{
  int k;

  for (k = 1; k <= get_n_Kids(abs(Fnum)); k++)
    set_bn_s1(Fnum, k, n, bn);
}


void no_mpoles(const int n)
{
  int j;

  printf("\nzeroing multipoles: %d\n", n);
  for (j = 0; j < n_elem; j++)
    if (elem[j].kind == Mpole)
      set_bn(elem[j].Fnum, elem[j].Knum, n, 0e0);
}


void get_S(void)
{
  int    j;
  double S;

  S = 0e0;
  for (j = 0; j <= n_elem; j++) {
    S += elem[j].L; elem[j].S = S;
  }
}


void get_twiss(const int i0, const int i1, const ss_vect<tps> &A)
{
  int          j, k;
  double       alpha1[2], beta1[2], eta1[2], etap1[2], dnu1[2], dnu2[2];
  ss_vect<tps> A1;

  // Include parameter dependence.
  danot_(2);

  for (k = 0; k < 2; k++)
    dnu1[k] = 0e0;
  A1 = A;
  for (j = i0; j <= i1; j++) {
    A1.propagate(j, j);
    elem_tps[j-1].A1 = get_A_CS(2, A1, dnu2);

    // Store linear optics for convenience.
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
    for (k = 0; k < 2; k++)
      dnu1[k] = dnu2[k];
  }
}


void get_twiss(const int i0, const int i1,
	       const double alpha[], const double beta[],
	       const double eta[], const double etap[])
{
  ss_vect<tps> A;

  // Include parameter dependence.
  danot_(2);

  A = get_A(alpha, beta, eta, etap); get_twiss(i0, i1, A);
}


void get_twiss(void)
{
  // Periodic.

  danot_(1);
  get_Map(); K = MapNorm(Map, g, A1, A0, Map_res, 1);
  // Include dispersion.
  A1[x_]  += A0[x_][delta_]*tps(0e0, delta_+1);
  A1[px_] += A0[px_][delta_]*tps(0e0, delta_+1);
  get_twiss(0+1, n_elem, A1);
}


void get_dnu(const int i0, const int i1, const double alpha0[],
	     const double beta0[], const double dp, double dnu[])
{
  int          k;
  double       m11, m12;
  ss_vect<tps> map;

  map.identity(); map[delta_] += dp;
  map.propagate(i0+1, i1);

  for (k = 0; k < 2; k++) {
    m11 = map[2*k][2*k]; m12 = map[2*k][2*k+1];
    dnu[k] = atan(m12/(beta0[k]*m11-alpha0[k]*m12))/(2e0*M_PI);
    if (m11 < 0e0) dnu[k] += (m12 >= 0)? 0.5e0 : -0.5e0;
    if (dnu[k] < 0e0) dnu[k] += 1e0;
  }
}


void get_dnu_dp(const int i0, const int i1, const double alpha0[],
		const double beta0[], const double dp, double dksi[])
{
  // To evaluate linear chromaticity the linear dispersion for the periodic
  // solution must be known.
  int    k;
  double dnu1[2], dnu0[2];

  get_dnu(loc[0], loc[1], alpha0, beta0, dp, dnu1);
  get_dnu(loc[0], loc[1], alpha0, beta0, -dp, dnu0);

  for (k = 0; k < 2; k++)
    dksi[k] = (dnu1[k]-dnu0[k])/(2e0*dp);
}


void prt_match(const param_type &b2_prms, const double *b2)
{
  int  k;
  FILE *outf;

  std::string file_name = "match.out";

  outf = file_write(file_name.c_str());

  k = 0;
  fprintf(outf, "Q01:  quadrupole, l = 0.234, k = %8.5f, N = Nquad"
	  ", Method = Meth;\n", b2[++k]);
  fprintf(outf, "Q03:  quadrupole, l = 0.434, k = %8.5f, N = Nquad"
	  ", Method = Meth;\n", b2[++k]);

  fprintf(outf, "\nEQ01: quadrupole, l = 0.234, k = %8.5f, N = Nquad"
	  ", Method = Meth;\n", b2[++k]);
  fprintf(outf, "EQ02: quadrupole, l = 0.234, k = %8.5f, N = Nquad"
	  ", Method = Meth;\n", b2[++k]);
  fprintf(outf, "Q02:  quadrupole, l = 0.234, k = %8.5f, N = Nquad"
	  ", Method = Meth;\n", b2[++k]);

  fprintf(outf, "EQ04: quadrupole, l = 0.234, k = %8.5f, N = Nquad"
	  ", Method = Meth;\n", b2[++k]);
  fprintf(outf, "EQ05: quadrupole, l = 0.234, k = %8.5f, N = Nquad"
	  ", Method = Meth;\n", b2[++k]);
  fprintf(outf, "EQ06: quadrupole, l = 0.234, k = %8.5f, N = Nquad"
	  ", Method = Meth;\n", b2[++k]);

  fprintf(outf, "\nD_Q01_L  = %8.5f;\n", b2[++k]);
  fprintf(outf, "D_Q03_L  = %8.5f;\n", b2[++k]);

  fprintf(outf, "\nD_EQ01_L = %8.5f;\n", b2[++k]);
  fprintf(outf, "D_EQ02_L = %8.5f;\n", b2[++k]);
  fprintf(outf, "D_Q02_L  = %8.5f;\n", b2[++k]);

  fprintf(outf, "D_EQ04_L = %8.5f;\n", b2[++k]);
  fprintf(outf, "D_EQ05_L = %8.5f;\n", b2[++k]);
  fprintf(outf, "D_EQ06_L = %8.5f;\n", b2[++k]);

  fprintf(outf, "\nD_B10_L  = %8.5f;\n", b2[++k]);

  // fprintf(outf, "\nU561: drift, L = %8.5f;\n", b2[18]);

  fclose(outf);
}


void prt_lin_opt(const int loc[])
{
  printf("\n      s    alpha_x  beta_x  eta_x  etap_x  alpha_y  beta_y\n");
  printf("  %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f\n",
	 elem[loc[0]].S,
	 elem[loc[0]].Alpha[X_], elem[loc[0]].Beta[X_],
	 elem[loc[0]].Eta[X_], elem[loc[0]].Etap[X_],
	 elem[loc[0]].Alpha[Y_], elem[loc[0]].Beta[Y_]);
  printf("  %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f\n",
	 elem[loc[1]].S,
	 elem[loc[1]].Alpha[X_], elem[loc[1]].Beta[X_],
	 elem[loc[1]].Eta[X_], elem[loc[1]].Etap[X_],
	 elem[loc[1]].Alpha[Y_], elem[loc[1]].Beta[Y_]);
  printf("  %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f\n",
	 elem[loc[2]].S,
	 elem[loc[2]].Alpha[X_], elem[loc[2]].Beta[X_],
	 elem[loc[2]].Eta[X_], elem[loc[2]].Etap[X_],
	 elem[loc[2]].Alpha[Y_], elem[loc[2]].Beta[Y_]);
  printf("  %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f\n",
	 elem[loc[3]].S,
	 elem[loc[3]].Alpha[X_], elem[loc[3]].Beta[X_],
	 elem[loc[3]].Eta[X_], elem[loc[3]].Etap[X_],
	 elem[loc[3]].Alpha[Y_], elem[loc[3]].Beta[Y_]);
  printf("  %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f\n",
	 elem[loc[4]].S,
	 elem[loc[4]].Alpha[X_], elem[loc[4]].Beta[X_],
	 elem[loc[4]].Eta[X_], elem[loc[4]].Etap[X_],
	 elem[loc[4]].Alpha[Y_], elem[loc[4]].Beta[Y_]);
}


double f_match(double *b2)
{
  static double chi2_ref = 1e30, chi2_prt = 1e30;

  int          i, loc1, loc2;
  double       chi2, dksi[2], L;
  ss_vect<tps> Ascr;

  const int n_prt = 50;

  b2_prms.set_prm(b2);

  Ascr = get_A(ic[0], ic[1], ic[2], ic[3]);
  get_twiss(loc[0]+1, loc[4]+1, Ascr);

  // get_dnu_dp(loc[0], loc[5], ic[0], ic[1], 1e-5, dksi);

  chi2 = 0e0;
  // Downstream of 10 degree dipole.
  chi2 += 1e12*sqr(elem[loc[1]].Eta[X_]);
  chi2 += 1e12*sqr(elem[loc[1]].Etap[X_]);
  // chi2 += 1e5*sqr(elem[loc[1]].Beta[Y_]);

  // Center of 1st straight.
  chi2 += 1e9*sqr(elem[loc[2]].Alpha[X_]);
  chi2 += 1e9*sqr(elem[loc[2]].Alpha[Y_]);
  chi2 += 1e6*sqr(elem[loc[2]].Beta[X_]-9.58);

  // Center of 2nd straight.
  chi2 += 1e9*sqr(elem[loc[3]].Alpha[X_]);
  chi2 += 1e9*sqr(elem[loc[3]].Alpha[Y_]);
  chi2 += 1e6*sqr(elem[loc[3]].Beta[X_]-8.0); 

  for (i = 1; i <= b2_prms.n_prm; i++) {
    loc1 = get_loc(b2_prms.Fnum[i-1], 1);
    L = elem[loc1].L;
    if (i == 3) {
      // Upstream of the quadrupole.
      chi2 += 1e-10*sqr(b2[i]*L*elem[loc1].Beta[X_]);
      chi2 += 1e6*sqr(b2[i]*L*elem[loc1].Beta[Y_]);
    } else if (i == 4) {
      chi2 += 1e-10*sqr(b2[i]*L*elem[loc1-1].Beta[X_]);
      chi2 += 1e-10*sqr(b2[i]*L*elem[loc1-1].Beta[Y_]);
    } else if (i == 5) {
      chi2 += 1e-10*sqr(b2[i]*L*elem[loc1].Beta[X_]);
      chi2 += 1e-10*sqr(b2[i]*L*elem[loc1].Beta[Y_]);
    } else if (i == 7) {
      chi2 += 1e-10*sqr(b2[i]*L*elem[loc1].Beta[X_]);
      chi2 += 1e-10*sqr(b2[i]*L*elem[loc1].Beta[Y_]);
    } else {
      chi2 += 1e-10*sqr(b2[i]*L*elem[loc1].Beta[X_]);
      chi2 += 1e-10*sqr(b2[i]*L*elem[loc1].Beta[Y_]);
    }
  }

  for (i = 1; i <= b2_prms.n_prm; i++) {
    if (b2[i] > b2_prms.bn_max[i-1])
      chi2 += 1e12*sqr(b2[i]-b2_prms.bn_max[i-1]);
    if (b2[i] < b2_prms.bn_min[i-1])
      chi2 += 1e12*sqr(b2[i]-b2_prms.bn_min[i-1]);
  }

  if (chi2 < chi2_ref) {
    n++;

    if (n % n_prt == 0) {
      printf("\n%3d chi2: %12.5e -> %12.5e\n", n, chi2_prt, chi2);
      printf("b2s:\n");
      b2_prms.prt_prm(b2);

      // printf("\ndksi: %8.5f %8.5f\n", dksi[X_], dksi[Y_]);

      // Downstream of 10 degree dipole.
      printf("\nDownstream of 10 degree dipole:\n");
      printf("eta_x   = %8.5f etap_x  = %8.5f\n",
	     elem[loc[1]].Eta[X_], elem[loc[1]].Etap[X_]);
      printf("beta_x  = %8.5f beta_y  = %8.5f\n",
	     elem[loc[1]].Beta[X_], elem[loc[1]].Beta[Y_]);
      // Center of 1st straight.
      printf("\nCenter of 1st straight:\n");
      printf("alpha_x = %8.5f alpha_y = %8.5f\n",
	     elem[loc[2]].Alpha[X_], elem[loc[2]].Alpha[Y_]);
      printf("beta_x  = %8.5f beta_y  = %8.5f\n",
	     elem[loc[2]].Beta[X_], elem[loc[2]].Beta[Y_]);

      // Center of 2nd straight.
      printf("\nCenter of 2nd straight:\n");
      printf("\nalpha_x = %8.5f alpha_y = %8.5f\n",
	     elem[loc[3]].Alpha[X_], elem[loc[3]].Alpha[Y_]);
      printf("beta_x  = %8.5f beta_y  = %8.5f\n",
	     elem[loc[3]].Beta[X_], elem[loc[3]].Beta[Y_]);

      loc1 = get_loc(get_Fnum("s_s_1"), 1);
      loc2 = get_loc(get_Fnum("s_s_1"), 2);
      printf("\nLength of 1st straight: %6.3f m\n", elem[loc2].S-elem[loc1].S);
      loc1 = get_loc(get_Fnum("s_s_2"), 1);
      loc2 = get_loc(get_Fnum("s_s_2"), 2);
      printf("Length of 2nd straight: %6.3f m\n", elem[loc2].S-elem[loc1].S);
      loc1 = get_loc(get_Fnum("s_s_3"), 1);
      loc2 = get_loc(get_Fnum("s_s_3"), 2);
      printf("Length of 3rd straight: %6.3f m\n", elem[loc2].S-elem[loc1].S);

      prt_match(b2_prms, b2);

      prt_mfile("flat_file.fit");
      prt_lat(loc[0], loc[4], "linlat1.out");
      prt_lat(loc[0], loc[4], "linlat.out", 10);

      chi2_prt = min(chi2, chi2_prt);
    }

    chi2_ref = min(chi2, chi2_ref);
  }

  return chi2;
}


void fit_match(param_type &b2_prms)
{
  int          n_b2, i, j, iter;
  double       *b2, *b2_lim, **xi, fret;
  ss_vect<tps> Ascr;

  n_b2 = b2_prms.n_prm;

  b2 = dvector(1, n_b2); b2_lim = dvector(1, n_b2);
  xi = dmatrix(1, n_b2, 1, n_b2);

  // Upstream of 20 degree dipole.
  loc[0] = get_loc(get_Fnum("sb"),  7) - 1;
  // Downstream of 10 degree dipole.
  loc[1] = get_loc(get_Fnum("b10"), 1) - 1;
  // Center of 1st straight.
  loc[2] = get_loc(get_Fnum("ef2"), 4) - 1;
  // Center of 2nd straight.
  loc[3] = get_loc(get_Fnum("ef2"), 16) - 1;
  // Downstream of 20 degree dipole.
  loc[4] = get_loc(get_Fnum("b20"), 5) - 1;

  Ascr = get_A(ic[0], ic[1], ic[2], ic[3]);
  get_twiss(loc[0]+1, loc[4]+1, Ascr);

  prt_lin_opt(loc);
  printf("\n%8.5f %8.5f\n",
	 elem[loc[5]].Nu[X_]-elem[loc[0]].Nu[X_],
	 elem[loc[5]].Nu[Y_]-elem[loc[0]].Nu[Y_]);

  prt_lat(loc[0], loc[4], "linlat1.out");
  prt_lat(loc[0], loc[4], "linlat.out", 10);

  b2_prms.ini_prm(b2, b2_lim);

  // Set initial directions (unit vectors).
  for (i = 1; i <= n_b2; i++)
    for (j = 1; j <= n_b2; j++)
      xi[i][j] = (i == j)? 1e0 : 0e0;

  n = 0;
  dpowell(b2, xi, n_b2, 1e-16, &iter, &fret, f_match);

  b2_prms.set_prm(b2);
  Ascr = get_A(ic[0], ic[1], ic[2], ic[3]);
  get_twiss(loc[0]+1, loc[4]+1, Ascr);

  prt_lat(loc[0], loc[4], "linlat1.out");
  prt_lat(loc[0], loc[4], "linlat.out", 10);

  free_dvector(b2, 1, n_b2);  free_dvector(b2_lim, 1, n_b2);
  free_dmatrix(xi, 1, n_b2, 1, n_b2);
}


int main(int argc, char *argv[])
{

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

  daeps_(1e-30);

  no_mpoles(Sext);
  no_mpoles(Oct);

  get_twiss();
  prt_lat("linlat1.out");
  prt_lat("linlat.out", 10);
  exit(0);
  
  if (true) {
    b2_prms.add_prm("q01",  2, -4.2, 4.2, 1.0);
    b2_prms.add_prm("q03",  2, -4.2, 4.2, 1.0);

    b2_prms.add_prm("eq01", 2, -4.2, 4.2, 1.0);
    b2_prms.add_prm("eq02", 2, -4.2, 4.2, 1.0);
    b2_prms.add_prm("q02",  2, -4.2, 4.2, 1.0);

    b2_prms.add_prm("eq04", 2, -4.2, 4.2, 1.0);
    b2_prms.add_prm("eq05", 2, -4.2, 4.2, 1.0);
    b2_prms.add_prm("eq06", 2, -4.2, 4.2, 1.0);

    b2_prms.add_prm("q01",  -2,  0.0,  0.05, 1.0);
    b2_prms.add_prm("q03",  -2,  0.0,  0.05, 1.0);

    b2_prms.add_prm("eq01", -2,  0.0,  0.05, 1.0);
    b2_prms.add_prm("eq02", -2,  0.0,  0.05, 1.0);
    b2_prms.add_prm("q02",  -2,  0.0,  0.05, 1.0);

    b2_prms.add_prm("eq04", -2, -0.05, 0.05, 1.0);
    b2_prms.add_prm("eq05", -2,  0.0,  0.05, 1.0);
    b2_prms.add_prm("eq06", -2,  0.0,  0.05, 1.0);

    b2_prms.add_prm("b10",  -2, -0.01, 0.01,  1.0);

    // U561 + U562: 2.14.
    // b2_prms.add_prm("u561", -1, 2.14, 2.14, 1.0);

    b2_prms.bn_tol = 1e-6; b2_prms.step = 1.0;

    fit_match(b2_prms);
  }
}
