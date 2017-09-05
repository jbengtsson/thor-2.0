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

int    loc[10], n;
double chi2 = 0e0, *f_lm, **A_lm;

const int n_prt = 8;

const double scl_eta = 1e4, scl_beta = 1e0, scl_alpha = 1e3;

double bn_internal(const double bn_bounded,
		   const double bn_min, const double bn_max);
double bn_bounded(const double bn_internal,
		  const double bn_min, const double bn_max);
double get_bn_s1(const int Fnum, const int Knum);
void set_bn_s1(const int Fnum, const double ds);
void set_s1_par(const int Fnum, const int j);
void clr_s1_par(const int Fnum);
void get_S(void);


struct param_type {
private:

public:
  int                 m_constr, n_prm;
  double              *bn;
  std::vector<double> bn_min, bn_max, bn_scl;
  std::vector<int>    Fnum, n;

  void add_prm(const std::string Fname, const int n,
	       const double bn_min, const double bn_max, const double bn_scl);
  void ini_prm(void);
  void set_prm_dep(const int k) const;
  void clr_prm_dep(const int k) const;
};


param_type   bn_prms;
int          n_iter, n_powell;
double       twoJ[2];
ss_vect<tps> Id_scl;

void param_type::add_prm(const std::string Fname, const int n,
			 const double bn_min, const double bn_max,
			 const double bn_scl)
{
  Fnum.push_back(get_Fnum(Fname.c_str()));
  this->n.push_back(n);
  this->bn_min.push_back(bn_scl*bn_min);
  this->bn_max.push_back(bn_scl*bn_max);
  this->bn_scl.push_back(bn_scl);
  n_prm = Fnum.size();
}


void param_type::ini_prm(void)
{
  int i;

  n_prm = Fnum.size();

  bn_prms.bn = dvector(1, n_prm);

  printf("\nInitial bn; (incl. scaling) (%d):\n", n_prm);
  for (i = 1; i <= n_prm; i++) {
    if (n[i-1] > 0)
      // Multipole.
      bn[i] = get_bn(Fnum[i-1], 1, n[i-1]);
    else if (n[i-1] == -1)
      // Drift.
      bn[i] = get_L(Fnum[i-1], 1);
    else if (n[i-1] == -2)
      // Location.
      bn[i] = get_bn_s1(Fnum[i-1], 1);
    printf(" %12.5e", bn[i]);
    if (i % n_prt == 0) printf("\n");
    // Bounded.
    bn[i] = bn_internal(bn[i]/bn_scl[i-1], bn_min[i-1], bn_max[i-1]);
  }
  if (n_prm % n_prt != 0) printf("\n");
}


void param_type::set_prm_dep(const int k) const
{

  if (n[k] > 0)
    set_bn_par(Fnum[k], n[k], 7);
  else if (n[k] == -1)
    set_L_par(Fnum[k], 7);
  else if (n[k] == -2)
    set_s1_par(Fnum[k], 7);
}


void param_type::clr_prm_dep(const int k) const
{

  if (n[k] > 0)
    clr_bn_par(Fnum[k], n[k]);
  else if (n[k] == -1)
    clr_L_par(Fnum[k]);
  else if (n[k] == -2)
    clr_s1_par(Fnum[k]);
}


void set_prm(const int Fnum, const int n, const double bn,
	     const double bn_min, const double bn_max)
{
  double bn_ext;

  // Bounded.
  bn_ext = bn_bounded(bn, bn_min, bn_max);
  if (n > 0)
    set_bn(Fnum, n, bn_ext);
  else if (n == -1)
    set_L(Fnum, bn_ext);
  else if (n == -2)
    set_bn_s1(Fnum, bn_ext);
}


double bn_internal(const double bn_bounded,
		   const double bn_min, const double bn_max)
{
  return asin((2e0*(bn_bounded-bn_min))/(bn_max-bn_min)-1e0);
}


double bn_bounded(const double bn_internal,
		  const double bn_min, const double bn_max)
{
  return bn_min + (sin(bn_internal)+1e0)*(bn_max-bn_min)/2e0;
}


double get_bn_s1(const int Fnum, const int Knum)
{
  int    loc;
  double ds;

  loc = get_loc(Fnum, Knum) - 1;
  switch (elem[loc-1].Name[1]) {
  case 'u':
    ds = elem[loc-1].L;
    break;
  case 'd':
    ds = -elem[loc-1].L;
    break;
  default:
    printf("/nget_bn_s1: configuration error %s (%d)\n",
	   elem[loc].Name, loc);
    exit(1);
    break;
  }

  return ds;
}


void set_bn_s1(const int Fnum, const int Knum, const double ds)
{
  char name[name_length];
  int  loc, loc_d;

  // Point to multipole.
  loc = get_loc(Fnum, Knum) - 1;

  if (elem[loc-1].Name[1] == 'u') {
    strcpy(name, elem[loc-1].Name); name[1] = 'd';
    loc_d = get_loc(get_Fnum(name), Knum) - 1;
    set_L(elem[loc-1].Fnum, Knum, ds);
    set_L(elem[loc_d].Fnum, Knum, -ds);
  } else if (elem[loc-1].Name[1] == 'd') {
    strcpy(name, elem[loc-1].Name); name[1] = 'u';
    loc_d = get_loc(get_Fnum(name), Knum) - 1;
    set_L(elem[loc-1].Fnum, Knum, -ds);
    set_L(elem[loc_d].Fnum, Knum, ds);
  } else if (elem[loc+1].Name[1] == 'd') {
    strcpy(name, elem[loc+1].Name); name[1] = 'u';
    loc_d = get_loc(get_Fnum(name), Knum) - 1;
    set_L(elem[loc+1].Fnum, Knum, -ds);
    set_L(elem[loc_d].Fnum, Knum, ds);
  } else if (elem[loc+1].Name[1] == 'u') {
    strcpy(name, elem[loc+1].Name); name[1] = 'd';
    loc_d = get_loc(get_Fnum(name), Knum) - 1;
    set_L(elem[loc+1].Fnum, Knum, ds);
    set_L(elem[loc_d].Fnum, Knum, -ds);
  } else {
    printf("\nset_bn_s1: configuration error %s (%d)\n",
	   elem[loc].Name, loc);
    exit(1);
  }
}


void set_bn_s1(const int Fnum, const double ds)
{
  int k;

  for (k = 1; k <= get_n_Kids(Fnum); k++)
    set_bn_s1(Fnum, k, ds);
}


void set_s1_par(const int Fnum, const int Knum, const int j)
{
  char   name[name_length];
  int    loc, loc_d;
  double L;

  const bool prt = false;

  // Point to multipole.
  loc = get_loc(Fnum, Knum) - 1;

  if (prt) printf("\nset_s1_par: %s, %d\n", elem[loc].Name, Knum);

  if (elem[loc-1].Name[1] == 'u') {
    strcpy(name, elem[loc-1].Name); name[1] = 'd';
    loc_d = get_loc(get_Fnum(name), Knum) - 1;
    L = elem_tps[loc-1].L.cst(); elem_tps[loc-1].L = tps(L, j);
    L = elem_tps[loc_d].L.cst(); elem_tps[loc_d].L = -tps(-L, j);
    if (prt)
      printf("set_s1_par: %s %s\n", elem_tps[loc-1].Name, elem[loc_d].Name);
  } else if (elem[loc-1].Name[1] == 'd') {
    strcpy(name, elem[loc-1].Name); name[1] = 'u';
    loc_d = get_loc(get_Fnum(name), Knum) - 1;
    L = elem_tps[loc-1].L.cst(); elem_tps[loc-1].L = -tps(-L, j);
    L = elem_tps[loc_d].L.cst(); elem_tps[loc_d].L = tps(L, j);
    if (prt)
      printf("set_s1_par: %s %s\n", elem_tps[loc-1].Name, elem[loc_d].Name);
  } else if (elem[loc+1].Name[1] == 'd') {
    strcpy(name, elem[loc+1].Name); name[1] = 'u';
    loc_d = get_loc(get_Fnum(name), Knum) - 1;
    L = elem_tps[loc+1].L.cst(); elem_tps[loc+1].L = -tps(-L, j);
    L = elem_tps[loc_d].L.cst(); elem_tps[loc_d].L = tps(L, j);
    if (prt)
      printf("set_s1_par: %s %s\n", elem_tps[loc+1].Name, elem[loc_d].Name);
  } else if (elem[loc+1].Name[1] == 'u') {
    strcpy(name, elem[loc+1].Name); name[1] = 'd';
    loc_d = get_loc(get_Fnum(name), Knum) - 1;
    L = elem_tps[loc+1].L.cst(); elem_tps[loc+1].L = tps(L, j);
    L = elem_tps[loc_d].L.cst(); elem_tps[loc_d].L = -tps(-L, j);
    if (prt)
      printf("set_s1_par: %s %s\n", elem_tps[loc+1].Name, elem[loc_d].Name);
  } else {
    printf("\nset_s1_par: configuration error %s (%d)\n",
	   elem[loc].Name, loc);
    exit(1);
  }
}


void set_s1_par(const int Fnum, const int j)
{
  int k;

  for (k = 1; k <= get_n_Kids(Fnum); k++)
    set_s1_par(Fnum, k, j);
}


void clr_s1_par(const int Fnum, const int Knum)
{
  char   name[name_length];
  int    loc, loc_d;
  double L;

  const bool  prt = false;

  // Point to multipole.
  loc = get_loc(Fnum, Knum) - 1;

  if (prt) printf("\nclr_s1_par: %s %d\n", elem[loc].Name, Knum);

  if (elem[loc-1].Name[1] == 'u') {
    strcpy(name, elem[loc-1].Name); name[1] = 'd';
    loc_d = get_loc(get_Fnum(name), Knum) - 1;
    L = elem_tps[loc-1].L.cst(); elem_tps[loc-1].L = L;
    L = elem_tps[loc_d].L.cst(); elem_tps[loc_d].L = L;
    if (prt)
      printf("%s %s\n", elem_tps[loc-1].Name, elem[loc_d].Name);
  } else if (elem[loc-1].Name[1] == 'd') {
    strcpy(name, elem[loc-1].Name); name[1] = 'u';
    loc_d = get_loc(get_Fnum(name), Knum) - 1;
    L = elem_tps[loc-1].L.cst(); elem_tps[loc-1].L = L;
    L = elem_tps[loc_d].L.cst(); elem_tps[loc_d].L = L;
    if (prt)
      printf("%s %s\n", elem_tps[loc-1].Name, elem[loc_d].Name);
  } else if (elem[loc+1].Name[1] == 'd') {
    strcpy(name, elem[loc+1].Name); name[1] = 'u';
    loc_d = get_loc(get_Fnum(name), Knum) - 1;
    L = elem_tps[loc+1].L.cst(); elem_tps[loc+1].L = L;
    L = elem_tps[loc_d].L.cst(); elem_tps[loc_d].L = L;
    if (prt)
      printf("%s %s\n", elem_tps[loc+1].Name, elem[loc_d].Name);
  } else if (elem[loc+1].Name[1] == 'u') {
    strcpy(name, elem[loc+1].Name); name[1] = 'd';
    loc_d = get_loc(get_Fnum(name), Knum) - 1;
    L = elem_tps[loc+1].L.cst(); elem_tps[loc+1].L = L;
    L = elem_tps[loc_d].L.cst(); elem_tps[loc_d].L = L;
    if (prt)
      printf("%s %s\n", elem_tps[loc+1].Name, elem[loc_d].Name);
  } else {
    printf("\nset_s1_par: configuration error %s (%d)\n",
	   elem[loc].Name, loc);
    exit(1);
  }
}


void clr_s1_par(const int Fnum)
{
  int k;

  for (k = 1; k <= get_n_Kids(Fnum); k++)
    clr_s1_par(Fnum, k);
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
    A1.propagate(j+1, j+1);
    elem_tps[j].A1 = get_A_CS(2, A1, dnu2);

    // Store linear optics for convenience.
    get_ab(A1, alpha1, beta1, dnu2, eta1, etap1);
    for (k = 0; k < 2; k++) {
      elem[j].Alpha[k] = alpha1[k]; elem[j].Beta[k] = beta1[k];
      elem[j].Eta[k] = eta1[k]; elem[j].Etap[k] = etap1[k];
    }
    // Assumes dnu < 360 degrees.
    for (k = 0; k < 2; k++) {
      if (j == 0)
	elem[j].Nu[k] = 0e0;
      else {
	elem[j].Nu[k] = floor(elem[j-1].Nu[k]) + dnu2[k];
	if ((dnu2[k] < dnu1[k]) && (elem[j].L >= 0e0)) elem[j].Nu[k] += 1e0;
      }
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
  get_twiss(0, n_elem, A1);
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


void prt_system(const int m, const int n_b2, double **A, double *b)
{
  int i, j;

  printf("\n Ax = b:\n");
  for (j = 1; j <= n_b2; j++)
    if (j == 1)
      printf("%11d", j);
    else
      printf("%11d", j);
  printf("\n");
  for (i = 1; i <= m; i++) {
    printf("%4d", i);
    for (j = 1; j <= n_b2; j++)
      printf("%11.3e", A[i][j]);
    printf("%11.3e\n", b[i]);
  }
}


void prt_b2(const param_type &b2_prms)
{
  int  k;
  FILE *outf;

  std::string file_name = "match.out";

  outf = file_write(file_name.c_str());

  k = 0;
  fprintf(outf, "Q01:  quadrupole, l = 0.234, k = %8.5f, N = Nquad"
	  ", Method = Meth;\n", bn_prms.bn[++k]);
  fprintf(outf, "Q03:  quadrupole, l = 0.434, k = %8.5f, N = Nquad"
	  ", Method = Meth;\n", bn_prms.bn[++k]);

  fprintf(outf, "\nEQ01: quadrupole, l = 0.234, k = %8.5f, N = Nquad"
	  ", Method = Meth;\n", bn_prms.bn[++k]);
  fprintf(outf, "EQ02: quadrupole, l = 0.234, k = %8.5f, N = Nquad"
	  ", Method = Meth;\n", bn_prms.bn[++k]);
  fprintf(outf, "Q02:  quadrupole, l = 0.234, k = %8.5f, N = Nquad"
	  ", Method = Meth;\n", bn_prms.bn[++k]);

  fprintf(outf, "EQ04: quadrupole, l = 0.234, k = %8.5f, N = Nquad"
	  ", Method = Meth;\n", bn_prms.bn[++k]);
  fprintf(outf, "EQ05: quadrupole, l = 0.234, k = %8.5f, N = Nquad"
	  ", Method = Meth;\n", bn_prms.bn[++k]);
  fprintf(outf, "EQ06: quadrupole, l = 0.234, k = %8.5f, N = Nquad"
	  ", Method = Meth;\n", bn_prms.bn[++k]);

  // fprintf(outf, "\nD_Q01_L  = %8.5f;\n", bn_prms.bn[++k]);
  // fprintf(outf, "D_Q03_L  = %8.5f;\n", bn_prms.bn[++k]);

  // fprintf(outf, "\nD_EQ01_L = %8.5f;\n", bn_prms.bn[++k]);
  // fprintf(outf, "D_EQ02_L = %8.5f;\n", bn_prms.bn[++k]);
  // fprintf(outf, "D_Q02_L  = %8.5f;\n", bn_prms.bn[++k]);

  // fprintf(outf, "D_EQ04_L = %8.5f;\n", bn_prms.bn[++k]);
  // fprintf(outf, "D_EQ05_L = %8.5f;\n", bn_prms.bn[++k]);
  // fprintf(outf, "D_EQ06_L = %8.5f;\n", bn_prms.bn[++k]);

  // fprintf(outf, "\nD_B10_L  = %8.5f;\n", bn_prms.bn[++k]);

  // fprintf(outf, "\nU561: drift, L = %8.5f;\n", b2[18]);

  fclose(outf);
}


void prt_lin_opt(void)
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


void prt_lev_marq(const int m, const int n)
{
  int          i, loc1, loc2;
  ss_vect<tps> Ascr;

  prt_system(m, n, A_lm, f_lm);

  prt_b2(bn_prms);

  printf("\n%d bn:\n", n_powell);
  for (i = 1; i <= bn_prms.n_prm; i++) {
    bn_prms.bn[i] = get_bn(bn_prms.Fnum[i-1], 1, bn_prms.n[i-1]);
    printf("%11.3e", bn_prms.bn[i]);
    if (i % n_prt == 0) printf("\n");
    // Bounded.
    bn_prms.bn[i] =
      bn_internal(bn_prms.bn[i]/bn_prms.bn_scl[i-1],
    		  bn_prms.bn_min[i-1], bn_prms.bn_max[i-1]);
  }
  if (bn_prms.n_prm % n_prt != 0) printf("\n");

  Ascr = get_A(ic[0], ic[1], ic[2], ic[3]);
  get_twiss(loc[0], loc[4], Ascr);

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

  loc1 = get_loc(get_Fnum("s_s_1"), 1) - 1;
  loc2 = get_loc(get_Fnum("s_s_1"), 2) - 1;
  printf("\nLength of 1st straight: %6.3f m\n", elem[loc2].S-elem[loc1].S);
  loc1 = get_loc(get_Fnum("s_s_2"), 1) - 1;
  loc2 = get_loc(get_Fnum("s_s_2"), 2) - 1;
  printf("Length of 2nd straight: %6.3f m\n", elem[loc2].S-elem[loc1].S);
  loc1 = get_loc(get_Fnum("s_s_3"), 1) - 1;
  loc2 = get_loc(get_Fnum("s_s_3"), 2) - 1;
  printf("Length of 3rd straight: %6.3f m\n", elem[loc2].S-elem[loc1].S);

  prt_mfile("flat_file.fit");
  prt_lat(loc[0], loc[4], "linlat1.out");
  prt_lat(loc[0], loc[4], "linlat.out", 10);
}


void get_f_grad(const int n_bn, double *f, double **A, double &chi2, int &m)
{
  int          i, j;
  tps          h_re, h_im, K_re, K_im, K_re_scl;
  ss_vect<tps> Ascr, AA_tp[3], A_disp[3];

  // printf("\n");
  for (i = 1; i <= n_bn; i++) {
    bn_prms.set_prm_dep(i-1);

    Ascr = get_A(ic[0], ic[1], ic[2], ic[3]);
    get_twiss(loc[0], loc[4], Ascr);

    for (j = 0; j < 3; j++) {
      AA_tp[j]  = elem_tps[loc[j+1]].A1*tp_S(2, elem_tps[loc[j+1]].A1);
      A_disp[j] = elem_tps[loc[j+1]].A1;
    }

    m = 0;
    A[++m][i] = get_a(scl_eta,    A_disp[0][x_],  0, 0, 0, 0, 1);
    A[++m][i] = get_a(scl_eta,    A_disp[0][px_], 0, 0, 0, 0, 1);
    A[++m][i] = get_a(scl_alpha, -AA_tp[1][x_],   0, 1, 0, 0, 0);
    A[++m][i] = get_a(scl_alpha, -AA_tp[1][y_],   0, 0, 0, 1, 0);
    A[++m][i] = get_a(scl_beta,   AA_tp[1][x_],   1, 0, 0, 0, 0);
    A[++m][i] = get_a(scl_alpha, -AA_tp[2][x_],   0, 1, 0, 0, 0);
    A[++m][i] = get_a(scl_alpha, -AA_tp[2][y_],   0, 0, 0, 1, 0);
    A[++m][i] = get_a(scl_beta,   AA_tp[2][x_],   1, 0, 0, 0, 0);

    for (j = 1; j <= m; j++)
      A[j][i] *= bn_prms.bn_scl[i-1];

    bn_prms.clr_prm_dep(i-1);
  }

  m = 0;
  f[++m] = get_b(scl_eta,    A_disp[0][x_],  0, 0, 0, 0, 1);
  f[++m] = get_b(scl_eta,    A_disp[0][px_], 0, 0, 0, 0, 1);
  f[++m] = get_b(scl_alpha, -AA_tp[1][x_],   0, 1, 0, 0, 0);
  f[++m] = get_b(scl_alpha, -AA_tp[1][y_],   0, 0, 0, 1, 0);
  f[++m] = get_b(scl_beta,   AA_tp[1][x_],   1, 0, 0, 0, 0) - scl_beta*9.58;
  f[++m] = get_b(scl_alpha, -AA_tp[2][x_],   0, 1, 0, 0, 0);
  f[++m] = get_b(scl_alpha, -AA_tp[2][y_],   0, 0, 0, 1, 0);
  f[++m] = get_b(scl_beta,   AA_tp[2][x_],   1, 0, 0, 0, 0) - scl_beta*8.0;

  chi2 = 0e0;
  for (j = 1; j <= m; j++)
    chi2 += sqr(f[j]);
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

    // Don't change best values.
    printf("\nget_f_der:\n");
    for (i = 1; i <= n; i++) {
      set_prm(bn_prms.Fnum[i-1], bn_prms.n[i-1], bn_prms.bn_scl[i-1]*bn[i],
	      bn_prms.bn_min[i-1], bn_prms.bn_max[i-1]);
      printf(" %12.5e", bn_prms.bn_scl[i-1]*bn_prms.bn[i]);
      if (i % n_prt == 0) printf("\n");
    }
    if (n % n_prt != 0) printf("\n");

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
  int          n_data, i, n, *ia;
  double       *x, *y, *sigma, **covar, **alpha, chisq, alambda, alambda0;
  ss_vect<tps> Ascr;

  const int n_bn = bn_prms.n_prm, n_iter = 500;

  n_data = 8;

  ia = ivector(1, n_bn);
  x = dvector(1, n_data); y = dvector(1, n_data); sigma = dvector(1, n_data);
  covar = dmatrix(1, n_bn, 1, n_bn); alpha = dmatrix(1, n_bn, 1, n_bn);
  f_lm = dvector(1, n_data); A_lm = dmatrix(1, n_data, 1, n_bn);

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
  get_twiss(loc[0], loc[4], Ascr);

  prt_lin_opt();

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
  } while (n < n_iter);

  alambda = 0e0;
  dmrqmin(x, y, sigma, n_data, bn_prms.bn, ia, n_bn,  covar, alpha, &chisq,
	  get_f_der, &alambda);

  free_dvector(f_lm, 1, n_data); free_dmatrix(A_lm, 1, n_data, 1, n_bn);
  free_ivector(ia, 1, n_bn);
  free_dvector(x, 1, n_data); free_dvector(y, 1, n_data);
  free_dvector(sigma, 1, n_data);
  free_dmatrix(covar, 1, n_bn, 1, n_bn); free_dmatrix(alpha, 1, n_bn, 1, n_bn);
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

  if (false) {
    get_twiss();
    prt_lat("linlat1.out");
    prt_lat("linlat.out", 10);
    exit(0);
  }
  
  bn_prms.add_prm("q01",  2, -5.0, 5.0, 1.0);
  bn_prms.add_prm("q03",  2, -5.0, 5.0, 1.0);

  bn_prms.add_prm("eq01", 2, -5.0, 5.0, 1.0);
  bn_prms.add_prm("eq02", 2, -5.0, 5.0, 1.0);
  bn_prms.add_prm("q02",  2, -5.0, 5.0, 1.0);

  bn_prms.add_prm("eq04", 2, -5.0, 5.0, 1.0);
  bn_prms.add_prm("eq05", 2, -5.0, 5.0, 1.0);
  bn_prms.add_prm("eq06", 2, -5.0, 5.0, 1.0);

  bn_prms.add_prm("q01",  -2,  0.0,  0.05, 1e-2);
  // bn_prms.add_prm("q03",  -2,  0.0,  0.05, 1e-2);

  // bn_prms.add_prm("eq01", -2,  0.0,  0.05, 1e-2);
  // bn_prms.add_prm("eq02", -2,  0.0,  0.05, 1e-2);
  // bn_prms.add_prm("q02",  -2,  0.0,  0.05, 1e-2);

  // bn_prms.add_prm("eq04", -2, -0.05, 0.05, 1e-2);
  // bn_prms.add_prm("eq05", -2,  0.0,  0.05, 1e-2);
  // bn_prms.add_prm("eq06", -2,  0.0,  0.05, 1e-2);

  // bn_prms.add_prm("b10",  -2, -0.01, 0.01, 1e-2);

  // U561 + U562: 2.14.
  // bn_prms.add_prm("u561", -1, 2.14, 2.14, 1.0);

  bn_prms.ini_prm();

  min_lev_marq();
}
