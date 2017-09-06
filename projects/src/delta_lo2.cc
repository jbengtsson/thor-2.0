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
// const double ic[][2] =
//   {{1.15199, -0.22236}, {0.65878, 5.53043 }, {0.03741, 0.0}, {-0.04304, 0.0}};
// Upstream of QD04.
// Periodic.
// const double ic[][2] =
//   {{4.10172, -2.96036}, {2.69522, 4.71532}, {0.12263, 0.0}, {-0.16484, 0.0}};
// Upstream of QF03.
// Periodic.
const double ic[][2] =
  {{-5.66627, 2.37465}, {7.07181, 2.86889}, {0.19967, 0.0}, {0.17750, 0.0}};

int    loc[10], n, n_strength;
double chi2 = 0e0, *f_lm, **A_lm;

const int n_prt = 8;

const double scl_eta = 1e4, scl_beta = 1e0, scl_alpha = 1e3, scl_bn = 1e2;

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
  std::vector<double> bn_min, bn_max, bn_scl;
  std::vector<int>    Fnum, n;

  void add_prm(const std::string Fname, const int n,
	       const double bn_min, const double bn_max, const double bn_scl);
  void ini_prm(double *bn);
  void set_prm(double *bn) const;
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


void param_type::ini_prm(double *bn)
{
  int i;

  n_prm = Fnum.size();

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


void param_type::set_prm(double *bn) const
{
  int    i;
  double bn_ext;

  printf("\nset_prm:\n");
  for (i = 1; i <= n_prm; i++) {
    // Bounded.
    bn_ext = bn_bounded(bn_scl[i-1]*bn[i], bn_min[i-1], bn_max[i-1]);
    if (n[i-1] > 0)
      set_bn(Fnum[i-1], n[i-1], bn_ext);
    else if (n[i-1] == -1)
      set_L(Fnum[i-1], bn_ext);
    else if (n[i-1] == -2)
      set_bn_s1(Fnum[i-1], bn_ext);
    printf(" %9.5f", bn_ext);
    if (i % n_prt == 0) printf("\n");
  }
  if (n_prm % n_prt != 0) printf("\n");
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

void get_s_loc(const int Fnum, const int Knum, int loc[])
{
  char name[name_length];
  
  // Point to multipole.
  loc[1] = get_loc(Fnum, Knum) - 1;
  if (elem[loc[1]-1].Name[1] == 'u') {
    loc[0] = loc[1] - 1;
    strcpy(name, elem[loc[1]-1].Name); name[1] = 'd';
    loc[2] = get_loc(get_Fnum(name), Knum) - 1;
  } else if (elem[loc[1]-1].Name[1] == 'd') {
    loc[2] = loc[1] - 1;
    strcpy(name, elem[loc[1]-1].Name); name[1] = 'u';
    loc[0] = get_loc(get_Fnum(name), Knum) - 1;
  } else if (elem[loc[1]+1].Name[1] == 'd') {
    loc[2] = loc[1] + 1;
    strcpy(name, elem[loc[1]+1].Name); name[1] = 'u';
    loc[0] = get_loc(get_Fnum(name), Knum) - 1;
  } else if (elem[loc[1]+1].Name[1] == 'u') {
    loc[0] = loc[1] + 1;
    strcpy(name, elem[loc[1]+1].Name); name[1] = 'd';
    loc[2] = get_loc(get_Fnum(name), Knum) - 1;
  } else {
    printf("\nget_s_loc: configuration error %s (%d)\n",
	   elem[loc[1]].Name, loc[1]);
    exit(1);
  }
}



double get_bn_s1(const int Fnum, const int Knum)
{
  
  int    loc[3];
  double ds;

  const bool prt = false;

  get_s_loc(Fnum, Knum, loc);
  ds = elem[loc[0]].L;
  if (prt)
    printf("\nget_bn_s1: %s %s(%d) %s %10.3e %10.3e\n",
	   elem_tps[loc[0]].Name, elem[loc[1]].Name, Knum, elem[loc[2]].Name,
	   elem[loc[0]].L, elem[loc[2]].L);
 
  return ds;
}


void set_bn_s1(const int Fnum, const int Knum, const double ds)
{
  int loc[3];

  const bool prt = false;

  get_s_loc(Fnum, Knum, loc);
  set_L(elem[loc[0]].Fnum, Knum, ds);
  set_L(elem[loc[2]].Fnum, Knum, -ds);
  if (prt)
    printf("\nset_bn_s1: %s %s(%d) %s %10.3e %10.3e\n",
	   elem[loc[0]].Name, elem[loc[1]].Name, Knum, elem[loc[2]].Name,
	   ds, -ds);
}


void set_bn_s1(const int Fnum, const double ds)
{
  int k;

  for (k = 1; k <= get_n_Kids(Fnum); k++)
    set_bn_s1(Fnum, k, ds);
}


void set_s1_par(const int Fnum, const int Knum, const int j)
{
  int    loc[3];
  double L;

  const bool prt = false;

  get_s_loc(Fnum, Knum, loc);
  L = elem_tps[loc[0]].L.cst(); elem_tps[loc[0]].L = tps(L, j);
  L = elem_tps[loc[2]].L.cst(); elem_tps[loc[2]].L = -tps(-L, j);
  if (prt)
    printf("\nset_s1_par: %s %s(%d) %s %10.3e %10.3e\n",
	   elem_tps[loc[0]].Name, elem[loc[1]].Name, Knum, elem[loc[2]].Name,
	   elem_tps[loc[0]].L.cst(), elem_tps[loc[2]].L.cst());
}


void set_s1_par(const int Fnum, const int j)
{
  int k;

  for (k = 1; k <= get_n_Kids(Fnum); k++)
    set_s1_par(Fnum, k, j);
}


void clr_s1_par(const int Fnum, const int Knum)
{
  int    loc[3];
  double L;

  const bool prt = false;

  get_s_loc(Fnum, Knum, loc);
  L = elem_tps[loc[0]].L.cst(); elem_tps[loc[0]].L = L;
  L = elem_tps[loc[2]].L.cst(); elem_tps[loc[2]].L = L;
  if (prt)
    printf("\nclr_s1_par: %s %s(%d) %s %10.3e %10.3e\n",
	   elem_tps[loc[0]].Name, elem[loc[1]].Name, Knum, elem[loc[2]].Name,
	   elem_tps[loc[0]].L.cst(), elem_tps[loc[2]].L.cst());
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


void get_der(const char *prm_name)
{
  // Numerical check of der for displacement.
  int          loc, loc0, loc1, Fnum;
  double       beta_x[3], s;
  ss_vect<tps> Ascr, AA_tp;

  const double eps = 1e-10;
  
  loc = get_loc(get_Fnum("ef2"), 16) - 1;

  loc0 = get_loc(get_Fnum("qf031"), 1) - 1;
  loc1 = get_loc(get_Fnum("b20"), 5) - 1;
  Ascr = get_A(ic[0], ic[1], ic[2], ic[3]);
  get_twiss(loc0, loc1, Ascr);

  AA_tp = elem_tps[loc].A1*tp_S(2, elem_tps[loc].A1);
  beta_x[1] = h_ijklm(AA_tp[x_], 1, 0, 0, 0, 0);

  Fnum = get_Fnum(prm_name);

  s = get_bn_s1(Fnum, 1);

  set_bn_s1(Fnum, s-eps);

  Ascr = get_A(ic[0], ic[1], ic[2], ic[3]);
  get_twiss(loc0, loc1, Ascr);

  AA_tp = elem_tps[loc].A1*tp_S(2, elem_tps[loc].A1);
  beta_x[0] = h_ijklm(AA_tp[x_], 1, 0, 0, 0, 0);

  set_bn_s1(Fnum, s+eps);

  Ascr = get_A(ic[0], ic[1], ic[2], ic[3]);
  get_twiss(loc0, loc1, Ascr);

  AA_tp = elem_tps[loc].A1*tp_S(2, elem_tps[loc].A1);
  beta_x[2] = h_ijklm(AA_tp[x_], 1, 0, 0, 0, 0);

  printf("der = %10.3e\n", (beta_x[2]-beta_x[0])/(2e0*eps));
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


void prt_b2(const param_type &b2_prms, double *b2)
{
  int  k;
  FILE *outf;

  std::string file_name = "match.out";

  outf = file_write(file_name.c_str());

  k = 1;
  fprintf(outf, "QF031: quadrupole, l = 0.217, k = %8.5f, N = Nquad"
	  ", Method = Meth;\n",
	  bn_bounded(b2[k], b2_prms.bn_min[k-1], b2_prms.bn_max[k-1]));
  k++;
  fprintf(outf, "QD041: quadrupole, l = 0.117, k = %8.5f, N = Nquad"
	  ", Method = Meth;\n",
	  bn_bounded(b2[k], b2_prms.bn_min[k-1], b2_prms.bn_max[k-1]));

  k++;
  fprintf(outf, "\nQ01:   quadrupole, l = 0.234, k = %8.5f, N = Nquad"
	  ", Method = Meth;\n",
	  bn_bounded(b2[k], b2_prms.bn_min[k-1], b2_prms.bn_max[k-1]));
  k++;
  fprintf(outf, "Q02:   quadrupole, l = 0.434, k = %8.5f, N = Nquad"
	  ", Method = Meth;\n",
	  bn_bounded(b2[k], b2_prms.bn_min[k-1], b2_prms.bn_max[k-1]));
  // k++;
  // fprintf(outf, "Q03:   quadrupole, l = 0.234, k = %8.5f, N = Nquad"
  // 	  ", Method = Meth;\n",
  // 	  bn_bounded(b2[k], b2_prms.bn_min[k-1], b2_prms.bn_max[k-1]));

  k++;
  fprintf(outf, "\nEQ01:  quadrupole, l = 0.234, k = %8.5f, N = Nquad"
	  ", Method = Meth;\n",
	  bn_bounded(b2[k], b2_prms.bn_min[k-1], b2_prms.bn_max[k-1]));
  k++;
  fprintf(outf, "EQ02:  quadrupole, l = 0.234, k = %8.5f, N = Nquad"
	  ", Method = Meth;\n",
	  bn_bounded(b2[k],
		     b2_prms.bn_min[k-1], b2_prms.bn_max[k-1]));

  k++;
  fprintf(outf, "\nEQ03:  quadrupole, l = 0.234, k = %8.5f, N = Nquad"
	  ", Method = Meth;\n",
	  bn_bounded(b2[k], b2_prms.bn_min[k-1], b2_prms.bn_max[k-1]));
  k++;
  fprintf(outf, "EQ04:  quadrupole, l = 0.234, k = %8.5f, N = Nquad"
	  ", Method = Meth;\n",
	  bn_bounded(b2[k], b2_prms.bn_min[k-1], b2_prms.bn_max[k-1]));
  k++;
  fprintf(outf, "EQ05:  quadrupole, l = 0.234, k = %8.5f, N = Nquad"
	  ", Method = Meth;\n",
	  bn_bounded(b2[k], b2_prms.bn_min[k-1], b2_prms.bn_max[k-1]));

  if (true) {
    k++;
    fprintf(outf, "\nD_Q01_L  = %8.5f;\n",
	    bn_bounded(b2[k], b2_prms.bn_min[k-1], b2_prms.bn_max[k-1]));
    k++;
    fprintf(outf, "D_Q02_L  = %8.5f;\n",
	    bn_bounded(b2[k], b2_prms.bn_min[k-1], b2_prms.bn_max[k-1]));
    // k++;
    // fprintf(outf, "D_Q03_L  = %8.5f;\n",
    // 	    bn_bounded(b2[k], b2_prms.bn_min[k-1], b2_prms.bn_max[k-1]));

    k++;
    fprintf(outf, "\nD_EQ01_L = %8.5f;\n",
	    bn_bounded(b2[k], b2_prms.bn_min[k-1], b2_prms.bn_max[k-1]));
    k++;
    fprintf(outf, "D_EQ02_L = %8.5f;\n",
	    bn_bounded(b2[k], b2_prms.bn_min[k-1], b2_prms.bn_max[k-1]));

    k++;
    fprintf(outf, "\nD_EQ03_L = %8.5f;\n",
	    bn_bounded(b2[k], b2_prms.bn_min[k-1], b2_prms.bn_max[k-1]));
    k++;
    fprintf(outf, "D_EQ04_L = %8.5f;\n",
	    bn_bounded(b2[k], b2_prms.bn_min[k-1], b2_prms.bn_max[k-1]));
    k++;
    fprintf(outf, "D_EQ05_L = %8.5f;\n",
	    bn_bounded(b2[k], b2_prms.bn_min[k-1], b2_prms.bn_max[k-1]));

    // k++;
    // fprintf(outf, "\nD_B10_L  = %8.5f;\n",
    // 	    bn_bounded(b2[k], b2_prms.bn_min[k-1], b2_prms.bn_max[k-1]));

    // k++;
    // fprintf(outf, "\nU561: drift, L = %8.5f;\n",
    // 	    bn_bounded(b2[18], b2_prms.bn_min[k-1], b2_prms.bn_max[k-1]));
  }

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


void prt_lev_marq(const int m, const int n, double *b2)
{
  int loc1, loc2;

  prt_system(m, n, A_lm, f_lm);

  printf("\n%d bn:\n", n_powell);
  prt_b2(bn_prms, b2);

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
  prt_lat(loc[0]+1, loc[4]-1, "linlat.out", 10);
}


void get_f_grad(const int n_bn, double *b2, double *f, double **A,
		double &chi2, int &m)
{
  int          i, j;
  tps          h_re, h_im, K_re, K_im, K_re_scl;
  ss_vect<tps> Ascr, AA_tp[3], A_disp[3];

  bn_prms.set_prm(b2);

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

    for (j = 1; j <= n_strength; j++)
      A[++m][i] = (i == j)? scl_bn*get_L(bn_prms.Fnum[j-1], 1) : 0e0;

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

  for (j = 1; j <= n_strength; j++)
    f[++m] = scl_bn*b2[j]*get_L(bn_prms.Fnum[j-1], 1);

  chi2 = 0e0;
  for (j = 1; j <= m; j++)
    chi2 += sqr(f[j]);
}


void get_f_der(double x, double *b2, double *yfit, double *dyda, int n)
{
  int        i, m1;
  static int m;
  double     chi2;

  m1 = (int)(x+0.5);

  if (m1 == 1) {
    n_powell++;
    get_f_grad(n, b2, f_lm, A_lm, chi2, m);
  }

  *yfit = f_lm[m1];

  for (i = 1; i <= n; i++)
    dyda[i] = A_lm[m1][i];
}


void min_lev_marq(void)
{
  int          n_data, i, n, *ia;
  double       *b2, *x, *y, *sigma, **covar, **alpha, chisq, alambda, alambda0;
  ss_vect<tps> Ascr;

  const int n_bn = bn_prms.n_prm, n_iter = 50;

  n_data = 8 + n_strength;

  b2 = dvector(1, n_bn);
  ia = ivector(1, n_bn);
  x = dvector(1, n_data); y = dvector(1, n_data); sigma = dvector(1, n_data);
  covar = dmatrix(1, n_bn, 1, n_bn); alpha = dmatrix(1, n_bn, 1, n_bn);
  f_lm = dvector(1, n_data); A_lm = dmatrix(1, n_data, 1, n_bn);

  // Upstream of QF03.
  loc[0] = get_loc(get_Fnum("qf031"), 1) - 1;
  // Upstream of QD04.
  // loc[0] = get_loc(get_Fnum("qd041"), 1) - 1;
  // Upstream of 20 degree dipole.
  // loc[0] = get_loc(get_Fnum("sb"), 7) - 1;
   // Upstream of 20 degree dipole.
  // loc[0] = get_loc(get_Fnum("sb"),  7) - 1;

  // Downstream of 10 degree dipole.
  loc[1] = get_loc(get_Fnum("b10"), 1) - 1;
  // Center of 1st straight.
  loc[2] = get_loc(get_Fnum("ef2"), 4) - 1;
  // Center of 2nd straight.
  loc[3] = get_loc(get_Fnum("ef2"), 16) - 1;

  // Downstream of 20 degree dipole.
  // loc[4] = get_loc(get_Fnum("b20"), 5) - 1;
  // Downstream of QF03.
  loc[4] = get_loc(get_Fnum("qf031"), 4) - 1;

  Ascr = get_A(ic[0], ic[1], ic[2], ic[3]);
  get_twiss(loc[0], loc[4], Ascr);

  prt_lin_opt();

  bn_prms.ini_prm(b2);

  get_f_grad(n_bn, b2, f_lm, A_lm, chi2, n_data);
  prt_system(n_data, n_bn, A_lm, f_lm);

  for (i = 1; i <= n_bn; i++)
    ia[i] = 1;

  for (i = 1; i <= n_data; i++) {
    sigma[i] = 1e0; x[i] = i; y[i] = 0e0;
  }

  alambda = -1e0; alambda0 = 1e-3;
  dmrqmin(x, y, sigma, n_data, b2, ia, n_bn, covar, alpha, &chisq,
	  get_f_der, &alambda);
  printf("\nalambda = %7.1e, chi2 = %9.3e\n", alambda, chisq);
  if (alambda < alambda0) prt_lev_marq(n_data, n_bn, b2);
  alambda0 = alambda;

  n = 0;
  do {
    n++;
    dmrqmin(x, y, sigma, n_data, b2, ia, n_bn,  covar, alpha, &chisq,
	    get_f_der, &alambda);
    printf("\nalambda = %7.1e, chi2 = %9.3e\n", alambda, chisq);
    if (alambda < alambda0) prt_lev_marq(n_data, n_bn, b2);
    alambda0 = alambda;
  } while (n < n_iter);

  alambda = 0e0;
  dmrqmin(x, y, sigma, n_data, b2, ia, n_bn,  covar, alpha, &chisq,
	  get_f_der, &alambda);

  free_dvector(b2, 1, n_bn);
  free_dvector(f_lm, 1, n_data);
  free_dmatrix(A_lm, 1, n_data, 1, n_bn);
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
  
  if (false) {
    printf("\n");
    get_der("q01");
    get_der("q03");
    get_der("eq01");
    get_der("eq02");
    // exit(0);
  }

  bn_prms.add_prm("qf031",  2, -4.2, 4.2, 1.0);
  bn_prms.add_prm("qd041",  2, -4.2, 4.2, 1.0);

  bn_prms.add_prm("q01",    2, -4.2, 4.2, 1.0);
  // bn_prms.add_prm("q02",    2, -4.2, 4.2, 1.0);
  bn_prms.add_prm("q03",    2, -4.2, 4.2, 1.0);

  bn_prms.add_prm("eq01",   2, -4.2, 4.2, 1.0);
  bn_prms.add_prm("eq02",   2, -4.2, 4.2, 1.0);

  bn_prms.add_prm("eq03",   2, -4.2, 4.2, 1.0);
  bn_prms.add_prm("eq04",   2, -4.2, 4.2, 1.0);
  bn_prms.add_prm("eq05",   2, -4.2, 4.2, 1.0);

  n_strength = 9;

  if (false) {
    bn_prms.add_prm("q01",  -2,  0.0,  0.05, 1e0);
    // bn_prms.add_prm("q02",  -2,  0.0,  0.05, 1e0);
    bn_prms.add_prm("q03",  -2,  0.0,  0.05, 1e0);

    bn_prms.add_prm("eq01", -2,  0.0,  0.05, 1e0);
    bn_prms.add_prm("eq02", -2,  0.0,  0.05, 1e0);

    bn_prms.add_prm("eq03", -2, -0.05, 0.05, 1e0);
    bn_prms.add_prm("eq04", -2,  0.0,  0.05, 1e0);
    bn_prms.add_prm("eq05", -2,  0.0,  0.05, 1e0);

    // bn_prms.add_prm("b10",  -2, -0.02, 0.02, 1e0);
  }

  // U561 + U562: 2.14.
  // bn_prms.add_prm("u561", -1, 2.14, 2.14, 1.0);

  min_lev_marq();
}
