// Optimise higher-order-achromat:
//
//   1. Optimise sextupoles - beware of singular values.
//   3. Optimise octupoles - for obtained sextupole strengths.
//   4. Optimise both sextupoles & octupoles.


#include <algorithm>

#include <assert.h>

#define NO 9

#include "thor_lib.h"

// User defined class for parameter dependence.
#include "param_type.h"

int
  no_tps   = NO,
  ndpt_tps = 5;


extern double b2_max;


const bool
  b_3_opt    = !false,
  b_4_opt    = !false,
  b_3_zero   = false,
  b_4_zero   = false;

const int
  max_iter  = 100,  
  svd_n_cut = 3;

const double
  A_max[]     = {6e-3, 3e-3},
  delta_max   = 6e-2,
  beta_inj[]  = {3.7, 3.9},

  bnL_scl[]   = {0e0, 0e0, 0e0,  1e0,  5e1,    1e4},
  bnL_min[]   = {0e0, 0e0, 0e0, -5e2, -5.0e4, -1.5e5},
  bnL_max[]   = {0e0, 0e0, 0e0,  5e2,  5.0e4,  1.5e5},

#if 0
  scl_h[]     = {1e-2, 1e-2},
  scl_ksi[]   = {0e0, 1e2, 1e0, 5e0, 1e0, 1e0},
  scl_a[]     = {1e0, 1e0, 1e0, 1e0},
  scl_K_sum[] = {1e0, 1e0},
#else
  scl_h[]     = {1e0, 1e0},
  scl_ksi[]   = {0e0, 1e2, 1e-15, 1e-15, 1e-15, 1e-15},
  scl_a[]     = {1e-15, 1e-15, 1e-15, 1e-15},
  scl_K_sum[] = {1e1, 1e1},
#endif

  step       = 0.15;


class Lie_gen_class {
private:
  std::string
    label;
  std::vector<int>
    index;
  double
    cst_scl,
    cst;
  std::vector<double>
    Jacobian;
public:
  void prt_Lie_gen(void) const;

  friend Lie_gen_class get_Lie_gen
  (const std::string label, const tps &h, const double scl, const int i,
   const int j, const int k, const int l, const int m);

  friend void get_system
  (const int m, const int n, const std::vector<Lie_gen_class> &Lie_gen,
   double **A, double *b);

  friend void get_Jacobian
  (const param_type &bns, const int k, const ss_vect<tps> &Id_scl,
   std::vector<Lie_gen_class> &Lie_gen);

  friend Lie_gen_class get_Lie_gen_sum
  (const std::string &name, const double scl, const int index[],
   const std::vector<Lie_gen_class> &Lie_gen);
};


//==============================================================================
// Function in tools.cc using globval variables.

ss_vect<tps> get_map(void)
{
  // 
  ss_vect<tps> M;

  M.identity();
  M.propagate(1, n_elem);
  return M;
}


void get_ab
(const ss_vect<tps> &A1, double alpha[], double beta[], const long int k)
{
  // The map A1 is a globval variable for the function in tools.cc.
  ss_vect<tps> a1, A1_A1tp;

  a1 = A1;
  a1.propagate(1, k);

  A1_A1tp = a1*tp_S(2, a1);

  alpha[X_] = -h_ijklm(A1_A1tp[x_], 0, 1, 0, 0, 0);
  alpha[Y_] = -h_ijklm(A1_A1tp[y_], 0, 0, 0, 1, 0);
  beta[X_]  =  h_ijklm(A1_A1tp[x_], 1, 0, 0, 0, 0);
  beta[Y_]  =  h_ijklm(A1_A1tp[y_], 0, 0, 1, 0, 0);
}

//==============================================================================

void Lie_gen_class::prt_Lie_gen(void) const
{
  printf(" %s", label.c_str());
  for (auto k = 0; k < (int)index.size(); k++)
    printf("%d", index[k]);
  printf(" %6.1e %10.3e", cst_scl, cst);
  for (auto k = 0; k < (int)Jacobian.size(); k++)
    printf(" %10.3e", Jacobian[k]);
  printf("\n");
}


Lie_gen_class get_Lie_gen
(const std::string label, const tps &h, const double scl, const int i,
 const int j, const int k, const int l, const int m)
{
  const std::vector<int> index = {i, j, k, l, m};

  Lie_gen_class     Lt;
  std::ostringstream str;

  Lt.label = label;
  Lt.index = index;
  Lt.cst_scl = scl;
  Lt.cst = Lt.cst_scl*h_ijklm(h, i, j, k, l, m);
  return Lt;
}


void prt_Lie_gen
(const std::string &str, const int i0, const int n,
 const std::vector<Lie_gen_class> &Lie_gen)
{
  printf("\n%s\n", str.c_str());
  for (auto k = i0; k < i0+n; k++)
    Lie_gen[k].prt_Lie_gen();
}


std::vector<Lie_gen_class> get_Lie_gen(const ss_vect<tps> &Id_scl)
{
  double                      nu[3], ksi[3];
  tps                         g, g_re, g_im, K_re, K_im;
  ss_vect<tps>                M, A1, A0, M_res;
  Lie_gen_class              h;
  std::vector<Lie_gen_class> Lie_gen;

  danot_(no_tps-1);
  M = get_map();
  danot_(no_tps);
  auto K = MapNorm(M, g, A1, A0, M_res, 1);
  auto nus = dHdJ(K);
  get_nu_ksi(nus, nu, ksi);

  CtoR(g*Id_scl, g_re, g_im);
  CtoR(K*Id_scl, K_re, K_im);

  if (false)
    std::cout << std::scientific << std::setprecision(3)
	      << "\nK:\n" << K_re << "\ng:\n" << g_im;

  if (false)
    printf("\n  nu  = [%6.3f, %6.3f]\n  ksi = [%6.3f, %6.3f]\n",
	   nu[X_], nu[Y_], ksi[X_], ksi[Y_]);

  Lie_gen.clear();

  Lie_gen.push_back(get_Lie_gen("k_", K_re, scl_ksi[1], 1, 1, 0, 0, 1));
  Lie_gen.push_back(get_Lie_gen("k_", K_re, scl_ksi[1], 0, 0, 1, 1, 1));

  Lie_gen.push_back(get_Lie_gen("g_", g_im, scl_h[0], 1, 0, 0, 0, 2));
  Lie_gen.push_back(get_Lie_gen("g_", g_im, scl_h[0], 2, 0, 0, 0, 1));
  Lie_gen.push_back(get_Lie_gen("g_", g_im, scl_h[0], 0, 0, 2, 0, 1));

  Lie_gen.push_back(get_Lie_gen("g_", g_im, scl_h[0], 1, 0, 1, 1, 0));
  Lie_gen.push_back(get_Lie_gen("g_", g_im, scl_h[0], 2, 1, 0, 0, 0));
  Lie_gen.push_back(get_Lie_gen("g_", g_im, scl_h[0], 3, 0, 0, 0, 0));
  Lie_gen.push_back(get_Lie_gen("g_", g_im, scl_h[0], 1, 0, 0, 2, 0));
  Lie_gen.push_back(get_Lie_gen("g_", g_im, scl_h[0], 1, 0, 2, 0, 0));

  Lie_gen.push_back(get_Lie_gen("g_", g_im, scl_h[1], 4, 0, 0, 0, 0));
  Lie_gen.push_back(get_Lie_gen("g_", g_im, scl_h[1], 3, 1, 0, 0, 0));
  Lie_gen.push_back(get_Lie_gen("g_", g_im, scl_h[1], 2, 0, 2, 0, 0));
  Lie_gen.push_back(get_Lie_gen("g_", g_im, scl_h[1], 1, 1, 2, 0, 0));
  Lie_gen.push_back(get_Lie_gen("g_", g_im, scl_h[1], 2, 0, 1, 1, 0));
  Lie_gen.push_back(get_Lie_gen("g_", g_im, scl_h[1], 0, 0, 3, 1, 0));
  Lie_gen.push_back(get_Lie_gen("g_", g_im, scl_h[1], 2, 0, 0, 2, 0));
  Lie_gen.push_back(get_Lie_gen("g_", g_im, scl_h[1], 0, 0, 4, 0, 0));

  Lie_gen.push_back(get_Lie_gen("k_", K_re, scl_a[0], 2, 2, 0, 0, 0));
  Lie_gen.push_back(get_Lie_gen("k_", K_re, scl_a[0], 1, 1, 1, 1, 0));
  Lie_gen.push_back(get_Lie_gen("k_", K_re, scl_a[0], 0, 0, 2, 2, 0));

  Lie_gen.push_back(get_Lie_gen("k_", K_re, scl_a[1], 3, 3, 0, 0, 0));
  Lie_gen.push_back(get_Lie_gen("k_", K_re, scl_a[1], 2, 2, 1, 1, 0));
  Lie_gen.push_back(get_Lie_gen("k_", K_re, scl_a[1], 1, 1, 2, 2, 0));
  Lie_gen.push_back(get_Lie_gen("k_", K_re, scl_a[1], 0, 0, 3, 3, 0));

  if (NO >= 9) {
    Lie_gen.push_back(get_Lie_gen("k_", K_re, scl_a[2], 4, 4, 0, 0, 0));
    Lie_gen.push_back(get_Lie_gen("k_", K_re, scl_a[2], 3, 3, 1, 1, 0));
    Lie_gen.push_back(get_Lie_gen("k_", K_re, scl_a[2], 2, 2, 2, 2, 0));
    Lie_gen.push_back(get_Lie_gen("k_", K_re, scl_a[2], 1, 1, 3, 3, 0));
    Lie_gen.push_back(get_Lie_gen("k_", K_re, scl_a[2], 0, 0, 4, 4, 0));
  }

  if (NO >= 11) {
    Lie_gen.push_back(get_Lie_gen("k_", K_re, scl_a[3], 5, 5, 0, 0, 0));
    Lie_gen.push_back(get_Lie_gen("k_", K_re, scl_a[3], 4, 4, 1, 1, 0));
    Lie_gen.push_back(get_Lie_gen("k_", K_re, scl_a[3], 3, 3, 2, 2, 0));
    Lie_gen.push_back(get_Lie_gen("k_", K_re, scl_a[3], 2, 2, 3, 3, 0));
    Lie_gen.push_back(get_Lie_gen("k_", K_re, scl_a[3], 1, 1, 4, 4, 0));
    Lie_gen.push_back(get_Lie_gen("k_", K_re, scl_a[3], 0, 0, 5, 5, 0));
  }

  Lie_gen.push_back(get_Lie_gen("k_", K_re, scl_ksi[2], 1, 1, 0, 0, 2));
  Lie_gen.push_back(get_Lie_gen("k_", K_re, scl_ksi[2], 0, 0, 1, 1, 2));

  Lie_gen.push_back(get_Lie_gen("k_", K_re, scl_ksi[3], 1, 1, 0, 0, 3));
  Lie_gen.push_back(get_Lie_gen("k_", K_re, scl_ksi[3], 0, 0, 1, 1, 3));

  Lie_gen.push_back(get_Lie_gen("k_", K_re, scl_ksi[4], 1, 1, 0, 0, 4));
  Lie_gen.push_back(get_Lie_gen("k_", K_re, scl_ksi[4], 0, 0, 1, 1, 4));

  if (NO >= 8) {
    Lie_gen.push_back(get_Lie_gen("k_", K_re, scl_ksi[5], 1, 1, 0, 0, 5));
    Lie_gen.push_back(get_Lie_gen("k_", K_re, scl_ksi[5], 0, 0, 1, 1, 5));
  }

  return Lie_gen;
}


inline double h_ijklm_p
(const tps &h, const std::vector<int> &ind)
{
  return h_ijklm_p(h, ind[x_], ind[px_], ind[y_], ind[py_], ind[delta_], 7);
}


void get_Jacobian
(const param_type &bns, const int k, const ss_vect<tps> &Id_scl,
 std::vector<Lie_gen_class> &Lie_gen)
{
  const bool
    prt     = false;
  const std::string
    name    = bns.name[k];
  const int
    // Driving terms.
    index[] = {2, 17},
    n       = bns.n[k];
  const double
    bn_scl  = (bns.L[k] == 0e0)? 1e0 : bns.bnL_scl[k]/bns.L[k];

  tps          g_re, g_im, K_re, K_im;
  ss_vect<tps> M, A1, A0, M_res;

  set_bn_par(bns.Fnum[k], n, 7);

  danot_(no_tps-1);
  M = get_map();
  danot_(no_tps);
  auto K = MapNorm(M, g, A1, A0, M_res, 1);
  CtoR(g*Id_scl, g_re, g_im);
  CtoR(K*Id_scl, K_re, K_im);

  if (prt)
    std::cout << std::scientific << std::setprecision(3)
	      << "\n" << std::setw(8) << name <<
      std::setw(4) << bns.Fnum[k] << ":\n"
	      << std::scientific << std::setprecision(3) << K_re << K_im;
 
  clr_bn_par(bns.Fnum[k], n);

  for (auto j = 0; j < (int)Lie_gen.size(); j++) {
    if ((j >= index[0]) and (j <= index[1]))
      Lie_gen[j].Jacobian.push_back
	(bn_scl*Lie_gen[j].cst_scl*h_ijklm_p(g_im, Lie_gen[j].index));
    else
      Lie_gen[j].Jacobian.push_back
	(bn_scl*Lie_gen[j].cst_scl*h_ijklm_p(K_re, Lie_gen[j].index));
  }
}


void get_Jacobian
(const ss_vect<tps> &Id_scl, const param_type &bns,
 std::vector<Lie_gen_class> &Lie_gen)
{
  for (auto k = 0; k < bns.n_prm; k++)
    get_Jacobian(bns, k, Id_scl, Lie_gen);
}


void get_system
(const int m, const int n, const std::vector<Lie_gen_class> &Lie_gen,
 double **A, double *b)
{
  const bool prt = false;

  for (auto j = 0; j < m; j++) {
    b[j+1] = -Lie_gen[j].cst;
    for (auto k = 0; k < n; k++)
      A[j+1][k+1] = Lie_gen[j].Jacobian[k];
  }

  if (prt) dmdump(stdout, (char *)"\nA:", A, m, n, (char *)" %10.3e");
}


Lie_gen_class get_Lie_gen_sum
(const std::string &name, const double scl, const int index[],
 const std::vector<Lie_gen_class> &Lie_gen)
{
  Lie_gen_class Lie_gen_sum;

  Lie_gen_sum.label = name;
  Lie_gen_sum.cst_scl = scl;
  Lie_gen_sum.cst = 0e0;
  for (auto k = 0; k < (int)Lie_gen[0].Jacobian.size(); k++)
    Lie_gen_sum.Jacobian.push_back(0e0);

  for (auto j = index[0]; j <= index[1]; j++) {
    if (Lie_gen[j].cst_scl == 0e0) {
      printf("\nget_Lie_gen_sum: cst_scl = 0 %s\n", Lie_gen[j].label.c_str());
      assert(false);
    }
    Lie_gen_sum.cst += Lie_gen_sum.cst_scl*Lie_gen[j].cst/Lie_gen[j].cst_scl;
    for (auto k = 0; k < (int)Lie_gen[j].Jacobian.size(); k++)
      Lie_gen_sum.Jacobian[k] +=
	Lie_gen_sum.cst_scl*Lie_gen[j].Jacobian[k]/Lie_gen[j].cst_scl;
  }

  return Lie_gen_sum;
}


void prt_system
(const param_type &bns, const std::vector<Lie_gen_class> &Lie_gen)
{
  printf("\n           scl.      cst.");
  for (auto k = 0; k < bns.n_prm; k++)
    printf("      %-5s", bns.name[k].c_str());
  printf("\n                         ");
  for (auto k = 0; k < bns.n_prm; k++)
    printf("       %1d   ", bns.n[k]);
  printf("\n                         ");
  for (auto k = 0; k < bns.n_prm; k++)
    printf("    %7.1e", bns.bnL_scl[k]);

  auto k = 0;
  prt_Lie_gen("Linear chromaticity:",         k, 2, Lie_gen);
  k += 2;
  prt_Lie_gen("3rd Order Chromatic terms:",   k, 3, Lie_gen);
  k += 3;
  prt_Lie_gen("3rd Order Geometric terms:",   k, 5, Lie_gen);
  k += 5;
  prt_Lie_gen("4th Order Geometric terms:",   k, 8, Lie_gen);
  k += 8;
  prt_Lie_gen("4th Order Anharmonic terms:",  k, 3, Lie_gen);
  k += 3;
  prt_Lie_gen("6th Order Anharmonic terms:",  k, 4, Lie_gen);
  k += 4;
  if (NO >= 9) {
    prt_Lie_gen("8th Order Anharmonic terms:",  k, 5, Lie_gen);
    k += 5;
  }
  if (NO >= 11) {
    prt_Lie_gen("10th Order Anharmonic terms:", k, 6, Lie_gen);
    k += 6;
  }
  prt_Lie_gen("2nd Order Chromaticity:",      k, 2, Lie_gen);
  k += 2;
  prt_Lie_gen("3rd Order Chromaticity:",      k, 2, Lie_gen);
  k += 2;
  prt_Lie_gen("4th Order Chromaticity:",      k, 2, Lie_gen);
  k += 2;
  if (NO >= 8) {
    prt_Lie_gen("5th Order Chromaticity:",      k,  2, Lie_gen);
    k += 2;
  }

  prt_Lie_gen("Sigma{K_dnu}:",                k, 1, Lie_gen);
  k += 1;
  prt_Lie_gen("Sigma{K_dxi}:",                k, 1, Lie_gen);
}


std::vector<Lie_gen_class> analyze
(const ss_vect<tps> &Id_scl, const param_type &bns)
{
  int
    index_dnu[] = {18, 24},
    index_dxi[] = {25, 30};

  if (NO == 8)
    index_dxi[1] += 2;
  else if (NO == 9) {
    index_dnu[1] += 5;
    index_dxi[0] += 5;
    index_dxi[1] += 5 + 2;
  } else if (NO == 11) {
    index_dnu[1] += 11;
    index_dxi[0] += 11;
    index_dxi[1] += 11 + 2;
  }

  auto Lie_gen = get_Lie_gen(Id_scl);
  get_Jacobian(Id_scl, bns, Lie_gen);

  Lie_gen.push_back(get_Lie_gen_sum
		    ("K_sum  ", scl_K_sum[0], index_dnu, Lie_gen));
  Lie_gen.push_back(get_Lie_gen_sum
		    ("K_sum  ", scl_K_sum[1], index_dxi, Lie_gen));

  prt_system(bns, Lie_gen);

  return Lie_gen;
}


std::vector<int> sort_sing_val(const int n, const double w[])
{
  const bool prt = false;

  std::vector<int> index;

  for (auto k = 0; k < n; k++)
    index.push_back(k+1);
  std::sort(index.begin(), index.end(), [&w](int a, int b) {
    return w[a] > w[b];
  });

  if (prt) {
    printf("\n");
    for (auto k = 0; k < n; k++)
      printf(" %d", index[k]);
    printf("\n");
  }

  return index;
}


void get_sing_val(const int n, double w[], const int svd_n_cut)
{
  const int n_prt = 8;

  std::vector<int> ind;

  ind = sort_sing_val(n, w);
  for (auto k = 1; k <= n; k++) {
    printf("  %9.3e", w[ind[k-1]]);
    if (k > n-svd_n_cut) {
      w[ind[k-1]] = 0e0;
      printf(" (zeroed)");
    }
    if (k % n_prt == 0) printf("\n");
  }
  if (n % n_prt != 0) printf("\n");
}


void set_Fam(param_type &bns, const int k, const double scl, const double *dbnL)
{
  double bnL_ext;

  bnL_ext =
    get_bnL(bns.Fnum[k], 1, bns.n[k]) + scl*bns.bnL_scl[k]*dbnL[k+1];
  bns.bnL[k] = bnL_internal(bnL_ext, bns.bnL_min[k], bns.bnL_max[k]);
  set_bnL(bns.Fnum[k], bns.n[k],
	  bnL_bounded(bns.bnL[k], bns.bnL_min[k], bns.bnL_max[k]));
}


void set_bnL(const double scl, const double *dbnL, param_type &bns)
{
  for (auto k = 0; k < bns.n_prm; k++)
    set_Fam(bns, k, scl, dbnL);
}


void prt_bend(FILE *outf, const int loc, const int n)
{
  const elem_type<double> *elemp = &elem[loc-1];

  fprintf(outf,
	  "%-8s: Multipole, L = %7.5f, Phi = %7.5f, Phi_1 = %7.5f"
	  ", Phi_2 = %7.5f,\n"
	  "          HOM = (%d, %12.5e, 0e0, %d, %12.5e, 0e0),\n"
	  "          N = nbend;\n",
	  elemp->Name, elemp->L,
	  elemp->L*elemp->mpole->h_bend*180e0/M_PI,
	  elemp->mpole->edge1, elemp->mpole->edge2,
	  Quad, get_bn(elem[loc-1].Fnum, elem[loc-1].Knum, Quad),
	  n, get_bn(elem[loc-1].Fnum, elem[loc-1].Knum, n));
}


void prt_quad(FILE *outf, const int loc, const int n)
{
  fprintf
    (outf,
     "%-8s: Multipole, L = %7.5f,\n          HOM = (%d, %12.5e, 0e0,"
     " %d, %12.5e, 0e0),\n          N = nquad;\n",
     elem[loc-1].Name, elem[loc-1].L, Quad,
     get_bn(elem[loc-1].Fnum, elem[loc-1].Knum, Quad), n,
     get_bn(elem[loc-1].Fnum, elem[loc-1].Knum, n));
}


void prt_single_mult(FILE *outf, const int loc, const int n)
{
  switch (n) {
  case Sext:
    fprintf
      (outf,
       "%-8s: Sextupole, L = %7.5f, B_3 = %12.5e, N = %d;\n",
       elem[loc-1].Name, elem[loc-1].L, elem[loc-1].mpole->bn[Sext-1],
       elem[loc-1].mpole->n_step);
    break;
  case Oct:
    fprintf
      (outf,
       "%-8s: Octupole, L = %7.5f, B_4 = %12.5e, N = %d;\n",
       elem[loc-1].Name, elem[loc-1].L, elem[loc-1].mpole->bn[Oct-1],
       elem[loc-1].mpole->n_step);
    break;
  default:
    printf("\nprt_single_mult - undefined multipole order: %d\n", n);
    break;
  }
}


int get_n_mpole(const int loc)
{
  int n_mpole = 0;

  for (auto k = 0; k < elem[loc-1].mpole->order; k++) {
    if ((elem[loc-1].mpole->bn[k] != 0e0) || (elem[loc-1].mpole->an[k] != 0e0))
      n_mpole++;
  }
  return n_mpole;
}


void prt_mult(FILE *outf, const int loc, const int n)
{
  std::string name;
  bool        first = true;
  int         n_step;
  double      L;

  if (get_n_mpole(loc) == 1)
    prt_single_mult(outf, loc, n);
  else {
    name = elem[loc-1].Name;
    L = elem[loc-1].L;
    fprintf
      (outf, "%-8s: multipole, l = %7.5f, hom = (", name.c_str(), L);
    for (auto k =  0; k < elem[loc-1].mpole->order; k++) {
      if ((elem[loc-1].mpole->bn[k] != 0e0)
	  || (elem[loc-1].mpole->an[k] != 0e0)) {
	if (first) {
	  fprintf(outf,
		  "\n            %d, %12.5e, %12.5e",
		  k+1, elem[loc-1].mpole->bn[k], elem[loc-1].mpole->an[k]);
	  first = false;
	} else
	  fprintf(outf,
		  ",\n            %d, %12.5e, %12.5e",
		  k+1, elem[loc-1].mpole->bn[k], elem[loc-1].mpole->an[k]);
      }
    }
    n_step = elem[loc-1].mpole->n_step;
    fprintf(outf, "), n = %d;\n", n_step);
  }
}


void prt_bn(const param_type &bns)
{
  const std::string file_name = "b4.out";

  FILE *outf;

  outf = file_write(file_name.c_str());

  fprintf(outf, "\n");
  for (auto k = 0; k < bns.n_prm; k++) {
    auto loc = (bns.Fnum[k] > 0)? get_loc(bns.Fnum[k], 1) : bns.locs[k][0];
    if (elem[loc-1].mpole->n_design == Dip)
      prt_bend(outf, loc, bns.n[k]);
    else if (elem[loc-1].mpole->n_design == Quad)
      prt_quad(outf, loc, bns.n[k]);
    else
      prt_mult(outf, loc, bns.n[k]);
  }

  fclose(outf);
}


void correct
(param_type &bns, const std::vector<Lie_gen_class> &Lie_gen,
 const int svd_n_cut,
 const double scl)
{
  const int
    m = Lie_gen.size(),
    n = bns.n_prm;

  double **A, **U, **V, *w, *b, *dbnL, *bnL_max, *bnL;

  printf("\nsvd:\n  m = %d n = %d\n", m, n);

  A = dmatrix(1, m, 1, n); U = dmatrix(1, m, 1, n); V = dmatrix(1, n, 1, n);
  w = dvector(1, n); b = dvector(1, m); dbnL = dvector(1, n);
  bnL_max = dvector(1, n); bnL = dvector(1, n);

  get_system(m, n, Lie_gen, A, b);

#if 1
  dmcopy(A, m, n, U);
  dsvdcmp(U, m, n, w, V);
  get_sing_val(n, w, svd_n_cut);

  dsvbksb(U, w, V, m, n, b, dbnL);
#else
  for (auto k = 0; k < n; k++) {
    auto Fnum = bns.Fnum[k];
    bnL_max[k+1] = (bns.L[k] == 0e0)?
      bns.bnL_max[k]/bns.bnL_scl[k] :
      bns.bnL_max[k]*bns.L[k]/bns.bnL_scl[k];
    bnL[k+1] = get_bnL(Fnum, 1, bns.n[k])/bns.bnL_scl[k];
  }

  SVD_lim(m, n, A, b, bnL_max, 1e-11, bnL, dbnL);

#endif

  set_bnL(scl, dbnL, bns);
  bns.print();
  prt_bn(bns);
  prt_mfile("flat_file.fit");

  free_dmatrix(A, 1, m, 1, n); free_dmatrix(U, 1, m, 1, n);
  free_dmatrix(V, 1, n, 1, n); free_dvector(w, 1, n); free_dvector(b, 1, m);
  free_dvector(dbnL, 1, n); free_dvector(bnL_max, 1, n);
  free_dvector(bnL, 1, n);
}


void no_mpoles(const int n)
{
  printf("\nzeroing multipoles: %d\n", n);
  for (auto j = 0; j < n_elem; j++)
    if (elem[j].kind == Mpole)
      set_bn(elem[j].Fnum, elem[j].Knum, n, 0e0);
}


void get_bns(param_type &bns)
{
  const int lat = 4;

  if (b_3_zero)
    no_mpoles(Sext);
  if (b_4_zero)
    no_mpoles(Oct);

  switch (lat) {
  case 1:
    if (b_3_opt) {
      bns.add_Fam("s1", Sext, bnL_min[Sext], bnL_max[Sext], bnL_scl[Sext]);
      bns.add_Fam("s2", Sext, bnL_min[Sext], bnL_max[Sext], bnL_scl[Sext]);
      bns.add_Fam("s3", Sext, bnL_min[Sext], bnL_max[Sext], bnL_scl[Sext]);
      bns.add_Fam("s4", Sext, bnL_min[Sext], bnL_max[Sext], bnL_scl[Sext]);
    }

    if (b_4_opt) {
      bns.add_Fam("o1",  Oct, bnL_min[Oct], bnL_max[Oct], bnL_scl[Oct]);
      bns.add_Fam("o2",  Oct, bnL_min[Oct], bnL_max[Oct], bnL_scl[Oct]);
      bns.add_Fam("o3",  Oct, bnL_min[Oct], bnL_max[Oct], bnL_scl[Oct]);
    }
    break;
  case 2:
    if (b_3_opt) {
      bns.add_Fam("s1_f1", Sext, bnL_min[Sext], bnL_max[Sext], bnL_scl[Sext]);
      bns.add_Fam("s2_f1", Sext, bnL_min[Sext], bnL_max[Sext], bnL_scl[Sext]);
      bns.add_Fam("s3_f1", Sext, bnL_min[Sext], bnL_max[Sext], bnL_scl[Sext]);
      bns.add_Fam("s4_f1", Sext, bnL_min[Sext], bnL_max[Sext], bnL_scl[Sext]);
    }

    if (b_4_opt) {
      bns.add_Fam("o1_f1",  Oct, bnL_min[Oct], bnL_max[Oct], bnL_scl[Oct]);
      bns.add_Fam("o2_f1",  Oct, bnL_min[Oct], bnL_max[Oct], bnL_scl[Oct]);
      bns.add_Fam("o3_f1",  Oct, bnL_min[Oct], bnL_max[Oct], bnL_scl[Oct]);
    }
    break;
  case 3:
    if (b_3_opt) {
      bns.add_Fam("s1_h1", Sext, bnL_min[Sext], bnL_max[Sext], bnL_scl[Sext]);
      bns.add_Fam("s2_h1", Sext, bnL_min[Sext], bnL_max[Sext], bnL_scl[Sext]);
      bns.add_Fam("s3_h1", Sext, bnL_min[Sext], bnL_max[Sext], bnL_scl[Sext]);
      bns.add_Fam("s4_h1", Sext, bnL_min[Sext], bnL_max[Sext], bnL_scl[Sext]);
    }

    if (b_4_opt) {
      bns.add_Fam("o1_h1",  Oct, bnL_min[Oct], bnL_max[Oct], bnL_scl[Oct]);
      bns.add_Fam("o2_h1",  Oct, bnL_min[Oct], bnL_max[Oct], bnL_scl[Oct]);
      bns.add_Fam("o3_h1",  Oct, bnL_min[Oct], bnL_max[Oct], bnL_scl[Oct]);
    }
    break;
  case 4:
    if (b_3_opt) {
      bns.add_Fam("s1_h2", Sext, bnL_min[Sext], bnL_max[Sext], bnL_scl[Sext]);
      bns.add_Fam("s2_h2", Sext, bnL_min[Sext], bnL_max[Sext], bnL_scl[Sext]);
      bns.add_Fam("s3_h2", Sext, bnL_min[Sext], bnL_max[Sext], bnL_scl[Sext]);
      bns.add_Fam("s4_h2", Sext, bnL_min[Sext], bnL_max[Sext], bnL_scl[Sext]);
      if (false)
	bns.add_Fam("s5_h2", Sext, bnL_min[Sext], bnL_max[Sext], bnL_scl[Sext]);
    }

    if (b_4_opt) {
      bns.add_Fam("o1_h2",  Oct, bnL_min[Oct], bnL_max[Oct], bnL_scl[Oct]);
      bns.add_Fam("o2_h2",  Oct, bnL_min[Oct], bnL_max[Oct], bnL_scl[Oct]);
      bns.add_Fam("o3_h2",  Oct, bnL_min[Oct], bnL_max[Oct], bnL_scl[Oct]);
    }
    break;
  case 5:
    if (b_3_opt) {
      bns.add_Fam("sfm", Sext, bnL_min[Sext], bnL_max[Sext], bnL_scl[Sext]);
      bns.add_Fam("sfi", Sext, bnL_min[Sext], bnL_max[Sext], bnL_scl[Sext]);
      bns.add_Fam("sdqd_1", Sext, bnL_min[Sext], bnL_max[Sext], bnL_scl[Sext]);
      bns.add_Fam("sdqd_2", Sext, bnL_min[Sext], bnL_max[Sext], bnL_scl[Sext]);
      bns.add_Fam("sdqd_3", Sext, bnL_min[Sext], bnL_max[Sext], bnL_scl[Sext]);
      bns.add_Fam("sdqd_4", Sext, bnL_min[Sext], bnL_max[Sext], bnL_scl[Sext]);
      bns.add_Fam("sdqd_5", Sext, bnL_min[Sext], bnL_max[Sext], bnL_scl[Sext]);
      bns.add_Fam("sdendq", Sext, bnL_min[Sext], bnL_max[Sext], bnL_scl[Sext]);
      bns.add_Fam("sfo", Sext, bnL_min[Sext], bnL_max[Sext], bnL_scl[Sext]);
    }

    if (b_4_opt) {
      bns.add_Fam("oxxo",  Oct, bnL_min[Oct], bnL_max[Oct], bnL_scl[Oct]);
      bns.add_Fam("oxyo",  Oct, bnL_min[Oct], bnL_max[Oct], bnL_scl[Oct]);
      bns.add_Fam("oyyo",  Oct, bnL_min[Oct], bnL_max[Oct], bnL_scl[Oct]);
    }
    break;
  }
}


void chk_lat(void)
{
  double       alpha[2], beta[2], nu[3], ksi[2];
  tps          g, K;
  ss_vect<tps> nus, M, A1, A0, M_res;

  danot_(2);
  M = get_map();
  danot_(3);
  K = MapNorm(M, g, A1, A0, M_res, 1);
  nus = dHdJ(K);
  get_nu_ksi(nus, nu, ksi);
  get_ab(A1, alpha, beta, 0);
  printf
    ("\n  alpha = [%6.3f, %6.3f]\n  beta  = [%6.3f, %6.3f]\n"
     "  nu    = [%6.3f, %6.3f]\n  ksi   = [%6.3f, %6.3f]\n",
     alpha[X_], alpha[Y_], beta[X_], beta[Y_], nu[X_], nu[Y_], ksi[X_],
     ksi[Y_]);
}


void set_state(void)
{
  rad_on         = false;
  H_exact        = false;
  totpath_on     = false;
  cavity_on      = false;
  quad_fringe_on = false;
  emittance_on   = false;
  IBS_on         = false;
}


int main(int argc, char *argv[])
{
  double                      twoJ;
  ss_vect<tps>                Id_scl;
  param_type                  bns;
  std::vector<Lie_gen_class> Lie_gen;

  set_state();

  rd_mfile(argv[1], elem);
  rd_mfile(argv[1], elem_tps);

  // Initialize the symplectic integrator after the energy has been defined.
  ini_si();

  // Disable from TPSALib & LieLib log messages
  idprset(-1);

  daeps_(1e-30);

  chk_lat();

  Id_scl.identity();
  for (auto k = 0; k < 2; k++) {
    twoJ = sqr(A_max[k])/beta_inj[k];
    Id_scl[2*k] *= sqrt(twoJ);
    Id_scl[2*k+1] *= sqrt(twoJ);
  }
  Id_scl[delta_] *= delta_max;

  get_bns(bns);
  bns.ini_prm();
  bns.print();

  printf("\n");
  for (auto k = 1; k <= max_iter; k++) {
    printf("\nk = %d:", k);
    Lie_gen = analyze(Id_scl, bns);
    correct(bns, Lie_gen, svd_n_cut, step);

    prt_mfile("flat_file.fit");
  }
  Lie_gen = analyze(Id_scl, bns);
}
