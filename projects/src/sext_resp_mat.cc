#include <assert.h>

#define NO 3

#include <assert.h>

#include "thor_lib.h"

#include "param_type.h"

int
  no_tps   = NO,
  ndpt_tps = 5;


tps get_h_k(const tps &h, const int k)
{
  // Take in Forest's F77 LieLib.
  // Get monomials of order k.
  long int no;
  tps      h_k;

  no = getno_();
  danot_(k-1);
  h_k = -h;
  danot_(k);
  h_k += h;
  danot_(no);
  return h_k;
}


ss_vect<tps> get_M_k(const ss_vect<tps> &x, const int k)
{
  // Taked in Forest's F77 LieLib.
  int          i;
  ss_vect<tps> map_k;

  for (i = 0; i < nv_tps; i++)
    map_k[i] = get_h_k(x[i], k);
  return map_k;
}


tps v_to_tps(const ss_vect<tps> &v, const tps &x)
{
  // Daflo in Forest's F77 LieLib.
  //   y = v * nabla * x
  int k;
  tps y;

  y = 0e0;
  for (k = 0; k < 2*nd_tps; k++)
    y += v[k]*Der(x, k+1);
  return y;
}


tps exp_v_to_tps(const ss_vect<tps> &v, const tps &x, const double eps,
	      const int n_max)
{
  // Expflo in Forest's F77 LieLib:
  //   y = exp(v*nabla) * x
  int    k;
  double eps1;
  tps    y_k, y;

  y_k = y = x;
  for (k = 1; k <= n_max; k++) {
    y_k = v_to_tps(v, y_k/k);
    y += y_k;
    eps1 = abs(y_k);
    if (eps1 < eps)
      break;
  }
  if (eps1 < eps)
    return y;
  else {
    printf("\n*** exp_v_to_tps: did not converge eps = %9.3e (eps = %9.3e)"
	   " n_max = %1d\n", eps1, eps, n_max);
    return NAN;
  }
}


tps exp_v_fac_to_tps(const ss_vect<tps> &v, const tps &x, const int k1,
		      const int k2, const double scl)
{
  // Facflo in Forest's F77 LieLib.
  //   y = exp(D_k1) * exp(D_k1+1) ...  * exp(D_k2) * x
  int          k;
  tps          y;
  ss_vect<tps> v_k;

  const int n_max = 100; 

  y = x;
  for (k = k1; k <= k2; k++) {
    v_k = scl*get_M_k(v, k);
    y = exp_v_to_tps(v_k, y, eps_tps, n_max);
  }
  return y;
}


ss_vect<tps> exp_v_to_map(const ss_vect<tps> &v, const ss_vect<tps> &map)
{
  const int n_max = 100; 

  ss_vect<tps> M;

  for (auto k = 0; k < 2*nd_tps; k++)
    M[k] = exp_v_to_tps(v, map[k], eps_tps, n_max);
  return map;
}


ss_vect<tps>M_to_M_fact(const ss_vect<tps> &map)
{
  // Flofac in Forest's F77 LieLib.
  // Factor map:
  //   M = M_2 ... * M_n
  int          k;
  tps          y;
  ss_vect<tps> map_lin, v_k, map_res, map_fact;

  std::cout << "\nM_to_M_fact:\n1:\n" << map[x_];
  map_lin = get_M_k(map, 1);
  std::cout << "\n2:\n" << map_lin;
  prt_lin_map(3, Inv(map_lin));
  std::cout << "\n3:\n" << Inv(map_lin);
  map_res = map*Inv(map_lin);
  std::cout << "\n4:\n" << map_res;

  map_fact.zero();
  for (k = 2; k <= no_tps; k++) {
    map_fact += get_M_k(map_res, k);
    v_k = get_M_k(-map_fact, k);
    map_res = exp_v_to_map(v_k, map_res);
  }
  std::cout << "\n5:\n" << map_fact[x_];
  return map_fact;
}


ss_vect<tps>M_to_M_fact2(const ss_vect<tps> &map)
{
  // Obsolete.
  // Flofac in Forest's F77 LieLib.
  // Factor map:
  //   M = M_2 ... * M_n
  int          j, k;
  ss_vect<tps> map_lin, map_res, map_fact;

  map_lin = get_M_k(map, 1);
  map_res = map*Inv(map_lin);
  map_fact.zero();
  for (k = 2; k <= no_tps; k++) {
    map_fact += get_M_k(map_res, k);
    for (j = 0; j < 2*nd_tps; j++)
      map_res[j] = exp_v_fac_to_tps(map_fact, map_res[j], k, k, -1e0);
  }
  return map_fact;
}


tps tps_fun
(const tps &a, std::function<double (const long int [])> fun)
{
  // Dacfu in Forest's F77 LieLib.
  // Multiplies mononials I_vec with function f(I_vec).
  char     name[name_len_for+1];
  int      k, n;
  long int ibuf1[bufsize], ibuf2[bufsize], jj[nv_tps];
  double   rbuf[bufsize];
  tps      b;

  a.exprt(rbuf, ibuf1, ibuf2, name);
  n = rbuf[0];
  for (k = 1; k <= n; k++) {
    dehash_(no_tps, nv_tps, ibuf1[k-1], ibuf2[k-1], jj);
    rbuf[k] *= fun(jj);
  }
  b.imprt(n, rbuf, ibuf1, ibuf2);
  return b;
}


double f_int_mon(const long int jj[])
{
  // Integrate monomials:
  //   scl = 1/(|I_vec|+1)
  int    k;
  double scl;

  scl = 0e0;
  for (k = 0; k < 2*nd_tps; k++)
    scl += jj[k];
  scl += 1e0;
  scl = 1e0/scl;
  return scl;
}


tps M_to_h(const ss_vect<tps> &map)
{
  // Intd in Forest's F77 LieLib.
  // E. Forest, M. Berz, J. Irwin "Normal Form Methods for Complicated
  // Periodic Systems: A Complete Solution Using Differential Algebra and Lie
  // Operators" Part. Accel. 24, 91-107 (1989):
  //   Eqs. (34)-(37).
  // Integrate monomials:
  //   M -> exp(:h:)
  int          k;
  tps          f_x, f_px, h;
  ss_vect<tps> Id;

  Id.identity();
  h = 0e0;
  for (k = 0; k < nd_tps; k++) {
    // Integrate monomials.
    f_x = tps_fun(map[2*k+1], f_int_mon)*Id[2*k];
    f_px = tps_fun(map[2*k], f_int_mon)*Id[2*k+1];
    h += f_x - f_px;
  }
  return h;
}


tps M_to_h_DF(const ss_vect<tps> &map)
{
  // Liefact in Forest's F77 LieLib.
  // A. Dragt, J. Finn "Lie Series and Invariant Functions for Analytic
  // Symplectic maps" J. Math. Phys. 17, 2215-2227 (1976).
  // Dragt-Finn factorization:
  //   M ->  M_lin * exp(:h_3:) * exp(:h_4:) ...  * exp(:h_n:)
  return M_to_h(M_to_M_fact(map));
}


tps get_h(ss_vect<tps> &map)
{
  tps          h;
  ss_vect<tps> map1, R;

  if (true) {
    // Dragt-Finn factorization
    h = LieFact_DF(Inv(A0*A1)*map*A0*A1, R);
    R[6] = tps(0e0, 7); h = h*R;
    return h;
  } else {
    // single Lie exponent
    danot_(1); map1 = map; danot_(no_tps);
    return LieFact(Inv(A0*A1)*map*Inv(map1)*A0*A1);
  }
}


void no_mpoles(const int n)
{
  int j;

  printf("\nzeroing multipoles: %d\n", n);
  for (j = 0; j < n_elem; j++)
    if (elem[j].kind == Mpole)
      set_bn(elem[j].Fnum, elem[j].Knum, n, 0e0);
}


void dragt_finn_fact()
{
  tps          h_1, h_2;
  ss_vect<tps> R;

  set_bn_par(get_Fnum("uq3"), Quad, 7);

  Map.identity();
  Map.propagate(1, n_elem);

  danot_(no_tps);

  prt_lin_map(3, Map);

  if (false)
    std::cout << std::scientific << std::setprecision(5)
	      << std::setw(13) << Map << "\n";

  h_1 = LieFact_DF(Map, R);
  R[6] = tps(0e0, 7); h_1 = h_1*R;

  h_2 = M_to_h_DF(Map);

  assert(false);
  std::cout << std::scientific << std::setprecision(5)
	    << "\nh:\n" << std::setw(13) << h_2 << "\n";
}


void print_map(const std::string &str, ss_vect<tps> M)
{
  std::cout << str << "\n";
  for (auto j = 0; j < ss_dim; j++) {
    for (auto k = 0; k < ss_dim; k++)
      std::cout << std::scientific << std::setprecision(5)
		<< std::setw(13) << M[j][k];
    std::cout << "\n";
  }
}


void test(void)
{
  ss_vect<tps> Id, M;

  Id.identity();

  M.zero();
  for (auto k = 0; k < 2*nd_tps; k++)
    M[k] += (k+1e0)*Id[k];
  M[x_] += 3.14e-3*Id[x_]*Id[6];

  print_map("\nId:", Id);
  print_map("\nM:", M);
  std::cout << std::scientific << std::setprecision(5)
	    << std::setw(13) << M[x_];
  std::cout << std::scientific << std::setprecision(5)
	    << std::setw(13) << (M*Id)[x_];
  std::cout << std::scientific << std::setprecision(5)
	    << std::setw(13) << (Id*M)[x_];
  std::cout << std::scientific << std::setprecision(5)
	    << std::setw(13) << (M*M)[x_];
}


int main(int argc, char *argv[])
{

  rad_on    = false; H_exact        = false; totpath_on   = false;
  cavity_on = false; quad_fringe_on = false; emittance_on = false;
  IBS_on    = false;

  // rd_mfile(argv[1], elem); rd_mfile(argv[1], elem_tps);

  // initialize the symplectic integrator after energy has been defined
  ini_si();

  // disable from TPSALib- and LieLib log messages
  idprset(-1);

  daeps_(1e-30);

  danot_(no_tps-1);

  if (false)
    no_mpoles(Sext);

  // dragt_finn_fact();

  danot_(no_tps);
  test();
}
