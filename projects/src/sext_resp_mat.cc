#include <eigen3/Eigen/Dense>
#include <complex>
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>
#include <boost/math/special_functions/sign.hpp>
#include <assert.h>

#define NO 5

#include <assert.h>

#include "thor_lib.h"

#include "param_type.h"

int
  no_tps   = NO,
  ndpt_tps = 5;


typedef struct MNFType
{
  tps
    K,              // Normalised generator.
    g;              // Generator for canonical transformation to Floquet space.
  ss_vect<tps>
    M,              // Poincaré map.
    M_res,          // Residual map.
    A0, A0_inv,     // Linear transformation to fixed point.
    A1, A1_inv,     // Linear transformation to Floquet space.
    A_nl, A_nl_inv, // Nonlinear transformation to Floquet space.
    R,              // Floquet space rotation.
    nus;            // Tune shift.
} MNFType;


ss_vect<tps> mat2map(const Eigen::MatrixXd &A)
{
  ss_vect<tps> Id, B;

  Id.identity();
  B.zero();
  for (auto j = 0; j < 2*nd_tps; j++)
    for (auto k = 0; k < 2*nd_tps; k++)
      B[j] += A(j, k)*Id[k];

  return B;
}


void print_map(const std::string &str, ss_vect<tps> M)
{
  const int n_dec = 5;

  std::cout << std::scientific << std::setprecision(n_dec)
	    << str << "cst:\n" << std::setw(n_dec+8) << M.cst()
	    << "\nlinear:\n";
  for (auto j = 0; j < 2*nd_tps; j++) {
    for (auto k = 0; k < 2*nd_tps; k++)
      std::cout << std::scientific << std::setprecision(n_dec)
		<< std::setw(n_dec+8) << M[j][k];
    std::cout << "\n";
  }
}


void print_int_vec(const std::string &str, const Eigen::VectorXi &v)
{
  std::cout << str;
  for (auto k = 0; k < v.size(); k++)
    std::cout << std::setw(1) << v(k) << "\n";
}


void print_vec(const std::string &str, const Eigen::VectorXd &v)
{
  const int n_dec = 5;

  std::cout << str;
  for (auto k = 0; k < v.size(); k++)
    std::cout << std::scientific << std::setprecision(n_dec)
	      << std::setw(n_dec+8) << v(k) << "\n";
}


void print_mat(const std::string &str, const Eigen::MatrixXd &M)
{
  const int n_dec = 5;

  std::cout << str;
  for (auto j = 0; j < M.col(0).size(); j++) {
    for (auto k = 0; k < M.row(0).size(); k++)
      std::cout << std::scientific << std::setprecision(n_dec)
		<< std::setw(n_dec+8) << M(j, k);
    std::cout << "\n";
  }
}


void print_complex_vec(const std::string &str, const Eigen::VectorXcd &v)
{
  const int n_dec = 5;

  std::cout << str;
  for (auto k = 0; k < v.size(); k++)
    std::cout << std::scientific << std::setprecision(n_dec)
	      << std::setw(n_dec+8) << v(k).real()
	      << ((v(k).imag() > 0e0)? "+i":"-i")
	      << std::setw(n_dec+6) << fabs(v(k).imag()) << "\n";
}


void print_complex_mat(const std::string &str, const Eigen::MatrixXcd &M)
{
  const int n_dec = 5;

  std::cout << str;
  for (auto j = 0; j < M.col(0).size(); j++) {
    for (auto k = 0; k < M.row(0).size(); k++)
      std::cout << std::scientific << std::setprecision(n_dec)
		<< std::setw(n_dec+8) << M(j, k).real()
		<< ((M(j, k).imag() > 0e0)? "+i":"-i")
		<< std::setw(n_dec+6) << fabs(M(j, k).imag());
    std::cout << "\n";
  }
}


double acos2(const double sin, const double cos)
{
  if (fabs(cos) > 1e0) {
    std::cout << std::scientific << std::setprecision(5)
	      << "\nacos2 - argument > 1: " << cos << "\n";
    return NAN;
  }
  auto phi = acos(cos);
  if (sin < 0e0)
    phi = 2e0*M_PI - phi;
  return phi;
}


void no_mpoles(const int n)
{
  printf("\nzeroing multipoles: %d\n", n);
  for (auto j = 0; j < n_elem; j++)
    if (elem[j].kind == Mpole)
      set_bn(elem[j].Fnum, elem[j].Knum, n, 0e0);
}


tps get_mns(const tps &a, const int no1, const int no2)
{
  tps  b;

  danot_(no1-1);
  b = -a;
  danot_(no2);
  b += a;
  danot_(no_tps);

  return b;
}


ss_vect<tps> get_mns(const ss_vect<tps> &x, const int no1, const int no2)
{
  ss_vect<tps> y;

  for (auto k = 0; k < nv_tps; k++)
    y[k] = get_mns(x[k], no1, no2);

  return y;
}


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
  ss_vect<tps> map_k;

  for (auto i = 0; i < nv_tps; i++)
    map_k[i] = get_h_k(x[i], k);
  return map_k;
}


#if 1
// Obsolete.

void CtoR_JB(const tps &a, tps &a_re, tps &a_im)
{
  tps          b, c;
  ss_vect<tps> Id, map;

  Id.identity();

  // b = p_k_cmplx_sgn_corr_fun(nd_tps, a);
  std::cout << "\np_k_cmplx_sgn_corr_fun() not implemented\n";
  assert(false);

  // q_k -> (q_k + p_k) / 2
  // p_k -> (q_k - p_k) / 2
  // Complex space:
  // q_k =   (h_q_k^+ + h_q_k^-) / 2
  // p_k = i (h_q_k^+ - h_q_k^-) / 2
  map.identity();
  for (auto k = 0; k < nd_tps; k++) {
    map[2*k]   = (Id[2*k]+Id[2*k+1])/2e0;
    map[2*k+1] = (Id[2*k]-Id[2*k+1])/2e0;
  }
  b = b*map;

  // q_k -> p_k
  // p_k -> q_k
  // Complex space:
  // i (q_k -/+ i p_k) = (i q_k +/- p_k)
  map.identity();
  for (auto k = 0; k < nd_tps; k++) {
    map[2*k]   = Id[2*k+1];
    map[2*k+1] = Id[2*k];
  }
  c = b*map;

  a_re = (b+c)/2e0;
  a_im = (b-c)/2e0;
}


tps RtoC_JB(tps &a_re, tps &a_im)
{
  int          k;
  tps          b;
  ss_vect<tps> Id, map;

  Id.identity();

  b = a_re + a_im;

  // q_k -> q_k + p_k
  // p_k -> q_k - p_k
  // Complex space:
  // h_q_k^+ = q_k - i h_p_k
  // h_q_k^- = q_k + i h_p_k
  map.identity();
  for (k = 0; k < nd_tps; k++) {
    map[2*k]   = Id[2*k] + Id[2*k+1];
    map[2*k+1] = Id[2*k] - Id[2*k+1];
  }
  b = b*map;
  // b = p_k_cmplx_sgn_corr_fun(nd_tps, b);
  std::cout << "\np_k_cmplx_sgn_corr_fun() not implemented\n";
  assert(false);
  return b;
}

#else

double f_p_k_cmplx_sgn_corr(const long int jj[])
{
  // Adjust the sign for the momenta for the oscillating planes.
  // Correct sign for complex vs. real momenta p_k.
  //   q_k =  (h_q_k^+ + h_q_k^-) / 2
  // i p_k = -(h_q_k^+ - h_q_k^-) / 2
  // Adjust the sign for the momenta for the oscillating planes.
  int ord, sgn = 0;

  // Compute the sum of exponents for the momenta for the oscillating planes:
  ord = 0;
  for (auto k = 0; k < nd_tps; k++)
    ord += jj[2*k+1];
  ord = (ord % 4);
  //  Sum_k c_ijkl x^i p_x^j y^k p_y^l
  //  j + l mod 4 = [0, 3: +1; 1, 2: -1]
  switch (ord) {
  case 0:
  case 3:
    sgn = 1;
    break;
  case 1:
  case 2:
    sgn = -1;
    break;
  default:
    printf("\n: undefined case %d\n", ord);
    break;
  }
  return sgn;
}


tps tps_compute_function
(const tps &a, std::function<double (const long int [])> fun)
{
  // Dacfu in Forest's LieLib.
  char     name[name_len_for+1];
  int      n;
  long int ibuf1[bufsize], ibuf2[bufsize], jj[nv_tps];
  double   rbuf[bufsize];
  tps      b;

  a.exprt(rbuf, ibuf1, ibuf2, name);
  n = rbuf[0];
  for (auto k = 1; k <= n; k++) {
    dehash_(no_tps, nv_tps, ibuf1[k-1], ibuf2[k-1], jj);
    rbuf[k] *= fun(jj);
  }
  b.imprt(n, rbuf, ibuf1, ibuf2);
  return b;
}


ctps p_k_cmplx_sgn_corr(const ctps &a)
{
  return
    ctps(tps_compute_function(a.real(), f_p_k_cmplx_sgn_corr),
	 tps_compute_function(a.imag(), f_p_k_cmplx_sgn_corr));
}


ctps CtoR(const ctps &a)
{
  // Cartesian to resonance basis:
  //   q_k =  (h_q_k^+ + h_q_k^-) / 2
  // i p_k = -(h_q_k^+ - h_q_k^-) / 2
  int          k;
  ss_vect<tps> Id, Zero;

  Id.identity();
  Zero.zero();
  map.identity();
  for (k = 0; k < nd_tps; k++) {
    map[2*k]   = (Id[2*k]+Id[2*k+1])/2e0;
    map[2*k+1] = (Id[2*k]-Id[2*k+1])/2e0;
  }
  return p_k_cmplx_sgn_corr(a)*css_vect(map, Zero);
}


ctps RtoC(const ctps &a)
{
  // Resonance to Cartesian basis.
  // h_q_k^+ = q_k - i h_p_k
  // h_q_k^- = q_k + i h_p_k
  int          k;
  ss_vect<tps> Id, Zero;

  Id.identity();
  Zero.zero();
  map.identity();
  for (k = 0; k < nd_tps; k++) {
    map[2*k]   = Id[2*k] + Id[2*k+1];
    map[2*k+1] = Id[2*k] - Id[2*k+1];
  }
  return p_k_cmplx_sgn_corr(a*css_vect(map, Zero));
}

#endif

ss_vect<tps> h_to_v(const tps &h)
{
  // Difd in Forest's F77 LieLib:
  // Compute vector flow operator from Lie operator :h:
  //   v = Omega * [del_x H, del_px H]^T
  ss_vect<tps> v;

  for (auto k = 0; k < nd_tps; k++) {
    v[2*k]   = -Der(h, 2*k+2);
    v[2*k+1] =  Der(h, 2*k+1);
  }
  return v;
}


tps v_to_tps(const ss_vect<tps> &v, const tps &x)
{
  // Daflo in Forest's F77 LieLib.
  //   y = v * nabla * x
  tps y = 0e0;

  for (auto k = 0; k < 2*nd_tps; k++)
    y += v[k]*Der(x, k+1);
  return y;
}


tps exp_v_to_tps
(const ss_vect<tps> &v, const tps &x, const double eps, const int n_max)
{
  // Expflo in Forest's F77 LieLib:
  //   y = exp(v*nabla) * x
  double eps1;
  tps    y_k, y;

  y_k = y = x;
  for (auto k = 1; k <= n_max; k++) {
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


tps exp_v_fac_to_tps
(const ss_vect<tps> &v, const tps &x, const int k1, const int k2,
 const double scl)
{
  // Facflo in Forest's F77 LieLib.
  //   y = exp(D_k1) * exp(D_k1+1) ...  * exp(D_k2) * x

  const int n_max = 100; 

  tps          y;
  ss_vect<tps> v_k;

  y = x;
  for (auto k = k1; k <= k2; k++) {
    v_k = scl*get_M_k(v, k);
    y = exp_v_to_tps(v_k, y, eps_tps, n_max);
  }
  return y;
}


ss_vect<tps> exp_v_fac_to_M(const ss_vect<tps> &v, const ss_vect<tps> &x,
			    const int k1, const int k2, const double scl)
{
  // Facflod in Forest's F77 LieLib.
  ss_vect<tps> M;

  for (auto k = 0; k < 2*nd_tps; k++)
    M[k] = exp_v_fac_to_tps(v, x[k], k1, k2, scl);
  return M;
}


ss_vect<tps> exp_v_to_map(const ss_vect<tps> &v, const ss_vect<tps> &map)
{
  const int n_max = 100; 

  ss_vect<tps> M;

  M = map;
  for (auto k = 0; k < 2*nd_tps; k++)
    M[k] = exp_v_to_tps(v, M[k], eps_tps, n_max);
  return M;
}

#if 1

ss_vect<tps>M_to_M_fact(const ss_vect<tps> &map)
{
  // Flofac in Forest's F77 LieLib.
  // Factor map:
  //   M = M_2 ... * M_n
  tps          y;
  ss_vect<tps> map_lin, map_k, map_res, map_fact;

  map_lin = get_M_k(map, 1);
  map_res = map*Inv(map_lin);

  map_fact.zero();
  for (auto k = 2; k <= no_tps; k++) {
    map_fact += get_M_k(map_res, k);
    map_k = get_M_k(-map_fact, k);
    map_res = exp_v_to_map(map_k, map_res);
  }

  return map_fact;
}

#else

ss_vect<tps>M_to_M_fact(const ss_vect<tps> &map)
{
  // Obsolete.
  // Flofac in Forest's F77 LieLib.
  // Factor map:
  //   M = M_2 ... * M_n
  ss_vect<tps> map_lin, map_res, map_fact;

  map_lin = get_M_k(map, 1);
  map_res = map*Inv(map_lin);

  map_fact.zero();
  for (auto k = 2; k <= no_tps; k++) {
    map_fact += get_M_k(map_res, k);

    for (auto j = 0; j < 2*nd_tps; j++)
      map_res[j] = exp_v_fac_to_tps(map_fact, map_res[j], k, k, -1e0);
  }

  return map_fact;
}

#endif

tps tps_fun
(const tps &a, std::function<double (const long int [])> fun)
{
  // Dacfu in Forest's F77 LieLib.
  // Multiplies mononials I_vec with function f(I_vec).
  char     name[name_len_for+1];
  int      n;
  long int ibuf1[bufsize], ibuf2[bufsize], jj[nv_tps];
  double   rbuf[bufsize];
  tps      b;

  a.exprt(rbuf, ibuf1, ibuf2, name);
  n = rbuf[0];
  for (auto k = 1; k <= n; k++) {
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
  double scl = 0e0;

  for (auto k = 0; k < 2*nd_tps; k++)
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
  tps          f_x, f_px, h;
  ss_vect<tps> Id;

  Id.identity();
  h = 0e0;
  for (auto k = 0; k < nd_tps; k++) {
    // Integrate monomials.
    f_x  = tps_fun(map[2*k], f_int_mon)*Id[2*k+1];
    f_px = tps_fun(map[2*k+1], f_int_mon)*Id[2*k];
    h -= f_x - f_px;
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


#if 1

ss_vect<tps> h_DF_to_M
(const tps &h_DF, const ss_vect<tps> &x, const int k1, const int k2)
{
  // Fexpo in Forest's F77 LieLib.
  // Compute map from Dragt-Finn factorisation:
  //   M = exp(:h_3:) * exp(:h_4:) ...  * exp(:h_n:) * X
  ss_vect<tps> v_DF;

  v_DF = h_to_v(h_DF);
  std::cout << v_DF;
  assert(false);
  return exp_v_fac_to_M(v_DF, x, k1, k2, 1e0);
}

#else

ss_vect<tps> h_DF_to_M
(const tps &h, const ss_vect<tps> &map, const int k1, const int k2)
{
  // Fexpo in Forest's LieLib.
  // Compute map from Dragt-Finn factorisation:
  //   exp(:h_3:) exp(:h_4:) ... exp(:h_no:)
  tps          h_k;
  ss_vect<tps> map1;

  map1.identity();
  for (auto k = k2; k >= k1; k--) {
    h_k = get_h_k(h, k);
    map1 = map1*LieExp(h_k, map);
  }
  return map1;
}

#endif

tps get_h(ss_vect<tps> &map)
{
  tps          h;
  ss_vect<tps> map1, R;

  if (true) {
    // Dragt-Finn factorization
    h = LieFact_DF(Inv(A0*A1)*map*A0*A1, R);
    R[6] = tps(0e0, 7);
    h = h*R;
    return h;
  } else {
    // single Lie exponent
    danot_(1); map1 = map; danot_(no_tps);
    return LieFact(Inv(A0*A1)*map*Inv(map1)*A0*A1);
  }
}


Eigen::MatrixXd get_lin_map(ss_vect<tps> &map)
{
  Eigen::MatrixXd M(6, 6);

  for (auto j = 0; j < 2*nd_tps; j++) {
    for (auto k = 0; k < 2*nd_tps; k++)
      M(j, k) = map[j][k];
  }
  return M;
}


Eigen::VectorXd compute_nu_symp(const Eigen::MatrixXd &M)
{
  // Eigenvalues for a 4x4 symplectic periodic matrix.

  const int
    n_dof = 2,
    n_dim = 2*n_dof;

  double
    x[n_dof];
  Eigen::VectorXd
    nu(n_dof);
  Eigen::MatrixXd
    Id = Eigen::MatrixXd::Identity(n_dim, n_dim);

  auto Pp1 = (M-Id).determinant();
  auto Pm1 = (M+Id).determinant();
  auto p = (Pp1-Pm1)/8e0;
  auto q = (Pp1+Pm1)/8e0 - 1e0;
  auto sgn =
    (M.block(0, 0, n_dof, n_dof).trace()
     > M.block(n_dof, n_dof, n_dof, n_dof).trace())?
    1 : -1;
  x[0] = -p/2e0 + sgn*sqrt(sqr(p/2e0)-q);
  x[1] = -p/2e0 - sgn*sqrt(sqr(p/2e0)-q);
  for (auto k = 0; k < n_dof; k++) {
    nu[k] = acos(x[k])/(2e0*M_PI);
    if (M(2*k, 2*k+1) < 0e0)
      nu[k] = 1e0 - nu[k];
  }
  std::cout << std::fixed << std::setprecision(5)
	    << "\ncompute_nu:\n  nu = [" << nu[X_] << ", " << nu[Y_] << "]\n";

  return nu;
}


int find_closest_nu(const double nu, const Eigen::VectorXcd &w)
{
  int    ind;
  double nu_k, diff;

  auto min = 1e30;
  for (auto k = 0; k < w.size(); k++) {
    nu_k = acos2(w(k).imag(), w(k).real())/(2e0*M_PI);
    diff = fabs(nu_k-nu);
    if (diff < min) {
      ind = k;
      min = diff;
    }
  }
  return ind;
}


Eigen::VectorXi sort_eigen_vec
(const Eigen::VectorXd &nu, const Eigen::VectorXcd &w)
{
  const int
    n_dof = 2,
    n_dim = 2*n_dof;

  Eigen::VectorXi order(n_dim);

  for (auto k = 0; k < nu.size(); k++) {
    order(2*k) = find_closest_nu(nu[k], w);
    order(2*k+1) = find_closest_nu(1e0-nu[k], w);
  }
  return order;
}


Eigen::MatrixXd compute_S(const int n_dof)
{
  // Remark: Beware of sign for longitudinal dim.
  const int n_dim = 2*n_dof;

  Eigen::MatrixXd S = Eigen::MatrixXd::Zero(n_dim, n_dim);
  for (auto k = 0; k < n_dof; k++) {
    S(2*k, 2*k+1) = 1e0;
    S(2*k+1, 2*k) = -1e0;
  }

  return S;
}


Eigen::VectorXd compute_dispersion(const Eigen::MatrixXd &M)
{
  const int
    n_dof = 2,
    n_dim = 2*n_dof;
  const Eigen::MatrixXd
    Id = Eigen::MatrixXd::Identity(n_dim, n_dim);

  auto D = M.col(delta_).segment(0, n_dim);
  auto M_4x4 = M.block(0, 0, n_dim, n_dim);

  return (Id-M_4x4).inverse()*D;
}


Eigen::MatrixXd compute_A0(const Eigen::MatrixXd &M)
{
  const int
    n_dof = 2;

  Eigen::MatrixXd
    A0 = Eigen::MatrixXd::Identity(6, 6);

  auto eta = compute_dispersion(M);

  if (n_dof == 2) {
    // Coasting beam - translate to momentum dependent fix point.
    for (auto k = 0; k < n_dof; k++) {
      A0(2*k, delta_)   =  eta(2*k);
      A0(2*k+1, delta_) =  eta(2*k+1);
      A0(ct_, 2*k)      =  eta(2*k+1);
      A0(ct_, 2*k+1)    = -eta(2*k);
    }
  }

  return A0;
}


Eigen::MatrixXd compute_A1(const int n_dof, Eigen::MatrixXcd &u_ord)
{
  const int
    n_dim = 2*n_dof;
  const std::complex<double>
    I = std::complex<double>(0e0, 1e0);

  Eigen::MatrixXd
    A1 = Eigen::MatrixXd::Identity(6, 6);
  Eigen::MatrixXcd
    u(n_dim, n_dim);

  auto S = compute_S(n_dof);

  // Normalise eigenvectors: A^T.omega.A = omega.
  for (auto k = 0; k < n_dof; k++) {
    auto z = u_ord.col(2*k).real().dot(S*u_ord.col(2*k).imag());
    auto sgn_im = boost::math::sign(z);
    auto scl = sqrt(fabs(z));
    auto sgn_vec = boost::math::sign(u_ord(2*k, 2*k).real());
    u.col(2*k) =
      sgn_vec*(u_ord.col(2*k).real()+sgn_im*u_ord.col(2*k).imag()*I)/scl;
    u.col(2*k+1) =
      sgn_vec*(u_ord.col(2*k+1).real()+sgn_im*u_ord.col(2*k+1).imag()*I)/scl;
  }
    
  u_ord = u;

  for (auto k = 0; k < n_dof; k++) {
    A1.block(0, 0, n_dim, n_dim).col(2*k)   = u.col(2*k).real();
    A1.block(0, 0, n_dim, n_dim).col(2*k+1) = u.col(2*k).imag();
  }

  return A1;
}


void compute_M_diag(MNFType &MNF)
{
  const int
    n_dof = 2,
    n_dim = 2*n_dof;

  Eigen::VectorXd
    nu_eig(n_dim),
    nu_eig_ord(n_dim);
  Eigen::VectorXcd
    w_ord(n_dim);
  Eigen::MatrixXcd
    u(n_dim, n_dim),
    u_ord(n_dim, n_dim);

  auto M = get_lin_map(MNF.M);

  // Check if stable.
  auto Tr_x = M.block(0, 0, n_dof, n_dof).trace();
  if (fabs(Tr_x) >= 2e0) {
    std::cout << std::scientific << std::setprecision(5)
	      << "\nEigenvalues - unstable in the horizontal plane:"
	      << " Tr{M_x} = " << Tr_x << "\n";
    return;
  }
  auto Tr_y = M.block(n_dof, n_dof, n_dof, n_dof).trace();
  if (fabs(Tr_y) >= 2e0) {
    std::cout << std::scientific << std::setprecision(5)
	      << "\nEigenvalues - unstable in the vertical plane:"
	      << " Tr{M_y} = " << Tr_y << "\n";
    return;
  }

  auto M_4x4 = M.block(0, 0, n_dim, n_dim);
  Eigen::EigenSolver<Eigen::Matrix<double, n_dim, n_dim> > s(M_4x4);

  auto nu_symp = compute_nu_symp(M_4x4);

  for (auto k = 0; k < n_dim; k++)
    nu_eig[k] =
      acos2(s.eigenvalues()(k).imag(), s.eigenvalues()(k).real())/(2e0*M_PI);

  auto order = sort_eigen_vec(nu_symp, s.eigenvalues());

  for (auto k = 0; k < n_dim; k++) {
    w_ord(k) = s.eigenvalues()(order(k));
    u_ord.col(k) = s.eigenvectors().col(order(k));
    nu_eig_ord(k) =
      acos2(w_ord(k).imag(), w_ord(k).real())/(2e0*M_PI);
  }

  auto A1_mat = compute_A1(n_dof, u_ord);
  auto A0_mat = compute_A0(M);
  auto R_mat  = (A0_mat*A1_mat).inverse()*M*A0_mat*A1_mat;

  MNF.A0 = mat2map(A0_mat);
  MNF.A1 = mat2map(A1_mat);
  MNF.R  = mat2map(R_mat);

  print_map("\nA0:\n",    MNF.A0);
  print_map("\nA1:\n",    MNF.A1);
  print_map("\nA0*A1:\n", MNF.A0*MNF.A1);
  print_map("\nR:\n",     MNF.R);
}


tps get_g(const tps nu_x, const tps nu_y, const tps &h)
{
  // Compute g = (1-R)^-1 * h 

  long int     jj1[ss_dim], jj2[ss_dim];
  double       re, im;
  tps          h_re, h_im, g_re, g_im, mn1, mn2, cotan;
  ss_vect<tps> Id;

  CtoR(h, h_re, h_im);

  for (auto k = 0; k < nv_tps; k++) {
    jj1[k] = jj2[k] = 0;
  }

  Id.identity();
  g_re = g_im = 0e0;
  for (auto i = 0; i <= no_tps; i++) {
    jj1[x_] = jj2[px_] = i;
    for (auto j = 0; j <= no_tps; j++) {
      jj1[px_] = jj2[x_] = j;
      for (auto k = 0; k <= no_tps; k++) {
	jj1[y_] = jj2[py_] = k;
	for (auto l = 0; l <= no_tps; l++) {
	  jj1[py_] = jj2[y_] = l;
	  if ((i+j+k+l <= no_tps) && ((i-j != 0) || (k-l != 0))) {
	    cotan = 1e0/tan(((i-j)*nu_x+(k-l)*nu_y)*M_PI);
	    mn1 =
	      pow(Id[x_], i)*pow(Id[px_], j)*pow(Id[y_], k)*pow(Id[py_], l);
	    mn2 =
	      pow(Id[x_], j)*pow(Id[px_], i)*pow(Id[y_], l)*pow(Id[py_], k);

	    for (auto m = 0; m <= no_tps-i-j-k-l; m++) {
	      jj1[delta_] = jj2[delta_] = m;
	      for (auto n = 0; n <= no_tps-i-j-k-l; n++) {
		jj1[ss_dim-1] = jj2[ss_dim-1] = n;
		re = h_re[jj1];
		im = h_im[jj1];
		// Compute g.
		g_re +=
		  (re-cotan*im)*(mn1+mn2)*pow(Id[delta_], m)
		  *pow(Id[ss_dim-1], n)/2e0;
		g_im +=
		  (im+cotan*re)*(mn1-mn2)*pow(Id[delta_], m)
		  *pow(Id[ss_dim-1], n)/2e0;
		h_re.pook(jj2, 0e0);
		h_im.pook(jj2, 0e0);
	      }
	    }
	  }
	}
      }
    }
  }

  return RtoC(g_re, g_im);
}


tps get_Ker(const tps &h)
{
  long int     jj[ss_dim];
  tps          h_Ke;
  ss_vect<tps> Id;

  for (auto k = 0; k < ss_dim; k++)
    jj[k] = 0;

  Id.identity();
  h_Ke = 0e0;
  for (auto i = 0; i <= no_tps; i++) {
    jj[x_] = jj[px_] = i;
    for (auto j = 0; j <= no_tps; j++) {
      jj[y_] = jj[py_] = j;
      for (auto k = 0; k <= no_tps; k++) {
	jj[delta_] = k;
	for (auto l = 0; l <= no_tps; l++) {
	  jj[ss_dim-1] = l;
	  if ((2*i+2*j+k+l <= no_tps) && ((i != 0) || (j != 0) || (k != 0))) {
	    h_Ke +=
	      h[jj]
	      *pow(Id[x_], i)*pow(Id[px_], i)
	      *pow(Id[y_], j)*pow(Id[py_], j)
	      *pow(Id[delta_], k)
	      *pow(Id[ss_dim-1], l);
	  }
	}
      }
    }
  }

  return h_Ke;
}


double get_coeff
(const tps &h, const int i_x, const int i_p_x, const int i_y, const int i_p_y,
 const int i_delta, const int i_ct, const int i_prm)
{
  const long int jj[] = {i_x, i_p_x, i_y, i_p_y, i_delta, i_ct, i_prm};

  return h[jj];
}


tps get_coeff_with_prm
(const tps &h, const int i_x, const int i_p_x, const int i_y, const int i_p_y,
 const int i_delta, const int i_ct, const int i_prm)
{
  tps coeff = 0e0;

  auto ord_max = no_tps - i_x - i_p_x - i_y - i_p_y - i_delta - i_ct - i_prm;
  for (auto ord = 0; ord < ord_max; ord++)
    coeff += get_coeff(h, i_x, i_p_x, i_y, i_p_y, i_delta, i_ct, ord);

  return coeff;
}



MNFType map_norm(const ss_vect<tps> &map)
{
  int           n;
  double        nu0[2];
  tps           hn, hn_re, hn_im, h_ke, gn, Kn;
  ss_vect<tps>  Id, A, nus, M_Fl, map2;
  MNFType       MNF;

  Id.identity();

  danot_(no_tps-1);

  MNF.M_res = MNF.M = map;

  danot_(no_tps);

#if 0
  // Find fixed point.
  GoFix(map, MNF.A0, MNF.A0_inv, no_tps);

  // Translate to fix point.
  map = MNF.A0_inv*map*MNF.A0;

  print_map("\nA0:", A0);
  print_map("\nM_map:", map);
#endif

  compute_M_diag(MNF);

  M_Fl = Inv(MNF.A0*MNF.A1)*MNF.M*MNF.A0*MNF.A1;

  print_map("\nM_Fl:\n", M_Fl);

  MNF.K = 0e0;
  for (auto k = 0; k < 2; k++) {
    nu0[k] = atan2(MNF.R[2*k][2*k+1], MNF.R[2*k][2*k]);
    if (nu0[k] < 0e0) nu0[k] += 2e0*M_PI;
    nu0[k] /= 2e0*M_PI;
    MNF.K -= M_PI*nu0[k]*(sqr(Id[2*k])+sqr(Id[2*k+1]));
  }
  std::cout << "\n";
  std::cout << std::fixed << std::setprecision(5)
       << "nu0 = (" << nu0[X_] << ", " << nu0[Y_] << ")" << "\n";

  // Coasting beam.
  MNF.K += h_ijklm(M_Fl[ct_], 0, 0, 0, 0, 1)*sqr(Id[delta_])/2e0;

  MNF.g = 0e0;
  for (auto k = 3; k <= no_tps; k++) {
    n = pow(2, k-3);

    map2 = M_Fl*Inv(MNF.R*FExpo(MNF.K, Id, 3, k-1, -1));
    hn = Intd(get_mns(map2, k-1, k-1), -1e0);
    gn = get_g(nu0[X_], nu0[Y_], hn);
    MNF.g += gn;
    CtoR(hn, hn_re, hn_im);
    Kn = RtoC(get_Ker(hn_re), get_Ker(hn_im));
    MNF.K += Kn;

    A = FExpo(gn, Id, k, k, -1);
    M_Fl = Inv(A)*M_Fl*A;
  }

  return MNF;
}


void test_map_norm(void)
{
  tps          g, g_re, g_im, K, k_re, k_im;
  ss_vect<tps> M, A0, A1, map_res;
  MNFType    MNF;

  if (true) {
    if (true)
      set_bn_par(get_Fnum("sf"), Sext, 7);
    else
      set_bn_par(get_Fnum("uq3"), Quad, 7);
  }

  danot_(no_tps-1);

  M.identity();
  M.propagate(1, n_elem);

  print_map("\nM:\n", M);

  if (false)
    std::cout << std::scientific << std::setprecision(5)
	      << std::setw(13) << M << "\n";

  danot_(no_tps);

  MNF = map_norm(M);
  CtoR(MNF.K, k_re, k_im);

  K = MapNorm(M, g, A1, A0, map_res, 1);

#if 0
  std::cout << "\nMNF.K:\n" << MNF.K;
  std::cout << "\nK:\n" << K;
#else
  std::cout << "\nMNF.K-K:\n" << MNF.K-K;
  std::cout << "\nMNF.g-g:\n" << MNF.g-g;
#endif

//  CtoR(MNF.K-K, hn_re, hn_im);
//  std::cout << hn_re;
}


void dragt_finn_fact(void)
{
  tps          h_1, h_2;
  ss_vect<tps> R;

  if (true) {
    if (true)
      set_bn_par(get_Fnum("sf"), Sext, 7);
    else
      set_bn_par(get_Fnum("uq3"), Quad, 7);
  }

  danot_(no_tps-1);

  Map.identity();
  Map.propagate(1, n_elem);

  danot_(no_tps);

  print_map("\nMap:", Map);

  if (false)
    std::cout << std::scientific << std::setprecision(5)
	      << std::setw(13) << Map << "\n";

  h_1 = LieFact_DF(Map, R);
  h_2 = M_to_h_DF(Map);

  std::cout << std::scientific << std::setprecision(5)
	    << "\nh:\n" << std::setw(13) << h_2 << "\n";
  std::cout << std::scientific << std::setprecision(5)
	    << "\nh_2-h_1:\n" << std::setw(13) << h_2-h_1 << "\n";
}

#if 0

void tst_ctor(MNFType &MNF)
{
  tps          K_re, K_im, K_re_JB, K_im_JB, g_re, g_im, g_re_JB, g_im_JB;
  ss_vect<tps> Id;
  ctps         cK, cg;

  Id.identity();

  CtoR(MNF.K, K_re, K_im);
  CtoR(MNF.g, g_re, g_im);
  cK = CtoR(ctps(MNF.K, 0e0));
  cg = CtoR(ctps(0e0, MNF.g));

  cout << "\n[K_re-K_re_JB, K_im-K_im_JB]:\n"
       << K_re-cK.real() << K_im-cK.imag();
  daeps_(1e-7);
  cout << "\n[g_re-g_re_JB, g_im-g_im_JB]:\n"
       << g_re-cg.real() << g_im-cg.imag();
  daeps_(eps_tps);

  cK = RtoC(ctps(K_re, K_im));
  cg = RtoC(ctps(g_re, g_im));

  cout << "\nRtoC_JB(2, cK)-K:\n" << cK.real()-MNF.K << cK.imag();

  daeps_(1e-6);
  cout << "\nRtoC_JB(2, cg)-g:\n" << 1e0*cg.real() << cg.imag()-MNF.g;
  daeps_(eps_tps);
}

#endif

void test_param_dep(void)
{
  ss_vect<tps> Id, M;

  Id.identity();

  M.zero();
  M[x_] = -1e0;
  for (auto k = 0; k < 2*nd_tps; k++)
    M[k] += (k+1e0)*Id[k];
  M[x_] += 3.14e-3*Id[x_]*Id[6];

  Id *= 2e0;
  print_map("\nId:", Id);
  print_map("\nM:", M);
  std::cout << std::scientific << std::setprecision(5)
	    << "\nM:\n" << std::setw(13) << M;
  std::cout << std::scientific << std::setprecision(5)
	    << "\n(M*Id)[x_]:\n" << std::setw(13) << (M*Id)[x_];
  std::cout << std::scientific << std::setprecision(5)
	    << "\n(Id*M)[x_]:\n" << std::setw(13) << (Id*M)[x_];
  std::cout << std::scientific << std::setprecision(5)
	    << "\n(M*M)[x_]:\n" << std::setw(13) << (M*M)[x_];
}


int main(int argc, char *argv[])
{

  rad_on    = false; H_exact        = false; totpath_on   = false;
  cavity_on = false; quad_fringe_on = false; emittance_on = false;
  IBS_on    = false;

  rd_mfile(argv[1], elem);
  rd_mfile(argv[1], elem_tps);

  // initialize the symplectic integrator after energy has been defined
  ini_si();

  // disable from TPSALib- and LieLib log messages
  idprset(-1);

  daeps_(1e-30);

  danot_(no_tps-1);

  if (false)
    no_mpoles(Sext);

 if (!false)
   test_map_norm();

  if (false)
    dragt_finn_fact();

  if (false) {
    danot_(no_tps);
    test_param_dep();
  }
}
