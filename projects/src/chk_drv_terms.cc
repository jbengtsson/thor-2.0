#define NO 4

#include<string>

#include "thor_lib.h"

int
  no_tps   = NO,
  ndpt_tps = 5;


inline long int* vec2arr(const std::vector<long int> &vec)
{
  static long int jj[ss_dim];
  for (auto k = 0; k < vec.size(); k++)
    jj[k] = vec[k];
  return jj;
 }


void prt_drv_term
(const std::vector<long int> &vec, const tps &h_re, const tps &h_im)
{
  std::string index;

  auto jj = vec2arr(vec);
  for (auto k = 0; k < 6; k++)
    index += '0' + jj[k];
  std::cout << std::scientific << std::setprecision(16)
	    << "h_" << index << " = [" << std::setw(16) << h_re[jj] << ", "
	    << std::setw(16) <<  h_im[jj] << "]\n";
}


void sxt_h1(const tps &h_re, const tps &h_im)
{
  // Linear chromaticity.
  std::cout << "\n";
  prt_drv_term({1, 1, 0, 0, 1, 0, 0}, h_re, h_im);
  prt_drv_term({0, 0, 1, 1, 1, 0, 0}, h_re, h_im);

  // First order chromatic terms.
  std::cout << "\n";
  prt_drv_term({2, 0, 0, 0, 1, 0, 0}, h_re, h_im);
  prt_drv_term({0, 0, 2, 0, 1, 0, 0}, h_re, h_im);
  prt_drv_term({1, 0, 0, 0, 2, 0, 0}, h_re, h_im);

  // Normal sextupoles.
  std::cout << "\n";
  prt_drv_term({2, 1, 0, 0, 0, 0, 0}, h_re, h_im);
  prt_drv_term({3, 0, 0, 0, 0, 0, 0}, h_re, h_im);
  prt_drv_term({1, 0, 1, 1, 0, 0, 0}, h_re, h_im);
  prt_drv_term({1, 0, 2, 0, 0, 0, 0}, h_re, h_im);
  prt_drv_term({1, 0, 0, 2, 0, 0, 0}, h_re, h_im);

  if (false) {
    // Skew sextupoles.
    std::cout << "\n";
    prt_drv_term({0, 0, 2, 1, 0, 0, 0}, h_re, h_im);
    prt_drv_term({0, 0, 3, 0, 0, 0, 0}, h_re, h_im);
    prt_drv_term({1, 1, 1, 0, 0, 0, 0}, h_re, h_im);
    prt_drv_term({0, 2, 1, 0, 0, 0, 0}, h_re, h_im);
    prt_drv_term({2, 0, 1, 0, 0, 0, 0}, h_re, h_im);
  }
}


void sxt_h2(const tps &h_re, const tps &h_im, const tps K_re)
{
  long int jj[ss_dim];

  for (int k = 0; k < ss_dim; k++)
    jj[k] = 0;

  jj[x_] = 4; jj[px_] = 0; jj[y_] = 0; jj[py_] = 0;
  std::cout << std::scientific << std::setprecision(16)
	    << "\nh_40000:" << std::setw(24) << h_re[jj]
	    << std::setw(24) << h_im[jj] << "\n";
  jj[x_] = 3; jj[px_] = 1; jj[y_] = 0; jj[py_] = 0;
  std::cout << std::scientific << std::setprecision(16)
	    << "h_31000:" << std::setw(24) << h_re[jj]
	    << std::setw(24) << h_im[jj] << "\n";
  jj[x_] = 2; jj[px_] = 0; jj[y_] = 1; jj[py_] = 1;
  std::cout << std::scientific << std::setprecision(16)
	    << "h_20110:" << std::setw(24) << h_re[jj]
	    << std::setw(24) << h_im[jj] << "\n";
  jj[x_] = 2; jj[px_] = 0; jj[y_] = 0; jj[py_] = 2;
  std::cout << std::scientific << std::setprecision(16)
	    << "h_20020:" << std::setw(24) << h_re[jj]
	    << std::setw(24) << h_im[jj] << "\n";
  jj[x_] = 2; jj[px_] = 0; jj[y_] = 2; jj[py_] = 0;
  std::cout << std::scientific << std::setprecision(16)
	    << "h_20200:" << std::setw(24) << h_re[jj]
	    << std::setw(24) << h_im[jj] << "\n";
  jj[x_] = 1; jj[px_] = 1; jj[y_] = 2; jj[py_] = 0;
  std::cout << std::scientific << std::setprecision(16)
	    << "h_11200:" << std::setw(24) << h_re[jj]
	    << std::setw(24) << h_im[jj] << "\n";
  jj[x_] = 0; jj[px_] = 0; jj[y_] = 3; jj[py_] = 1;
  std::cout << std::scientific << std::setprecision(16)
	    << "h_00310:" << std::setw(24) << h_re[jj]
	    << std::setw(24) << h_im[jj] << "\n";
  jj[x_] = 0; jj[px_] = 0; jj[y_] = 4; jj[py_] = 0;
  std::cout << std::scientific << std::setprecision(16)
	    << "h_00400:" << std::setw(24) << h_re[jj]
	    << std::setw(24) << h_im[jj] << "\n";

  jj[x_] = 2; jj[px_] = 2; jj[y_] = 0; jj[py_] = 0;
  std::cout << std::scientific << std::setprecision(16)
	    << "\nh_22000:" << std::setw(24) << h_re[jj]
	    << std::setw(24) << h_im[jj]
	    << "\n";
  jj[x_] = 1; jj[px_] = 1; jj[y_] = 1; jj[py_] = 1;
  std::cout << std::scientific << std::setprecision(16)
	    << "h_11110:" << std::setw(24) << h_re[jj]
	    << std::setw(24) << h_im[jj] << "\n";
  jj[x_] = 0; jj[px_] = 0; jj[y_] = 2; jj[py_] = 2;
  std::cout << std::scientific << std::setprecision(16)
	    << "h_00220:" << std::setw(24) << h_re[jj]
	    << std::setw(24) << h_im[jj] << "\n";

  jj[x_] = 2; jj[px_] = 2; jj[y_] = 0; jj[py_] = 0;
  std::cout << std::scientific << std::setprecision(16) << "\na_xx ="
	    << std::setw(24) << K_re[jj]/M_PI << "\n";
  jj[x_] = 1; jj[px_] = 1; jj[y_] = 1; jj[py_] = 1;
  std::cout << std::scientific << std::setprecision(16) << "a_xy ="
	    << std::setw(24) << K_re[jj]/M_PI << "\n";
  jj[x_] = 0; jj[px_] = 0; jj[y_] = 2; jj[py_] = 2;
  std::cout << std::scientific << std::setprecision(16) << "a_yy ="
	    << std::setw(24) << K_re[jj]/M_PI << "\n";
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
  int           k;
  ss_vect<tps>  y;

  for (k = 0; k < nv_tps; k++)
    y[k] = get_mns(x[k], no1, no2);

  return y;
}


tps LieFact_JB(const ss_vect<tps> &map)
{
  /* Dragt-Finn factorization:

       M = exp(:h_no:)...exp(:h_4:)exp(:h_3:)R

     Input is map in Floquet space.                                           */

  int          k;
  tps          h, hn;
  ss_vect<tps> Id, A0_inv, R, Fn, map1;

  Id.identity();

  danot_(1);
  R = map;
  danot_(no_tps);
  map1 = map*Inv(R);

  h = 0e0;
  for (k = 3; k <= no_tps; k++) {
    Fn = Taked(map1, k-1);
    hn = Intd(Fn, -1e0);
    h += hn;
    map1 = map1*FExpo(-hn, Id, k, k, -1);
  }

  return h;
}


tps get_h_local(const ss_vect<tps> &map, const bool dragt_finn)
{
  ss_vect<tps>  map1, R;

  if (dragt_finn)
    // Dragt-Finn factorization.
#if 1
    return LieFact_DF(map, R);
#else
    return LieFact_JB(map);
#endif
  else {
    // Single Lie exponent.
    danot_(1);
    R = map;
    danot_(no_tps);
#if 1
    return LieFact(map*Inv(R));
#else
    return LieFact(map);
#endif
  }
}


void compute_drv_terms(void)
{
  tps          K, K_re, K_im, h, h_re, h_im;
  ss_vect<tps> Id_scl, M_Fl;

  Id_scl.identity();
  // Id_scl[x_] *= sqrt(2e0*Jx); Id_scl[px_] *= sqrt(2e0*Jx);
  // Id_scl[y_] *= sqrt(2e0*Jy); Id_scl[py_] *= sqrt(2e0*Jy);
  // Id_scl[delta_] *= delta;

  danot_(no_tps-1);

  get_Map();

  prt_lin_map(3, Map);

  danot_(no_tps);

  K = MapNorm(Map, g, A1, A0, Map_res, 1);
  CtoR(K, K_re, K_im);

  printf("\nA0");
  prt_lin_map(3, A0);

  printf("\nA1");
  prt_lin_map(3, A1);

  M_Fl = Inv(A0*A1)*Map*A0*A1;

  h = get_h_local(M_Fl, false);
  h = h*Id_scl;
  CtoR(h, h_re, h_im);

  sxt_h1(h_re, h_im);
  sxt_h2(h_re, h_im, K_re);
}


void set_lat_state(void)
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

  set_lat_state();

  rd_mfile(argv[1], elem); rd_mfile(argv[1], elem_tps);

  //Iinitialize the symplectic integrator after energy has been defined.
  ini_si();

  // Disable from TPSALib & LieLib log messages.
  idprset(-1);

  compute_drv_terms();
}
