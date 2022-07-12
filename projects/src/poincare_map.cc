#define NO 4

#include "thor_lib.h"

int
  no_tps   = NO,
  ndpt_tps = 5;


ss_vect<tps> get_map_Fl(const ss_vect<tps> &map)
{
  const int no_delta = 1;

  K = MapNorm(Map, g, A1, A0, Map_res, no_delta);
  return Inv(A0*A1)*map*A0*A1;
}


void Dragt_Finn_Fact(const ss_vect<tps> &map, ss_vect<tps> &M, tps &h)
{
  h = LieFact_DF(map, M);
}


tps get_k_2(const ss_vect<tps> &R)
{
  int          k;
  double       mu[2];
  tps          k_2, k_2_re, k_2_im;
  ss_vect<tps> Id;

  Id.identity();

  k_2_re = 0e0;
  for (k = 0; k < 2; k++) {
    mu[k] = atan2(R[2*k][2*k+1], R[2*k][2*k]);
    if (mu[k] < 0e0) mu[k] += 2e0*M_PI;
    k_2_re -= mu[k]*Id[2*k]*Id[2*k+1]/2e0;
  }
  printf("\nget_k_2: nu = [%5.3f, %5.3f]\n",
	 mu[X_]/(2e0*M_PI), mu[Y_]/(2e0*M_PI));
  k_2_im = 0e0;
  k_2 = RtoC(k_2_re, k_2_im);
  return k_2;
}


void analyse(void)
{
  tps          h, h_re, h_im, k_2;
  ss_vect<tps> Id, map_Fl, map_Fl_2, R;

  Id.identity();

  danot_(no_tps-1);
  get_Map();
  danot_(no_tps);

  printf("\nM:\n");
  prt_lin_map(3, Map);

  map_Fl = get_map_Fl(Map);

  Dragt_Finn_Fact(map_Fl, R, h);

  k_2 = get_k_2(R);
  std::cout << std::scientific << std::setprecision(3) << "\nk_2:\n"
	    << std::setw(11) << k_2 << "\n";

  std::cout << std::scientific << std::setprecision(3) << "\nh:\n"
	    << std::setw(11) << h << "\n";

  h = BCH(k_2, h, 7);
  std::cout << std::scientific << std::setprecision(3) << "\nh:\n"
	    << std::setw(11) << h << "\n";
}


int main()
{
  danot_(no_tps-1);

  rad_on    = false; H_exact        = false; totpath_on   = false;
  cavity_on = false; quad_fringe_on = false; emittance_on = false;
  IBS_on    = false;


  rd_mfile("flat_file.dat", elem); rd_mfile("flat_file.dat", elem_tps);

  // initialize the symplectic integrator after energy has been defined
  ini_si();

  // disable from TPSALib- and LieLib log messages
  idprset(-1);

  analyse();
}
