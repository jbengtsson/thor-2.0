#define NO 3

#include "thor_lib.h"

int  no_tps   = NO,
     ndpt_tps = 5;

extern tps           K, g;
extern ss_vect<tps>  Map, A0, A1, Map_res;


void prt_lin_map(const int n_DOF, const ss_vect<tps> &map)
{
  int  i, j, i1 = 0, j1 = 0;

  cout << endl;
  for (i = 0; i < 2*n_DOF; i++) {
    for (j = 0; j < 2*n_DOF; j++) {
      switch (i) {
      case 0 ... 3: i1 = i;
	break;
      case 4: i1 = ct_;
	break;
      case 5: i1 = delta_;
      }
      switch (j) {
      case 0 ... 3: j1 = j;
	break;
      case 4: j1 = ct_;
	break;
      case 5: j1 = delta_;
      }
      if (true)
	cout << scientific << setprecision(6)
	     << setw(14) << map[i1][j1];
      else
	cout << scientific << setprecision(16)
	     << setw(24) << map[i1][j1];
    }
    cout << endl;
  }
}


int main(int argc, char *argv[])
{
  double        nu[2], ksi[2];
  tps           K_re, K_im, g_re, g_im;
  ss_vect<tps>  nus;

  danot_(no_tps-1);

  rad_on    = false; H_exact        = false; totpath_on   = false;
  cavity_on = false; quad_fringe_on = false; emittance_on = false;
  IBS_on    = false;

  string home_dir = "/home/johan/projects/src/";

  rd_mfile((home_dir+"flat_file.dat").c_str(), elem);
  rd_mfile((home_dir+"flat_file.dat").c_str(), elem_tps);

  // initialize the symplectic integrator after energy has been defined
  ini_si();

  // disable from TPSALib- and LieLib log messages
  idprset(-1);

  danot_(no_tps-1);
  get_Map();
  danot_(no_tps);

  prt_lin_map(3, Map);

  K = MapNorm(Map, g, A1, A0, Map_res, no_tps);
  CtoR(K, K_re, K_im); CtoR(g, g_re, g_im);
  nus = dHdJ(K); get_nu_ksi(nus, nu, ksi);

  cout << endl;
  prt_nu(nus);

  cout << endl;
  cout << scientific << setprecision(6)
       << "R_56:    " << setw(14) << h_ijklm(Map[ct_], 0, 0, 0, 0, 1) << endl;
  cout << scientific << setprecision(6)
       << "T_566:   " << setw(14) << h_ijklm(Map[ct_], 0, 0, 0, 0, 2) << endl;
  cout << scientific << setprecision(6) << setw(13) << K << endl;
  cout << scientific << setprecision(6) << setw(13) << g_im << endl;
}
