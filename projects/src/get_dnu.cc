#define NO 5

#include "thor_lib.h"

int no_tps   = NO,
    ndpt_tps = 5;


const double
  // beta_inj[] = {8.7, 2.1},
  // A_max[]    = {5e-3, 1.5e-3},
  // delta_max  = 3e-2,
  beta_inj[] = {10.7, 6.5},
  A_max[]    = {3.5e-3, 1.5e-3},
  delta_max  = 2e-2,
  // beta_inj[] = {4.1, 1.8},
  // A_max[]    = {4e-3, 2.5e-3},
  // delta_max  = 4e-2,
  twoJ[]     = {sqr(A_max[X_])/beta_inj[X_], sqr(A_max[Y_])/beta_inj[Y_]};

const char home_dir[] = "/home/bengtsson";

double       alpha_rad_x, nu0[2];
ss_vect<tps> A_inv, nus;


void get_dnu(const ss_vect<tps> &A, double dnu[])
{   
  int  k;

  for (k = 0; k <= 1; k++) {
    dnu[k] = atan2(A[2*k][2*k+1], A[2*k][2*k])/(2.0*M_PI);
    if (dnu[k] < 0.0) dnu[k] += 1.0;
  }
} 


void get_ab(const ss_vect<tps> A1, tps ab[])
{
  /* parameter dependance for the periodic solution:


        -1       T           |  0  I |        T   | beta   -alpha |
       A   = -S A  S,    S = |       |,    A A  = |               |
                             | -I  0 |            | -alpha  gamma |

       alpha_x = -h_01000(A_Atp[x_]), alpha_y = -h_00010(A_Atp[y_]),
       beta_x  =  h_10000(A_Atp[x_]), beta_y  =  h_00100(A_Atp[y_])

       eta_x   =  h_00001(A[x_]),     eta_y   =  h_00001(A[y_]),
       eta'_x  =  h_00001(A[px_]),    eta'_y  =  h_00001(A[py_])

  */

  ss_vect<tps>  A1_A1tp;

  A1_A1tp = A1*tp_S(2, A1); ab[X_] = A1_A1tp[x_]; ab[Y_] = A1_A1tp[y_];
}


void calc_twiss(void)
{
  long int  j;
  int       k;
  double    dnu0[2], dnu[2];
  tps       ab[2];

  elem[0].Nu[X_] = 0.0; elem[0].Nu[Y_] = 0.0;
  for (j = 0; j < n_elem; j++) {
    for (k = 0; k <= 1; k++) {
      elem[j].Eta[k]  = h_ijklm(elem_tps[j].A1[2*k],  0, 0, 0, 0, 1);
      elem[j].Etap[k] = h_ijklm(elem_tps[j].A1[2*k+1], 0, 0, 0, 0, 1);
    }

    get_ab(elem_tps[j].A1, ab);

    elem[j].Alpha[X_] = -h_ijklm(ab[X_], 0, 1, 0, 0, 0);
    elem[j].Alpha[Y_] = -h_ijklm(ab[Y_], 0, 0, 0, 1, 0);
    elem[j].Beta[X_]  =  h_ijklm(ab[X_], 1, 0, 0, 0, 0);
    elem[j].Beta[Y_]  =  h_ijklm(ab[Y_], 0, 0, 1, 0, 0);

    if (j > 0) {
      get_dnu(elem_tps[j-1].A1, dnu0); get_dnu(elem_tps[j].A1, dnu);
      for (k = 0; k <= 1; k++) {
	dnu[k] -= dnu0[k];
	if ((dnu[k] < 0.0) && (elem[j].L > 0.0)) dnu[k] += 1.0;
	elem[j].Nu[k] = elem[j-1].Nu[k] + dnu[k];
      }

      elem[j].S = elem[j-1].S + elem[j].L;
    }
  }
}


void get_twiss(const long int i1, const long int i2)
{
  // computes linear dispersion for delta = 0
  ss_vect<tps>  A1_prm;

  emittance_on = true; rad_on = false;

  danot_(no_tps-1); Map.identity(); Map.propagate(i1, i2);

  danot_(no_tps); K = MapNorm(Map, g, A1, A0, Map_res, 1);

  nus_ = dHdJ(K); get_nu_ksi(nus_, nu_, ksi_);

  // compute A and diffusion coeffs
  A1_prm = LieExp(g, A0*A1); A1_prm.propagate(i1, i2);

  // radiation damping is constant
  eps_[X_] = -D_[X_]/(4.0*alpha_rad_x);
}


void get_map_n(const int n)
{
  ss_vect<tps> map2, map4;

  get_Map();

  switch (n) {
  case 1:
    break;
  case 2:
    Map = Map*Map;
    break;
  case 3:
    map2 = Map*Map; Map = map2*Map;
    break;
  case 4:
    map2 = Map*Map; Map = map2*map2;
    break;
  case 5:
    map2 = Map*Map; Map = map2*map2*Map;
    break;
  case 6:
    map2 = Map*Map; Map = map2*map2*map2;
    break;
  case 7:
    map2 = Map*Map; map4 = map2*map2; Map = map4*map2*Map;
    break;
  case 8:
    map2 = Map*Map; map4 = map2*map2; Map = map4*map4;
    break;
  case 9:
    map2 = Map*Map; map4 = map2*map2; Map = map4*map4*Map;
    break;
  case 10:
    map2 = Map*Map; map4 = map2*map2; Map = map4*map4*map2;
    break;
  case 11:
    map2 = Map*Map; map4 = map2*map2; Map = map4*map4*map2*Map;
    break;
  case 12:
    map2 = Map*Map; map4 = map2*map2; Map = map4*map4*map4;
    break;
  case 13:
    map2 = Map*Map; map4 = map2*map2; Map = map4*map4*map4*Map;
    break;
  case 14:
    map2 = Map*Map; map4 = map2*map2; Map = map4*map4*map4*map2;
    break;
  case 15:
    map2 = Map*Map; map4 = map2*map2; Map = map4*map4*map4*map2*Map;
    break;
  default:
    std::cout << "get_map_n: n not defined " << n << std::endl;
    exit(1);
    break;
  }
}


void get_map_normal_form()
{

  danot_(no_tps);

  K = MapNorm(Map, g, A1, A0, Map_res, no_tps);
}


void get_A(void)
{
  int          j;
  long int     jj[ss_dim];
  tps          gn;
  ss_vect<tps> Id, A;

  Id.identity(); A = A1;
  for (j = no_tps; j >= 3; j--) {
    gn = Take(g, j); A = A*LieExp(gn, Id);
  }

  for (j = 0; j < ss_dim; j++)
    jj[j] = (j < 4)? 1 : 0;
  A_inv = PInv(A, jj);
}


void get_twoJ(const ss_vect<double> &ps, double twoJ[])
{
  int             j;
  ss_vect<double> z;
  ss_vect<tps>    Id;

  z = (A_inv*ps).cst();

  for (j = 0; j < 2; j++)
    twoJ[j] = sqr(z[2*j]) + sqr(z[2*j+1]);
}


void get_dnu(const double Ax_max, const double Ay_max, const double delta_max)
{
  char            str[max_str];
  int             i;
  double          twoJ[2], nux, nuy;
  ss_vect<double> ps;
  ss_vect<tps>    Id;
  std::ifstream   inf;
  std::ofstream   outf;

  const int    n_ampl = 25, n_delta = 20;
  const double A_min = 1e-6;

  if (false) {
    sprintf(str, "%s%s", home_dir, "/Thor-2.0/thor/wrk");
    file_rd(inf, strcat(str, "/nus.dat"));
    inf >> nus[3] >> nus[4];
    inf.close();
  }

//  sprintf(str, "%s%s", home_dir, "/projects/src/");
//  file_wr(outf, strcat(str, "dnu_dAx_pert.out"));
  file_wr(outf, "dnu_dAx_pert.out");
  Id.zero(); ps.zero();
  for (i = -n_ampl; i <= n_ampl; i++) {
    ps[x_] = i*Ax_max/n_ampl;
    if (ps[x_] == 0.0) ps[x_] = A_min;
    ps[y_] = A_min;
    get_twoJ(ps, twoJ);
    Id[x_] = sqrt(twoJ[X_]); Id[px_] = sqrt(twoJ[X_]);
    Id[y_] = sqrt(twoJ[Y_]); Id[py_] = sqrt(twoJ[Y_]);
    nux = (nus[3]*Id).cst(); nuy = (nus[4]*Id).cst();

    outf << std::scientific << std::setprecision(3)
	 << std::setw(12) << 1e3*ps[x_] << std::setw(12) << 1e3*ps[y_]
	 << std::fixed << std::setprecision(5)
	 << std::setw(9) << nux << std::setw(9) << nuy << std::endl;
  }
  outf.close();

//  sprintf(str, "%s%s", home_dir, "/projects/src/");
//  file_wr(outf, strcat(str, "dnu_dAy_pert.out"));
  file_wr(outf, "dnu_dAy_pert.out");
  Id.zero(); ps.zero();
  for (i = -n_ampl; i <= n_ampl; i++) {
    ps[x_] = A_min;
    ps[y_] = i*Ay_max/n_ampl;
    if (ps[y_] == 0.0) ps[y_] = A_min;
//    get_twoJ(2, ps, A1, twoJ);
    get_twoJ(ps, twoJ);
    Id[x_] = sqrt(twoJ[X_]); Id[px_] = sqrt(twoJ[X_]);
    Id[y_] = sqrt(twoJ[Y_]); Id[py_] = sqrt(twoJ[Y_]);
    nux = (nus[3]*Id).cst(); nuy = (nus[4]*Id).cst();

    outf << std::scientific << std::setprecision(3)
	 << std::setw(12) << 1e3*ps[x_] << std::setw(12) << 1e3*ps[y_]
	 << std::fixed << std::setprecision(6)
	 << std::setw(10) << nux << std::setw(10) << nuy << std::endl;
  }
  outf.close();

//  sprintf(str, "%s%s", home_dir, "/projects/src/");
//  file_wr(outf, strcat(str, "chrom2_pert.out"));
  file_wr(outf, "chrom2_pert.out");
  Id.zero(); ps.zero();
  for (i = -n_delta; i <= n_delta; i++) {
    ps[delta_] = i*delta_max/n_delta; Id[delta_] = ps[delta_];

    nux = (nus[3]*Id).cst(); nuy = (nus[4]*Id).cst();

    outf << std::scientific << std::setprecision(3)
	 << std::setw(12) << 1e2*ps[delta_]
	 << std::fixed << std::setprecision(5)
	 << std::setw(9) << nux << std::setw(9) << nuy << std::endl;
  }
  outf.close();
}


void get_dnu_x_delta(const double Ax_max, const double Ay,
		     const double delta_max)
{
  int             i, j;
  double          twoJ[2], nux, nuy;
  ss_vect<double> ps;
  ss_vect<tps>    Id;
  std::ifstream   inf;
  std::ofstream   outf;

  const int n_ampl = 10, n_delta = 5;

  file_wr(outf, "dnu_dA_x_delta_pert.out");
  Id.zero();
  ps.zero(); ps[y_] = Ay;
  for (i = -n_delta; i <= n_delta; i++) {
    ps[delta_] = i*delta_max/n_delta;
    for (j = -n_ampl; j <= n_ampl; j++) {
      ps[x_] = j*Ax_max/n_ampl;
      get_twoJ(ps, twoJ);
      Id[x_] = sqrt(twoJ[X_]); Id[px_] = sqrt(twoJ[X_]);
      Id[y_] = sqrt(twoJ[Y_]); Id[py_] = sqrt(twoJ[Y_]);
      Id[delta_] = ps[delta_];
      nux = (nus[3]*Id).cst(); nuy = (nus[4]*Id).cst();

      outf << std::scientific << std::setprecision(3)
	   << std::setw(12) << 1e3*ps[x_] << std::setw(12) << ps[delta_]
	   << std::fixed << std::setprecision(5)
	   << " " << std::setw(8) << nux << " " << std::setw(8) << nuy
	   << std::endl;
    }
    outf << std::endl;
  }
  outf.close();
}


void get_dnu2(const double Ax_max, const double Ay_max, const double delta)
{
  char            str[max_str];
  int             i, j;
  double          twoJ[2], nux, nuy;
  ss_vect<double> ps;
  ss_vect<tps>    Id;
  std::ifstream   inf;
  std::ofstream   outf;

  const int n_ampl = 10;

  if (false) {
    sprintf(str, "%s%s", home_dir, "/Thor-2.0/thor/wrk");
    file_rd(inf, strcat(str, "/nus.dat"));
    inf >> nus[3] >> nus[4];
    inf.close();
  }

  file_wr(outf, "dnu_dAxy_pert.out");
  Id.zero(); Id[delta_] = delta;
  ps.zero();
  for (i = -n_ampl; i <= n_ampl; i++) {
    for (j = -n_ampl; j <= n_ampl; j++) {
      ps[x_] = i*Ax_max/n_ampl; ps[y_] = j*Ay_max/n_ampl;
      get_twoJ(ps, twoJ);
      Id[x_] = sqrt(twoJ[X_]); Id[px_] = sqrt(twoJ[X_]);
      Id[y_] = sqrt(twoJ[Y_]); Id[py_] = sqrt(twoJ[Y_]);
      nux = (nus[3]*Id).cst(); nuy = (nus[4]*Id).cst();

      outf << std::scientific << std::setprecision(3)
	   << std::setw(12) << 1e3*ps[x_] << std::setw(12) << 1e3*ps[y_]
	   << std::fixed << std::setprecision(5)
	   << " " << std::setw(8) << nux << " " << std::setw(8) << nuy
	   << std::endl;
    }
    outf << std::endl;
  }
  outf.close();
}


void wtf()
{
  long int jj[ss_dim];
  int      k;
  tps      g_re, g_im, K_re, K_im;
  
  get_twiss(1, n_elem); calc_twiss();

  get_Map();

  danot_(no_tps);
  get_map_normal_form(); nus = dHdJ(K);

  CtoR(K, K_re, K_im); CtoR(g, g_re, g_im);

  std::cout << std::scientific << std::setprecision(5) << K_re;
  for (k = 0; k < ss_dim; k++)
    jj[k] = 0;
  jj[x_] = 1; jj[px_] = 1; jj[delta_] = 1;
  printf("\nTune Shift Terms:\n  h_11001 = %12.5e\n", -2e0*K_re[jj]);
  jj[x_] = 0; jj[px_] = 0; jj[y_] = 1; jj[py_] = 1;
  printf("  h_00111 = %12.5e\n", -2e0*K_re[jj]);
  jj[x_] = 2; jj[px_] = 2; jj[y_] = 0; jj[py_] = 0; jj[delta_] = 0;
  printf("  h_22000 = %12.5e\n", -4e0*K_re[jj]);
  jj[x_] = 0; jj[px_] = 0; jj[y_] = 2; jj[py_] = 2;
  printf("  h_00220 = %12.5e\n", -4e0*K_re[jj]);
  jj[x_] = 1; jj[px_] = 1; jj[y_] = 1; jj[py_] = 1;
  printf("  h_11110 = %12.5e\n", -4e0*K_re[jj]);
}


void prt_drv_terms(const tps &h_re, const tps &h_im, const tps &K_re,
		   const ss_vect<tps> &nus)
{
  printf("\nLinear chromaticity:\n");
  printf("  ksi1    = [%10.3e, %10.3e]\n",
	 h_ijklm(nus[3], 0, 0, 0, 0, 1), h_ijklm(nus[4], 0, 0, 0, 0, 1));

  printf("\nFirst order chromatic terms:\n");
  printf("\n  h_20001 = [%10.3e, %10.3e]\n",
	 h_ijklm(h_re, 2, 0, 0, 0, 1), h_ijklm(h_im, 2, 0, 0, 0, 1));
  printf("  h_00201 = [%10.3e, %10.3e]\n",
	 h_ijklm(h_re, 0, 0, 2, 0, 1), h_ijklm(h_im, 0, 0, 2, 0, 1));
  printf("  h_10002 = [%10.3e, %10.3e]\n",
	 h_ijklm(h_re, 1, 0, 0, 0, 2), h_ijklm(h_im, 1, 0, 0, 0, 2));

  printf("\nFirst order geometric terms:\n");
  printf("  h_30000 = [%10.3e, %10.3e]\n",
	 h_ijklm(h_re, 3, 0, 0, 0, 0), h_ijklm(h_im, 3, 0, 0, 0, 0));
  printf("  h_21000 = [%10.3e, %10.3e]\n",
	 h_ijklm(h_re, 2, 1, 0, 0, 0), h_ijklm(h_im, 2, 1, 0, 0, 0));
  printf("  h_10200 = [%10.3e, %10.3e]\n",
	 h_ijklm(h_re, 1, 0, 2, 0, 0), h_ijklm(h_im, 1, 0, 2, 0, 0));
  printf("  h_10110 = [%10.3e, %10.3e]\n",
	 h_ijklm(h_re, 1, 0, 1, 1, 0), h_ijklm(h_im, 1, 0, 1, 1, 0));
  printf("  h_10020 = [%10.3e, %10.3e]\n",
	 h_ijklm(h_re, 1, 0, 0, 2, 0), h_ijklm(h_im, 1, 0, 0, 2, 0));
 
  printf("\nSecond order geometric terms:\n");
  printf("  h_40000 = [%10.3e, %10.3e]\n",
	 h_ijklm(h_re, 4, 0, 0, 0, 0), h_ijklm(h_im, 4, 0, 0, 0, 0));
  printf("  h_31000 = [%10.3e, %10.3e]\n",
	 h_ijklm(h_re, 3, 1, 0, 0, 0), h_ijklm(h_im, 3, 1, 0, 0, 0));
  printf("  h_22000 = [%10.3e, %10.3e]\n",
	 h_ijklm(h_re, 2, 2, 0, 0, 0), h_ijklm(h_im, 2, 2, 0, 0, 0));
  printf("  h_20200 = [%10.3e, %10.3e]\n",
	 h_ijklm(h_re, 2, 0, 2, 0, 0), h_ijklm(h_im, 2, 0, 2, 0, 0));
  printf("  h_11200 = [%10.3e, %10.3e]\n",
	 h_ijklm(h_re, 1, 1, 2, 0, 0), h_ijklm(h_im, 1, 1, 2, 0, 0));
  printf("  h_20110 = [%10.3e, %10.3e]\n",
	 h_ijklm(h_re, 2, 0, 1, 1, 0), h_ijklm(h_im, 2, 0, 1, 1, 0));
  printf("  h_11110 = [%10.3e, %10.3e]\n",
	 h_ijklm(h_re, 1, 1, 1, 1, 0), h_ijklm(h_im, 1, 1, 1, 1, 0));
  printf("  h_00310 = [%10.3e, %10.3e]\n",
	 h_ijklm(h_re, 0, 0, 3, 1, 0), h_ijklm(h_im, 0, 0, 3, 1, 0));
  printf("  h_20020 = [%10.3e, %10.3e]\n",
	 h_ijklm(h_re, 2, 0, 0, 2, 0), h_ijklm(h_im, 2, 0, 0, 2, 0));
  printf("  h_00400 = [%10.3e, %10.3e]\n",
	 h_ijklm(h_re, 0, 0, 4, 0, 0), h_ijklm(h_im, 0, 0, 4, 0, 0));
  printf("  h_00220 = [%10.3e, %10.3e]\n",
	 h_ijklm(h_re, 0, 0, 2, 2, 0), h_ijklm(h_im, 0, 0, 2, 2, 0));

  printf("\nSecond order anharmonic terms:\n");
  printf("  k_22000 = %10.3e\n", h_ijklm(K_re, 2, 2, 0, 0, 0));
  printf("  k_11110 = %10.3e\n", h_ijklm(K_re, 1, 1, 1, 1, 0));
  printf("  k_00220 = %10.3e\n", h_ijklm(K_re, 0, 0, 2, 2, 0));

  printf("\nSecond order chromaticity:\n");
  printf("  ksi2    = [%10.3e, %10.3e]\n",
	 h_ijklm(nus[3], 0, 0, 0, 0, 2), h_ijklm(nus[4], 0, 0, 0, 0, 2));

  printf("\nSecond order cross terms:\n");
  printf("  k_22001 = %10.3e\n", h_ijklm(K_re, 2, 2, 0, 0, 1));
  printf("  k_11111 = %10.3e\n", h_ijklm(K_re, 1, 1, 1, 1, 1));
  printf("  k_00221 = %10.3e\n", h_ijklm(K_re, 0, 0, 2, 2, 1));

  printf("\nThird order chromaticity:\n");
  printf("  ksi3    = [%10.3e, %10.3e]\n",
	 h_ijklm(nus[3], 0, 0, 0, 0, 3), h_ijklm(nus[4], 0, 0, 0, 0, 3));
}


int main(int argc, char *argv[])
{
  tps           h_re, h_im, g_re, g_im, K_re, K_im;
  ss_vect<tps>  Id_scl;
  std::ofstream outf;

  rad_on    = false; H_exact        = false; totpath_on   = false;
  cavity_on = false; quad_fringe_on = false; emittance_on = false;
  IBS_on    = false;

  rd_mfile(argv[1], elem); rd_mfile(argv[1], elem_tps);

  // Initialize the symplectic integrator after energy has been defined.
  ini_si();

  // Disable log messages from TPSALib and LieLib.
  idprset(-1);

  if (false) {
    wtf();
    exit(0);
  }

  danot_(1);

  if (!false) {
    // get_twiss(1, n_elem); calc_twiss();
    prt_lat("linlat.out");
  }

  danot_(no_tps-1);

//  get_map_n(3);
  get_Map();

  danot_(no_tps);
  get_map_normal_form(); nus = dHdJ(K);
  CtoR(K, K_re, K_im); CtoR(g, g_re, g_im);

  Id_scl.identity();
  Id_scl[x_] *= sqrt(twoJ[X_]); Id_scl[px_] *= sqrt(twoJ[X_]);
  Id_scl[y_] *= sqrt(twoJ[Y_]); Id_scl[py_] *= sqrt(twoJ[Y_]);
  Id_scl[delta_] *= delta_max;

  if (true) {
    CtoR(get_h(), h_re, h_im);

    Id_scl.identity();
    prt_drv_terms(h_re*Id_scl, h_im*Id_scl, K_re*Id_scl, nus*Id_scl);

    outf.open("h.dat", std::ios::out);
    outf << h_re << h_im;
    outf.close();

    outf.open("K.dat", std::ios::out);
    outf << K_re;
    outf.close();

    outf.open("g.dat", std::ios::out);
    outf << g_re << g_im;
    outf.close();

//    cout << N*nus[3]*Id_scl << N*nus[4]*Id_scl;

    outf.open("nus.dat", std::ios::out);
    outf << nus[3] << nus[4];
    outf.close();
  }

  daeps_(eps_tps);

  // Note, nus are in Floquet space.
  get_A();

  // get_dnu(A_max[X_], A_max[Y_], delta_max);

  get_dnu_x_delta(A_max[X_], 0.1e-3, delta_max);

  // get_dnu2(Ax, Ay, 0.0);
}
