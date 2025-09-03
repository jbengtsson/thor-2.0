#define NO 4

#include <map>

#include "thor_lib.h"

int no_tps   = NO,
    ndpt_tps = 5;


class mml_type {
  using KeyType   = std::string;
  using ValueType = std::vector<std::string>;

  typedef std::map<KeyType, ValueType> DictType;

private:
  DictType dict;

public:
  mml_type(void);

  void set_bn_fam_par(const std::string &Fam, const int n);
  void clr_bn_fam_par(const std::string &Fam, const int n);
  void prt_dict(void);
};

void mml_type::set_bn_fam_par(const std::string &Fam, const int n)
{
  const bool prt = false;

  if (prt)
    std::cout << "\nset_bn_fam_par - " << Fam << ":\n";
  for (int k = 0; k < dict[Fam].size(); k++) {
    int Fnum = get_Fnum(dict[Fam][k].c_str());
    if (prt)
      std::cout << "  " << std::left << std::setw(10) << dict[Fam][k]
		<< std::right << " " << Fnum << "\n";
    set_bn_par(Fnum, n, 7);
  }
}

void mml_type::clr_bn_fam_par(const std::string &Fam, const int n)
{
  const bool prt = false;

  if (prt)
    std::cout << "\nclr_bn_fam_par - " << Fam << ":\n";
  for (int k = 0; k < dict[Fam].size(); k++) {
    int Fnum = get_Fnum(dict[Fam][k].c_str());
    if (prt)
      std::cout << "  " << std::left << std::setw(10) << dict[Fam][k]
		<< std::right << " " << Fnum << "\n";
    clr_bn_par(Fnum, n);
  }
}

mml_type::mml_type(void)
{
  dict["s1pr"] = {
    "s1md1r", "s1mt1r", "s1md2r", "s1mt2r", "s1md3r", "s1mt3r", "s1md4r",
    "s1mt4r", "s1md5r", "s1mt5r", "s1md6r", "s1mt6r", "s1md7r", "s1mt7r",
    "s1md8r", "s1mt8r"
  };
  dict["s2pdr"] = {
    "s2m1d1r", "s2m2d1r", "s2m1d2r", "s2m2d2r", "s2m1d3r", "s2m2d3r",
    "s2m1d4r", "s2m2d4r", "s2m1d5r", "s2m2d5r", "s2m1d6r", "s2m2d6r",
    "s2m1d7r", "s2m2d7r", "s2m1d8r",
    "s2m2d8r"
  };
  dict["s2ptr"] = {
    "s2m1t1r", "s2m2t1r", "s2m1t2r", "s2m2t2r", "s2m1t3r", "s2m2t3r",
    "s2m1t4r", "s2m2t4r", "s2m1t5r", "s2m2t5r", "s2m1t6r", "s2m2t6r",
    "s2m1t7r", "s2m2t7r", "s2m1t8r",
    "s2m2t8r"
  };
  dict["s3pdr"] = {
    "s3m1d2r", "s3m2d2r", "s3m1d3r", "s3m2d3r", "s3m1d4r", "s3m2d4r",
    "s3m1d5r", "s3m2d5r", "s3m1d6r", "s3m2d6r", "s3m1d7r", "s3m2d7r",
    "s3m1d8r", "s3m2d8r"
  };
  dict["s3ptr"] = {
    "s3m1t1r", "s3m2t1r", "s3m1t2r", "s3m2t2r", "s3m1t3r", "s3m2t3r",
    "s3m1t4r", "s3m2t4r", "s3m1t5r", "s3m2t5r", "s3m1t7r", "s3m2t7r",
    "s3m1t8r", "s3m2t8r"
  };
  dict["s4pdr"] = {
    "s4m1d2r", "s4m2d2r", "s4m1d3r", "s4m2d3r", "s4m1d4r", "s4m2d4r",
    "s4m1d5r", "s4m2d5r", "s4m1d6r", "s4m2d6r", "s4m1d7r", "s4m2d7r",
    "s4m1d8r", "s4m2d8r"};
  dict["s4ptr"] = {
    "s4m1t1r", "s4m2t1r", "s4m1t2r", "s4m2t2r", "s4m1t3r", "s4m2t3r",
    "s4m1t4r", "s4m2t4r", "s4m1t5r", "s4m2t5r", "s4m1t7r", "s4m2t7r",
    "s4m1t8r", "s4m2t8r"};
  dict["s3pd1r"] = {
    "s3m1d1r", "s3m2d1r"};
  dict["s4pd1r"] = {
    "s4m1d1r", "s4m2d1r"};
  dict["s3p1t6r"] = {
    "s3m1t6r"};
  dict["s3p2t6r"] = {
    "s3m2t6r"};
  dict["s4p1t6r"] = {
    "s4m1t6r"};
  dict["s4p2t6r"] = {
    "s4m2t6r"
  };
}

void mml_type::prt_dict(void)
{
  const int n_prt = 5;

  for(const auto &elem : dict) {
    std::cout << "\n" << elem.first << ":\n";
    auto n = elem.second.size();
    for (int k = 0; k < n; k++) {
      if (k == 0)
	std::cout << "  ";
      else if (k % n_prt == 0)
	std::cout << "\n  ";
      std::cout << std::left << std::setw(10) << elem.second[k];
    }
    if (n % n_prt != 0)
      std::cout << "\n";
  }
}


void lat_glob_prop(const bool cod)
{
  double       nu[3], ksi[2], alpha[2], beta[2];
  ss_vect<tps> nus;

  if (!cod)
    get_Map();
  else
    get_COD(10, 1e-10, 0.0, true);
  K = MapNorm(Map, g, A1, A0, Map_res, 1);
  nus = dHdJ(K);
  get_nu_ksi(nus, nu, ksi);
  get_ab(alpha, beta, 0);

  prt_lin_map(3, Map);

  printf("\nalpha = [%5.3f, %5.3f]\n", alpha[X_], alpha[Y_]);
  printf("beta  = [%5.3f, %5.3f]\n", beta[X_], beta[Y_]);
  printf("nu    = [%18.16f, %18.16f]\n", nu[X_], nu[Y_]);
  printf("ksi   = [%8.6f, %8.6f]\n", ksi[X_], ksi[Y_]);
}


void get_twiss
(const double alpha[], const double beta[], const double eta[],
 const double etap[])
{
  int      j, k;
  double   alpha1[2], beta1[2], eta1[2], etap1[2], dnu1[2], dnu2[2];
  std::ofstream outf;

  const std::string file_name = "linlat.out";

  // Crucial; to only store linear part of A.
  danot_(1);

  file_wr(outf, file_name.c_str());

  for (k = 0; k < 2; k++)
    dnu1[k] = 0e0;
  A1 = get_A(alpha, beta, eta, etap);
  for (j = 1; j <= n_elem; j++) {
    A1.propagate(j, j);
    elem_tps[j-1].A1 = get_A_CS(2, A1, dnu2);
    get_ab(A1, alpha1, beta1, dnu2, eta1, etap1);
    for (k = 0; k < 2; k++) {
      elem[j-1].Alpha[k] = alpha1[k]; elem[j-1].Beta[k] = beta1[k];
      elem[j-1].Eta[k] = eta1[k]; elem[j-1].Etap[k] = etap1[k];
    }

    // Assumes dnu < 1.0.
    for (k = 0; k < 2; k++) {
      elem[j-1].Nu[k] = floor(elem[j-2].Nu[k]) + dnu2[k];
      if ((dnu2[k] < dnu1[k]) && (elem[j-1].L >= 0e0))
	elem[j-1].Nu[k] += 1e0;
    }
    for (k = 0; k < 2; k++)
      dnu1[k] = dnu2[k];

    outf << std::setprecision(5) << std::fixed
	 << std::setw(4) << j-1 << " "
	 << std::setw(15) << std::left << elem[j-1].Name << std::right
	 << std::setw(10) << elem[j-1].S
	 << std::setw(10) << elem[j-1].Alpha[X_]
	 << std::setw(9) << elem[j-1].Beta[X_]
	 << std::setw(9) << elem[j-1].Nu[X_]
	 << std::setw(10) << elem[j-1].Eta[X_]
	 << std::setw(10) << elem[j-1].Etap[X_]
	 << std::setw(10) << elem[j-1].Alpha[Y_]
	 << std::setw(9) << elem[j-1].Beta[Y_]
	 << std::setw(9) << elem[j-1].Nu[Y_]
	 << std::setw(10) << elem[j-1].Eta[Y_]
	 << std::setw(10) << elem[j-1].Etap[Y_]
	 << std::endl;
  }

  outf.close();
}


tps get_h_local(const ss_vect<tps> &map, const bool dragt_finn)
{
  ss_vect<tps>  map1, R;

  if (dragt_finn)
    // Dragt-Finn factorization.
    return LieFact_DF(map, R);
  else {
    // Single Lie exponent.
    danot_(1);
    map1 = map;
    danot_(no_tps);
    return LieFact(map*Inv(map1));
  }
}


void prt_drv_terms(const tps &h)
{
  tps h_re, h_im;

  CtoR(h, h_re, h_im);

  printf("\nh_11001 %10.3e %10.3e\n",
	 h_ijklm_p(h_re, 1, 1, 0, 0, 1, 7),
	 h_ijklm_p(h_im, 1, 1, 0, 0, 1, 7));
  printf("h_00111 %10.3e %10.3e\n",
	 h_ijklm_p(h_re, 0, 0, 1, 1, 1, 7),
	 h_ijklm_p(h_im, 0, 0, 1, 1, 1, 7));

  printf("\nh_20001 %10.3e %10.3e\n",
	 h_ijklm_p(h_re, 2, 0, 0, 0, 1, 7),
	 h_ijklm_p(h_im, 2, 0, 0, 0, 1, 7));
  printf("h_00201 %10.3e %10.3e\n",
	 h_ijklm_p(h_re, 0, 0, 2, 0, 1, 7),
	 h_ijklm_p(h_im, 0, 0, 2, 0, 1, 7));

  printf("\nh_21000 %10.3e %10.3e\n",
	 h_ijklm_p(h_re, 2, 1, 0, 0, 0, 7),
	 h_ijklm_p(h_im, 2, 1, 0, 0, 0, 7));
  printf("h_10110 %10.3e %10.3e\n",
	 h_ijklm_p(h_re, 1, 0, 1, 1, 0, 7),
	 h_ijklm_p(h_im, 1, 0, 1, 1, 0, 7));
  printf("h_30000 %10.3e %10.3e\n",
	 h_ijklm_p(h_re, 3, 0, 0, 0, 0, 7),
	 h_ijklm_p(h_im, 3, 0, 0, 0, 0, 7));
  printf("h_10200 %10.3e %10.3e\n",
	 h_ijklm_p(h_re, 1, 0, 2, 0, 0, 7),
	 h_ijklm_p(h_im, 1, 0, 2, 0, 0, 7));
  printf("h_10020 %10.3e %10.3e\n",
	 h_ijklm_p(h_re, 1, 0, 0, 2, 0, 7),
	 h_ijklm_p(h_im, 1, 0, 0, 2, 0, 7));
}


void compute_drv_terms(const std::string &pwr_sup_name)
{
  const bool prt = false;

  tps      K, h;
  mml_type mml;

  mml.set_bn_fam_par(pwr_sup_name, 3);

  danot_(no_tps-1);
  get_Map();
  prt_lin_map(3, Map);
  danot_(no_tps);
  K = MapNorm(Map, g, A1, A0, Map_res, 1);

  h = get_h_local(Map, false);
  if (prt)
    std::cout << h;

  prt_drv_terms(h);

  mml.clr_bn_fam_par(pwr_sup_name, 3);
}


void compute_sext_resp_mat(void)
{
  const std::string pwr_sup_name = "s1pr";

  compute_drv_terms(pwr_sup_name);
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
  mml_type mml;

  danot_(no_tps-1);

  set_lat_state();

  rd_mfile(argv[1], elem);
  rd_mfile(argv[1], elem_tps);

  // initialize the symplectic integrator after energy has been defined
  ini_si();

  // disable from TPSALib- and LieLib log messages
  idprset(-1);

  lat_glob_prop(false);

  if (false) {
    const double
      alpha[] = {0.04355,  -0.02641},
      beta[]  = {17.02439,  3.79429},
      eta[]   = {-0.01098,  0.00052},
      etap[]  = {0.0,       0.0};

    get_twiss(alpha, beta, eta, etap);
  }

  if (false)
    mml.prt_dict();

  if (false) {
    mml.set_bn_fam_par("s1pr", 3);
    mml.clr_bn_fam_par("s1pr", 3);
  }

  if (!false)
    compute_sext_resp_mat();
}
