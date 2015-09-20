const double  max_ampl[] = { 100e-3, 100e-3 };

const int     n_alphac = 3;

const int     n_prm_max  = 40;
const double  scl_alpha_mp = 1e-2, scl_beta_mp = 1e-2;

extern double           b2_max, ds_max, scl_ds;
extern double           nu_[], ksi_[], rad_[], part_numb_[], tau_[];
extern tps              K, g, eps_[];
extern ss_vect<double>  fixed_point;
extern ss_vect<tps>     Map, A0, A0_inv, A1, A1_inv, Map_res, nus_;

FILE* file_read(const char file_name[]);

FILE* file_write(const char file_name[]);

void file_rd(ifstream &inf, const char file_name[]);

void file_wr(ofstream &outf, const char file_name[]);

void set_to_cout(ofstream &fp_out);

void prt_lat(const char *file_name);

void get_matrix(const ss_vect<tps> &Map, float **M);

void copy_mat(const int n, float **A, float **B);

double Det(const int n, float **A);

void prt_mat(const int n, const double **M);

int sign(const double x);

void prt_vec(const int n, const double *x);

void vec_cp(const int n, double *x, double *y);

void file_rd(ifstream &inf, const char file_name[]);

void file_wr(ofstream &outf, const char file_name[]);

void set_to_cout(ofstream &fp_out);

void get_matrix(const ss_vect<tps> &Map, float **M);

void copy_mat(const int n, float **A, float **B);

double Det(const int n, float **A);

void prt_mat(const int n, const double **M);

int sign(const double x);

void prt_vec(const int n, const double *x);

void vec_cp(const int n, double *x, double *y);

double scl_prod(const int n, double *x, double *y);

double vec_abs(const int n, double *x);

void vec_add(const int n, double *x, double *y, double *z);

void vec_sub(const int n, double *x, double *y, double *z);

void lin_trans(const int m, const int n, double **A, double *x, double *y);

void mat_cp(const int m, const int n, double **A, double **B);

void mat_mul(const int m, const int n, double **A, double **B, double **C);

void mat_tp(const int m, const int n, double **A, double **B);

void prt_mat(const int m, const int n, double **A);

void SVD_lim(const int m, const int n, double **A, double *b,
	     const double *corr_max, const double s_cut, double *corr0,
	     double *dcorr);

void SVD(const int m, const int n, double **A, double *b, double *dcorr);

void get_A_inv(const int m, const int n, float **U, float *w, float **V,
	       float **A_inv);

double GetAngle(const double x, const double y);

tps do_PBs(const tps &a, const tps &b, const int n, const int jj[]);

tps BCH(const tps &a, const tps &b, const int n_ord);

int get_Fnum(const char name[]);

char* get_Name(const int Fnum);

int get_n_Kids(const int Fnum);

long int get_loc(const int Fnum, const int Knum);

double get_L(const int Fnum, const int Knum);

double get_bn(const int Fnum, const int Knum, const int n);

double get_bnL(const int Fnum, const int Knum, const int n);

void set_bn(const int Fnum, const int Knum, const int n, const double bn);

void set_bn(const int Fnum, const int n, const double bn);

void set_bnL(const int Fnum, const int Knum, const int n, const double bnL);

void set_bnL(const int Fnum, const int n, const double bnL);

void set_dbn(const int Fnum, const int Knum, const int n, const double dbn);

void set_dbn(const int Fnum, const int n, const double dbn);

void set_dbnL(const int Fnum, const int Knum, const int n, const double dbnL);

void set_dbnL(const int Fnum, const int n, const double dbnL);

void set_L(const int Fnum, const int Knum, const double L);

void set_L(const int Fnum, const double L);

void set_dL(const int Fnum, const int Knum, const double dL);

void set_dL(const int Fnum, const double dL);

void set_bn_par(const int Fnum, const int Knum, const int n, const int j);

void set_bn_par(const int Fnum, const int n, const int j);

void clr_bn_par(const int Fnum, const int Knum, const int n);

void clr_bn_par(const int Fnum, const int n);

void set_s_par(const int Fnum, const int Knum, const int j);

void set_s_par(const int Fnum, const int j);

void clr_s_par(const int Fnum, const int Knum);

void clr_s_par(const int Fnum);

double get_BoBrho(const int Fnum, const int Knum);

void set_BoBrho(const int Fnum, const int Knum, const double BoBrho);

void get_cav(const int Fnum, const int Knum,
	     int &h_rf, double &V_rf, double &f_rf);

void get_Map(void);

void get_Map(const ss_vect<double> &fixed_point);

ss_vect<tps> get_Map(const int k);

bool get_COD(const int i_max, const double eps, const double delta,
	     const bool prt);

void get_A1(const double alpha_x, const double beta_x,
	    const double alpha_y, const double beta_y);

tps get_h(void);

tps get_H(void);

// transform from phase space to Floquet space
inline ss_vect<double> get_FS(const ss_vect<double> ps)
{
  return (A1_inv*ps).cst();
}

ss_vect<double> get_J_phi(const ss_vect<double> ps);

ss_vect<tps> get_A_inv();

ss_vect<tps> get_S(const int n_DOF);

ss_vect<tps> tp_S(const int n_DOF, const ss_vect<tps> &A);

bool is_h_ijklm(const int i1, const int j1, const int k1,
		const int l1, const int m1,
		const int i2, const int j2, const int k2,
		const int l2, const int m2);

double h_ijklm(const tps &h, const int i, const int j, const int k,
	       const int l, const int m);

double h_ijklm_scl(const tps &h, const int i, const int j, const int k,
		   const int l, const int m,
		   const double nu_x, const double nu_y);

double h_ijklm_p(const tps &h, const int i, const int j, const int k,
		 const int l, const int m, const int p);

double h_ijklm_p_scl(const tps &h, const int i, const int j, const int k,
		     const int l, const int m, const int p,
		     const double nu_x, const double nu_y);

void get_nu_ksi(const ss_vect<tps> &nus, double nu[], double ksi[]);

void prt_nu(const ss_vect<tps> &nus);

void get_k_J(const tps &K,
	     double k_J2[], double k_J3[], double k_J4[], double k_delta2[]);

void prt_dnu(const tps &K);

void get_alphac(double alphac[]);

void prt_alphac();

double H_long(const double phi, const double delta,
	      const int h_rf, const double V_rf, const double phi0,
	      const double alphac[]);

void prt_H_long(const int n, const double phi_max, const double delta_max,
		const double U0);

double arccos_(double x);

void get_nu_symp(const ss_vect<tps> &Map, double nu[]);

inline bool Check_Ampl(const ss_vect<double> &x);

inline bool Check_Ampl(const ss_vect<tps> &x);

bool track(const double x, const double px, const double y, const double py,
	   const double delta, const long int n, const double f_rf,
	   const bool prt);

bool map_track(const double x, const double px,
	       const double y, const double py,
	       const double delta, const long int n,
	       const double f_rf, const bool prt);

void get_r_stable(double &r, const double phi, const double delta,
		  const long int n, const double eps, bool map);

ss_vect<tps> get_eta(const long int k);

void get_ab(double alpha[], double beta[], const long int k);

ss_vect<tps> get_A1_CS(const ss_vect<tps> &A1, tps dnu[]);

void get_ab(tps ab[], tps nu[], const long int loc);

void get_ab(tps ab[], tps dnu[], const long int k1, const long int k2);

void get_nu(double nu[], const long int k);

void bend_cal(const int Fnum1);

double get_bn_s(const int Fnum, const int Knum, const int n);

double get_bnL_s(const int Fnum, const int Knum, const int n);

void set_bn_s(const int Fnum, const int Knum, const int n, const double bn);

void set_bn_s(const int Fnum, const int n, const double bn);

void set_bnL_s(const int Fnum, const int Knum, const int n, const double bnL);

void set_bnL_s(const int Fnum, const int n, const double bnL);

void set_dbn_s(const int Fnum, const int Knum, const int n, const double dbn);

void set_dbn_s(const int Fnum, const int n, const double dbn);

void fit_alpha(const double alpha0_x, const double beta0_x,
	       const double alpha0_y, const double beta0_y,
	       const long int k1, const long int k2,
	       const int n_b2, const int b2s[],
	       const double eps, const bool prt);

void fit_beta(const double alpha0_x, const double beta0_x,
	      const double alpha0_y, const double beta0_y,
	      const double beta1_x, const double beta1_y,
	      const long int k1, const long int k2,
	      const int n_b2, const int b2s[],
	      const double eps, const bool prt);

void fit_alpha_beta(const double alpha0_x, const double beta0_x,
		    const double alpha0_y, const double beta0_y,
		    const double beta1_x, const double beta1_y,
		    const long int k1, const long int k2,
		    const int n_b2, const int b2s[],
		    const double eps, const bool prt);

void fit_tune(const double nu_x, const double nu_y,
	      const int n_b2, const int b2s[],
	      const double eps, const bool prt);

void fit_tune(const double nu_x, const double nu_y,
	      const double beta1_x, const double beta1_y, const int k1,
	      const double beta2_x, const double beta2_y, const int k2,
	      const double beta3_x, const double beta3_y, const int k3,
	      const int n_b2, const int b2s[],
	      const double eps, const bool prt);

void fit_tune_SLS(const double nu_x, const double nu_y,
		  const double beta1_x, const double beta1_y, const int k1,
		  const double beta2_x, const double beta2_y, const int k2,
		  const double beta3_x, const double beta3_y, const int k3,
		  const double beta4_x, const double beta4_y, const int k4,
		  const int n_b2, const int b2s[],
		  const double eps, const bool prt);

double get_diff(const double alpha1[], const double beta1[],
		const double dnu12[], const long int locs[],
		tps ab2[], tps nu[], const bool prt);

void fit_femto_SLS(const double dnu_x, const double dnu_y,
		   const double alpha1_x, const double alpha1_y,
		   const double beta1_x, const double beta1_y, const int k1,
		   const int k2,
		   const int n_b2, const int b2s[],
		   const double eps, const bool prt);

void fit_dnu_SLS(const double nu_x,
		 const double beta1_x, const double beta1_y, const int k1,
		 const int k2,
		 const int n_b2, const int b2s[],
		 const double eps, const bool prt);

void fit_chrom(const double ksi_x, const double ksi_y,
	       const int n_b3, const int b3s[], const bool prt);

void fit_chrom1(const double ksi_x, const double ksi_y,
	        const int n_b3, const int b3s[], const double eps,
	        const bool prt);

void fit_disp(const double eta1_x, const int k1,
	      const double eta2_x, const int k2,
	      const int n_b2, const int b2s[],
	      const double eps, const bool prt);

void fit_DBA(const double nu_x, const double nu_y, const long int k,
	     const int n_b2, const int b2s[], const double eps,
	     const bool prt);

double get_dynap(const double r, const double delta, const int n,
		 const double eps, const int n_pts,
		 double x_min[], double x_max[], bool map);

void get_rad(char *cav_name, double I_beam);

void get_emittance(void);

