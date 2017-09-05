/* Author:       Johan Bengtsson

   Definitions:  Tools for lattice studies.  */


double b2_max = 3.0, ds_max = 0.15,
       scl_ds = 0.1; // Obsolete.

double          nu_[3], ksi_[3], rad_[3], part_numb_[3], tau_[3];
tps             K, g, eps_[3];
ss_vect<double> fixed_point;
ss_vect<tps>    Map, A0, A0_inv, A1, A1_inv, Map_res, nus_;


FILE* file_read(const char file_name[])
{
  FILE *fp;
  
  fp = fopen(file_name, "r");
  if (fp == NULL) {
    printf("File not found: %s\n", file_name);
    exit(-1);
  } else
    return(fp);
}


FILE* file_write(const char file_name[])
{
  FILE *fp;

  fp = fopen(file_name, "w");
  if (fp == NULL) {
    printf("Could not create file: %s\n", file_name);
    exit(-1);
  } else
    return(fp);
}


void file_rd(std::ifstream &inf, const char file_name[])
{

  inf.open(file_name, std::ios::in);
  if (!inf.is_open()) {
    printf("File not found: %s\n", file_name);
    exit(-1);
  }
}


void file_wr(std::ofstream &outf, const char file_name[])
{

  outf.open(file_name, std::ios::out);
  if (!outf.is_open()) {
    printf("Could not create file: %s\n", file_name);
    exit(-1);
  }
}


void set_to_cout(std::ofstream &fp_out)
{
  // redirect to cout
  fp_out.copyfmt(std::cout); fp_out.clear(std::cout.rdstate());
  fp_out.std::basic_ios<char>::rdbuf(std::cout.rdbuf());
}


void prt_lin_map(const int n_DOF, const ss_vect<tps> &map)
{
  int i, j;

  std::cout << std::endl;
  for (i = 0; i < 2*n_DOF; i++) {
    for (j = 0; j < 2*n_DOF; j++) {
      if (true)
	std::cout << std::scientific << std::setprecision(6)
		  << std::setw(14) << map[i][j];
      else
	std::cout << std::scientific << std::setprecision(16)
		  << std::setw(24) << map[i][j];
    }
    std::cout << std::endl;
  }
}


double get_code(elem_type<double> &elem)
{
  double code;

  switch (elem.kind) {
  case Drift:
    code = 0e0;
    break;
  case Mpole:
    if (elem.mpole->h_bend != 0e0)
      code = 0.5e0;
    else if (elem.mpole->bn[Quad-1] != 0)
      code = sgn(elem.mpole->bn[Quad-1]);
    else if (elem.mpole->bn[Sext-1] != 0)
      code = 1.5*sgn(elem.mpole->bn[Sext-1]);
    else
      code = 0e0;
    break;
  default:
    code = 0e0;
    break;
  }
  return code;
}

void prt_lat(const int i0, const int i1, const char *file_name)
{
  long int i;
  FILE     *outf;

  outf = file_write(file_name);
  fprintf(outf, "#        name           s   code"
	        "  alphax  betax   nux   etax   etapx");
  fprintf(outf, "  alphay  betay   nuy   etay   etapy\n");
  fprintf(outf, "#                      [m]"
	        "                 [m]           [m]");
  fprintf(outf, "                   [m]           [m]\n");
  fprintf(outf, "#\n");

  for (i = i0; i <= i1; i++) {
    fprintf(outf, "%4ld %15s %9.5f %4.1f"
	    " %9.5f %8.5f %8.5f %8.5f %8.5f"
	    " %9.5f %8.5f %8.5f %8.5f %8.5f\n",
	    i, elem[i].Name, elem[i].S, get_code(elem[i]),
	    elem[i].Alpha[X_], elem[i].Beta[X_], elem[i].Nu[X_],
	    elem[i].Eta[X_], elem[i].Etap[X_],
	    elem[i].Alpha[Y_], elem[i].Beta[Y_], elem[i].Nu[Y_],
	    elem[i].Eta[Y_], elem[i].Etap[Y_]);
  }

  fclose(outf);
}

void prt_lat(const char *file_name)
{
  prt_lat(0, n_elem, file_name);
}

void get_dnu(const int n, const ss_vect<tps> &A, double dnu[])
{
  int k;

  for (k = 0; k < n; k++) {
    dnu[k] = atan2(A[2*k][2*k+1], A[2*k][2*k])/(2.0*M_PI);
    if (dnu[k] < 0.0) dnu[k] += 1.0;
  }
}


void get_ab(const ss_vect<tps> &A,
	    double alpha[], double beta[], double dnu[],
	    double eta[], double etap[])
{
  int          k;
  ss_vect<tps> A_Atp;

  A_Atp = A*tp_S(2, A);

  for (k = 0; k <= 1; k++) {
    eta[k] = A[2*k][delta_]; etap[k] = A[2*k+1][delta_];

    alpha[k] = -A_Atp[2*k][2*k+1]; beta[k] = A_Atp[2*k][2*k];
  }

  get_dnu(2, A, dnu);
}


ss_vect<tps> get_A(const double alpha[], const double beta[],
		   const double eta[], const double etap[])
{
  int          k;
  ss_vect<tps> A, Id;

  Id.identity();

  A.identity();
  for (k = 0; k < 2; k++) {
    A[2*k]  = sqrt(beta[k])*Id[2*k];
    A[2*k+1] = -alpha[k]/sqrt(beta[k])*Id[2*k] + 1.0/sqrt(beta[k])*Id[2*k+1];

    A[2*k] += eta[k]*Id[delta_]; A[2*k+1] += etap[k]*Id[delta_];
  }

  return A;
}


ss_vect<tps> get_A_CS(const int n, const ss_vect<tps> &A, double dnu[])
{
  int          k;
  double       c, s;
  ss_vect<tps> Id, R;

  Id.identity(); R.identity(); get_dnu(n, A, dnu);
  for (k = 0; k < n; k++) {
    c = cos(2.0*M_PI*dnu[k]); s = sin(2.0*M_PI*dnu[k]);
    R[2*k] = c*Id[2*k] - s*Id[2*k+1]; R[2*k+1] = s*Id[2*k] + c*Id[2*k+1];
  }

  return A*R;
}


void get_twoJ(const int n_DOF, const ss_vect<double> &ps,
	      const ss_vect<tps> &A, double twoJ[])
{
  int             j, k, jj[ss_dim];
  ss_vect<double> z;
  ss_vect<tps>    A1;

  for (j = 0; j < ss_dim; j++)
    jj[j] = (j < 2*n_DOF)? 1 : 0;

  // Get linear part.
  A1.zero();
  for (j = 0; j < ps_dim; j++)
    for (k = 0; k < ps_dim; k++)
      A1[j] += A[j][k]*tps(0e0, k+1);

  z = (PInv(A1, jj)*ps).cst();

  for (j = 0; j < n_DOF; j++)
    twoJ[j] = sqr(z[2*j]) + sqr(z[2*j+1]);
}


void get_twoJ(const int n_DOF, const ss_vect<tps> &ps,
	      const ss_vect<tps> &A, tps twoJ[])
{
  int          j, k, jj[ss_dim];
  ss_vect<tps> A1, z;

  for (j = 0; j < ss_dim; j++)
    jj[j] = ((j < 2*n_DOF) || (j == ss_dim-1))? 1 : 0;

  // Get linear part.
  A1.zero();
  for (j = 0; j < ps_dim; j++)
    for (k = 0; k < ps_dim; k++)
      A1[j] += A[j][k]*tps(0e0, k+1);
  // Include parameter dependance.
  for (j = 0; j < 2; j++) {
    A1[j]   += h_ijklm_p(A[j],   1, 0, 0, 0, 0, 7)*tps(0e0, 1)*tps(0e0, 7);
    A1[j]   += h_ijklm_p(A[j],   0, 1, 0, 0, 0, 7)*tps(0e0, 2)*tps(0e0, 7);
    A1[j+2] += h_ijklm_p(A[j+2], 0, 0, 1, 0, 0, 7)*tps(0e0, 3)*tps(0e0, 7);
    A1[j+2] += h_ijklm_p(A[j+2], 0, 0, 0, 1, 0, 7)*tps(0e0, 4)*tps(0e0, 7);
  }

  z = PInv(A1, jj)*ps;
 
  for (j = 0; j < n_DOF; j++)
    twoJ[j] = sqr(z[2*j]) + sqr(z[2*j+1]);
}


void prt_lat(const int i0, const int i1, const char *fname, const int n)
{
  long int           i;
  int                j, k;
  double             s, h;
  double             alpha[2], beta[2], nu[2], dnu[2], eta[2], etap[2], dnu1[2];
  mpole_type<double> *mpole;
  ss_vect<double>    eta_Fl;
  ss_vect<tps>       A_CS;
  FILE               *outf;

  const double  c1 = 1e0/(2e0*(2e0-pow(2e0, 1e0/3e0))), c2 = 0.5e0-c1;
  const double  d1 = 2e0*c1, d2 = 1e0-2e0*d1;

  outf = file_write(fname);
  fprintf(outf, "#        name           s   code"
	        "  alphax  betax   nux   etax   etapx");
  fprintf(outf, "  alphay  betay   nuy   etay   etapy\n");
  fprintf(outf, "#                      [m]"
	        "                 [m]           [m]");
  fprintf(outf, "                   [m]           [m]\n");
  fprintf(outf, "#\n");

  for (i = i0; i <= i1; i++) {
    if ((i != 0) &&
	((elem[i].kind == Drift) ||
	 ((elem[i].kind == Mpole) && (elem[i].L != 0e0)))) {
      mpole = elem[i].mpole;

      for (k = 0; k < 2; k++) {
	alpha[k] = elem[i-1].Alpha[k]; beta[k] = elem[i-1].Beta[k];
	nu[k] = elem[i-1].Nu[k];
	eta[k] = elem[i-1].Eta[k]; etap[k] = elem[i-1].Etap[k];
      }

      A1 = get_A(alpha, beta, eta, etap);
      s = elem[i].S - elem[i].L; h = elem[i].L/n;
      for (j = 1; j <= n; j++) {
	s += h;

	if (elem[i].kind == Drift)
	  drift_pass(h, A1);
	else if (elem[i].kind == Mpole) {
	  if ((j == 1) && (mpole->h_bend != 0e0))
	    bend_HE_fringe(mpole->h_bend, mpole->edge1, mpole->gap, A1);

	  drift_pass(c1*h, A1);
	  thin_kick(Quad, mpole->an, mpole->bn, d1*h, mpole->h_bend,
		    mpole->h_bend, true, A1);
	  drift_pass(c2*h, A1);
	  thin_kick(Quad, mpole->an, mpole->bn, d2*h, mpole->h_bend,
		    mpole->h_bend, true, A1);
	  drift_pass(c2*h, A1);
	  thin_kick(Quad, mpole->an, mpole->bn, d1*h, mpole->h_bend,
		    mpole->h_bend, true, A1);
	  drift_pass(c1*h, A1);

	  if ((j == n) && (mpole->h_bend != 0e0))
	    bend_HE_fringe(mpole->h_bend, mpole->edge2, mpole->gap, A1);
	}

	get_ab(A1, alpha, beta, dnu, eta, etap);

	if(elem[i].L < 0e0)
	  for (k = 0; k < 2; k++)
	    dnu[k] -= 1e0;

	A_CS = get_A_CS(2, A1, dnu1);

	eta_Fl.zero();
	for (k = 0; k < 2; k++) {
	  eta_Fl[2*k] = eta[k]; eta_Fl[2*k+1] = etap[k];
	}

	fprintf(outf, "%4ld %15s %6.2f %4.1f"
		" %7.3f %6.3f %6.3f %6.3f %6.3f"
		" %7.3f %6.3f %6.3f %6.3f %6.3f %10.3e %10.3e\n",
		i, elem[i].Name, s, get_code(elem[i]),
		alpha[X_], beta[X_], nu[X_]+dnu[X_], eta[X_], etap[X_],
		alpha[Y_], beta[Y_], nu[Y_]+dnu[Y_], eta[Y_], etap[Y_],
		eta_Fl[x_], eta_Fl[px_]);
      }
    } else {
      A1 = get_A(elem[i].Alpha, elem[i].Beta, elem[i].Eta, elem[i].Etap);

      eta_Fl.zero();
      for (k = 0; k < 2; k++) {
	eta_Fl[2*k] = elem[i].Eta[k]; eta_Fl[2*k+1] = elem[i].Etap[k];
      }
      eta_Fl = (Inv(A1)*eta_Fl).cst();

      fprintf(outf, "%4ld %15s %6.2f %4.1f"
	      " %7.3f %6.3f %6.3f %6.3f %6.3f"
	      " %7.3f %6.3f %6.3f %6.3f %6.3f %10.3e %10.3e\n",
	      i, elem[i].Name, elem[i].S, get_code(elem[i]),
	      elem[i].Alpha[X_], elem[i].Beta[X_], elem[i].Nu[X_],
	      elem[i].Eta[X_], elem[i].Etap[X_],
	      elem[i].Alpha[Y_], elem[i].Beta[Y_], elem[i].Nu[Y_],
	      elem[i].Eta[Y_], elem[i].Etap[Y_],
	      eta_Fl[x_], eta_Fl[px_]);
    }
  }

  fclose(outf);
}


void prt_lat(const char *fname, const int n)
{
  prt_lat(0, n_elem, fname, n);
}


void get_matrix(const ss_vect<tps> &Map, float **M)
{
  int i, j, jj[ss_dim];

  for (i = 0; i < ss_dim; i++)
    jj[i] = 0;

  for (i = 0; i < ps_dim; i++)
    for (j = 0; j < ps_dim; j++) {
      jj[j] = 1; M[i+1][j+1] = Map[i][jj]; jj[j] = 0;
    }
}


void copy_mat(const int n, float **A, float **B)
{
  int i, j;

  for (i = 1; i <= n; i++)
    for (j = 1; j <= n; j++)
      B[i][j] = A[i][j];
}


double Det(const int n, float **A)
{
  int   j, *indx;
  float d, **A1;

  indx = ivector(1, n); A1 = matrix(1, n, 1, n);

  copy_mat(n, A, A1); ludcmp(A1, n, indx, &d);
  for (j = 1; j <= n; j++)
    d *= A1[j][j];

  free_ivector(indx, 1, n); free_matrix(A1, 1, n, 1, n);

  return d;
}


void prt_mat(const int n, const double **M)
{
  int i, j;

  std::cout << std::endl;
  std::cout << std::scientific << std::setprecision(5);
  for (i = 1; i <= n; i++) {
    for (j = 1; j <= n; j++)
      std::cout << " " << std::setw(12) << M[i][j];
    std::cout << std::endl;
  }
}


int sign(const double x)
{
  if (x >= 0.0)
    return 1;
  else
    return -1;
}


void prt_vec(const int n, const double *x)
{
  int j;

  const int n_prt = 9;

  for (j = 1; j <= n; j++) {
    std::cout << std::scientific << std::setprecision(3)
	      << std::setw(11) << x[j];
    if (j % n_prt == 0) std::cout << std::endl;
  }
  if (n % n_prt != 0) std::cout << std::endl;
}


void vec_cp(const int n, double *x, double *y)
{
  int j;

  for (j = 1; j <= n; j++)
    y[j] = x[j];
}


double scl_prod(const int n, double *x, double *y)
{
  int    k;
  double xy;

  xy = 0.0;
  for (k = 1; k <= n; k++)
    xy += x[k]*y[k];
  
  return xy;
}


double vec_abs(const int n, double *x) { return sqrt(scl_prod(n, x, x)); }


void vec_add(const int n, double *x, double *y, double *z)
{
  int i;

  for (i = 1; i <= n; i++)
    z[i] = x[i] + y[i];
}


void vec_sub(const int n, double *x, double *y, double *z)
{
  int i;

  for (i = 1; i <= n; i++)
    z[i] = x[i] - y[i];
}


void lin_trans(const int m, const int n, double **A, double *x, double *y)
{
  int j, k;

  for (j = 1; j <= m; j++) {
    y[j] = 0.0;
    for (k = 1; k <= n; k++)
      y[j] += A[j][k]*x[k];
  }
}


void mat_cp(const int m, const int n, double **A, double **B)
{
  int i, j;

  for (i = 1; i <= m; i++)
    for (j = 1; j <= n; j++)
      B[i][j] = A[i][j];
}


void mat_mul(const int m, const int n, double **A, double **B, double **C)
{
  int i, j, k;

  for (i = 1; i <= m; i++)
    for (j = 1; j <= n; j++) {
      C[i][j] = 0.0;
      for (k = 1; k <= n; k++)
	C[i][j] += A[i][k]*B[k][j];
    }
}


void mat_tp(const int m, const int n, double **A, double **B)
{
  int i, j;

  for (i = 1; i <= m; i++)
    for (j = 1; j <= n; j++)
      B[j][i] = A[i][j];
}


void prt_mat(const int m, const int n, double **A)
{
  int i, j;

  for (i = 1; i <= m; i++) {
    for (j = 1; j <= n; j++)
      std::cout << std::scientific << std::setprecision(3)
	   << std::setw(12) << A[i][j];
    std::cout << std::endl;
  }
}


void SVD_lim(const int m, const int n, double **A, double *b,
	     const double corr_max[], const double s_cut, double *corr0,
	     double *dcorr)
{
  int    i, j, k, n_sing, max_ind = 0, sgn, max_sgn = 0, n1, mpn;
  double **U, **U_tp, **V, **AtU_inv, **n0_tp, n_norm;
  double **U_m, **V_m, **U_m_tp, **N_m;
  double *w, *w_m, *dUb, *Ub0, *b_m, *p_m, *b_ext;
  double **n_m_tp, Delta_m, **d0, max_dist, u, **C, **C_tp, s_max, v;
  double **A_ext;

  const bool   prt = false;
  const double eps = 1e-10;

  mpn = m + n;

  if (prt) {
    std::cout << std::endl;
    std::cout << "m = " << m << ", n = " << n << ", mpn = " << mpn << std::endl;
  }

  dUb = dvector(1, n); Ub0 = dvector(1, n);
  b_m = dvector(1, mpn);
  U = dmatrix(1, mpn, 1, n); w = dvector(1, n); V = dmatrix(1, n, 1, n);
  U_tp = dmatrix(1, n, 1, mpn);
  A_ext = dmatrix(1, mpn, 1, n); b_ext = dvector(1, mpn);
  AtU_inv = dmatrix(1, n, 1, n);
  n0_tp = dmatrix(1, n, 1, n); d0 = dmatrix(1, n, 1, 2);

  N_m = dmatrix(1, n, 1, n); U_m = dmatrix(1, n, 1, n);
  w_m = dvector(1, n); V_m = dmatrix(1, n, 1, n);
  U_m_tp = dmatrix(1, n, 1, n); n_m_tp = dmatrix(1, n, 1, n);
  p_m = dvector(1, n); 

  C = dmatrix(1, n, 1, n); C_tp = dmatrix(1, n, 1, n);

  n1 = m;

  if (prt) {
    std::cout << std::endl;
    std::cout << "A:" << std::endl;
    prt_mat(m, n, A);
  }

  for (i = 1; i <= m; i++) {
    b_ext[i] = b[i];
    for (j = 1; j <= n; j++) {
      U[i][j] = A[i][j]; A_ext[i][j] = A[i][j];
    }
  }

  // SVD of A
  dsvdcmp(U, m, n, w, V);
  
  mat_tp(m, n, U, U_tp);

  if (prt) {
    std::cout << std::endl;
    std::cout << "V:" << std::endl;
    prt_mat(n, n, V);
  }

  // transform b into U-basis
  lin_trans(n, m, U_tp, b, dUb);

  if (prt) {
    std::cout << std::endl;
    std::cout << "dUb:" << std::endl;
    prt_vec(n, dUb);
  }

  s_max = -1e30;
  for (i = 1; i <= n; i++)
    s_max = max(w[i], s_max);
  
  std::cout << std::endl << "singular values:" << std::endl;
  n_sing = 0;
  for (i = 1; i <= n; i++) {
    std::cout << std::scientific << std::setprecision(3)
	      << std::setw(10) << w[i];
//    if (w[i]/s_max < s_cut) {
    if (w[i] < s_cut) {
      w[i] = 0.0;
      std::cout << " (zeroed)";
      // if A is singular, extend with null space
      n_sing++; n1++;
      for (j = 1; j <= n; j++)
	A_ext[n1][j] = V[j][i];
      b_ext[n1] = 0.0;
   }
    if (i % 8 == 0) std::cout << std::endl;
  }
  if (n % 8 != 0) std::cout << std::endl;

  if (prt) {
    std::cout << std::endl;
    std::cout << "n_sing = " << n_sing << std::endl;
    std::cout << std::endl;
    std::cout << "A_ext:" << std::endl;
    prt_mat(n1, n, A_ext);
  }

  if (n_sing > 0) {
    for (i = 1; i <= n1; i++)
      for (j = 1; j <= n; j++) {
	U[i][j] = A_ext[i][j];
      }

    // SVD of A
    dsvdcmp(U, n1, n, w, V);
  
    mat_tp(n1, n, U, U_tp);

    if (prt) {
      std::cout << std::endl;
      std::cout << "V:" << std::endl;
      prt_mat(n, n, V);
    }

    // transform b into U-basis
    lin_trans(n, n1, U_tp, b_ext, dUb);

    if (prt) {
      std::cout << std::endl;
      std::cout << "dUb:" << std::endl;
      prt_vec(n, dUb);
    }

    s_max = -1e30;
    for (i = 1; i <= n; i++)
      s_max = max(w[i], s_max);
  
    std::cout << std::endl;
    std::cout << "singular values:" << std::endl;
    for (i = 1; i <= n; i++) {
      std::cout << std::scientific << std::setprecision(3)
		<< std::setw(10) << w[i];
//      if (fabs(w[i])/s_max < s_cut) {
      if (fabs(w[i]) < s_cut) {
	w[i] = 0.0;
	std::cout << " (zeroed, SVD_lim)" << std::endl;
	exit(0);
      }
      if (i % 8 == 0) std::cout << std::endl;
    }
    if (n % 8 != 0) std::cout << std::endl;
    
    n_sing = 0;
  }

  // compute corrector configurations
  for (i = 1; i <= n; i++)
    for (j = 1; j <= n; j++)
      C_tp[j][i] = w[i]*V[j][i];

  if (prt) {
    std::cout << std::endl;
    std::cout << "C^T:" << std::endl;
    prt_mat(n, n, C_tp);
  }

  mat_tp(n, n, C_tp, C); lin_trans(n, n, C, corr0, Ub0);

  // compute the normal vectors for the limiting planes j
  for (i = 1; i <= n; i++)
    for (j = 1; j <= n; j++)
      if (w[i] != 0.0)
	AtU_inv[i][j] = V[j][i]/w[i];
      else
	AtU_inv[i][j] = 0.0;
  
  if (prt) {
    std::cout << std::endl;
    std::cout << "AtU^-1:" << std::endl;
    prt_mat(n, n, AtU_inv);
  }

  mat_tp(n, n, AtU_inv, n0_tp);
  // normalize
  for (j = 1; j <= n; j++) {
    n_norm = vec_abs(n, n0_tp[j]);
    for (k = 1; k <= n; k++)
      n0_tp[j][k] /= n_norm;
  }

  if (prt) {
    std::cout << std::endl;
    std::cout << "n0^T:" << std::endl;
    prt_mat(n, n, n0_tp);
  }

  if (prt) {
    std::cout << std::endl;
    std::cout << "d0:" << std::endl;
    std::cout << std::endl;
  }
  for (j = 1; j <= n; j++) {
    u = corr_max[j]*scl_prod(n, C_tp[j], n0_tp[j]);
    d0[j][1] = -u; d0[j][2] = u;
    if (prt)
      std::cout << std::scientific << std::setprecision(3)
	   << std::setw(10) << d0[j][1] << " - "
		<< std::setw(10) << d0[j][2] << std::endl;
  }

  max_dist = 1.0;
//  while ((n_sing < m) && (max_dist > eps)) {
  while ((n_sing < n) && (max_dist > eps)) {
    if (n_sing > 0) {
      if (prt) {
	std::cout << std::endl;
	std::cout << "N_m:" << std::endl;
	prt_mat(n, n_sing, N_m);
      }

      // remove null space from orthogonal basis
      for (i = 1; i <= n; i++)
	for (j = 1; j <= n_sing; j++)
	  U_m[i][j] = N_m[i][j];
      
      dsvdcmp(U_m, n, n_sing, w_m, V_m); mat_tp(n, n_sing, U_m, U_m_tp);

      if (prt) {
	std::cout << std::endl;
	std::cout << "singular values:" << std::endl;
	if (i % 8 == 0) std::cout << std::endl;
	for (i = 1; i <= n; i++) {
	  std::cout << std::scientific << std::setprecision(3)
		    << std::setw(10) << w[i];
	  if (n % 8 != 0) std::cout << std::endl;
	}
	std::cout << std::endl;

	std::cout << std::endl;
	std::cout << "U_m:" << std::endl;
	prt_mat(n, n_sing, U_m);

	std::cout << std::endl;
	std::cout << "V_m:" << std::endl;
	prt_mat(n_sing, n_sing, V_m);
      }

      for (j = 1; j <= n; j++) {
	lin_trans(n_sing, n, U_m_tp, n0_tp[j], p_m);
	lin_trans(n, n_sing, U_m, p_m, n_m_tp[j]);
	vec_sub(n, n0_tp[j], n_m_tp[j], n_m_tp[j]);
	// normalize
	n_norm = vec_abs(n, n_m_tp[j]);
	for (k = 1; k <= n; k++)
	  n_m_tp[j][k] /= n_norm;
      }
    } else
      mat_cp(n, n, n0_tp, n_m_tp);

    if (prt) {
      std::cout << std::endl;
      std::cout << "n_m^T:" << std::endl;
      prt_mat(n, n, n_m_tp);
    }

    if (prt) {
      std::cout << std::endl;
      std::cout << "Delta" << std::endl;
      std::cout << std::endl;
    }

    max_dist = 0.0;
    for (j = 1; j <= n; j++) {
      u = scl_prod(n, dUb, n0_tp[j]) + scl_prod(n, Ub0, n0_tp[j]);
      v = scl_prod(n, n_m_tp[j], n0_tp[j]);
      if (u >= 0.0) {
	sgn = 1; Delta_m = (u-d0[j][2])/v;
      } else {
	sgn = -1; Delta_m = -(u-d0[j][1])/v;
      }
      if (Delta_m > max_dist) {
	max_ind = j; max_dist = Delta_m; max_sgn = sgn;
      }
      if (prt) {
	std::cout << std::scientific << std::setprecision(3)
	     << "n0 = " << std::setw(11) << d0[j][2]
	     << ", sign = " << std::setw(2) << sgn
	     << ", b.n0 = " << std::setw(11) << u
	     << ", Delta_m = " << std::setw(11) << Delta_m
	     << ", n0.n_m = " << std::setw(11) << v
	     << ", max_dist = " << max_dist << std::endl;
      }
    }

    if (max_dist > eps) {
      n_sing++;
      for (k = 1; k <= n; k++) {
	dUb[k] -= max_sgn*max_dist*n_m_tp[max_ind][k];
	N_m[k][n_sing] = n0_tp[max_ind][k];
      }
    }

    if (prt) {
      std::cout << std::endl;
      std::cout << "dUb:" << std::endl;
      prt_vec(n, dUb);
    }
  }

  if (prt) {
    std::cout << std::endl;
    std::cout << "n_sing = " << n_sing << std::endl;
  }

  if (prt) {
    std::cout << std::endl;
    std::cout << "b:" << std::endl;
    prt_vec(m, b);
  }

  // transform dUb back to original basis
  lin_trans(n1, n, U, dUb, b_m);

  if (prt) {
    std::cout << std::endl;
    std::cout << "b_m:" << std::endl;
    prt_vec(n1, b_m);
  }

  dsvbksb(U, w, V, n1, n, b_m, dcorr);

  std::cout << "dcorr.:" << std::endl;
  prt_vec(n, dcorr);

  free_dvector(dUb, 1, n); free_dvector(Ub0, 1, n);
  free_dvector(b_m, 1, mpn);
  free_dmatrix(U, 1, mpn, 1, n); free_dvector(w, 1, n);
  free_dmatrix(V, 1, n, 1, n);
  free_dmatrix(U_tp, 1, n, 1, mpn);
  free_dmatrix(A_ext, 1, mpn, 1, n); free_dvector(b_ext, 1, mpn);
  free_dmatrix(AtU_inv, 1, n, 1, n);
  free_dmatrix(n0_tp, 1, n, 1, n); free_dmatrix(d0, 1, n, 1, 2);

  free_dmatrix(N_m, 1, n, 1, n); free_dmatrix(U_m, 1, n, 1, n);
  free_dvector(w_m, 1, n); free_dmatrix(V_m, 1, n, 1, n);
  free_dmatrix(U_m_tp, 1, n, 1, n); free_dmatrix(n_m_tp, 1, n, 1, n);
  free_dvector(p_m, 1, n); 

  free_dmatrix(C, 1, n, 1, n); free_dmatrix(C_tp, 1, n, 1, n);
}


void SVD(const int m, const int n, double **A, double *b, double *dcorr)
{
  int    i, j;
  double *w, **U, **V;

  U = dmatrix(1, m, 1, n); w = dvector(1, n); V = dmatrix(1, n, 1, n);

  for (i = 1; i <= m; i++)
    for (j = 1; j <= n; j++)
      U[i][j] = A[i][j];

  dsvdcmp(U, m, n, w, V); dsvbksb(U, w, V, m, n, b, dcorr);
  
  std::cout << std::endl;
  std::cout << "corr.:" << std::endl;
  prt_vec(n, dcorr);

  free_dmatrix(U, 1, m, 1, n); free_dvector(w, 1, n);
  free_dmatrix(V, 1, n, 1, n);
}


void get_A_inv(const int m, const int n, float **U, float *w, float **V,
	       float **A_inv)
{
  int   j, k, l;
  float **V1;

  V1 = matrix(1, n, 1, n);

  for (j = 1; j <= n; j++)
    for (k = 1; k <= n; k++)
      V1[j][k] = V[j][k];

  for (j = 1; j <= n; j++)
    for (k = 1; k <= n; k++) {
      if (w[k] != 0.0)
	V1[j][k] /= w[k];
      else
	V1[j][k] = 0.0;
    }

  for (j = 1; j <= n; j++) {
    for (k = 1; k <= m; k++) {
      A_inv[j][k] = 0.0;
      for (l = 1; l <= n; l++)
	A_inv[j][k] += V1[j][l]*U[k][l];
    }
  }

  std::cout << std::endl;
  std::cout << "A^-1:" << std::endl;
  std::cout << std::endl;
  for (j = 1; j <= n; j++) {
    for (k = 1; k <= m; k++)
      std::cout << std::scientific << std::setprecision(2)
	   << std::setw(10) << A_inv[j][k];
    std::cout << std::endl;
  }

  free_matrix(V1, 1, n, 1, n);
}


double GetAngle(const double x, const double y)
{
  double z;

  if (x != 0.0)
    z = atan(y/x);
  else
    z = sign(y)*pi/2;
  if (x < 0.0) z += pi;
  return z;
}


tps do_PBs(const tps &a, const tps &b, const int n, const int jj[])
{
  int  k, jj1[100];
  tps  c;

  for (k = 0; k < n; k++)
    jj1[k] = jj[k];
  if (jj1[n-1] == 1) {
    c = b; jj1[n-1]--;
  } else if ((jj1[n-1] == 0) && (jj1[n-2] == 1)) {
    c = a; jj1[n-2]--;
  } else
    std::cout << "do_PBs: inconsistent parameters" << std::endl;
  for (k = n; k >= 1; k--)
    while (jj1[k-1] > 0) {
      if (k % 2 == 0)
	c = PB(b, c);
      else
	c = PB(a, c);
      jj1[k-1]--;
    }
  return c;
}


tps BCH(const tps &a, const tps &b, const int n_ord)
/* Evaluation of the Campbell-Baker-Hausdorff formula based on
   L.Corwin & F.P. Greenleaf "Representation of Nilpotent Lie Groups
   and their Applications, Part 1: Basic Theory and Examples".         */
{
  bool skip;
  int  i, j, m, twom, n, jj[100], m_ord, sgn, coeff;
  tps  c;

  c = 0.0;
  for (n = 1; n <= n_ord; n++) {
//    cout << endl;
    sgn = -1;
    for (m = 1; m <= n; m++) {
      sgn = -sgn; twom = 2*m;
      for (i = 0; i < twom; i++)
	jj[i] = 0;
      for (i = 0; i < (int)(pow((double)(n/m+1+1), (double)twom)+0.5); i++) {
	skip = false;
	for (j = 0; j < m; j++)
	  if (jj[2*j]+jj[2*j+1] == 0) {
	    skip = true; break;
	  }
	m_ord = 0; coeff = 1;
	for (j = 0; j < twom; j++) {
	  m_ord += jj[j]; coeff *= fact(jj[j]);
	}
	coeff *= sgn*n*m;
	if (! skip && (m_ord == n)) {
//	  std::cout << setw(2) << m_ord << ":" << std::setw(5) << coeff << ",";
//	  for (j = 0; j < twom; j++)
//	    std::cout << std::setw(2) << jj[j];
//	  std::cout << std::endl;
	  if ((jj[twom-1] == 1) || ((jj[twom-1] == 0) && (jj[twom-2] == 1))) {
	    // avoid call if zero
	    c = c + do_PBs(a, b, twom, jj)/coeff;
	  }
	}
	jj[0]++;
	for (j = 0; j < twom; j++)
	  if (jj[j] == n/m+1+1) {
	    jj[j] = 0; jj[j+1]++;
	  }
      }
    }
  }
  return c;
}

int get_Fnum(const char name[])
{
  int i, Fnum = 0;

  for (i = 0; i < max_Family; i++)
    if (strcmp(Families[i].Name, name) == 0) {
      Fnum = i + 1;
      break;
    }
  if (Fnum == 0) {
    std::cout << "get_Fnum: " << name << " not defined" << std::endl;
    exit(0);
  }
  return Fnum;
}


char* get_Name(const int Fnum)
{
  return Families[Fnum-1].Name;
}


int get_n_Kids(const int Fnum) { return Families[Fnum-1].n_Kids; }


long int get_loc(const int Fnum, const int Knum)
{
  int k = 0;

  if (Fnum <= max_Family)
    if (Knum <= max_Kids)
      k = Families[Fnum-1].Kids[Knum-1] + 1;
    else {
      std::cout << "Knum (" << Knum << ") > max_Kids (" << max_Kids
	   << ") for Fnum = " << Fnum << std::endl;
      exit(0);
    }
  else {
    std::cout << "Fnum (" << Fnum << ") > max_Family (" << max_Family << ")"
	 << std::endl;
    exit(0);
  }
  return k;
}


double get_L(const int Fnum, const int Knum)
{
  return elem[get_loc(Fnum, Knum)-1].L;
}


double get_bn(const int Fnum, const int Knum, const int n)
{
  return elem[get_loc(Fnum, Knum)-1].mpole->bn[n-1];
}


double get_bnL(const int Fnum, const int Knum, const int n)
{
  int    k;
  double bnL = 0.0;

  k = get_loc(Fnum, Knum) - 1;
  if (elem[k].L == 0.0)
    bnL = elem[k].mpole->bn[n-1];
  else
    bnL = elem[k].mpole->bn[n-1]*elem[k].L;
  return bnL;
}


void set_bn(const int Fnum, const int Knum, const int n, const double bn)
{
  int k;

  k = get_loc(Fnum, Knum) - 1;
  elem[k].mpole->bn[n-1] = bn; elem_tps[k].mpole->bn[n-1] = bn;
  if (n > elem[k].mpole->order) {
    elem[k].mpole->order = n; elem_tps[k].mpole->order = n;
  }
}


void set_bn(const int Fnum, const int n, const double bn)
{
  int j;

  for (j = 1; j <= get_n_Kids(Fnum); j++)
    set_bn(Fnum, j, n, bn);
}


void set_bnL(const int Fnum, const int Knum, const int n, const double bnL)
{
  int k;

  k = get_loc(Fnum, Knum) - 1;
  if (elem[k].L == 0.0) {
    elem[k].mpole->bn[n-1] = bnL; elem_tps[k].mpole->bn[n-1] = bnL;
  } else {
    elem[k].mpole->bn[n-1] = bnL/elem[k].L;
    elem_tps[k].mpole->bn[n-1] = bnL/elem[k].L;
  }
  if (n > elem[k].mpole->order) {
    elem[k].mpole->order = n; elem_tps[k].mpole->order = n;
  }
}


void set_bnL(const int Fnum, const int n, const double bnL)
{
  int j;

  for (j = 1; j <= get_n_Kids(Fnum); j++)
    set_bnL(Fnum, j, n, bnL);
}


void set_dbn(const int Fnum, const int Knum, const int n, const double dbn)
{
  int k;

  k = get_loc(Fnum, Knum) - 1;
  elem[k].mpole->bn[n-1] += dbn; elem_tps[k].mpole->bn[n-1] += dbn;
  if (n > elem[k].mpole->order) {
    elem[k].mpole->order = n; elem_tps[k].mpole->order = n;
  }
}


void set_dbn(const int Fnum, const int n, const double dbn)
{
  int j;

  for (j = 1; j <= get_n_Kids(Fnum); j++)
    set_dbn(Fnum, j, n, dbn);
}


void set_dbnL(const int Fnum, const int Knum, const int n, const double dbnL)
{
  int k;

  k = get_loc(Fnum, Knum) - 1;
  if (elem[k].L == 0.0) {
    elem[k].mpole->bn[n-1] += dbnL; elem_tps[k].mpole->bn[n-1] += dbnL;
  } else {
    elem[k].mpole->bn[n-1] += dbnL/elem[k].L;
    elem_tps[k].mpole->bn[n-1] += dbnL/elem[k].L;
  }
  if (n > elem[k].mpole->order) {
    elem[k].mpole->order = n; elem_tps[k].mpole->order = n;
  }
}


void set_dbnL(const int Fnum, const int n, const double dbnL)
{
  int j;

  for (j = 1; j <= get_n_Kids(Fnum); j++)
    set_dbnL(Fnum, j, n, dbnL);
}


void set_L(const int Fnum, const int Knum, const double L)
{
  int  k;

  k = get_loc(Fnum, Knum) - 1;
  elem[k].L = L; elem_tps[k].L = L;
}


void set_L(const int Fnum, const double L)
{
  int j;

  for (j = 1; j <= get_n_Kids(Fnum); j++)
    set_L(Fnum, j, L);
}


void set_dL(const int Fnum, const int Knum, const double dL)
{
  int k;

  k = get_loc(Fnum, Knum) - 1;
  elem[k].L += dL; elem_tps[k].L += dL;
}


void set_dL(const int Fnum, const double dL)
{
  int j;

  for (j = 1; j <= get_n_Kids(Fnum); j++)
    set_dL(Fnum, j, dL);
}


void set_bn_par(const int Fnum, const int Knum, const int n, const int j)
{
  // Set parameter dependence.
  int    k;
  double bn;

  k = get_loc(Fnum, Knum) - 1;
  bn = elem_tps[k].mpole->bn[n-1].cst();
  elem_tps[k].mpole->bn[n-1] = tps(bn, j);
  if (n > elem_tps[k].mpole->order) elem_tps[k].mpole->order = n;
}


void set_bn_par(const int Fnum, const int n, const int j)
{
  // Set parameter dependence.
  int  k;

  for (k = 1; k <= get_n_Kids(Fnum); k++)
    set_bn_par(Fnum, k, n, j);
}


void clr_bn_par(const int Fnum, const int Knum, const int n)
{
  // Clear parameter dependence.
  int    k;
  double bn;

  k = get_loc(Fnum, Knum) - 1;
  bn = elem_tps[k].mpole->bn[n-1].cst(); elem_tps[k].mpole->bn[n-1] = bn;
  // clear order
}


void clr_bn_par(const int Fnum, const int n)
{
  // Clear parameter dependence.
  int k;

  for (k = 1; k <= get_n_Kids(Fnum); k++)
    clr_bn_par(Fnum, k, n);
}


void set_L_par(const int Fnum, const int Knum, const int j)
{
  // Set parameter dependence.
  int    k;
  double L;

  k = get_loc(Fnum, Knum) - 1;
  L = elem_tps[k].L.cst(); elem_tps[k].L = tps(L, j);
}


void set_L_par(const int Fnum, const int j)
{
  // Set parameter dependence.
  int k;

  for (k = 1; k <= get_n_Kids(Fnum); k++)
    set_L_par(Fnum, k, j);
}


void clr_L_par(const int Fnum, const int Knum)
{
  // Clear parameter dependence.
  int    k;
  double L;

  k = get_loc(Fnum, Knum) - 1;
  L = elem_tps[k].L.cst(); elem_tps[k].L = L;
}


void clr_L_par(const int Fnum)
{
  // Clear parameter dependence.
  int k;

  for (k = 1; k <= get_n_Kids(Fnum); k++)
    clr_L_par(Fnum, k);
}


void set_s_par(const int Fnum, const int Knum, const int j)
{
  // Set s-dependence.

  long int k;
  double   L;

  const bool prt = false;

  // Point to multipole.
  k = get_loc(Fnum, Knum) - 1;

  if (prt)
    std::cout << "set_s_par: " << elem[k].Name << ", "  << std::setw(2) << Knum
	      << std::endl;

  switch (elem[k-1].Name[1]) {
  case 'u':
    if (elem[k+1].Name[1] == 'd') {
      L = elem_tps[k-1].L.cst(); elem_tps[k-1].L = tps(L, j);
      L = elem_tps[k+1].L.cst(); elem_tps[k+1].L = -tps(-L, j);
      printf("\n%s %s\n", elem_tps[k-1].Name, elem[k+1].Name);
    } else {
      std::cout << "set_s_par: configuration error " << elem[k+1].Name
	   << " (" << k+2 << ")" << std::endl;
      exit(1);
    }
    break;
  case 'd':
    if (elem[k+1].Name[1] == 'u') {
      L = elem_tps[k-1].L.cst(); elem_tps[k-1].L = -tps(-L, j);
      L = elem_tps[k+1].L.cst(); elem_tps[k+1].L = tps(L, j);
    } else {
      std::cout << "set_s_par: configuration error " << elem[k+1].Name
	   << " (" << k+2 << ")" << std::endl;
      exit(1);
    }
    break;
  default:
    std::cout << "set_s_par: configuration error " << elem[k-1].Name
	 << " (" << k << ")" << std::endl;
    exit(1);
    break;
  }
}


void set_s_par(const int Fnum, const int j)
{
  // Set s-dependence.
  int k;

  for (k = 1; k <= get_n_Kids(Fnum); k++)
    set_s_par(Fnum, k, j);
}


void clr_s_par(const int Fnum, const int Knum)
{
  // Clear s-dependence.
  int    k;
  double L;

  const bool  prt = false;

  // point to multipole
  k = get_loc(Fnum, Knum) - 1;

  if (prt)
    std::cout << "clr_s_par: " << elem[k].Name << ", "  << std::setw(2) << Knum
	      << std::endl;

  L = elem_tps[k-1].L.cst(); elem_tps[k-1].L = L;
  L = elem_tps[k+1].L.cst(); elem_tps[k+1].L = L;
}


void clr_s_par(const int Fnum)
{
  // Clear s-dependence.
  int k;

  for (k = 1; k <= get_n_Kids(Fnum); k++)
    clr_s_par(Fnum, k);
}


double get_BoBrho(const int Fnum, const int Knum)
{
  return elem[get_loc(Fnum, Knum)-1].wiggler->BoBrhoV[0];
}


void set_BoBrho(const int Fnum, const int Knum, const double BoBrho)
{
  int k;

  k = get_loc(Fnum, Knum) - 1;
  elem[k].wiggler->BoBrhoV[0] = BoBrho;
  elem_tps[k].wiggler->BoBrhoV[0] = BoBrho;
}


void get_cav(const int Fnum, const int Knum,
	     int &h_rf, double &V_rf, double &f_rf)
{
  int k;

  k = get_loc(Fnum, Knum) - 1;
  h_rf = elem[k].cavity->h_rf; V_rf = elem[k].cavity->V_rf;
  f_rf = elem[k].cavity->f_rf;
}


void get_Map(void)
{
  // Compute one-turn map.
  Map.identity(); Map.propagate(1, n_elem);
}


void get_Map(const ss_vect<double> &fixed_point)
{
  // Compute one-turn map.
  Map.identity();
  Map += fixed_point; Map.propagate(1, n_elem); Map -= fixed_point;
}


ss_vect<tps> get_Map(const int k)
{
  ss_vect<tps> Mk;

  // Propagate one-turn map to location k.
  Mk.identity(); Mk.propagate(1, k); Mk = Mk*Map*Inv(Mk);

  return Mk;
}


bool get_COD(const int i_max, const double eps, const double delta,
	     const bool prt)
{
  bool            cod;
  int             n_dim, n_iter, j, jj[ss_dim];
  double          dx_abs = 0.0;
  ss_vect<double> x1, dx;
  ss_vect<tps>    I, dx0;

  // danot_(1);

  n_dim = (cavity_on)? 6 : 4;

  fixed_point.zero(); 

  if (n_dim == 4) { fixed_point[delta_] = delta; fixed_point[ct_] = 0.0; }

  for (j = 0; j < ss_dim; j++)
    jj[j] = (j < n_dim)? 1 : 0;

  n_iter = 0; I.identity();
  for (n_iter = 1; n_iter <= i_max; n_iter++) {
    Map.identity(); Map += fixed_point; Map.propagate(1, n_elem);

    x1 = Map.cst(); dx = fixed_point - x1; dx0 = PInv(Map-I-x1, jj)*dx;

    dx_abs = 0.0;
    for (j = 0; j < n_dim; j++)
      dx_abs += sqr(dx0[j].cst());
    dx_abs = sqrt(dx_abs); fixed_point += dx0.cst();

    if (prt) {
      std::cout << std::endl;
      std::cout << std::scientific << std::setprecision(3)
	   << "get_COD: " << n_iter << ", Err = " << dx_abs 
	   << " (" << eps << ")" << std::endl;
      std::cout << std::scientific << std::setprecision(3)
	   << "dx:" << std::setw(11) << dx0.cst() << std::endl;
      if (true)
	std::cout << std::scientific << std::setprecision(5)
	     << "x:" << std::setw(13) << fixed_point << std::endl;
      else
	std::cout << std::scientific << std::setprecision(16)
	     << "x:" << std::setw(24) << fixed_point << std::endl;
    }

    if (dx_abs < eps) break;
  }
  
  cod = dx_abs < eps;

  // danot_(no);

  if (cod)
    get_Map(fixed_point);
  else {
    std::cout << "get_COD failed" << std::endl;
    exit(1);
  }

  return cod;
}


void get_A1(const double alpha[], const double beta[],
	    const double eta[], const double etap[])
{
  int          k;
  ss_vect<tps> Id;

  Map.zero(); A0.zero(); A1.zero(); g = 0.0; Id.identity();

  for (k = 0; k < 2; k++) {
    A1[2*k]  = sqrt(beta[k])*Id[2*k];
    A1[2*k+1] = -alpha[k]/sqrt(beta[k])*Id[2*k] + 1.0/sqrt(beta[k])*Id[2*k+1];
    A1[2*k] += eta[k]*Id[delta_]; A1[2*k+1] += etap[k]*Id[delta_];
  }
}


void get_A1(const double alphax, const double betax,
	    const double alphay, const double betay)
{
  int    k;
  double alpha[2], beta[2], eta[2], etap[2];

  for (k = 0; k < 2; k++) {
    eta[k] = 0e0; etap[k] = 0e0;
  }

  alpha[X_] = alphax; beta[X_] = betax; alpha[Y_] = alphay; beta[Y_] = betay;

  get_A1(alpha, beta, eta, etap);
}


tps get_h(void)
{
  // Parallel transport nonlinear kick to start of lattice.
  tps          h;
  ss_vect<tps> Map1, R;

  if (true) {
    // Dragt-Finn factorization
    h = LieFact_DF(Inv(A0*A1)*Map*A0*A1, R);
    R[6] = tps(0e0, 7); h = h*R;
    return h;
  } else {
    // single Lie exponent
    danot_(1); Map1 = Map; danot_(no_tps);
    return LieFact(Inv(A0*A1)*Map*Inv(Map1)*A0*A1);
  }
}


tps get_H(void)
{
  int          i;
  tps          H, gn;
  ss_vect<tps> Id, Mn;

  const bool prt = false;

  // Construct generator.
  // K is in Dragt-Finn form but the generators commute.
  Id.identity(); H = K;
  for (i = no_tps; i >= 3; i--) {
    gn = Take(g, i); H = H*LieExp(-gn, Id);
  }

  if (prt) {
    std::cout << H;
    std::cout << (H-H*Inv(A0*A1)*Map*A0*A1);
  }

  if (false) {
    // Normalize map (=> Map_res)
    for (i = 3; i <= no_tps; i++) {
      gn = Take(g, i); Mn = LieExp(gn, Id); Map = Inv(Mn)*Map*Mn;
    }
  }

  return H;
}


// transform from Floquet space to action-angle variables
ss_vect<double> get_J_phi(const ss_vect<double> ps)
{
  ss_vect<double>  ps1;

  ps1[0] = sqr(ps[x_]) + sqr(ps[px_]); ps1[1] = atan2(ps[px_], ps[x_]);
  ps1[2] = sqr(ps[y_]) + sqr(ps[py_]); ps1[3] = atan2(ps[py_], ps[y_]);

  return ps1;
}


ss_vect<tps> get_A_inv()
{
  int          i;
  ss_vect<tps> Id, A_inv;

  Id.identity(); A_inv = Inv(A1);
  for (i = 3; i <= no_tps; i++)
    A_inv = LieExp(Take(-g, i), Id)*A_inv;

  return A_inv;
}


ss_vect<tps> get_S(const int n_DOF)
{
  int          j;
  ss_vect<tps> S;

  S.zero();
  for (j = 0; j < n_DOF; j++) {
    S[2*j] = tps(0.0, 2*j+2); S[2*j+1] = -tps(0.0, 2*j+1);
  }

  return S;
}


ss_vect<tps> tp_S(const int n_DOF, const ss_vect<tps> &A)
{
  int          j, jj[ss_dim];
  ss_vect<tps> S;

  for (j = 1; j <= ss_dim; j++)
    jj[j-1] = (j <= 2*n_DOF)? 1 : 0;

  S = get_S(n_DOF);

  return -S*PInv(A, jj)*S;
}


bool is_h_ijklm(const int i1, const int j1, const int k1,
		const int l1, const int m1,
		const int i2, const int j2, const int k2,
		const int l2, const int m2)
{
  if ((i1 == i2) && (j1 == j2) && (k1 == k2) && (l1 == l2) && (m1 == m2))
    return true;
  else
    return false;
}


double h_ijklm(const tps &h, const int i, const int j, const int k,
	       const int l, const int m)
{
  int i1, jj[ss_dim];

  for (i1 = 0; i1 < ss_dim; i1++)
    jj[i1] = 0;
  jj[x_] = i; jj[px_] = j; jj[y_] = k; jj[py_] = l; jj[delta_] = m;
  return h[jj];
}


double h_ijklm_scl(const tps &h, const int i, const int j, const int k,
		   const int l, const int m,
		   const double nu_x, const double nu_y)
{
  int    i1, jj[ss_dim];
  double scl;

  for (i1 = 0; i1 < ss_dim; i1++)
    jj[i1] = 0;
  jj[x_] = i; jj[px_] = j; jj[y_] = k; jj[py_] = l; jj[delta_] = m;

  scl = 1.0;
  if ((i != j) || (k != l))
    scl /= max(sin(M_PI*((i-j)*nu_x+(k-l)*nu_y)), 0.1);
  else if ((i != 0) && (k != 0))
    // i == j, k == l
    scl /= max(sin(2.0*M_PI*nu_x)*sin(2.0*M_PI*nu_y), 0.1);
  else if (i != 0)
    // i == j, k == l
    scl /= max(sin(2.0*M_PI*nu_x), 0.1);
  else if (k != 0)
    // i == j, k == l
    scl /= max(sin(2.0*M_PI*nu_y), 0.1);

  return scl*h[jj];
}


double h_ijklm_p(const tps &h, const int i, const int j, const int k,
		 const int l, const int m, const int p)
{
  int  i1, jj[ss_dim];

  for (i1 = 0; i1 < ss_dim; i1++)
    jj[i1] = 0;
  jj[x_] = i; jj[px_] = j; jj[y_] = k; jj[py_] = l; jj[delta_] = m;
  jj[p-1] = 1;
  return h[jj];
}


double h_ijklm_p_scl(const tps &h, const int i, const int j, const int k,
		     const int l, const int m, const int p,
		     const double nu_x, const double nu_y)
{
  int    i1, jj[ss_dim];
  double scl;

  for (i1 = 0; i1 < ss_dim; i1++)
    jj[i1] = 0;
  jj[x_] = i; jj[px_] = j; jj[y_] = k; jj[py_] = l; jj[delta_] = m;
  jj[p-1] = 1;

  if ((i != j) || (k != l))
    scl = 1.0/max(sin(M_PI*((i-j)*nu_x+(k-l)*nu_y)), 0.1);
  else
    scl = 1.0;

  return scl*h[jj];
}


void get_nu_ksi(const ss_vect<tps> &nus, double nu[], double ksi[])
{
  int i, jj[ss_dim];

  nu[X_] = nus[3].cst(); nu[Y_] = nus[4].cst(); nu[Z_] = nus[5].cst();

  for (i = 0; i < ss_dim; i++)
    jj[i] = 0;

  jj[delta_] = 1; ksi[X_] = nus[3][jj]; ksi[Y_] = nus[4][jj];
}


void prt_nu(const ss_vect<tps> &nus)
{
  double nu[3], ksi[2];

  const int precision = 16, width = 23;

  get_nu_ksi(nus, nu, ksi);
  std::cout << std::fixed << std::setprecision(precision)
       << " nu_x = " << nu[0] << ",  nu_y = " << nu[1];
  if (!cavity_on)
   std::cout << std::scientific << std::setprecision(precision)
	     << ", alphac = " << std::setw(width)
	     << h_ijklm(Map[ct_], 0, 0, 0, 0, 1)/elem[n_elem-1].S
	     << std::endl;
  else
    std::cout << std::scientific << std::setprecision(precision) << ", nu_s = "
	 << std::setw(width) << nu[2] << std::endl;
  std::cout << "ksi_x = " << std::setw(width) << ksi[0]
       << ", ksi_y = " << std::setw(width) << ksi[1] << std::endl;
}


void get_k_J(const tps &K,
	     double k_J2[], double k_J3[], double k_J4[], double k_delta2[])
{
  // k_J^2
  k_J2[0] = h_ijklm(K, 2, 2, 0, 0, 0);
  k_J2[1] = h_ijklm(K, 1, 1, 1, 1, 0);
  k_J2[2] = h_ijklm(K, 0, 0, 2, 2, 0);
 
  // k_J^3
  k_J3[0] = h_ijklm(K, 3, 3, 0, 0, 0);
  k_J3[1] = h_ijklm(K, 2, 2, 1, 1, 0);
  k_J3[2] = h_ijklm(K, 1, 1, 2, 2, 0);
  k_J3[3] = h_ijklm(K, 0, 0, 3, 3, 0);

  // k_J^4
  k_J4[0] = h_ijklm(K, 4, 4, 0, 0, 0);
  k_J4[1] = h_ijklm(K, 3, 3, 1, 1, 0);
  k_J4[2] = h_ijklm(K, 2, 2, 2, 2, 0);
  k_J4[3] = h_ijklm(K, 1, 1, 3, 3, 0);
  k_J4[4] = h_ijklm(K, 0, 0, 4, 4, 0);

  // k_delta^2
  k_delta2[0] = h_ijklm(K, 1, 1, 0, 0, 1);
  k_delta2[1] = h_ijklm(K, 0, 0, 1, 1, 1);
}


void prt_dnu(const tps &K)
{
  double k_J2[3], k_J3[4], k_J4[5], k_delta2[2];

  get_k_J(K, k_J2, k_J3, k_J4, k_delta2);
  std::cout << "K_J:         "
       << std::scientific << std::setprecision(2)
       << "a_11000 = " << std::setw(9) << k_J2[0]
       << ", a_11000 = " << std::setw(9) << k_J2[1]
       << ", a_00110 = " << std::setw(9) << k_J2[2] << std::endl;
/*  cout << "K_J^2:     "
       << std::scientific << std::setprecision(2)
       << "a_22000 = " << std::setw(9) << k_J3[0]
       << ", a_11110 = " << std::setw(9) << k_J3[1]
       << ", a_11110 = " << std::setw(9) << k_J3[2]
       << ", a_00220 = " << std::setw(9) << k_J3[3] << std::endl;
  cout << "K_delta^2: "
       << std::scientific << std::setprecision(2)
       << "a_11002 = " << std::setw(9) << k_delta2[0]
       << ", a_00112 = " << std::setw(9) << k_delta2[1] << std::endl;*/
}


void get_alphac(double alphac[])
{
  /* alphac = d_M[ct]/d_delta */
  int i;

  for (i = 0; i < n_alphac; i++)
    alphac[i] = h_ijklm(Map[ct_], 0, 0, 0, 0, i+1)/elem[n_elem-1].S;
}


void prt_alphac()
{
  double alphac[n_alphac], po2, q, pm;

  get_alphac(alphac);
  std::cout << std::endl;
  std::cout << std::scientific << std::setprecision(3)
       << "alphac = " << std::setw(10) << alphac[0] << " + "
        << std::setw(10) << alphac[1] << "*delta + "
        << std::setw(10) << alphac[2] << "*delta^2" << std::endl;

  if (n_alphac >= 2) {
    std::cout << std::endl;
    std::cout << std::fixed << std::setprecision(2)
	 << "Fixed points to O(3): " << 0
	 << ", " << -1e2*alphac[0]/alphac[1] <<"%" << std::endl;
  }

  if (n_alphac >= 3) {
    po2 = alphac[1]/(2.0*alphac[2]); q = alphac[0]/alphac[2];
    pm = sqrt(sqr(po2)-q);
    std::cout << std::endl;
    std::cout << std::fixed << std::setprecision(2)
	 << "Fixed points to O(4): " << 0
	 << ", " << -1e2*(po2+pm) << "%, "
	 << -1e2*(po2-pm) << "%" << std::endl;
  }
}


double H_long(const double phi, const double delta,
	      const int h_rf, const double V_rf, const double phi0,
	      const double alphac[])
{
  int    i;
  double H;

  H = V_rf/(E0*1e9)*(cos(phi+phi0)+phi*sin(phi0));
  for (i = 2; i <= n_alphac+1; i++)
    H += 2.0*pi*h_rf*alphac[i-2]*pow(delta, (double)i)/i;
  return H;
}


void prt_H_long(const int n, const double phi_max, const double delta_max,
		const double U0, const bool neg_alphac)
{
  int            i, j, h_rf;
  double         V_rf, f_rf, alphac[n_alphac], phi, delta, H, delta_rf, phi0;
  std::ofstream  os;

  os.open("H_long.dat", std::ios::out);

  get_alphac(alphac); get_cav(get_Fnum("cav"), 1, h_rf, V_rf, f_rf);

  phi0 = - fabs(asin(U0/V_rf));
  if (neg_alphac) phi0 += pi;

  delta_rf = sqrt(-V_rf*cos(pi+phi0)*(2e0-(pi-2e0*(pi+phi0))*tan(pi+phi0))
	     /(alphac[0]*pi*h_rf*E0*1e9));
  std::cout << std::endl;
  std::cout << std::fixed << std::setprecision(1)
	    << "U0               = " << 1e-3*U0 << " keV" << std::endl;
  if (!neg_alphac) 
    std::cout << std::fixed << std::setprecision(2)
	      << "phi0             = " << fabs(phi0)*180e0/pi-180e0
	      << " deg" << std::endl;
  else
    std::cout << std::fixed << std::setprecision(2)
	      << "phi0             = 180 - " << fabs(phi0)*180e0/pi-180e0
	      << " deg" << std::endl;
  std::cout << std::endl;
  std::cout << std::fixed << std::setprecision(2)
	    << "RF bucket height = " << 1e2*delta_rf << "%" << std::endl;

  for (i = -n; i <= n ; i++) {
    for (j = -n; j <= n ; j++) {
      phi = i*phi_max/n; delta = j*delta_max/n;
      H = H_long(phi, delta, h_rf, V_rf, pi+phi0, alphac);
      os << std::fixed
	 << std::setprecision(2)<< std::setw(8) <<  (phi0+phi)*180e0/pi
	 << std::setprecision(5) << std::setw(10) << 1e2*delta
	 << std::scientific << std::setw(13) << H << std::endl;
    }
    os << std::endl;
  }
  os.close();
}


double arccos_(double x)
{
  double y;

  if (x == 0e0)
    y = 1.0;
  else if (x >= 0.0)
    y = atan(sqrt(1.0/(x*x)-1.0));
  else
    y = (pi-atan(sqrt(1.0/(x*x)-1.0)));

  return y;
}


void get_nu_symp(const ss_vect<tps> &Map, double nu[])
{
  /* Get nu from a symplectic matrix */
  int    i, j;
  double sgn, detp, detm, b, c, th, tv, b2mc;
  float  **M;

  const int n_dim = 4;

  M = matrix(1, ps_dim, 1, ps_dim);

  get_matrix(Map, M);

  for (i = 1; i <= n_dim; i++)
    M[i][i] -= 1.0;
  detp = Det(n_dim, M);
  for (i = 1; i <= n_dim; i++)
    M[i][i] += 2.0;
  detm = Det(n_dim, M);
  for (i = 1; i <= n_dim; i++)
    M[i][i] -= 1.0;
  b = (detp-detm)/16.0; c = (detp+detm)/8.0-1.0;
  th = (M[1][1]+M[2][2])/2.0; tv = (M[3][3]+M[4][4])/2.0;
  b2mc = b*b - c;
  if (b2mc < 0e0) {
    nu[X_] = -1.0; nu[Y_] = -1.0;
    std::cout << "*** get_nu_symp: unstable in tune" << std::endl;
    return;
  }
  sgn = (th > tv)? 1.0 : -1.0;

  nu[X_] = arccos_(sgn*sqrt(b2mc)-b)/(2.0*pi);
  nu[Y_] = arccos_(-b-sgn*sqrt(b2mc))/(2.0*pi);

  for (i = 0; i <= 1; i++) {
    j = (i+1)*2 - 1;
    if (M[j][j+1] < 0.0) nu[i] = 1.0 - nu[i];
  }
  std::cout << std::fixed << std::setprecision(5)
       << "nu_x = " << nu[X_] << ", nu_y = " << nu[Y_] << std::endl;

  free_matrix(M, 1, ps_dim, 1, ps_dim);
}


inline bool Check_Ampl(const ss_vect<double> &x)
{
  bool lost;

  // check for inf and NaN
  if (std::isinf(x[x_]) || std::isnan(x[x_]) || std::isinf(x[y_])
      || std::isnan(x[y_]))
    lost = true;
  else
    lost = (fabs(x[x_]) > max_ampl[X_]) || (fabs(x[y_]) > max_ampl[Y_]);
  return lost;
}


inline bool Check_Ampl(const ss_vect<tps> &x)
{
  bool lost;

  // check for inf and NaN
  if (std::isinf(x[x_].cst()) || std::isnan(x[x_].cst()) ||
      std::isinf(x[y_].cst()) || std::isnan(x[y_].cst()))
    lost = true;
  else
    lost = (fabs(x[x_].cst()) > max_ampl[X_]) ||
           (fabs(x[y_].cst()) > max_ampl[Y_]);
  return lost;
}


bool track(const double x, const double px, const double y, const double py,
	   const double delta, const long int n, const double f_rf,
	   const bool prt)
{
  long int        i;
  ss_vect<double> ps;
  std::ofstream   os;

  danot_(1);

  ps[x_] = x; ps[px_] = px; ps[y_] = y; ps[py_] = py;
  ps[delta_] = delta; ps[ct_] = 0.0;

  if (prt) {
    os.open("track.out", std::ios::out);
    os << "# Tracking with Thor" << std::endl;
    os << "#" << std::endl;
    os << "#  n       x           p_x          y            p_y  "
       << "       delta         cdt" << std::endl;
    if (f_rf == 0.0) {
      os << "#         [mm]        [mrad]       [mm]         [mrad]"
	 << "                    [mm]" << std::endl;
      os << std::scientific << std::setprecision(16)
	 << std::setw(4) << 0
	 << std::setw(24) << 1e3*ps[x_] << std::setw(24) << 1e3*ps[px_]
	 << std::setw(24) << 1e3*ps[y_] << std::setw(24) << 1e3*ps[py_]
	 << std::setw(24) << 1e2*ps[delta_] 
	 << std::setw(24) << 1e3*ps[ct_] << std::endl;
    } else {
      os << "#         [mm]        [mrad]       [mm]         [mrad]"
	 << "                    [deg]" << std::endl;
      os << std::scientific << std::setprecision(16)
	 << std::setw(4) << 0
	 << std::setw(24) << 1e3*ps[x_] << std::setw(24) << 1e3*ps[px_]
	 << std::setw(24) << 1e3*ps[y_] << std::setw(24) << 1e3*ps[py_]
	 << std::setw(24) << 1e2*ps[delta_] 
	 << std::setw(24) << 2.0*f_rf*180.0*ps[ct_]/clight << std::endl;
    }
    os << "#" << std::endl;
  }

  for (i = 1; i <= n; i++) {
    if (ps.propagate(1, n_elem)) {
      if (prt) {
	if (f_rf == 0.0)
	  os << std::scientific << std::setprecision(16)
	     << std::setw(4) << i
	     << std::setw(24) << 1e3*ps[x_] << std::setw(24) << 1e3*ps[px_]
	     << std::setw(24) << 1e3*ps[y_] << std::setw(24) << 1e3*ps[py_]
	     << std::setw(24) << 1e2*ps[delta_] 
	     << std::setw(24) << 1e3*ps[ct_] << std::endl;
	else
	  os << std::scientific << std::setprecision(16)
	     << std::setw(4) << i
	     << std::setw(24) << 1e3*ps[x_] << std::setw(24) << 1e3*ps[px_]
	     << std::setw(24) << 1e3*ps[y_] << std::setw(24) << 1e3*ps[py_]
	     << std::setw(24) << 1e2*ps[delta_] 
	     << std::setw(24) << 2.0*f_rf*180.0*ps[ct_]/clight << std::endl;
      }
    } else
      return false;
  }
  if (prt) os.close();

  return true;
}


bool map_track(const double x, const double px,
	       const double y, const double py,
	       const double delta, const long int n,
	       const double f_rf, const bool prt)
{
  bool          lost;
  long int      i;
  ss_vect<tps>  ps;
  std::ofstream os;

  danot_(1);

  ps[x_] = x; ps[px_] = px; ps[y_] = y; ps[py_] = py;
  ps[delta_] = delta; ps[ct_] = 0.0;
  lost = Check_Ampl(ps);
  if (lost) return false;

  if (prt) {
    os.open("track.out", std::ios::out);
    os << "# Tracking with Thor" << std::endl;
    os << "#" << std::endl;
    os << "#  n       x           p_x          y            p_y  "
       << "       delta         cdt" << std::endl;
    if (f_rf == 0.0) {
      os << "#         [mm]        [mrad]       [mm]         [mrad]"
	 << "                    [mm]" << std::endl;
      os << std::endl;
      os << std::scientific << std::setprecision(16)
	 << std::setw(4) << 0
	 << std::setw(24) << 1e3*ps[x_].cst()
	 << std::setw(24) << 1e3*ps[px_].cst()
	 << std::setw(24) << 1e3*ps[y_].cst()
	 << std::setw(24) << 1e3*ps[py_].cst()
	 << std::setw(24) << 1e2*ps[delta_].cst() 
	 << std::setw(24) << 1e3*ps[ct_].cst() << std::endl;
    } else {
      os << "#         [mm]        [mrad]       [mm]         [mrad]"
	 << "                    [deg]" << std::endl;
      os << std::endl;
      os << std::scientific << std::setprecision(16)
	 << std::setw(4) << 0
	 << std::setw(24) << 1e3*ps[x_].cst()
	 << std::setw(24) << 1e3*ps[px_].cst()
	 << std::setw(24) << 1e3*ps[y_].cst()
	 << std::setw(24) << 1e3*ps[py_].cst()
	 << std::setw(24) << 1e2*ps[delta_].cst() 
	 << std::setw(24) << 2.0*f_rf*180.0*ps[ct_].cst()/clight << std::endl;
    }
  }

  for (i = 1; i <= n; i++) {
    ps = Map*ps; lost = Check_Ampl(ps);
    if (lost) return false;
    if (prt) {
      if (f_rf == 0.0)
	os << std::scientific << std::setprecision(16)
	   << std::setw(4) << 0
	   << std::setw(24) << 1e3*ps[x_].cst()
	   << std::setw(24) << 1e3*ps[px_].cst()
	   << std::setw(24) << 1e3*ps[y_].cst()
	   <<std:: setw(24) << 1e3*ps[py_].cst()
	   << std::setw(24) << 1e2*ps[delta_].cst() 
	   << std::setw(24) << 1e3*ps[ct_].cst() << std::endl;
      else
	os << std::scientific << std::setprecision(16)
	   << std::setw(4) << 0
	   << std::setw(24) << 1e3*ps[x_].cst()
	   << std::setw(24) << 1e3*ps[px_].cst()
	   << std::setw(24) << 1e3*ps[y_].cst()
	   << std::setw(24) << 1e3*ps[py_].cst()
	   << std::setw(24) << 1e2*ps[delta_].cst() 
	   << std::setw(24) << 2.0*f_rf*180.0*ps[ct_].cst()/clight << std::endl;
    }
  }
  if (prt) os.close();

  return true;
}


void get_r_stable(double &r, const double phi, const double delta,
		  const long int n, const double eps, bool map)
{
  /* Binary search for dynamical aperture. */
  bool   lost = false;
  double r_min = 0.0, r_max = r;

  while (!lost ) {
    if (!map)
      lost = ! track(r_max*cos(phi), 0.0, r_max*sin(phi), 0.0, delta,
		     n, 0, false);
    else
      lost = ! map_track(r_max*cos(phi), 0.0, r_max*sin(phi), 0.0, delta,
			 n, 0, false);
    r_max *= 2.0;
  }
  while (r_max-r_min >= eps) {
    r = r_min + (r_max-r_min)/2.0;
    if (!map)
      lost = !track(r*cos(phi), 0.0, r*sin(phi), 0.0, delta, n, 0, false);
    else
      lost = !map_track(r*cos(phi), 0.0, r*sin(phi), 0.0, delta, n, 0, false);
    if (!lost)
      r_min = r;
    else
      r_max = r;
  }
  r = r_min + (r_max-r_min)/2.0;
}


ss_vect<tps> get_eta(const long int k)
{
  // parameter dependance for the periodic solution
  int          i;
  ss_vect<tps> eta;

  eta.zero();
  for (i = 1; i < ps_dim; i++)
    eta[i-1] = tps(0.0, i);

  eta = LieExp(g, A0*eta);

//  eta += fixed_point;
  eta.propagate(1, k);
//  eta -= eta.cst();

  return eta;
}


void get_ab(double alpha[], double beta[], const long int k)
{
  ss_vect<tps> a1, A1_A1tp;

  a1 = A1;
//  a1 += fixed_point;
  a1.propagate(1, k);
//  a1 -= a1.cst();

  A1_A1tp = a1*tp_S(2, a1);

  alpha[X_] = -h_ijklm(A1_A1tp[x_], 0, 1, 0, 0, 0);
  alpha[Y_] = -h_ijklm(A1_A1tp[y_], 0, 0, 0, 1, 0);
  beta[X_]  =  h_ijklm(A1_A1tp[x_], 1, 0, 0, 0, 0);
  beta[Y_]  =  h_ijklm(A1_A1tp[y_], 0, 0, 1, 0, 0);
}


ss_vect<tps> get_A1_CS(const ss_vect<tps> &A1, tps dnu[])
{
  int          k;
  double       c, s;
  tps          a11, a12;
  ss_vect<tps> Id, R;

  Id.identity(); R.identity();
  for (k = 0; k <= 1; k++) {
    a11 = Der(A1[2*k], 2*k+1); a12 = Der(A1[2*k], 2*k+2);
    dnu[k] = atan2(a12, a11)/(2.0*M_PI);
    if (dnu[k] < 0.0) dnu[k] += 1.0;

    c = cos(2.0*M_PI*dnu[k].cst()); s = sin(2.0*M_PI*dnu[k].cst());
    R[2*k] = c*Id[2*k] - s*Id[2*k+1]; R[2*k+1] = s*Id[2*k] + c*Id[2*k+1];
  }

  return A1*R;
}


void get_ab(tps ab[], tps nu[], const long int loc)
{
  /* parameter dependance for the periodic solution:


        -1       T           |  0  I |        T   | beta   -alpha |
       A   = -S A  S,    S = |       |,    A A  = |               |
                             | -I  0 |            | -alpha  gamma |

       alpha_x = -h_01000(A_Atp[x_]), alpha_y = -h_00010(A_Atp[y_]),
       beta_x  =  h_10000(A_Atp[x_]), beta_y  =  h_00100(A_Atp[y_])

  */

  ss_vect<tps> A1_prm, A1_A1tp;

  A1_prm = LieExp(g, A1);

//  A1_prm += fixed_point;
  A1_prm.propagate(1, loc);
//  A1_prm -= A1_prm.cst();

  get_A1_CS(A1_prm, nu);

  A1_A1tp = A1_prm*tp_S(2, A1_prm);

  ab[X_] = A1_A1tp[x_]; ab[Y_] = A1_A1tp[y_];
}


void get_ab(tps ab[], tps dnu[], const long int k1, const long int k2)
{
  ss_vect<tps> A1_prm, A1_A1tp;

  A1_prm = LieExp(g, A1);

//  A1_prm += fixed_point;
  A1_prm.propagate(k1, k2);
// A1_prm -= A1_prm.cst();

  get_A1_CS(A1_prm, dnu);

  A1_A1tp = A1_prm*tp_S(2, A1_prm);

  ab[X_] = A1_A1tp[x_]; ab[Y_] = A1_A1tp[y_];
}


void get_nu(double nu[], const long int k)
{
  double       a11, a12;
  ss_vect<tps> a1;

  a1 = A1;

//  a1 += fixed_point;
  a1.propagate(1, k);
//  a1 -= a1.cst();

  a11 = h_ijklm(a1[x_], 1, 0, 0, 0, 0); a12 = h_ijklm(a1[x_], 0, 1, 0, 0, 0);
  nu[X_] = atan2(a12, a11)/(2.0*pi);
  if (nu[X_] < 0.0) nu[X_] += 1.0;

  a11 = h_ijklm(a1[y_], 0, 0, 1, 0, 0); a12 = h_ijklm(a1[y_], 0, 0, 0, 1, 0);
  nu[Y_] = atan2(a12, a11)/(2.0*pi);
  if (nu[Y_] < 0.0) nu[Y_] += 1.0;
}


void bend_cal(const int Fnum1)
{
  /* Assumes a thin dipole has been inserted up- and down stream
     of the dipole. */
  int          k, n, Fnum0, Fnum2;
  double       b0L[2], A[2][2], A_inv[2][2], det;
  ss_vect<tps> x;

  const int    n_max = 100;
  const double eps   = 1e-15;

  std::cout << std::endl;
  std::cout << "bend_cal" << std::endl;

  k = get_loc(Fnum1, 1); Fnum0 = elem[k-2].Fnum; Fnum2 = elem[k].Fnum;

  n = 0; b0L[0] = 0.0; b0L[1] = 0.0;
  do {
    n++;

    set_bn_par(Fnum0, 1, Dip, 7);
    x.zero(); x.propagate(get_loc(Fnum0, 1), get_loc(Fnum2, 1));
    clr_bn_par(Fnum0, 1, Dip);
    A[0][0] = h_ijklm_p(x[x_], 0, 0, 0, 0, 0, 7);
    A[1][0] = h_ijklm_p(x[px_], 0, 0, 0, 0, 0, 7);

    set_bn_par(Fnum2, 1, Dip, 7);
    x.zero(); x.propagate(get_loc(Fnum0, 1), get_loc(Fnum2, 1));
    clr_bn_par(Fnum2, 1, Dip);
    A[0][1] = h_ijklm_p(x[x_], 0, 0, 0, 0, 0, 7);
    A[1][1] = h_ijklm_p(x[px_], 0, 0, 0, 0, 0, 7);

    det = A[0][0]*A[1][1] - A[0][1]*A[1][0];

    A_inv[0][0] =  A[1][1]/det; A_inv[1][1] =  A[0][0]/det;
    A_inv[0][1] = -A[0][1]/det; A_inv[1][0] = -A[1][0]/det;

    b0L[0] -= A_inv[0][0]*x[x_].cst() + A_inv[0][1]*x[px_].cst();
    b0L[1] -= A_inv[1][0]*x[x_].cst() + A_inv[1][1]*x[px_].cst();

    set_bn(Fnum0, 1, Dip, b0L[0]); set_bn(Fnum2, 1, Dip, b0L[1]);
  } while (((fabs(x[x_].cst()) > eps) || (fabs(x[px_].cst()) > eps))
	   && (n <= n_max));

  std::cout << std::endl;
  std::cout << std::scientific << std::setprecision(1)
       << "bend_cal: x = " << x[x_].cst() << ", p_x = " << x[px_].cst()
       << std::setprecision(3)
       << ", b_1,entrance = " << b0L[0]
       << ", b_1,exit = " << b0L[1] << std::endl;

  if (n <= n_max) {
    for (k = 1; k <= get_n_Kids(Fnum0); k++)
      set_bn(Fnum0, k, Dip, b0L[0]);
    for (k = 1; k <= get_n_Kids(Fnum2); k++)
      set_bn(Fnum2, k, Dip, b0L[1]);
  }
}


double get_bn_s(const int Fnum, const int Knum, const int n)
{
  long int k;
  double   bn;

  if (Fnum > 0)
    bn = get_bn(Fnum, Knum, n);
  else {
    k = get_loc(abs(Fnum), Knum) - 1;

    switch (elem[k-1].Name[1]) {
    case 'u':
      bn = elem[k-1].L;
      break;
    case 'd':
      bn = elem[k+1].L;
      break;
    default:
      std::cout << "get_bn_s: configuration error " << elem[k-1].Name
	   << " (" << k << ")" << std::endl;
      exit(1);
      break;
    }
  }

  return bn;
}


double get_bnL_s(const int Fnum, const int Knum, const int n)
{
  long int k;
  double   bnL;

  if (Fnum > 0)
    bnL = get_bnL(Fnum, Knum, n);
  else {
    k = get_loc(abs(Fnum), Knum) - 1;

    switch (elem[k-1].Name[1]) {
    case 'u':
      bnL = elem[k-1].L;
      break;
    case 'd':
      bnL = elem[k+1].L;
      break;
    default:
      std::cout << "get_bnL_s: configuration error " << elem[k-1].Name
	   << " (" << k << ")" << std::endl;
      exit(1);
      break;
    }
  }

  return bnL;
}


void set_dbn_s(const int Fnum, const int Knum, const int n, const double dbn)
{
  long int k;

  if (Fnum > 0)
    set_dbn(Fnum, Knum, n, dbn);
  else {
    // point to multipole
    k = get_loc(abs(Fnum), Knum) - 1;

    switch (elem[k-1].Name[1]) {
    case 'u':
      if (elem[k+1].Name[1] == 'd') {
	set_dL(elem[k-1].Fnum, elem[k-1].Knum, dbn);
	set_dL(elem[k+1].Fnum, elem[k+1].Knum, -dbn);
      } else {
	std::cout << "set_dbn_s: configuration error " << elem[k+1].Name
	     << " (" << k+2 << ")" << std::endl;
	exit(1);
      }
      break;
    case 'd':
      if (elem[k+1].Name[1] == 'u') {
	set_dL(elem[k-1].Fnum, elem[k-1].Knum, -dbn);
	set_dL(elem[k+1].Fnum, elem[k+1].Knum, dbn);
      } else {
	std::cout << "set_dbn_s: configuration error " << elem[k+1].Name
	     << " (" << k+2 << ")" << std::endl;
	exit(1);
      }
      break;
    default:
      std::cout << "set_dbn_s: configuration error " << elem[k-1].Name
	   << " (" << k << ")" << std::endl;
      exit(1);
      break;
    }
  }
}


void set_dbn_s(const int Fnum, const int n, const double dbn)
{
  int k;

  for (k = 1; k <= get_n_Kids(abs(Fnum)); k++)
    set_dbn_s(Fnum, k, n, dbn);
}


void set_bn_s(const int Fnum, const int Knum, const int n, const double dbn)
{
  long int k;

  if (Fnum > 0)
    set_bn(Fnum, Knum, n, dbn);
  else {
    // point to multipole
    k = get_loc(abs(Fnum), Knum) - 1;

    switch (elem[k-1].Name[1]) {
    case 'u':
      if (elem[k+1].Name[1] == 'd') {
	set_L(elem[k-1].Fnum, elem[k-1].Knum, dbn);
	set_L(elem[k+1].Fnum, elem[k+1].Knum, -dbn);
      } else {
	std::cout << "set_dbn_s: configuration error " << elem[k+1].Name
	     << " (" << k+2 << ")" << std::endl;
	exit(1);
      }
      break;
    case 'd':
      if (elem[k+1].Name[1] == 'u') {
	set_L(elem[k-1].Fnum, elem[k-1].Knum, -dbn);
	set_L(elem[k+1].Fnum, elem[k+1].Knum, dbn);
      } else {
	std::cout << "set_dbn_s: configuration error " << elem[k+1].Name
	     << " (" << k+2 << ")" << std::endl;
	exit(1);
      }
      break;
    default:
      std::cout << "set_dbn_s: configuration error " << elem[k-1].Name
	   << " (" << k << ")" << std::endl;
      exit(1);
      break;
    }
  }
}


void set_bn_s(const int Fnum, const int n, const double dbn)
{
  int k;

  for (k = 1; k <= get_n_Kids(abs(Fnum)); k++)
    set_bn_s(Fnum, k, n, dbn);
}


void fit_alpha(const double alpha0_x, const double beta0_x,
	       const double alpha0_y, const double beta0_y,
	       const long int k1, const long int k2,
	       const int n_b2, const int b2s[],
	       const double eps, const bool prt)
{
  int           i, j, k;
  double        **A, *b, *b2_lim, *b2, *db2, step;
  double        dalpha1[2];
  tps           ab[2], nu[2];
  std::ofstream quad_out;

  const bool   debug = true;
  const int    m     = 2;
  const double s_cut = 1e-5, step0 = 0.7;

  b = dvector(1, m); b2_lim = dvector(1, n_b2);
  b2 = dvector(1, n_b2); db2 = dvector(1, n_b2);
  A = dmatrix(1, m, 1, n_b2);

  for (i = 1; i <= n_b2; i++) {
    if (b2s[i-1] > 0) {
      b2_lim[i] = b2_max; b2[i] = get_bn(b2s[i-1], 1, Quad);
    } else {
      b2_lim[i] = ds_max; b2[i] = get_L(abs(b2s[i-1]), 1);
    }
  }

  std::cout << std::endl;
  std::cout << std::fixed << std::setprecision(5)
       << "fit_alpha: " << k1 << " - " << k2
       << ", alpha0_x = " << std::setw(8) << alpha0_x
       << ", alpha0_y = " << std::setw(8) << alpha0_y
       << ", beta0_x = " << std::setw(8) << beta0_x
       << ", beta0_y = " << std::setw(8) << beta0_y << std::endl;

  danot_(3);

  get_A1(alpha0_x, beta0_x, alpha0_y, beta0_y); get_ab(ab, nu, k1+1, k2);

  dalpha1[X_] = -h_ijklm(ab[X_], 0, 1, 0, 0, 0);
  dalpha1[Y_] = -h_ijklm(ab[Y_], 0, 0, 0, 1, 0);

   while ((fabs(dalpha1[X_]) > eps) || (fabs(dalpha1[Y_]) > eps)) {
    step = step0;
    for (i = 1; i <= n_b2; i++) {
      for (j = 1; j <= get_n_Kids(abs(b2s[i-1])); j++)
	if (b2s[i-1] > 0)
	  set_bn_par(b2s[i-1], j, Quad, 7);
	else
	  set_s_par(abs(b2s[i-1]), j, 7);

      get_A1(alpha0_x, beta0_x, alpha0_y, beta0_y); get_ab(ab, nu, k1+1, k2);

      A[1][i] = h_ijklm_p(ab[X_], 0, 1, 0, 0, 0, 7);
      A[2][i] = h_ijklm_p(ab[Y_], 0, 0, 0, 1, 0, 7);

      if (b2s[i-1] < 0)
	for (j = 1; j <= m; j++)
	  A[j][i] *= scl_ds;

      for (j = 1; j <= get_n_Kids(abs(b2s[i-1])); j++)
	if (b2s[i-1] > 0)
	  clr_bn_par(b2s[i-1], j, Quad);
	else
	  clr_s_par(abs(b2s[i-1]), j);
    }

    b[1] = dalpha1[X_]; b[2] = dalpha1[Y_];

    SVD_lim(m, n_b2, A, b, b2_lim, s_cut, b2, db2);

    for (i = 1; i <= n_b2; i++) {
      set_dbn_s(b2s[i-1], Quad, step*db2[i]);
      b2[i] = get_bn_s(b2s[i-1], 1, Quad);
    }

    get_A1(alpha0_x, beta0_x, alpha0_y, beta0_y); get_ab(ab, nu, k1+1, k2);

    dalpha1[X_] = -h_ijklm(ab[X_], 0, 1, 0, 0, 0);
    dalpha1[Y_] = -h_ijklm(ab[Y_], 0, 0, 0, 1, 0);

    if (debug) {
      std::cout << std::endl;
      std::cout << " Ax = b:" << std::endl;
      std::cout << std::endl;
      for (i = 1; i <= m; i++) {
	for (j = 1; j <= n_b2; j++)
	  std::cout << std::scientific << std::setprecision(3)
		    << std::setw(11) << A[i][j];
	std::cout << std::scientific << std::setprecision(3)
		  << std::setw(11) << b[i] << std::endl;
      }
   
      std::cout << std::endl;
      std::cout << std::fixed << std::setprecision(5)
	   << "dalpha1_x = " << std::setw(8) << dalpha1[X_]
	   << ", dalpha1_y = " << std::setw(8) << dalpha1[Y_] << std::endl;
    }
  }

  if (prt) {
    quad_out.open("fit_alpha.dat", std::ios::out);
    quad_out << std::endl;
    quad_out << "n = 1:" << std::endl;
    for (i = 1; i <= n_b2; i++)
      for (j = 1; j <= get_n_Kids(abs(b2s[i-1])); j++)
	if (b2s[i-1] > 0)
	  quad_out << std::fixed << std::setprecision(7) 
		   << std::setw(6) << get_Name(b2s[i-1])
		   << "(" << j << ") = "
		   << std::setw(11) << get_bn(b2s[i-1], j, Quad)
		   << std::setw(2) << Quad << std::endl;
	else {
	  k = get_loc(abs(b2s[i-1]), j) - 1;
	  quad_out << std::fixed << std::setprecision(7) 
		   << std::setw(6) << get_Name(elem[k-1].Fnum)
		   << "(" << j << ") = "
		   << std::setw(11) << get_L(elem[k-1].Fnum, j)
		   << std::setw(3) << -Quad << std::endl;
	  quad_out << std::fixed << std::setprecision(7) 
		   << std::setw(6) << get_Name(elem[k+1].Fnum)
		   << "(" << j << ") = "
		   << std::setw(11) << get_L(elem[k+1].Fnum, j)
		   << std::setw(3) << -Quad <<  std::endl;
	}
    quad_out.close();
  }

  free_dvector(b, 1, m);
  free_dvector(b2_lim, 1, n_b2); free_dvector(b2, 1, n_b2);
  free_dvector(db2, 1, n_b2);
  free_dmatrix(A, 1, m, 1, n_b2);
}


void fit_beta(const double alpha0_x, const double beta0_x,
	      const double alpha0_y, const double beta0_y,
	      const double beta1_x, const double beta1_y,
	      const long int k1, const long int k2,
	      const int n_b2, const int b2s[],
	      const double eps, const bool prt)
{
  int           i, j, k;
  double        **A, *b, *b2_lim, *b2, *db2, step;
  double        dalpha1[2], dbeta1[2];
  tps           ab[2], nu[2];
  std::ofstream quad_out;

  const bool   debug = true;
  const int    m     = 2;
  const double s_cut = 1e-5, step0 = 0.7;

  b = dvector(1, m); b2_lim = dvector(1, n_b2);
  b2 = dvector(1, n_b2); db2 = dvector(1, n_b2);
  A = dmatrix(1, m, 1, n_b2);

  for (i = 1; i <= n_b2; i++) {
    if (b2s[i-1] > 0) {
      b2_lim[i] = b2_max; b2[i] = get_bn(b2s[i-1], 1, Quad);
    } else {
      b2_lim[i] = ds_max/scl_ds; b2[i] = get_L(abs(b2s[i-1]), 1);
    }
  }

  std::cout << std::endl;
  std::cout << std::fixed << std::setprecision(5)
       << "fit_beta: " << k1 << " - " << k2
       << ", alpha0_x = " << std::setw(8) << alpha0_x
       << ", alpha0_y = " << std::setw(8) << alpha0_y
       << ", beta0_x = " << std::setw(8) << beta0_x
       << ", beta0_y = " << std::setw(8) << beta0_y << std::endl;
  std::cout << "beta1_x = " << std::setw(8) << beta1_x
       << ", beta1_y = " << std::setw(8) << beta1_y << std::endl;

  danot_(3);

  get_A1(alpha0_x, beta0_x, alpha0_y, beta0_y); get_ab(ab, nu, k1+1, k2);

  dalpha1[X_] = -h_ijklm(ab[X_], 0, 1, 0, 0, 0);
  dalpha1[Y_] = -h_ijklm(ab[Y_], 0, 0, 0, 1, 0);
  dbeta1[X_]  =  h_ijklm(ab[X_], 1, 0, 0, 0, 0) - beta1_x;
  dbeta1[Y_]  =  h_ijklm(ab[Y_], 0, 0, 1, 0, 0) - beta1_y;

   while ((fabs(dbeta1[X_]) > eps) || (fabs(dbeta1[Y_]) > eps)) {
    step = step0;
    for (i = 1; i <= n_b2; i++) {
      for (j = 1; j <= get_n_Kids(abs(b2s[i-1])); j++)
	if (b2s[i-1] > 0)
	  set_bn_par(b2s[i-1], j, Quad, 7);
	else
	  set_s_par(abs(b2s[i-1]), j, 7);

      get_A1(alpha0_x, beta0_x, alpha0_y, beta0_y); get_ab(ab, nu, k1+1, k2);

      A[1][i] = -h_ijklm_p(ab[X_], 1, 0, 0, 0, 0, 7);
      A[2][i] = -h_ijklm_p(ab[Y_], 0, 0, 1, 0, 0, 7);

      if (b2s[i-1] < 0)
	for (j = 1; j <= m; j++)
	  A[j][i] *= scl_ds;

      for (j = 1; j <= get_n_Kids(abs(b2s[i-1])); j++)
	if (b2s[i-1] > 0)
	  clr_bn_par(b2s[i-1], j, Quad);
	else
	  clr_s_par(abs(b2s[i-1]), j);
    }

    b[1] = dbeta1[X_]; b[2] = dbeta1[Y_];

    SVD_lim(m, n_b2, A, b, b2_lim, s_cut, b2, db2);

    for (i = 1; i <= n_b2; i++) {
      set_dbn_s(b2s[i-1], Quad, step*db2[i]);
      b2[i] = get_bn_s(b2s[i-1], 1, Quad);
    }

    get_A1(alpha0_x, beta0_x, alpha0_y, beta0_y); get_ab(ab, nu, k1+1, k2);

    dalpha1[X_] = -h_ijklm(ab[X_], 0, 1, 0, 0, 0);
    dalpha1[Y_] = -h_ijklm(ab[Y_], 0, 0, 0, 1, 0);
    dbeta1[X_]  =  h_ijklm(ab[X_], 1, 0, 0, 0, 0) - beta1_x;
    dbeta1[Y_]  =  h_ijklm(ab[Y_], 0, 0, 1, 0, 0) - beta1_y;

    if (debug) {
      std::cout << std::endl;
      std::cout << " Ax = b:" << std::endl;
      std::cout << std::endl;
      for (i = 1; i <= m; i++) {
	for (j = 1; j <= n_b2; j++)
	  std::cout << std::scientific << std::setprecision(3)
		    << std::setw(11) << A[i][j];
	std::cout << std::scientific << std::setprecision(3)
		  << std::setw(11) << b[i] << std::endl;
      }
   
      std::cout << std::endl;
      std::cout << std::fixed << std::setprecision(5)
	   << "dalpha1_x = " << std::setw(8) << dalpha1[X_]
	   << ", dalpha1_y = " << std::setw(8) << dalpha1[Y_]
	   << ", dbeta1_x = " << std::setw(8) << dbeta1[X_]
	   << ", dbeta1_y = " << std::setw(8) << dbeta1[Y_] << std::endl;
    }
  }

  if (prt) {
    quad_out.open("fit_beta.dat", std::ios::out);
    quad_out << std::endl;
    quad_out << "n = 1:" << std::endl;
    for (i = 1; i <= n_b2; i++)
      for (j = 1; j <= get_n_Kids(abs(b2s[i-1])); j++)
	if (b2s[i-1] > 0)
	  quad_out << std::fixed << std::setprecision(7) 
		   << std::setw(6) << get_Name(b2s[i-1]) << "(" << j << ") = "
		   << std::setw(11) << get_bn(b2s[i-1], j, Quad)
		   << std::setw(2) << Quad << std::endl;
	else {
	  k = get_loc(abs(b2s[i-1]), j) - 1;
	  quad_out << std::fixed << std::setprecision(7) 
		   << std::setw(6) << get_Name(elem[k-1].Fnum)
		   << "(" << j << ") = "
		   << std::setw(11) << get_L(elem[k-1].Fnum, j)
		   << std::setw(3) << -Quad << std::endl;
	  quad_out << std::fixed << std::setprecision(7) 
		   << std::setw(6) << get_Name(elem[k+1].Fnum)
		   << "(" << j << ") = "
		   << std::setw(11) << get_L(elem[k+1].Fnum, j)
		   << std::setw(3) << -Quad << std::endl;
	}
    quad_out.close();
  }

  free_dvector(b, 1, m);
  free_dvector(b2_lim, 1, n_b2); free_dvector(b2, 1, n_b2);
  free_dvector(db2, 1, n_b2);
  free_dmatrix(A, 1, m, 1, n_b2);
}


void fit_alpha_beta(const double alpha0_x, const double beta0_x,
		    const double alpha0_y, const double beta0_y,
		    const double beta1_x, const double beta1_y,
		    const long int k1, const long int k2,
		    const int n_b2, const int b2s[],
		    const double eps, const bool prt)
{
  int           i, j, k;
  double        **A, *b, *b2_lim, *b2, *db2, step;
  double        dalpha1[2], dbeta1[2];
  tps           ab[2], nu[2];
  std::ofstream quad_out;

  const bool   debug = true;
  const int    m     = 4;
  const double s_cut = 1e-5, step0 = 0.5, scl_beta = scl_beta_mp*1e-3;

  b = dvector(1, m); b2_lim = dvector(1, n_b2);
  b2 = dvector(1, n_b2); db2 = dvector(1, n_b2);
  A = dmatrix(1, m, 1, n_b2);

  for (i = 1; i <= n_b2; i++) {
    if (b2s[i-1] > 0) {
      b2_lim[i] = b2_max; b2[i] = get_bn(b2s[i-1], 1, Quad);
    } else {
      b2_lim[i] = ds_max; b2[i] = get_L(abs(b2s[i-1]), 1);
    }
  }

  std::cout << std::endl;
  std::cout << std::fixed << std::setprecision(5)
       << "fit_alpha: " << k1 << " - " << k2
       << ", alpha0_x = " << std::setw(8) << alpha0_x
       << ", alpha0_y = " << std::setw(8) << alpha0_y
       << ", beta0_x = " << std::setw(8) << beta0_x
       << ", beta0_y = " << std::setw(8) << beta0_y << std::endl;
  std::cout << "beta1_x = " << std::setw(8) << beta1_x
       << ", beta1_y = " << std::setw(8) << beta1_y << std::endl;

  danot_(3);

  get_A1(alpha0_x, beta0_x, alpha0_y, beta0_y); get_ab(ab, nu, k1+1, k2);

  dalpha1[X_] = -h_ijklm(ab[X_], 0, 1, 0, 0, 0);
  dalpha1[Y_] = -h_ijklm(ab[Y_], 0, 0, 0, 1, 0);
  dbeta1[X_]  =  h_ijklm(ab[X_], 1, 0, 0, 0, 0) - beta1_x;
  dbeta1[Y_]  =  h_ijklm(ab[Y_], 0, 0, 1, 0, 0) - beta1_y;

   while ((fabs(dalpha1[X_]) > eps) || (fabs(dalpha1[Y_]) > eps)) {
    step = step0;
    for (i = 1; i <= n_b2; i++) {
      for (j = 1; j <= get_n_Kids(abs(b2s[i-1])); j++)
	if (b2s[i-1] > 0)
	  set_bn_par(b2s[i-1], j, Quad, 7);
	else
	  set_s_par(abs(b2s[i-1]), j, 7);

      get_A1(alpha0_x, beta0_x, alpha0_y, beta0_y); get_ab(ab, nu, k1+1, k2);

      A[1][i] =  scl_alpha_mp*h_ijklm_p(ab[X_], 0, 1, 0, 0, 0, 7);
      A[2][i] =  scl_alpha_mp*h_ijklm_p(ab[Y_], 0, 0, 0, 1, 0, 7);
      A[3][i] = -scl_beta*h_ijklm_p(ab[X_], 1, 0, 0, 0, 0, 7);
      A[4][i] = -scl_beta*h_ijklm_p(ab[Y_], 0, 0, 1, 0, 0, 7);

      if (b2s[i-1] < 0)
	for (j = 1; j <= m; j++)
	  A[j][i] *= scl_ds;

      for (j = 1; j <= get_n_Kids(abs(b2s[i-1])); j++)
	if (b2s[i-1] > 0)
	  clr_bn_par(b2s[i-1], j, Quad);
	else
	  clr_s_par(abs(b2s[i-1]), j);
    }

    b[1] = scl_alpha_mp*dalpha1[X_]; b[2] = scl_alpha_mp*dalpha1[Y_];
    b[3] = scl_beta*dbeta1[X_];   b[4] = scl_beta*dbeta1[Y_];

    SVD_lim(m, n_b2, A, b, b2_lim, s_cut, b2, db2);

    for (i = 1; i <= n_b2; i++) {
      set_dbn_s(b2s[i-1], Quad, step*db2[i]);
      b2[i] = get_bn_s(b2s[i-1], 1, Quad);
    }

    get_A1(alpha0_x, beta0_x, alpha0_y, beta0_y); get_ab(ab, nu, k1+1, k2);

    dalpha1[X_] = -h_ijklm(ab[X_], 0, 1, 0, 0, 0);
    dalpha1[Y_] = -h_ijklm(ab[Y_], 0, 0, 0, 1, 0);
    dbeta1[X_]  =  h_ijklm(ab[X_], 1, 0, 0, 0, 0) - beta1_x;
    dbeta1[Y_]  =  h_ijklm(ab[Y_], 0, 0, 1, 0, 0) - beta1_y;

    if (debug) {
      std::cout << std::endl;
      std::cout << " Ax = b:" << std::endl;
      std::cout << std::endl;
      for (i = 1; i <= m; i++) {
	for (j = 1; j <= n_b2; j++)
	  std::cout << std::scientific << std::setprecision(3)
		    << std::setw(11) << A[i][j];
	std::cout << std::scientific << std::setprecision(3)
		  << std::setw(11) << b[i] << std::endl;
      }
   
      std::cout << std::endl;
      std::cout << std::fixed <<std:: setprecision(5)
	   << "dalpha1_x = " << std::setw(8) << dalpha1[X_]
	   << ", dalpha1_y = " << std::setw(8) << dalpha1[Y_]
	   << ", dbeta1_x = " << std::setw(8) << dbeta1[X_]
	   << ", dbeta1_y = " << std::setw(8) << dbeta1[Y_] << std::endl;
    }
  }

  if (prt) {
    quad_out.open("fit_alpha.dat", std::ios::out);
    quad_out << std::endl;
    quad_out << "n = 1:" << std::endl;
    for (i = 1; i <= n_b2; i++)
      for (j = 1; j <= get_n_Kids(abs(b2s[i-1])); j++)
	if (b2s[i-1] > 0)
	  quad_out << std::fixed << std::setprecision(7) 
		   << std::setw(6) << get_Name(b2s[i-1]) << "(" << j << ") = "
		   << std::setw(11) << get_bn(b2s[i-1], j, Quad)
		   << std::setw(2) << Quad << std::endl;
	else {
	  k = get_loc(abs(b2s[i-1]), j) - 1;
	  quad_out << std::fixed << std::setprecision(7) 
		   << std::setw(6) << get_Name(elem[k-1].Fnum)
		   << "(" << j << ") = "
		   << std::setw(11) << get_L(elem[k-1].Fnum, j)
		   << std::setw(3) << -Quad << std::endl;
	  quad_out << std::fixed << std::setprecision(7) 
		   << std::setw(6) << get_Name(elem[k+1].Fnum)
		   << "(" << j << ") = "
		   << std::setw(11) << get_L(elem[k+1].Fnum, j)
		   << std::setw(3) << -Quad << std::endl;
	}
    quad_out.close();
  }

  free_dvector(b, 1, m);
  free_dvector(b2_lim, 1, n_b2); free_dvector(b2, 1, n_b2);
  free_dvector(db2, 1, n_b2);
  free_dmatrix(A, 1, m, 1, n_b2);
}


void fit_tune(const double nu_x, const double nu_y,
	      const int n_b2, const int b2s[],
	      const double eps, const bool prt)
{
  // Periodic solution: [nu_x, nu_y, beta_x, beta_y]

  int           i, j;
  double        **A, *b, *b2_lim, *b2, *db2, step;
  double        nu_fract[2], dnu[2];
  ss_vect<tps>  nus, dnus;
  std::ofstream quad_out;

  const bool    debug = true;
  const int     m     = 2;
  const double  s_cut = 1e-7, step0 = 0.7;

  b = dvector(1, m); b2_lim = dvector(1, n_b2);
  b2 = dvector(1, n_b2); db2 = dvector(1, n_b2);
  A = dmatrix(1, m, 1, n_b2);

  nu_fract[X_] = fract(nu_x); nu_fract[Y_] = fract(nu_y);
  std::cout << std::endl;
  std::cout << std::fixed << std::setprecision(5)
       << "fit_tune: nu_x = " << nu_fract[X_] << ", nu_y = " << nu_fract[Y_]
       << std::endl;

  for (i = 1; i <= n_b2; i++) {
    if (b2s[i-1] > 0) {
      b2_lim[i] = b2_max; b2[i] = get_bn(b2s[i-1], 1, Quad);
    } else {
      b2_lim[i] = ds_max; b2[i] = get_L(abs(b2s[i-1]), 1);
    }
  }

  danot_(3);

  get_Map(); K = MapNorm(Map, g, A1, A0, Map_res, 1); nus = dHdJ(K);

  dnu[X_] = nus[3].cst() - nu_fract[X_];
  dnu[Y_] = nus[4].cst() - nu_fract[Y_];
  
  while ((fabs(dnu[X_]) > eps) || (fabs(dnu[Y_]) > eps)) {
    step = step0;
    for (i = 1; i <= n_b2; i++) {
      for (j = 1; j <= get_n_Kids(abs(b2s[i-1])); j++)
	if (b2s[i-1] > 0)
	  set_bn_par(b2s[i-1], j, Quad, 7);
	else
	  set_s_par(abs(b2s[i-1]), j, 7);

      get_Map(); K = MapNorm(Map, g, A1, A0, Map_res, 1); nus = dHdJ(K);

      A[1][i] = -h_ijklm_p(nus[3], 0, 0, 0, 0, 0, 7);
      A[2][i] = -h_ijklm_p(nus[4], 0, 0, 0, 0, 0, 7);

      if (b2s[i-1] < 0)
	for (j = 1; j <= m; j++)
	  A[j][i] *= scl_ds;

      for (j = 1; j <= get_n_Kids(abs(b2s[i-1])); j++)
	if (b2s[i-1] > 0)
	  clr_bn_par(b2s[i-1], j, Quad);
	else
	  clr_s_par(abs(b2s[i-1]), j);
    }

    b[1] = dnu[X_]; b[2] = dnu[Y_];

    SVD_lim(m, n_b2, A, b, b2_lim, s_cut, b2, db2);

    for (i = 1; i <= n_b2; i++) {
      set_dbn_s(b2s[i-1], Quad, step*db2[i]);
      b2[i] = get_bn_s(b2s[i-1], 1, Quad);
    }

    get_Map(); K = MapNorm(Map, g, A1, A0, Map_res, 1);

    while (!stable) {
      // roll back
      for (i = 1; i <= n_b2; i++) {
	set_dbn_s(b2s[i-1], Quad, -step*db2[i]);
	b2[i] = get_bn_s(b2s[i-1], 1, Quad);
      }

      step /= 2.0;
      std::cout << std::endl;
      std::cout << std::scientific << std::setprecision(3)
		<< "step = " << step << std::endl;

      for (i = 1; i <= n_b2; i++) {
	set_dbn_s(b2s[i-1], Quad, step*db2[i]);
	b2[i] = get_bn_s(b2s[i-1], 1, Quad);
      }
	
      get_Map(); K = MapNorm(Map, g, A1, A0, Map_res, 1);
    }
      
    nus = dHdJ(K);

    dnu[X_] = nus[3].cst() - nu_fract[X_];
    dnu[Y_] = nus[4].cst() - nu_fract[Y_];
    
    if (debug) {
      std::cout << std::endl;
      std::cout << " Ax = b:" << std::endl;
      std::cout << std::endl;
      for (i = 1; i <= m; i++) {
	for (j = 1; j <= n_b2; j++)
	  std::cout << std::scientific << std::setprecision(3)
		    << std::setw(11) << A[i][j];
	std::cout << std::scientific << std::setprecision(3)
		  << std::setw(11) << b[i] << std::endl;
      }
	
      std::cout << std::endl;
      std::cout << std::fixed << std::setprecision(5)
	   << "dnu_x = " << dnu[X_] << ", dnu_y = " << dnu[Y_]
	   << std::endl;
    }
   
    std::cout << std::endl;
    std::cout << std::fixed << std::setprecision(5)
	      << "nu_x = " << nus[3].cst() << ", nu_y = " << nus[4].cst()
	      << std::endl;
  }

  if (prt) {
    quad_out.open("fit_tune.dat", std::ios::out);
    quad_out << std::endl;
    quad_out << "n = 1:" << std::endl;
    for (i = 1; i <= n_b2; i++)
      for (j = 1; j <= get_n_Kids(abs(b2s[i-1])); j++)
	if (b2s[i-1] > 0)
	  quad_out << std::fixed << std::setprecision(7) 
		   << std::setw(6) << get_Name(b2s[i-1]) << "(" << j << ") = "
		   << std::setw(11) << get_bnL(b2s[i-1], j, Quad)
		   << std::setw(2) << Quad << std::endl;
	else {
	  quad_out << std::fixed << std::setprecision(7) 
		   << std::setw(6) << get_Name(abs(b2s[i-1])) << "(" << j << ") = "
		   << std::setw(11) << get_L(abs(b2s[i-1]), j)
		   << std::setw(3) << -Quad << std::endl;
	}
    quad_out.close();
  }

  free_dvector(b, 1, m);
  free_dvector(b2_lim, 1, n_b2); free_dvector(b2, 1, n_b2);
  free_dvector(db2, 1, n_b2);
  free_dmatrix(A, 1, m, 1, n_b2);
}


void fit_tune(const double nu_x, const double nu_y,
	      const double beta1_x, const double beta1_y, const int k1,
	      const double beta2_x, const double beta2_y, const int k2,
	      const double beta3_x, const double beta3_y, const int k3,
	      const int n_b2, const int b2s[],
	      const double eps, const bool prt)
{
  // Periodic solution: [nu_x, nu_y, beta_x, beta_y]

  int           i, j;
  long int      k;
  double        **A, *b, *b2_lim, *b2, *db2, step;
  double        nu_fract[2], dnu[2];
  double        dalpha1[2], dbeta1[2], dbeta2[2], dbeta3[2], db2_max;
  tps           ab1[2], ab2[2], ab3[2], dnu1[2], dnu2[2], dnu3[2];
  ss_vect<tps>  nus, dnus, Mk;
  std::ofstream quad_out;

  const bool   debug  = true;
  const int    m      = 10;
  const double s_cut  = 1e-4, step0 = 0.7;
  const double scl_beta_mp_y = scl_beta_mp*1e-3;
  const double scl_beta_ss = scl_beta_mp*1e-3, scl_beta_ls = scl_beta_mp*1e-3;
  const double db2_tol = 1e-4;

  b = dvector(1, m); b2_lim = dvector(1, n_b2);
  b2 = dvector(1, n_b2); db2 = dvector(1, n_b2);
  A = dmatrix(1, m, 1, n_b2);

  nu_fract[X_] = fract(nu_x); nu_fract[Y_] = fract(nu_y);
  std::cout << std::endl;
  std::cout << std::fixed << std::setprecision(5)
       << "fit_tune: nu_x = " << nu_fract[X_] << ", nu_y = " << nu_fract[Y_]
       << std::endl;
  std::cout << std::fixed << std::setprecision(5)
       << "beta1_x = " << std::setw(8) << beta1_x
       << ", beta1_y = " << std::setw(8) << beta1_y << std::endl;
  std::cout << std::fixed << std::setprecision(5)
       << "beta2_x = " << std::setw(8) << beta2_x
       << ", beta2_y = " << std::setw(8) << beta2_y << std::endl;
  std::cout << std::fixed << std::setprecision(5)
       << "beta3_x = " << std::setw(8) << beta3_x
       << ", beta3_y = " << std::setw(8) << beta3_y << std::endl;

  for (i = 1; i <= n_b2; i++) {
    if (b2s[i-1] > 0) {
      b2_lim[i] = b2_max; b2[i] = get_bn(b2s[i-1], 1, Quad);
    } else {
      b2_lim[i] = ds_max;
      k = get_loc(abs(b2s[i-1]), 1) - 1; b2[i] = get_L(elem[k-1].Fnum, 1);
    }
  }

  danot_(3);

  get_Map(); K = MapNorm(Map, g, A1, A0, Map_res, 1); nus = dHdJ(K);
  get_ab(ab1, dnu1, k1); get_ab(ab2, dnu2, k2); get_ab(ab3, dnu3, k3);

  dnu[X_]     =  nus[3].cst() - nu_fract[X_];
  dnu[Y_]     =  nus[4].cst() - nu_fract[Y_];
  dalpha1[X_] = -h_ijklm(ab1[X_], 0, 1, 0, 0, 0);
  dalpha1[Y_] = -h_ijklm(ab1[Y_], 0, 0, 0, 1, 0);
  dbeta1[X_]  =  h_ijklm(ab1[X_], 1, 0, 0, 0, 0) - beta1_x;
  dbeta1[Y_]  =  h_ijklm(ab1[Y_], 0, 0, 1, 0, 0) - beta1_y;
  dbeta2[X_]  =  h_ijklm(ab2[X_], 1, 0, 0, 0, 0) - beta2_x;
  dbeta2[Y_]  =  h_ijklm(ab2[Y_], 0, 0, 1, 0, 0) - beta2_y;
  dbeta3[X_]  =  h_ijklm(ab3[X_], 1, 0, 0, 0, 0) - beta3_x;
  dbeta3[Y_]  =  h_ijklm(ab3[Y_], 0, 0, 1, 0, 0) - beta3_y;

  db2_max = 1e30;

/*  while ((fabs(dnu[X_]) > eps) || (fabs(dnu[Y_]) > eps) ||
	 (scl_alpha_mp*fabs(dalpha1[X_]) > eps) ||
	 (scl_alpha_mp*fabs(dalpha1[Y_]) > eps) ||
	 (scl_beta_mp*fabs(dbeta1[X_]) > eps) ||
	 (scl_beta_mp_y*fabs(dbeta1[Y_]) > eps) ||
	 (scl_beta_ss*fabs(dbeta2[X_]) > eps) ||
	 (scl_beta_ss*fabs(dbeta2[Y_]) > eps) ||
	 (scl_beta_ls*fabs(dbeta3[X_]) > eps) ||
	 (scl_beta_ls*fabs(dbeta3[Y_]) > eps)) {*/
  while (((fabs(dnu[X_]) > eps) || (fabs(dnu[Y_]) > eps)) &&
	 (db2_max > db2_tol)) {
    step = step0;
    for (i = 1; i <= n_b2; i++) {
      for (j = 1; j <= get_n_Kids(abs(b2s[i-1])); j++)
	if (b2s[i-1] > 0)
	  set_bn_par(b2s[i-1], j, Quad, 7);
	else
	  set_s_par(abs(b2s[i-1]), j, 7);

      get_Map(); K = MapNorm(Map, g, A1, A0, Map_res, 1); nus = dHdJ(K);
      get_ab(ab1, dnu1, k1); get_ab(ab2, dnu2, k2); get_ab(ab3, dnu3, k3);

      A[1][i]  = -h_ijklm_p(nus[3], 0, 0, 0, 0, 0, 7);
      A[2][i]  = -h_ijklm_p(nus[4], 0, 0, 0, 0, 0, 7);
      A[3][i]  =  scl_alpha_mp*h_ijklm_p(ab1[X_], 0, 1, 0, 0, 0, 7);
      A[4][i]  =  scl_alpha_mp*h_ijklm_p(ab1[Y_], 0, 0, 0, 1, 0, 7);
      A[5][i]  = -scl_beta_mp*h_ijklm_p(ab1[X_], 1, 0, 0, 0, 0, 7);
      A[6][i]  = -scl_beta_mp_y*h_ijklm_p(ab1[Y_], 0, 0, 1, 0, 0, 7);
      A[7][i]  = -scl_beta_ss*h_ijklm_p(ab2[X_], 1, 0, 0, 0, 0, 7);
      A[8][i]  = -scl_beta_ss*h_ijklm_p(ab2[Y_], 0, 0, 1, 0, 0, 7);
      A[9][i]  = -scl_beta_ls*h_ijklm_p(ab3[X_], 1, 0, 0, 0, 0, 7);
      A[10][i] = -scl_beta_ls*h_ijklm_p(ab3[Y_], 0, 0, 1, 0, 0, 7);

      if (b2s[i-1] < 0)
	for (j = 1; j <= m; j++)
	  A[j][i] *= scl_ds;

      for (j = 1; j <= get_n_Kids(abs(b2s[i-1])); j++)
	if (b2s[i-1] > 0)
	  clr_bn_par(b2s[i-1], j, Quad);
	else
	  clr_s_par(abs(b2s[i-1]), j);
    }

    b[1] = dnu[X_];                  b[2] = dnu[Y_];
    b[3] = scl_alpha_mp*dalpha1[X_]; b[4] = scl_alpha_mp*dalpha1[Y_];
    b[5] = scl_beta_mp*dbeta1[X_];   b[6] = scl_beta_mp_y*dbeta1[Y_];
    b[7] = scl_beta_ss*dbeta2[X_];   b[8] = scl_beta_ss*dbeta2[Y_];
    b[9] = scl_beta_ls*dbeta3[X_];   b[10] = scl_beta_ls*dbeta3[Y_];

    SVD_lim(m, n_b2, A, b, b2_lim, s_cut, b2, db2);

    db2_max = 0.0;
    for (i = 1; i <= n_b2; i++) {
      set_dbn_s(b2s[i-1], Quad, step*db2[i]);
      b2[i] = get_bn_s(b2s[i-1], 1, Quad);
      db2_max = max(step*db2[i], db2_max);
    }

    get_Map(); K = MapNorm(Map, g, A1, A0, Map_res, 1);

    while (!stable) {
      // roll back
      for (i = 1; i <= n_b2; i++) {
	set_dbn_s(b2s[i-1], Quad, -step*db2[i]);
	b2[i] = get_bn_s(b2s[i-1], 1, Quad);
      }

      step /= 2.0;
      std::cout << std::endl;
      std::cout << std::scientific << std::setprecision(3)
		<< "step = " << step << std::endl;

      for (i = 1; i <= n_b2; i++) {
	set_dbn_s(b2s[i-1], Quad, step*db2[i]);
	b2[i] = get_bn_s(b2s[i-1], 1, Quad);
      }
	
      get_Map(); K = MapNorm(Map, g, A1, A0, Map_res, 1);
    }
      
    nus = dHdJ(K);
    get_ab(ab1, dnu1, k1); get_ab(ab2, dnu2, k2); get_ab(ab3, dnu3, k3);

    dnu[X_]     =  nus[3].cst() - nu_fract[X_];
    dnu[Y_]     =  nus[4].cst() - nu_fract[Y_];
    dalpha1[X_] = -h_ijklm(ab1[X_], 0, 1, 0, 0, 0);
    dalpha1[Y_] = -h_ijklm(ab1[Y_], 0, 0, 0, 1, 0);
    dbeta1[X_]  =  h_ijklm(ab1[X_], 1, 0, 0, 0, 0) - beta1_x;
    dbeta1[Y_]  =  h_ijklm(ab1[Y_], 0, 0, 1, 0, 0) - beta1_y;
    dbeta2[X_]  =  h_ijklm(ab2[X_], 1, 0, 0, 0, 0) - beta2_x;
    dbeta2[Y_]  =  h_ijklm(ab2[Y_], 0, 0, 1, 0, 0) - beta2_y;
    dbeta3[X_]  =  h_ijklm(ab3[X_], 1, 0, 0, 0, 0) - beta3_x;
    dbeta3[Y_]  =  h_ijklm(ab3[Y_], 0, 0, 1, 0, 0) - beta3_y;
    
    if (debug) {
      std::cout << std::endl;
      std::cout << " Ax = b:" << std::endl;
      std::cout << std::endl;
      for (j = 1; j <= n_b2; j++)
	if (j == 1)
	  std::cout << std::setw(8) << j;
	else
	  std::cout << std::setw(11) << j;
      std::cout << std::endl;
      for (i = 1; i <= m; i++) {
	std::cout << std::setw(2) << i;
	for (j = 1; j <= n_b2; j++)
	  std::cout << std::scientific << std::setprecision(3)
		    << std::setw(11) << A[i][j];
	std::cout << std::scientific << std::setprecision(3)
		  << std::setw(11) << b[i] << std::endl;
      }
	
      std::cout << std::endl;
      std::cout << std::fixed << std::setprecision(5)
	   << "dnu_x = " << dnu[X_] << ", dnu_y = " << dnu[Y_] << std::endl;
      std::cout << std::fixed << std::setprecision(5)
	   << "dbeta1_x = " << std::setw(8) << dbeta1[X_]
	   << ", dbeta1_y = " << std::setw(8) << dbeta1[Y_]
	   << ", dalpha1_x = " << std::setw(8) << dalpha1[X_]
	   << ", dalpha1_y = " << std::setw(8) << dalpha1[Y_] << std::endl;
      std::cout << std::fixed << std::setprecision(5)
	   << "dbeta2_x = " << std::setw(8) << dbeta2[X_]
	   << ", dbeta2_y = " << std::setw(8) << dbeta2[Y_] << std::endl;
      std::cout << std::fixed << std::setprecision(5)
	   << "dbeta3_x = " << std::setw(8) << dbeta3[X_]
	   << ", dbeta3_y = " << std::setw(8) << dbeta3[Y_] << std::endl;
    }
   
    std::cout << std::endl;
    std::cout << std::fixed << std::setprecision(5)
	 << "nu_x = " << nus[3].cst() << ", nu_y = " << nus[4].cst()
	      << std::endl;
  }

  if (prt) {
    quad_out.open("fit_tune.dat", std::ios::out);
    quad_out << std::endl;
    quad_out << "n = 1:" << std::endl;
    for (i = 1; i <= n_b2; i++)
      for (j = 1; j <= get_n_Kids(abs(b2s[i-1])); j++)
	if (b2s[i-1] > 0)
	  quad_out << std::fixed << std::setprecision(7) 
		   << std::setw(6) << get_Name(b2s[i-1])
		   << "(" << j << ") = "
		   << std::setw(11) << get_bnL(b2s[i-1], j, Quad)
		   << std::setw(2) << Quad << std::endl;
	else {
	  k = get_loc(abs(b2s[i-1]), j) - 1;
	  quad_out << std::fixed << std::setprecision(7) 
		   << std::setw(6) << get_Name(elem[k-1].Fnum)
		   << "(" << j << ") = "
		   << std::setw(11) << get_L(elem[k-1].Fnum, j)
		   << std::setw(3) << -Quad << std::endl;
	  quad_out << std::fixed << std::setprecision(7) 
		   << std::setw(6) << get_Name(elem[k+1].Fnum)
		   << "(" << j << ") = "
		   << std::setw(11) << get_L(elem[k+1].Fnum, j)
		   << std::setw(3) << -Quad << std::endl;
	}
    quad_out.close();
  }

  free_dvector(b, 1, m);
  free_dvector(b2_lim, 1, n_b2); free_dvector(b2, 1, n_b2);
  free_dvector(db2, 1, n_b2);
  free_dmatrix(A, 1, m, 1, n_b2);
}


void fit_tune_SLS(const double nu_x, const double nu_y,
		  const double beta1_x, const double beta1_y, const int k1,
		  const double beta2_x, const double beta2_y, const int k2,
		  const double beta3_x, const double beta3_y, const int k3,
		  const double beta4_x, const double beta4_y, const int k4,
		  const int n_b2, const int b2s[],
		  const double eps, const bool prt)
{
  // Periodic solution: [nu_x, nu_y, beta_x, beta_y]

  int           i, j;
  long int      k;
  double        **A, *b, *b2_lim, *b2, *db2, step;
  double        nu_fract[2], dnu[2];
  double        dalpha1[2], dbeta1[2], dalpha2[2], dbeta2[2], dbeta3[2];
  double        dbeta4[2];
  tps           ab1[2], ab2[2], ab3[2], ab4[2];
  tps           dnu1[2], dnu2[2], dnu3[2], dnu4[2];
  ss_vect<tps>  nus, dnus, Mk;
  std::ofstream quad_out;

  const bool   debug  = true;
  const int    m      = 14;
  const double s_cut  = 1e-5, step0 = 0.7;
  const double scl_beta_mp_y = scl_beta_mp*1e-3;
  const double scl_beta_ss = scl_beta_mp*1e-3, scl_beta_ls = scl_beta_mp*1e-3;
  const double scl_beta_ms = scl_beta_mp*1e-3;

  b = dvector(1, m); b2_lim = dvector(1, n_b2);
  b2 = dvector(1, n_b2); db2 = dvector(1, n_b2);
  A = dmatrix(1, m, 1, n_b2);

  nu_fract[X_] = fract(nu_x); nu_fract[Y_] = fract(nu_y);
  std::cout << std::endl;
  std::cout << std::fixed << std::setprecision(5)
       << "fit_tune_SLS: nu_x = " << nu_fract[X_]
       << ", nu_y = " << nu_fract[Y_] << std::endl;
  std::cout << std::fixed << std::setprecision(5)
       << "beta1_x = " << std::setw(8) << beta1_x
       << ", beta1_y = " << std::setw(8) << beta1_y << std::endl;
  std::cout << std::fixed << std::setprecision(5)
       << "beta2_x = " << std::setw(8) << beta2_x
       << ", beta2_y = " << std::setw(8) << beta2_y << std::endl;
  std::cout << std::fixed << std::setprecision(5)
       << "beta3_x = " << std::setw(8) << beta3_x
       << ", beta3_y = " << std::setw(8) << beta3_y << std::endl;
  std::cout << std::fixed << std::setprecision(5)
       << "beta4_x = " << std::setw(8) << beta4_x
       << ", beta4_y = " << std::setw(8) << beta4_y << std::endl;

  for (i = 1; i <= n_b2; i++) {
    if (b2s[i-1] > 0) {
      b2_lim[i] = b2_max; b2[i] = get_bn(b2s[i-1], 1, Quad);
    } else {
      b2_lim[i] = ds_max;
      k = get_loc(abs(b2s[i-1]), 1) - 1; b2[i] = get_L(elem[k-1].Fnum, 1);
    }
  }

  danot_(3);

  get_Map(); K = MapNorm(Map, g, A1, A0, Map_res, 1); nus = dHdJ(K);
  get_ab(ab1, dnu1, k1); get_ab(ab2, dnu2, k2); get_ab(ab3, dnu3, k3);
  get_ab(ab4, dnu4, k4);

  dnu[X_]     =  nus[3].cst() - nu_fract[X_];
  dnu[Y_]     =  nus[4].cst() - nu_fract[Y_];
  dalpha1[X_] = -h_ijklm(ab1[X_], 0, 1, 0, 0, 0);
  dalpha1[Y_] = -h_ijklm(ab1[Y_], 0, 0, 0, 1, 0);
  dbeta1[X_]  =  h_ijklm(ab1[X_], 1, 0, 0, 0, 0) - beta1_x;
  dbeta1[Y_]  =  h_ijklm(ab1[Y_], 0, 0, 1, 0, 0) - beta1_y;
  dalpha2[X_] = -h_ijklm(ab2[X_], 0, 1, 0, 0, 0);
  dalpha2[Y_] = -h_ijklm(ab2[Y_], 0, 0, 0, 1, 0);
  dbeta2[X_]  =  h_ijklm(ab2[X_], 1, 0, 0, 0, 0) - beta2_x;
  dbeta2[Y_]  =  h_ijklm(ab2[Y_], 0, 0, 1, 0, 0) - beta2_y;
  dbeta3[X_]  =  h_ijklm(ab3[X_], 1, 0, 0, 0, 0) - beta3_x;
  dbeta3[Y_]  =  h_ijklm(ab3[Y_], 0, 0, 1, 0, 0) - beta3_y;
  dbeta4[X_]  =  h_ijklm(ab4[X_], 1, 0, 0, 0, 0) - beta4_x;
  dbeta4[Y_]  =  h_ijklm(ab4[Y_], 0, 0, 1, 0, 0) - beta4_y;

  while ((fabs(dnu[X_]) > eps) || (fabs(dnu[Y_]) > eps) ||
	 (scl_alpha_mp*fabs(dalpha1[X_]) > eps) ||
	 (scl_alpha_mp*fabs(dalpha1[Y_]) > eps) ||
	 (scl_beta_mp*fabs(dbeta1[X_]) > eps) ||
	 (scl_beta_mp_y*fabs(dbeta1[Y_]) > eps) ||
	 (scl_alpha_mp*fabs(dalpha2[X_]) > eps) ||
	 (scl_alpha_mp*fabs(dalpha2[Y_]) > eps) ||
	 (scl_beta_ss*fabs(dbeta2[X_]) > eps) ||
	 (scl_beta_ss*fabs(dbeta2[Y_]) > eps) ||
	 (scl_beta_ls*fabs(dbeta3[X_]) > eps) ||
	 (scl_beta_ls*fabs(dbeta3[Y_]) > eps) ||
	 (scl_beta_ms*fabs(dbeta4[X_]) > eps) ||
	 (scl_beta_ms*fabs(dbeta4[Y_]) > eps)) {
    step = step0;
    for (i = 1; i <= n_b2; i++) {
      for (j = 1; j <= get_n_Kids(abs(b2s[i-1])); j++)
	if (b2s[i-1] > 0)
	  set_bn_par(b2s[i-1], j, Quad, 7);
	else
	  set_s_par(abs(b2s[i-1]), j, 7);

      get_Map(); K = MapNorm(Map, g, A1, A0, Map_res, 1); nus = dHdJ(K);
      get_ab(ab1, dnu1, k1); get_ab(ab2, dnu2, k2); get_ab(ab3, dnu3, k3);
      get_ab(ab4, dnu4, k4);

      A[1][i]  = -h_ijklm_p(nus[3], 0, 0, 0, 0, 0, 7);
      A[2][i]  = -h_ijklm_p(nus[4], 0, 0, 0, 0, 0, 7);
      A[3][i]  =  scl_alpha_mp*h_ijklm_p(ab1[X_], 0, 1, 0, 0, 0, 7);
      A[4][i]  =  scl_alpha_mp*h_ijklm_p(ab1[Y_], 0, 0, 0, 1, 0, 7);
      A[5][i]  = -scl_beta_mp*h_ijklm_p(ab1[X_], 1, 0, 0, 0, 0, 7);
      A[6][i]  = -scl_beta_mp_y*h_ijklm_p(ab1[Y_], 0, 0, 1, 0, 0, 7);
      A[7][i]  =  scl_alpha_mp*h_ijklm_p(ab2[X_], 0, 1, 0, 0, 0, 7);
      A[8][i]  =  scl_alpha_mp*h_ijklm_p(ab2[Y_], 0, 0, 0, 1, 0, 7);
      A[9][i]  = -scl_beta_ss*h_ijklm_p(ab2[X_], 1, 0, 0, 0, 0, 7);
      A[10][i] = -scl_beta_ss*h_ijklm_p(ab2[Y_], 0, 0, 1, 0, 0, 7);
      A[11][i] = -scl_beta_ls*h_ijklm_p(ab3[X_], 1, 0, 0, 0, 0, 7);
      A[12][i] = -scl_beta_ls*h_ijklm_p(ab3[Y_], 0, 0, 1, 0, 0, 7);
      A[13][i] = -scl_beta_ms*h_ijklm_p(ab4[X_], 1, 0, 0, 0, 0, 7);
      A[14][i] = -scl_beta_ms*h_ijklm_p(ab4[Y_], 0, 0, 1, 0, 0, 7);

      if (b2s[i-1] < 0)
	for (j = 1; j <= m; j++)
	  A[j][i] *= scl_ds;

      for (j = 1; j <= get_n_Kids(abs(b2s[i-1])); j++)
	if (b2s[i-1] > 0)
	  clr_bn_par(b2s[i-1], j, Quad);
	else
	  clr_s_par(abs(b2s[i-1]), j);
    }

    b[1]  = dnu[X_];                  b[2]  = dnu[Y_];
    b[3]  = scl_alpha_mp*dalpha1[X_]; b[4]  = scl_alpha_mp*dalpha1[Y_];
    b[5]  = scl_beta_mp*dbeta1[X_];   b[6]  = scl_beta_mp_y*dbeta1[Y_];
    b[7]  = scl_alpha_mp*dalpha2[X_]; b[8]  = scl_alpha_mp*dalpha2[Y_];
    b[9]  = scl_beta_ss*dbeta2[X_];   b[10] = scl_beta_ss*dbeta2[Y_];
    b[11] = scl_beta_ls*dbeta3[X_];   b[12] = scl_beta_ls*dbeta3[Y_];
    b[13] = scl_beta_ms*dbeta4[X_];   b[14] = scl_beta_ms*dbeta4[Y_];

    SVD_lim(m, n_b2, A, b, b2_lim, s_cut, b2, db2);

    for (i = 1; i <= n_b2; i++) {
      set_dbn_s(b2s[i-1], Quad, step*db2[i]);
      b2[i] = get_bn_s(b2s[i-1], 1, Quad);
    }

    get_Map(); K = MapNorm(Map, g, A1, A0, Map_res, 1);

    while (!stable) {
      // roll back
      for (i = 1; i <= n_b2; i++) {
	set_dbn_s(b2s[i-1], Quad, -step*db2[i]);
	b2[i] = get_bn_s(b2s[i-1], 1, Quad);
      }

      step /= 2.0;
      std::cout << std::endl;
      std::cout << std::scientific << std::setprecision(3)
		<< "step = " << step << std::endl;

      for (i = 1; i <= n_b2; i++) {
	set_dbn_s(b2s[i-1], Quad, step*db2[i]);
	b2[i] = get_bn_s(b2s[i-1], 1, Quad);
      }
	
      get_Map(); K = MapNorm(Map, g, A1, A0, Map_res, 1);
    }
      
    nus = dHdJ(K);
    get_ab(ab1, dnu1, k1); get_ab(ab2, dnu2, k2); get_ab(ab3, dnu3, k3);
    get_ab(ab4, dnu4, k4);

    dnu[X_]     =  nus[3].cst() - nu_fract[X_];
    dnu[Y_]     =  nus[4].cst() - nu_fract[Y_];
    dalpha1[X_] = -h_ijklm(ab1[X_], 0, 1, 0, 0, 0);
    dalpha1[Y_] = -h_ijklm(ab1[Y_], 0, 0, 0, 1, 0);
    dbeta1[X_]  =  h_ijklm(ab1[X_], 1, 0, 0, 0, 0) - beta1_x;
    dbeta1[Y_]  =  h_ijklm(ab1[Y_], 0, 0, 1, 0, 0) - beta1_y;
    dalpha2[X_] = -h_ijklm(ab2[X_], 0, 1, 0, 0, 0);
    dalpha2[Y_] = -h_ijklm(ab2[Y_], 0, 0, 0, 1, 0);
    dbeta2[X_]  =  h_ijklm(ab2[X_], 1, 0, 0, 0, 0) - beta2_x;
    dbeta2[Y_]  =  h_ijklm(ab2[Y_], 0, 0, 1, 0, 0) - beta2_y;
    dbeta3[X_]  =  h_ijklm(ab3[X_], 1, 0, 0, 0, 0) - beta3_x;
    dbeta3[Y_]  =  h_ijklm(ab3[Y_], 0, 0, 1, 0, 0) - beta3_y;
    dbeta4[X_]  =  h_ijklm(ab4[X_], 1, 0, 0, 0, 0) - beta4_x;
    dbeta4[Y_]  =  h_ijklm(ab4[Y_], 0, 0, 1, 0, 0) - beta4_y;
    
    if (debug) {
      std::cout << std::endl;
      std::cout << " Ax = b:" << std::endl;
      std::cout << std::endl;
      for (j = 1; j <= n_b2; j++)
	if (j == 1)
	  std::cout << std::setw(8) << j;
	else
	  std::cout << std::setw(11) << j;
      std::cout << std::endl;
      for (i = 1; i <= m; i++) {
	std::cout << std::setw(2) << i;
	for (j = 1; j <= n_b2; j++)
	  std::cout << std::scientific << std::setprecision(3)
		    << std::setw(11) << A[i][j];
	std::cout << std::scientific << std::setprecision(3)
		  << std::setw(11) << b[i] << std::endl;
      }
	
      std::cout << std::endl;
      std::cout << std::fixed <<std:: setprecision(5)
	   << "dnu_x = " << dnu[X_] << ", dnu_y = " << dnu[Y_] << std::endl;
      std::cout << std::fixed << std::setprecision(5)
	   << "dbeta1_x = " << std::setw(8) << dbeta1[X_]
	   << ", dbeta1_y = " << std::setw(8) << dbeta1[Y_]
	   << ", dalpha1_x = " << std::setw(8) << dalpha1[X_]
	   << ", dalpha1_y = " << std::setw(8) << dalpha1[Y_] << std::endl;
      std::cout << std::fixed << std::setprecision(5)
		<< "dbeta2_x = " << std::setw(8) << dbeta2[X_]
		<< ", dbeta2_y = " << std::setw(8) << dbeta2[Y_]
		<< ", dalpha2_x = " << std::setw(8) << dalpha2[X_]
		<< ", dalpha2_y = " << std::setw(8) << dalpha2[Y_] << std::endl;
      std::cout << std::fixed << std::setprecision(5)
	   << "dbeta3_x = " << std::setw(8) << dbeta3[X_]
	   << ", dbeta3_y = " << std::setw(8) << dbeta3[Y_] << std::endl;
      std::cout << std::fixed << std::setprecision(5)
	   << "dbeta4_x = " << std::setw(8) << dbeta4[X_]
	   << ", dbeta4_y = " << std::setw(8) << dbeta4[Y_] << std::endl;
    }
   
    std::cout << std::endl;
    std::cout << std::fixed << std::setprecision(5)
	      << "nu_x = " << nus[3].cst() << ", nu_y = " << nus[4].cst()
	      << std::endl;
  }

  if (prt) {
    quad_out.open("fit_tune.dat", std::ios::out);
    quad_out << std::endl;
    quad_out << "n = 1:" << std::endl;
    for (i = 1; i <= n_b2; i++)
      for (j = 1; j <= get_n_Kids(abs(b2s[i-1])); j++)
	if (b2s[i-1] > 0)
	  quad_out << std::fixed << std::setprecision(7) 
		   << std::setw(6) << get_Name(b2s[i-1]) << "(" << j << ") = "
		   << std::setw(11) << get_bnL(b2s[i-1], j, Quad)
		   << std::setw(2) << Quad << std::endl;
	else {
	  k = get_loc(abs(b2s[i-1]), j) - 1;
	  quad_out << std::fixed << std::setprecision(7) 
		   << std::setw(6) << get_Name(elem[k-1].Fnum)
		   << "(" << j << ") = "
		   << std::setw(11) << get_L(elem[k-1].Fnum, j)
		   << std::setw(3) << -Quad << std::endl;
	  quad_out << std::fixed << std::setprecision(7) 
		   << std::setw(6) << get_Name(elem[k+1].Fnum)
		   << "(" << j << ") = "
		   << std::setw(11) << get_L(elem[k+1].Fnum, j)
		   << std::setw(3) << -Quad << std::endl;
	}
    quad_out.close();
  }

  free_dvector(b, 1, m);
  free_dvector(b2_lim, 1, n_b2); free_dvector(b2, 1, n_b2);
  free_dvector(db2, 1, n_b2);
  free_dmatrix(A, 1, m, 1, n_b2);
}


double get_diff(const double alpha1[], const double beta1[],
		const double dnu12[], const long int locs[],
		tps ab2[], tps nu[], const bool prt)
{
  double dalpha[2], dbeta[2], dnu[2], diff;
  tps    ab1[2];


  get_A1(alpha1[X_], beta1[X_], alpha1[Y_], beta1[Y_]); get_ab(ab1, nu, 0);
  get_ab(ab2, nu, locs[0]+1, locs[1]);

  // alpha_2 = - alpha_1
  dalpha[X_] =
    - h_ijklm(ab2[X_], 0, 1, 0, 0, 0) - h_ijklm(ab1[X_], 0, 1, 0, 0, 0);
  dalpha[Y_] =
    - h_ijklm(ab2[Y_], 0, 0, 0, 1, 0) - h_ijklm(ab1[Y_], 0, 0, 0, 1, 0);
  dbeta[X_] =
    h_ijklm(ab2[X_], 1, 0, 0, 0, 0) - h_ijklm(ab1[X_], 1, 0, 0, 0, 0);
  dbeta[Y_] =
    h_ijklm(ab2[Y_], 0, 0, 1, 0, 0) - h_ijklm(ab1[Y_], 0, 0, 1, 0, 0);
  dnu[X_] = nu[X_].cst() - dnu12[X_];
  dnu[Y_] = nu[Y_].cst() - dnu12[Y_];

  if (prt) {
    std::cout << std::endl;
    std::cout << std::fixed << std::setprecision(5)
	 << "dnu_x = " << dnu[X_] << ", dnu_y = " << dnu[Y_] << std::endl;
    std::cout << std::fixed << std::setprecision(5)
	 << "dbeta_x = " << std::setw(8) << dalpha[X_]
	 << ", dbeta_y = " << std::setw(8) << dalpha[Y_]
	 << ", dalpha_x = " << std::setw(8) << dbeta[X_]
	 << ", dalpha_y = " << std::setw(8) << dbeta[Y_]
	 << std::endl;
  }

  diff = scl_alpha_mp*fabs(dalpha[X_]) + scl_alpha_mp*fabs(dalpha[Y_])
         + scl_beta_mp*fabs(dbeta[X_]) + scl_beta_mp*fabs(dbeta[Y_])
         + fabs(dnu[X_]) + fabs(dnu[Y_]);

  return diff;
}


void fit_femto_SLS(const double dnu_x, const double dnu_y,
		   const double alpha1_x, const double alpha1_y,
		   const double beta1_x, const double beta1_y, const int k1,
		   const int k2,
		   const int n_b2, const int b2s[],
		   const double eps, const bool prt)
{
  // Periodic solution: [nu_x, nu_y, beta_x, beta_y]

  int           i, j;
  long int      k, locs[2];
  double        **A, *b, *b2_lim, *b2, *db2, step;
  double        alpha1[2], beta1[2], dnu12[2], diff;
  tps           ab2[2], nu[2];
  std::ofstream quad_out;

  const bool   debug  = true;
  const int    m      = 6;
  const double s_cut  = 1e-5, step0 = 0.7;

  b = dvector(1, m); b2_lim = dvector(1, n_b2);
  b2 = dvector(1, n_b2); db2 = dvector(1, n_b2);
  A = dmatrix(1, m, 1, n_b2);

  std::cout << std::endl;
  std::cout << std::fixed << std::setprecision(5)
	    << "fit_femto_SLS: dnu_x = " << dnu_x << ", dnu_y = " << dnu_y
	    << std::endl;
  std::cout << std::fixed << std::setprecision(5)
	    << "beta1_x = " << std::setw(8) << beta1_x
	    << ", beta1_y = " << std::setw(8) << beta1_y << std::endl;

  for (i = 1; i <= n_b2; i++) {
    if (b2s[i-1] > 0) {
      b2_lim[i] = b2_max; b2[i] = get_bn(b2s[i-1], 1, Quad);
    } else {
      b2_lim[i] = ds_max;
      k = get_loc(abs(b2s[i-1]), 1) - 1; b2[i] = get_L(elem[k-1].Fnum, 1);
    }
  }

  danot_(3);

  alpha1[X_] = alpha1_x; alpha1[Y_] = alpha1_y;
  beta1[X_] = beta1_x; beta1[Y_] = beta1_y;
  dnu12[X_] = dnu_x; dnu12[Y_] = dnu_y; locs[0] = k1, locs[1] = k2;

  diff = get_diff(alpha1, beta1, dnu12, locs, ab2, nu, true);

  while (diff > eps) {
    step = step0;
    for (i = 1; i <= n_b2; i++) {
      for (j = 1; j <= get_n_Kids(abs(b2s[i-1])); j++)
	if (b2s[i-1] > 0)
	  set_bn_par(b2s[i-1], j, Quad, 7);
	else
	  set_s_par(abs(b2s[i-1]), j, 7);

      diff = get_diff(alpha1, beta1, dnu12, locs, ab2, nu, false);

      A[1][i]  = -h_ijklm_p(nu[X_], 0, 0, 0, 0, 0, 7);
      A[2][i]  = -h_ijklm_p(nu[Y_], 0, 0, 0, 0, 0, 7);
      A[3][i]  =  scl_alpha_mp*h_ijklm_p(ab2[X_], 0, 1, 0, 0, 0, 7);
      A[4][i]  =  scl_alpha_mp*h_ijklm_p(ab2[Y_], 0, 0, 0, 1, 0, 7);
      A[5][i]  = -scl_beta_mp*h_ijklm_p(ab2[X_], 1, 0, 0, 0, 0, 7);
      A[6][i]  = -scl_beta_mp*h_ijklm_p(ab2[Y_], 0, 0, 1, 0, 0, 7);

      if (b2s[i-1] < 0)
	for (j = 1; j <= m; j++)
	  A[j][i] *= scl_ds;

      for (j = 1; j <= get_n_Kids(abs(b2s[i-1])); j++)
	if (b2s[i-1] > 0)
	  clr_bn_par(b2s[i-1], j, Quad);
	else
	  clr_s_par(abs(b2s[i-1]), j);
    }

    b[1]  =  nu[X_].cst() - dnu_x;
    b[2]  =  nu[Y_].cst() - dnu_y;
    // alpha_2 = - alpha_1
    b[3]  = scl_alpha_mp*(-h_ijklm(ab2[X_], 0, 1, 0, 0, 0)+alpha1[X_]);
    b[4]  = scl_alpha_mp*(-h_ijklm(ab2[Y_], 0, 0, 0, 1, 0)+alpha1[Y_]);
    b[5]  = scl_beta_mp*(h_ijklm(ab2[X_], 1, 0, 0, 0, 0)-beta1[X_]);
    b[6]  = scl_beta_mp*(h_ijklm(ab2[Y_], 0, 0, 1, 0, 0)-beta1[Y_]);

    SVD_lim(m, n_b2, A, b, b2_lim, s_cut, b2, db2);

    for (i = 1; i <= n_b2; i++) {
      set_dbn_s(b2s[i-1], Quad, step*db2[i]);
      b2[i] = get_bn_s(b2s[i-1], 1, Quad);
    }

    diff = get_diff(alpha1, beta1, dnu12, locs, ab2, nu, true);

    if (debug) {
      std::cout << std::endl;
      std::cout << " Ax = b:" << std::endl;
      std::cout << std::endl;
      for (j = 1; j <= n_b2; j++)
	if (j == 1)
	  std::cout << std::setw(8) << j;
	else
	  std::cout << std::setw(11) << j;
      std::cout << std::endl;
      for (i = 1; i <= m; i++) {
	std::cout << std::setw(2) << i;
	for (j = 1; j <= n_b2; j++)
	  std::cout << std::scientific << std::setprecision(3)
		    << std::setw(11) << A[i][j];
	std::cout << std::scientific << std::setprecision(3)
		  << std::setw(11) << b[i] << std::endl;
      }
    }
  }

  if (prt) {
    quad_out.open("fit_femto.dat", std::ios::out);
    quad_out << std::endl;
    quad_out << "n = 1:" << std::endl;
    for (i = 1; i <= n_b2; i++)
      for (j = 1; j <= get_n_Kids(abs(b2s[i-1])); j++)
	if (b2s[i-1] > 0)
	  quad_out << std::fixed << std::setprecision(7) 
		   << std::setw(6) << get_Name(b2s[i-1]) << "(" << j << ") = "
		   << std::setw(11) << get_bnL(b2s[i-1], j, Quad)
		   << std::setw(2) << Quad << std::endl;
	else {
	  k = get_loc(abs(b2s[i-1]), j) - 1;
	  quad_out << std::fixed << std::setprecision(7) 
		   << std::setw(6) << get_Name(elem[k-1].Fnum)
		   << "(" << j << ") = "
		   << std::setw(11) << get_L(elem[k-1].Fnum, j)
		   << std::setw(3) << -Quad << std::endl;
	  quad_out << std::fixed << std::setprecision(7) 
		   << std::setw(6) << get_Name(elem[k+1].Fnum)
		   << "(" << j << ") = "
		   << std::setw(11) << get_L(elem[k+1].Fnum, j)
		   << std::setw(3) << -Quad << std::endl;
	}
    quad_out.close();
  }

  free_dvector(b, 1, m);
  free_dvector(b2_lim, 1, n_b2); free_dvector(b2, 1, n_b2);
  free_dvector(db2, 1, n_b2);
  free_dmatrix(A, 1, m, 1, n_b2);
}


void fit_dnu_SLS(const double nu_x,
		 const double beta1_x, const double beta1_y, const int k1,
		 const int k2,
		 const int n_b2, const int b2s[],
		 const double eps, const bool prt)
{
  // Periodic solution: [nu_x, nu_y, beta_x, beta_y]

  int           i, j;
  long int      k;
  double        **A, *b, *w, **U, **V, *b2_lim, *b2, *db2, step, **A_inv;
  double        nu_fract[2], dnu[2];
  double        dalpha1[2], dbeta1[2], dalpha2[2], dbeta2[2], beta2[2];
  tps           ab1[2], ab2[2], dnu1[2], dnu2[2];
  ss_vect<tps>  nus;
  std::ofstream quad_out, A_out;

  const bool   debug  = true;
  const int    m      = 7;
  const double s_cut  = 1e-5, step0 = 0.7;

  b = dvector(1, m); w = dvector(1, n_b2); b2_lim = dvector(1, n_b2);
  b2 = dvector(1, n_b2); db2 = dvector(1, n_b2);
  A = dmatrix(1, m, 1, n_b2); U = dmatrix(1, m, 1, n_b2);
  V = dmatrix(1, n_b2, 1, n_b2); A_inv = dmatrix(1, n_b2, 1, m);
  nu_fract[X_] = fract(nu_x);
  std::cout << std::endl;
  std::cout << std::fixed << std::setprecision(5)
       << "fit_dnu_SLS: nu_x = " << nu_fract[X_] << std::endl;
  std::cout << std::fixed << std::setprecision(5)
       << "beta1_x = " << std::setw(8) << beta1_x
       << ", beta1_y = " << std::setw(8) << beta1_y << std::endl;

  for (i = 1; i <= n_b2; i++) {
    if (b2s[i-1] > 0) {
      b2_lim[i] = b2_max; b2[i] = get_bn(b2s[i-1], 1, Quad);
    } else {
      b2_lim[i] = ds_max;
      k = get_loc(abs(b2s[i-1]), 1) - 1; b2[i] = get_L(elem[k-1].Fnum, 1);
    }
  }

  danot_(3);

  get_Map(); K = MapNorm(Map, g, A1, A0, Map_res, 1); nus = dHdJ(K);
  get_ab(ab1, dnu1, k1); get_ab(ab2, dnu2, k2);

  dnu[X_]     =  nus[3].cst() - nu_fract[X_]; nu_fract[Y_] = nus[4].cst();
  dalpha1[X_] = -h_ijklm(ab1[X_], 0, 1, 0, 0, 0);
  dalpha1[Y_] = -h_ijklm(ab1[Y_], 0, 0, 0, 1, 0);
  dbeta1[X_]  =  h_ijklm(ab1[X_], 1, 0, 0, 0, 0) - beta1_x;
  dbeta1[Y_]  =  h_ijklm(ab1[Y_], 0, 0, 1, 0, 0) - beta1_y;
  dalpha2[X_] = -h_ijklm(ab2[X_], 0, 1, 0, 0, 0);
  dalpha2[Y_] = -h_ijklm(ab2[Y_], 0, 0, 0, 1, 0);
  beta2[X_]   =  h_ijklm(ab2[X_], 1, 0, 0, 0, 0);
  beta2[Y_]   =  h_ijklm(ab2[Y_], 0, 0, 1, 0, 0);

  while ((fabs(dnu[X_]) > eps) ||
	 (scl_alpha_mp*fabs(dalpha1[X_]) > eps) ||
	 (scl_alpha_mp*fabs(dalpha1[Y_]) > eps) ||
	 (scl_beta_mp*fabs(dbeta1[X_]) > eps) ||
	 (scl_beta_mp*fabs(dbeta1[Y_]) > eps) ||
	 (scl_alpha_mp*fabs(dalpha2[X_]) > eps) ||
	 (scl_alpha_mp*fabs(dalpha2[Y_]) > eps)) {
    step = step0;
    for (i = 1; i <= n_b2; i++) {
      for (j = 1; j <= get_n_Kids(abs(b2s[i-1])); j++)
	if (b2s[i-1] > 0)
	  set_bn_par(b2s[i-1], j, Quad, 7);
	else
	  set_s_par(abs(b2s[i-1]), j, 7);

      get_Map(); K = MapNorm(Map, g, A1, A0, Map_res, 1); nus = dHdJ(K);
      get_ab(ab1, dnu1, k1); get_ab(ab2, dnu2, k2);

      A[1][i]  = -h_ijklm_p(nus[3], 0, 0, 0, 0, 0, 7);
      A[2][i]  =  scl_alpha_mp*h_ijklm_p(ab1[X_], 0, 1, 0, 0, 0, 7);
      A[3][i]  =  scl_alpha_mp*h_ijklm_p(ab1[Y_], 0, 0, 0, 1, 0, 7);
      A[4][i]  = -scl_beta_mp*h_ijklm_p(ab1[X_], 1, 0, 0, 0, 0, 7);
      A[5][i]  = -scl_beta_mp*h_ijklm_p(ab1[Y_], 0, 0, 1, 0, 0, 7);
      A[6][i]  =  scl_alpha_mp*h_ijklm_p(ab2[X_], 0, 1, 0, 0, 0, 7);
      A[7][i]  =  scl_alpha_mp*h_ijklm_p(ab2[Y_], 0, 0, 0, 1, 0, 7);

      if (b2s[i-1] < 0)
	for (j = 1; j <= m; j++)
	  A[j][i] *= scl_ds;

      for (j = 1; j <= get_n_Kids(abs(b2s[i-1])); j++)
	if (b2s[i-1] > 0)
	  clr_bn_par(b2s[i-1], j, Quad);
	else
	  clr_s_par(abs(b2s[i-1]), j);
    }

    b[1]  = dnu[X_];
    b[2]  = scl_alpha_mp*dalpha1[X_]; b[3]  = scl_alpha_mp*dalpha1[Y_];
    b[4]  = scl_beta_mp*dbeta1[X_];   b[5]  = scl_beta_mp*dbeta1[Y_];
    b[6]  = scl_alpha_mp*dalpha2[X_]; b[7]  = scl_alpha_mp*dalpha2[Y_];

    SVD_lim(m, n_b2, A, b, b2_lim, s_cut, b2, db2);

    for (i = 1; i <= n_b2; i++) {
      set_dbn_s(b2s[i-1], Quad, step*db2[i]);
      b2[i] = get_bn_s(b2s[i-1], 1, Quad);
    }

    get_Map(); K = MapNorm(Map, g, A1, A0, Map_res, 1);

    while (!stable) {
      // roll back
      for (i = 1; i <= n_b2; i++) {
	set_dbn_s(b2s[i-1], Quad, -step*db2[i]);
	b2[i] = get_bn_s(b2s[i-1], 1, Quad);
      }

      step /= 2.0;
      std::cout << std::endl;
      std::cout << std::scientific << std::setprecision(3)
		<< "step = " << step << std::endl;

      for (i = 1; i <= n_b2; i++) {
	set_dbn_s(b2s[i-1], Quad, step*db2[i]);
	b2[i] = get_bn_s(b2s[i-1], 1, Quad);
      }
	
      get_Map(); K = MapNorm(Map, g, A1, A0, Map_res, 1);
    }
      
    nus = dHdJ(K);
    get_ab(ab1, dnu1, k1); get_ab(ab2, dnu2, k2);

    dnu[X_]     =  nus[3].cst() - nu_fract[X_];
    dnu[Y_]     =  nus[4].cst() - nu_fract[Y_];
    dalpha1[X_] = -h_ijklm(ab1[X_], 0, 1, 0, 0, 0);
    dalpha1[Y_] = -h_ijklm(ab1[Y_], 0, 0, 0, 1, 0);
    dbeta1[X_]  =  h_ijklm(ab1[X_], 1, 0, 0, 0, 0) - beta1_x;
    dbeta1[Y_]  =  h_ijklm(ab1[Y_], 0, 0, 1, 0, 0) - beta1_y;
    dalpha2[X_] = -h_ijklm(ab2[X_], 0, 1, 0, 0, 0);
    dalpha2[Y_] = -h_ijklm(ab2[Y_], 0, 0, 0, 1, 0);
    dbeta2[X_]  =  h_ijklm(ab2[X_], 1, 0, 0, 0, 0) - beta2[X_];
    dbeta2[Y_]  =  h_ijklm(ab2[Y_], 0, 0, 1, 0, 0) - beta2[Y_];
    
    A_out.open("A.dat", std::ios::out);

    for (i = 1; i <= m; i++) {
      for (j = 1; j <= n_b2; j++)
	A_out << std::scientific << std::setprecision(3)
	      << std::setw(11) << A[i][j];
      A_out << std::endl;
    }

    A_out.close();

    for (i = 1; i <= m; i++)
      for (j = 1; j <= n_b2; j++)
        U[i][j] = A[i][j];

    dsvdcmp(U, m, n_b2, w, V);

    for (i = 1; i <= n_b2; i++)
      for (j = 1; j <= m; j++) {
        A_inv[i][j] = 0.0;
        for (k = 1; k <= n_b2; k++)
          if (w[k] != 0.0) 
            A_inv[i][j] += V[i][k]/w[k]*U[j][k];
      }

    A_out.open("A_inv.dat", std::ios::out);

    for (i = 1; i <= n_b2; i++) {
      for (j = 1; j <= m; j++)
	A_out << std::scientific << std::setprecision(3)
	      << std::setw(11) << A_inv[i][j];
      A_out << std::endl;
    }

    A_out.close();

//    exit(0);

    if (debug) {
      std::cout << std::endl;
      std::cout << " Ax = b:" << std::endl;
      std::cout << std::endl;
      for (j = 1; j <= n_b2; j++)
	if (j == 1)
	  std::cout << std::setw(8) << j;
	else
	  std::cout << std::setw(11) << j;
      std::cout << std::endl;
      for (i = 1; i <= m; i++) {
	std::cout << std::setw(2) << i;
	for (j = 1; j <= n_b2; j++)
	  std::cout << std::scientific << std::setprecision(3)
		    << std::setw(11) << A[i][j];
	std::cout << std::scientific << std::setprecision(3)
		  << std::setw(11) << b[i] << std::endl;
      }

      std::cout << std::endl;
      std::cout << std::fixed << std::setprecision(5)
	   << "dnu_x = " << dnu[X_] << ", dnu_y = " << dnu[Y_] << std::endl;
      std::cout << std::fixed << std::setprecision(5)
	   << "dbeta1_x = " << std::setw(8) << dbeta1[X_]
	   << ", dbeta1_y = " << std::setw(8) << dbeta1[Y_]
	   << ", dalpha1_x = " << std::setw(8) << dalpha1[X_]
	   << ", dalpha1_y = " << std::setw(8) << dalpha1[Y_] << std::endl;
      std::cout << std::fixed << std::setprecision(5)
	   << "dbeta2_x = " << std::setw(8) << dbeta2[X_]
	   << ", dbeta2_y = " << std::setw(8) << dbeta2[Y_]
	   << ", dalpha2_x = " << std::setw(8) << dalpha2[X_]
	   << ", dalpha2_y = " << std::setw(8) << dalpha2[Y_] << std::endl;
    }
   
    std::cout << std::endl;
    std::cout << std::fixed << std::setprecision(5)
	 << "nu_x = " << nus[3].cst() << ", nu_y = " << nus[4].cst()
	      << std::endl;
  }

  if (prt) {
    quad_out.open("fit_dnu.dat", std::ios::out);
    quad_out << std::endl;
    quad_out << "n = 1:" << std::endl;
    for (i = 1; i <= n_b2; i++)
      for (j = 1; j <= get_n_Kids(abs(b2s[i-1])); j++)
	if (b2s[i-1] > 0)
	  quad_out << std::fixed << std::setprecision(7) 
		   << std::setw(6) << get_Name(b2s[i-1]) << "(" << j << ") = "
		   << std::setw(11) << get_bnL(b2s[i-1], j, Quad)
		   << std::setw(2) << Quad << std::endl;
	else {
	  k = get_loc(abs(b2s[i-1]), j) - 1;
	  quad_out << std::fixed << std::setprecision(7) 
		   << std::setw(6) << get_Name(elem[k-1].Fnum)
		   << "(" << j << ") = "
		   << std::setw(11) << get_L(elem[k-1].Fnum, j)
		   << std::setw(3) << -Quad << std::endl;
	  quad_out << std::fixed << std::setprecision(7) 
		   << std::setw(6) << get_Name(elem[k+1].Fnum)
		   << "(" << j << ") = "
		   << std::setw(11) << get_L(elem[k+1].Fnum, j)
		   << std::setw(3) << -Quad << std::endl;
	}
    quad_out.close();
  }

  free_dvector(b, 1, m); free_dvector(w, 1, n_b2);
  free_dvector(b2_lim, 1, n_b2); free_dvector(b2, 1, n_b2);
  free_dvector(db2, 1, n_b2);
  free_dmatrix(A, 1, m, 1, n_b2); free_dmatrix(U, 1, m, 1, n_b2);
  free_dmatrix(V, 1, n_b2, 1, n_b2); free_dmatrix(A_inv, 1, n_b2, 1, m);
}


void fit_chrom(const double ksi_x, const double ksi_y,
	       const int n_b3, const int b3s[], const bool prt)
{
  int           i, j;
  float         **A, *b, *w, **U, **V, *db3;
  double        ksi[2], dksi[2], s_max;
  ss_vect<tps>  nus;
  std::ofstream sext_out;

  const bool   debug = false;
  const int    m     = 2; // (ksi_x, ksi_y)
  const double s_cut = 1e-7;

  b = vector(1, m); w = vector(1, n_b3); db3 = vector(1, n_b3);
  A = matrix(1, m, 1, n_b3); U = matrix(1, m, 1, n_b3);
  V = matrix(1, n_b3, 1, n_b3);

  std::cout << std::endl;
  std::cout << std::fixed << std::setprecision(5)
       << "fit_chrom: ksi_x = " << std::setw(8) << ksi_x
       << ", ksi_y = " << std::setw(8) << ksi_y << std::endl;

  danot_(3);
  get_Map();
  danot_(4);
  // lieinit_(4, ss_dim, nd_tps, ndpt_tps, iref_tps, 0);
  K = MapNorm(Map, g, A1, A0, Map_res, 1); nus = dHdJ(K);
  ksi[X_] = h_ijklm(nus[3], 0, 0, 0, 0, 1);
  ksi[Y_] = h_ijklm(nus[4], 0, 0, 0, 0, 1);
  dksi[X_] = ksi[X_] - ksi_x; dksi[Y_] = ksi[Y_] - ksi_y;

  std::cout << std::fixed << std::setprecision(5)
       << "ksi_x = " << std::setw(8) << ksi[X_]
       << ", ksi_y = " << std::setw(8) << ksi[Y_] << std::endl;

  for (i = 1; i <= n_b3; i++) {
    for (j = 1; j <= get_n_Kids(b3s[i-1]); j++)
      set_bn_par(b3s[i-1], j, Sext, 7);

    danot_(3);
    get_Map();
    danot_(4);
    K = MapNorm(Map, g, A1, A0, Map_res, 1); nus = dHdJ(K);

    A[1][i] = -h_ijklm_p(nus[3], 0, 0, 0, 0, 1, 7);
    A[2][i] = -h_ijklm_p(nus[4], 0, 0, 0, 0, 1, 7);

    for (j = 1; j <= get_n_Kids(b3s[i-1]); j++)
      clr_bn_par(b3s[i-1], j, Sext);
  }

  b[1] = dksi[X_]; b[2] = dksi[Y_];

  if (debug) {
    std::cout << std::endl;
    std::cout << " Ax = b:" << std::endl;
    std::cout << std::endl;
    for (i = 1; i <= m; i++) {
      for (j = 1; j <= n_b3; j++)
	std::cout << std::scientific << std::setprecision(3)
		  << std::setw(11) << A[i][j];
      std::cout << std::scientific << std::setprecision(3)
		<< std::setw(11) << b[i] << std::endl;
    }
  }
    
  for (i = 1; i <= m; i++)
    for (j = 1; j <= n_b3; j++)
      U[i][j] = A[i][j];

  svdcmp(U, m, n_b3, w, V);

  s_max = -1e30;
  for (i = 1; i <= n_b3; i++)
    s_max = max(w[i], s_max);
  
  std::cout << std::endl;
  std::cout << "singular values:" << std::endl;
  for (i = 1; i <= n_b3; i++) {
    std::cout << std::scientific << std::setprecision(3)
	      << std::setw(10) << w[i];
    if (w[i]/s_max < s_cut) {
      w[i] = 0.0;
      std::cout << " (zeroed)";
    }
    std::cout << std::endl;
  }

  svbksb(U, w, V, m, n_b3, b, db3);

  if (debug) {
    std::cout << std::endl;
    std::cout << "b3s:" << std::endl;
    std::cout << std::endl;
  }
  for (i = 1; i <= n_b3; i++)
    for (j = 1; j <= get_n_Kids(b3s[i-1]); j++) {
      set_dbn(b3s[i-1], j, Sext, db3[i]);
      if (debug) {
	std::cout << std::fixed << std::setprecision(5)
	     << std::setw(8) << get_bn(b3s[i-1], j, Sext) << std::endl;
      }
    }
  if (debug) std::cout << std::endl;

  danot_(3);
  get_Map();
  danot_(4);
  K = MapNorm(Map, g, A1, A0, Map_res, 1); nus = dHdJ(K);
  ksi[X_] = h_ijklm(nus[3], 0, 0, 0, 0, 1);
  ksi[Y_] = h_ijklm(nus[4], 0, 0, 0, 0, 1);
  dksi[X_] = ksi[X_] - ksi_x; dksi[Y_] = ksi[Y_] - ksi_y;
  std::cout << std::endl;
  std::cout << std::fixed << std::setprecision(5)
       << "ksi_x = " << std::setw(8) << ksi[X_]
       << ", ksi_y = " << std::setw(8) << ksi[Y_] << std::endl;

  if (prt) {
    sext_out.open("fit_chrom.dat", std::ios::out);
    sext_out << std::endl;
    sext_out << "n = 1:" << std::endl;
    for (i = 1; i <= n_b3; i++)
      for (j = 1; j <= get_n_Kids(b3s[i-1]); j++)
	sext_out << std::fixed << std::setprecision(7) 
		 << std::setw(6) << get_Name(b3s[i-1]) << "(" << j << ") = "
		 << std::setw(11) << get_bnL(b3s[i-1], j, Sext)
		 << std::setw(2) << Sext << std::endl;
    sext_out.close();
  }

  free_vector(b, 1, m); free_vector(w, 1, n_b3);
  free_vector(db3, 1, n_b3);
  free_matrix(A, 1, m, 1, n_b3); free_matrix(U, 1, m, 1, n_b3);
  free_matrix(V, 1, n_b3, 1, n_b3);
}


void fit_chrom1(const double ksi_x, const double ksi_y,
	        const int n_b3, const int b3s[], const double eps,
	        const bool prt)
{
  int           i, j;
  double        **A, *b, *db3, *b3, *b3_max;
  double        ksi[2], dksi[2], h_10002, h_20001, h_00201;
  tps           h, h_re, h_im;
  ss_vect<tps>  nus;
  std::ofstream sext_out;

  const bool   debug = false;
  const int    m     = 5; // (ksi_x, ksi_y)
  const double scl_beta2 = 1e-3, scl_eta2 = 1e-3, b3L_max = 3.5;
  const double s_cut = 1e-10;

  b = dvector(1, m); db3 = dvector(1, n_b3);
  b3 = dvector(1, n_b3); b3_max = dvector(1, n_b3);
  A = dmatrix(1, m, 1, n_b3);

  std::cout << std::endl;
  std::cout << std::fixed << std::setprecision(5)
       << "fit_chrom: ksi_x = " << std::setw(8) << ksi_x
       << ", ksi_y = " << std::setw(8) << ksi_y << std::endl;

  danot_(4);

  for (i = 1; i <= n_b3; i++)
    b3_max[i] = b3L_max;

  get_Map(); K = MapNorm(Map, g, A1, A0, Map_res, 1); nus = dHdJ(K);
  h = get_h(); CtoR(h, h_re, h_im);
  ksi[X_] = h_ijklm(nus[3], 0, 0, 0, 0, 1);
  ksi[Y_] = h_ijklm(nus[4], 0, 0, 0, 0, 1);
  h_10002 = scl_eta2*h_ijklm(h_re, 1, 0, 0, 0, 2);
  h_20001 = scl_beta2*h_ijklm(h_re, 2, 0, 0, 0, 1);
  h_00201 = scl_beta2*h_ijklm(h_re, 0, 0, 2, 0, 1);
  dksi[X_] = ksi[X_] - ksi_x; dksi[Y_] = ksi[Y_] - ksi_y;

  std::cout << std::fixed << std::setprecision(5)
       << "ksi_x = " << std::setw(8) << ksi[X_]
       << ", ksi_y = " << std::setw(8) << ksi[Y_]
       << ", h_10002 = " << std::setw(8) << h_10002
       << ", h_20001 = " << std::setw(8) << h_20001
       << ", h_00201 = " << std::setw(8) << h_00201 << std::endl;

  while ((fabs(dksi[X_]) > eps) || (fabs(dksi[Y_]) > eps)) {
    for (i = 1; i <= n_b3; i++) {
      for (j = 1; j <= get_n_Kids(b3s[i-1]); j++)
	set_bn_par(b3s[i-1], j, Sext, 7);
      b3[i] = get_bnL(b3s[i-1], 1, Sext);

      get_Map(); K = MapNorm(Map, g, A1, A0, Map_res, 1); nus = dHdJ(K);
      h = get_h(); CtoR(h, h_re, h_im);

      A[1][i] = -h_ijklm_p(nus[3], 0, 0, 0, 0, 1, 7);
      A[2][i] = -h_ijklm_p(nus[4], 0, 0, 0, 0, 1, 7);
      A[3][i] = -scl_eta2*h_ijklm_p(h_re, 1, 0, 0, 0, 2, 7);
      A[4][i] = -scl_eta2*h_ijklm_p(h_re, 2, 0, 0, 0, 1, 7);
      A[5][i] = -scl_eta2*h_ijklm_p(h_re, 0, 0, 2, 0, 1, 7);

      for (j = 1; j <= get_n_Kids(b3s[i-1]); j++)
	clr_bn_par(b3s[i-1], j, Sext);
    }

    b[1] = dksi[X_]; b[2] = dksi[Y_];
    b[3] = h_10002; b[4] = h_20001; b[5] = h_00201;

    if (debug) {
      std::cout << std::endl;
      std::cout << " Ax = b:" << std::endl;
      std::cout << std::endl;
      for (i = 1; i <= m; i++) {
	for (j = 1; j <= n_b3; j++)
	  std::cout << std::scientific << std::setprecision(3)
		    << std::setw(11) << A[i][j];
	std::cout << std::scientific << std::setprecision(3)
		  << std::setw(11) << b[i] << std::endl;
      }
    }
    
    SVD_lim(m, n_b3, A, b, b3_max, s_cut, b3, db3);

    if (debug) {
      std::cout << std::endl;
      std::cout << "b3s:" << std::endl;
      std::cout << std::endl;
    }
    for (i = 1; i <= n_b3; i++)
      for (j = 1; j <= get_n_Kids(b3s[i-1]); j++) {
	set_dbn(b3s[i-1], j, Sext, db3[i]);
	if (debug) {
	  std::cout << std::fixed << std::setprecision(5)
	       << std::setw(8) << get_bn(b3s[i-1], j, Sext) << std::endl;
	}
      }
    if (debug) std::cout << std::endl;

    get_Map(); K = MapNorm(Map, g, A1, A0, Map_res, 1); nus = dHdJ(K);
    h = get_h(); CtoR(h, h_re, h_im);
    ksi[X_] = h_ijklm(nus[3], 0, 0, 0, 0, 1);
    ksi[Y_] = h_ijklm(nus[4], 0, 0, 0, 0, 1);
    h_10002 = scl_eta2*h_ijklm(h_re, 1, 0, 0, 0, 2);
    h_20001 = scl_beta2*h_ijklm(h_re, 2, 0, 0, 0, 1);
    h_00201 = scl_beta2*h_ijklm(h_re, 0, 0, 2, 0, 1);
    dksi[X_] = ksi[X_] - ksi_x; dksi[Y_] = ksi[Y_] - ksi_y;
    std::cout << std::endl;
    std::cout << std::fixed << std::setprecision(5)
	 << "ksi_x = " << std::setw(8) << ksi[X_]
	 << ", ksi_y = " << std::setw(8) << ksi[Y_]
	 << ", h_10002 = " << std::setw(8) << h_10002
	 << ", h_20001 = " << std::setw(8) << h_20001
	 << ", h_00201 = " << std::setw(8) << h_00201 << std::endl;
  }


  if (prt) {
    sext_out.open("fit_chrom.dat", std::ios::out);
    sext_out << std::endl;
    sext_out << "n = 1:" << std::endl;
    for (i = 1; i <= n_b3; i++)
      for (j = 1; j <= get_n_Kids(b3s[i-1]); j++)
	sext_out << std::fixed << std::setprecision(7) 
		 << std::setw(6) << get_Name(b3s[i-1]) << "(" << j << ") = "
		 << std::setw(11) << get_bnL(b3s[i-1], j, Sext)
		 << std::setw(2) << Sext << std::endl;
    sext_out.close();
  }

  free_dvector(b, 1, m);
  free_dvector(db3, 1, n_b3);
  free_dvector(b3, 1, n_b3); free_dvector(b3_max, 1, n_b3);
  free_dmatrix(A, 1, m, 1, n_b3);
}


void fit_disp(const double eta1_x, const int k1,
	      const double eta2_x, const int k2,
	      const int n_b2, const int b2s[],
	      const double eps, const bool prt)
{
  int           i, j;
  float         **A, *b, *w, **U, **V, *b2, *db2;
  double        s_max, deta1, deta2;
  tps           eta1, eta2;
  ss_vect<tps>  Mk;
  std::ofstream quad_out;

  const bool   debug = true;
  const int    m     = 2; // periodic solution: eta_x
  const double s_cut = 1e-7, step = 0.7;

  b = vector(1, m); w = vector(1, n_b2);
  b2 = vector(1, n_b2); db2 = vector(1, n_b2);
  A = matrix(1, m, 1, n_b2); U = matrix(1, m, 1, n_b2);
  V = matrix(1, n_b2, 1, n_b2);

  std::cout << std::endl;
  std::cout << std::fixed << std::setprecision(5)
	    << "fit_disp: eta1_x = " << eta1_x << ", eta2_x = " << eta2_x
	    << std::endl;

  danot_(2);

  get_Map();
  // compute fixed point to second order to includ parameter dependance
  Mk = get_Map(k1); K = MapNorm(Mk, g, A1, A0, Map_res, 2);
  eta1 = get_eta(0)[x_];
  Mk = get_Map(k2); K = MapNorm(Mk, g, A1, A0, Map_res, 2);
  eta2 = get_eta(0)[x_];

  deta1 = h_ijklm(eta1, 0, 0, 0, 0, 1) - eta1_x;
  deta2 = h_ijklm(eta2, 0, 0, 0, 0, 1) - eta2_x;

  std::cout << std::fixed << std::setprecision(5)
       << "deta1 = " << deta1 << ", deta2 = " << deta2 << std::endl;

  while ((fabs(deta1) > eps) || (fabs(deta2) > eps)) {
    for (i = 1; i <= n_b2; i++) {
      for (j = 1; j <= get_n_Kids(abs(b2s[i-1])); j++)
	if (b2s[i-1] > 0)
	  set_bn_par(b2s[i-1], j, Quad, 7);
	else
	  set_s_par(abs(b2s[i-1]), j, 7);

      get_Map();
      // compute fixed point to second order to includ parameter dependance
      Mk = get_Map(k1); K = MapNorm(Mk, g, A1, A0, Map_res, 2);
      eta1 = get_eta(0)[x_];
      Mk = get_Map(k2); K = MapNorm(Mk, g, A1, A0, Map_res, 2);
      eta2 = get_eta(0)[x_];

      A[1][i] = -h_ijklm_p(eta1, 0, 0, 0, 0, 1, 7);
      A[2][i] = -h_ijklm_p(eta2, 0, 0, 0, 0, 1, 7);
 
      if (b2s[i-1] < 0)
	for (j = 1; j <= m; j++)
	  A[j][i] *= scl_ds;

      for (j = 1; j <= get_n_Kids(abs(b2s[i-1])); j++)
	if (b2s[i-1] > 0)
	  clr_bn_par(b2s[i-1], j, Quad);
	else
	  clr_s_par(abs(b2s[i-1]), j);
    }

    b[1] = h_ijklm(eta1, 0, 0, 0, 0, 1) - eta1_x;
    b[2] = h_ijklm(eta2, 0, 0, 0, 0, 1) - eta2_x;

    if (debug) {
      std::cout << std::endl;
      std::cout << " Ax = b:" << std::endl;
      std::cout << std::endl;
      for (i = 1; i <= m; i++) {
	for (j = 1; j <= n_b2; j++)
	  std::cout << std::scientific << std::setprecision(3)
		    << std::setw(11) << A[i][j];
	std::cout << std::scientific << std::setprecision(3)
		  << std::setw(11) << b[i] << std::endl;
      }
    }
    
    for (i = 1; i <= m; i++)
      for (j = 1; j <= n_b2; j++)
	U[i][j] = A[i][j];

    svdcmp(U, m, n_b2, w, V);

    s_max = -1e30;
    for (i = 1; i <= n_b2; i++)
      s_max = max(w[i], s_max);
  
    if (debug) {
      std::cout << std::endl;
      std::cout << "singular values:" << std::endl;
      for (i = 1; i <= n_b2; i++) {
	std::cout << std::scientific << std::setprecision(3)
		  << std::setw(10) << w[i];
	if (w[i]/s_max < s_cut) {
	  w[i] = 0.0;
	  std::cout << " (zeroed)";
	}
	std::cout << std::endl;
      }
    }

    svbksb(U, w, V, m, n_b2, b, db2);

    for (i = 1; i <= n_b2; i++) {
      set_dbn_s(b2s[i-1], Quad, step*db2[i]);
      b2[i] = get_bn_s(b2s[i-1], 1, Quad);
    }
    
    get_Map();
    Mk = get_Map(k1); K = MapNorm(Mk, g, A1, A0, Map_res, 2);
    eta1 = get_eta(0)[x_];
    Mk = get_Map(k2); K = MapNorm(Mk, g, A1, A0, Map_res, 2);
    eta2 = get_eta(0)[x_];

    deta1 = h_ijklm(eta1, 0, 0, 0, 0, 1) - eta1_x;
    deta2 = h_ijklm(eta2, 0, 0, 0, 0, 1) - eta2_x;

    std::cout << std::fixed << std::setprecision(5)
	      << "deta1 = " << deta1 << ", deta2 = " << deta2 << std::endl;
  }

  if (prt) {
    quad_out.open("fit_disp.dat", std::ios::out);
    quad_out << std::endl;
    quad_out << "n = 1:" << std::endl;
    for (i = 1; i <= n_b2; i++)
      for (j = 1; j <= get_n_Kids(abs(b2s[i-1])); j++)
	if (b2s[i-1] > 0)
	  quad_out << std::fixed << std::setprecision(7) 
		   << std::setw(6) << get_Name(b2s[i-1]) << "(" << j << ") = "
		   << std::setw(11) << get_bnL(b2s[i-1], j, Quad)
		   << std::setw(2) << Quad << std::endl;
	else {
	  quad_out << std::fixed << std::setprecision(7) 
		   << std::setw(6) << get_Name(abs(b2s[i-1]))
		   << "(" << j << ") = "
		   << std::setw(11) << get_L(abs(b2s[i-1]), j)
		   << std::setw(3) << -Quad << std::endl;
	}
    quad_out.close();
  }

  free_vector(b, 1, m); free_vector(w, 1, n_b2);
  free_vector(b2, 1, n_b2); free_vector(db2, 1, n_b2);
  free_matrix(A, 1, m, 1, n_b2); free_matrix(U, 1, m, 1, n_b2);
  free_matrix(V, 1, n_b2, 1, n_b2);
}


void fit_DBA(const double nu_x, const double nu_y, const long int k,
	     const int n_b2, const int b2s[], const double eps, const bool prt)
{
  // Periodic solution: [nu_x, nu_y, beta_x, beta_y]

  int           i, j;
  double        **A, *b, *b2_lim, *b2, *db2, step;
  double        nu_fract[2], dnu[2], dalpha1[2], dalpha2[2];
  tps           ab1[2], ab2[2], dnu1[2];
  ss_vect<tps>  nus, dnus, eta;
  std::ofstream quad_out;

  const bool   debug = true;
  const int    m     = 6;
  const double s_cut = 1e-7, step0 = 0.7, scl_eta = 1e3, scl_etap = 1e3;

  b = dvector(1, m); b2_lim = dvector(1, n_b2);
  b2 = dvector(1, n_b2); db2 = dvector(1, n_b2);
  A = dmatrix(1, m, 1, n_b2);

  nu_fract[X_] = fract(nu_x); nu_fract[Y_] = fract(nu_y);
  std::cout << std::endl;
  std::cout << std::fixed << std::setprecision(5)
       << "fit_DBA: nu_x = " << nu_fract[X_] << ", nu_y = " << nu_fract[Y_]
       << std::endl;

  for (i = 1; i <= n_b2; i++) {
    if (b2s[i-1] > 0) {
      b2_lim[i] = b2_max; b2[i] = get_bn(b2s[i-1], 1, Quad);
    } else {
      b2_lim[i] = ds_max; b2[i] = get_L(abs(b2s[i-1]), 1);
    }
  }

  danot_(3);

  get_Map(); K = MapNorm(Map, g, A1, A0, Map_res, 1); nus = dHdJ(K);

  dnu[X_]     =  nus[3].cst() - nu_fract[X_];
  dnu[Y_]     =  nus[4].cst() - nu_fract[Y_];
  dalpha1[X_] = -h_ijklm(ab1[X_], 0, 1, 0, 0, 0);
  dalpha1[Y_] = -h_ijklm(ab1[Y_], 0, 0, 0, 1, 0);
  dalpha2[X_] = -h_ijklm(ab2[X_], 0, 1, 0, 0, 0);
  dalpha2[Y_] = -h_ijklm(ab2[Y_], 0, 0, 0, 1, 0);

  while ((fabs(dnu[X_]) > eps) || (fabs(dnu[Y_]) > eps)) {
    step = step0;
    for (i = 1; i <= n_b2; i++) {
      for (j = 1; j <= get_n_Kids(abs(b2s[i-1])); j++)
	if (b2s[i-1] > 0)
	  set_bn_par(b2s[i-1], j, Quad, 7);
	else
	  set_s_par(abs(b2s[i-1]), j, 7);

      get_Map(); K = MapNorm(Map, g, A1, A0, Map_res, 1);
      nus = dHdJ(K);
      eta = get_eta(k); get_ab(ab1, dnu1, k);

      // eta_x
      A[1][i] = -scl_eta*h_ijklm_p(eta[x_], 1, 0, 0, 0, 1, 7);
      // eta'_x
      A[2][i] = -scl_etap*h_ijklm_p(eta[px_], 0, 1, 0, 0, 1, 7);
      // alpha_x,y
      A[3][i]  =  scl_alpha_mp*h_ijklm_p(ab1[X_], 0, 1, 0, 0, 0, 7);
      A[4][i]  =  scl_alpha_mp*h_ijklm_p(ab1[Y_], 0, 0, 0, 1, 0, 7);
      // beta_x,y
      A[5][i]  = -scl_beta_mp*h_ijklm_p(ab1[X_], 1, 0, 0, 0, 0, 7);
      A[6][i]  = -scl_beta_mp*h_ijklm_p(ab1[Y_], 0, 0, 1, 0, 0, 7);

      if (b2s[i-1] < 0)
	for (j = 1; j <= m; j++)
	  A[j][i] *= scl_ds;

      for (j = 1; j <= get_n_Kids(abs(b2s[i-1])); j++)
	if (b2s[i-1] > 0)
	  clr_bn_par(b2s[i-1], j, Quad);
	else
	  clr_s_par(abs(b2s[i-1]), j);
    }

    b[1]  = dnu[X_];                  b[2]  = dnu[Y_];
    b[3]  = scl_alpha_mp*dalpha1[X_]; b[4]  = scl_alpha_mp*dalpha1[Y_];
    b[5]  = scl_alpha_mp*dalpha2[X_]; b[6]  = scl_alpha_mp*dalpha2[Y_];

    SVD_lim(m, n_b2, A, b, b2_lim, s_cut, b2, db2);

    for (i = 1; i <= n_b2; i++) {
      set_dbn_s(b2s[i-1], Quad, step*db2[i]);
      b2[i] = get_bn_s(b2s[i-1], 1, Quad);
    }

    get_Map(); K = MapNorm(Map, g, A1, A0, Map_res, 1);

    while (!stable) {
      // roll back
      for (i = 1; i <= n_b2; i++) {
	set_dbn_s(b2s[i-1], Quad, -step*db2[i]);
	b2[i] = get_bn_s(b2s[i-1], 1, Quad);
      }

      step /= 2.0;
      std::cout << std::endl;
      std::cout << std::scientific << std::setprecision(3)
	<< "step = " << step << std::endl;

      for (i = 1; i <= n_b2; i++) {
	set_dbn_s(b2s[i-1], Quad, step*db2[i]);
	b2[i] = get_bn_s(b2s[i-1], 1, Quad);
      }
	
      get_Map(); K = MapNorm(Map, g, A1, A0, Map_res, 1);
    }
      
    nus = dHdJ(K);

    dnu[X_] = nus[3].cst() - nu_fract[X_];
    dnu[Y_] = nus[4].cst() - nu_fract[Y_];
    
    if (debug) {
      std::cout << std::endl;
      std::cout << " Ax = b:" << std::endl;
      std::cout << std::endl;
      for (i = 1; i <= m; i++) {
	for (j = 1; j <= n_b2; j++)
	  std::cout << std::scientific << std::setprecision(3)
		    << std::setw(11) << A[i][j];
	std::cout << std::scientific << std::setprecision(3)
		  << std::setw(11) << b[i] << std::endl;
      }
	
      std::cout << std::endl;
      std::cout << std::fixed << std::setprecision(5)
		<< "dnu_x = " << dnu[X_] << ", dnu_y = " << dnu[Y_]
		<< std::endl;
    }
   
    std::cout << std::endl;
    std::cout << std::fixed << std::setprecision(5)
	      << "nu_x = " << nus[3].cst() << ", nu_y = " << nus[4].cst()
	      << std::endl;
  }

  if (prt) {
    quad_out.open("fit_tune.dat", std::ios::out);
    quad_out << std::endl;
    quad_out << "n = 1:" << std::endl;
    for (i = 1; i <= n_b2; i++)
      for (j = 1; j <= get_n_Kids(abs(b2s[i-1])); j++)
	if (b2s[i-1] > 0)
	  quad_out << std::fixed << std::setprecision(7) 
		   << std::setw(6) << get_Name(b2s[i-1])
		   << "(" << j << ") = "
		   << std::setw(11) << get_bnL(b2s[i-1], j, Quad)
		   << std::setw(2) << Quad << std::endl;
	else {
	  quad_out << std::fixed << std::setprecision(7) 
		   << std::setw(6) << get_Name(abs(b2s[i-1]))
		   << "(" << j << ") = "
		   << std::setw(11) << get_L(abs(b2s[i-1]), j)
		   << std::setw(3) << -Quad << std::endl;
	}
    quad_out.close();
  }

  free_dvector(b, 1, m);
  free_dvector(b2_lim, 1, n_b2); free_dvector(b2, 1, n_b2);
  free_dvector(db2, 1, n_b2);
  free_dmatrix(A, 1, m, 1, n_b2);
}


double get_dynap(const double r, const double delta, const int n,
		 const double eps, const int n_pts,
		 double x_min[], double x_max[], bool map)
{
  /* Determine the dynamical aperture by tracking.
     Assumes mid-plane symmetry.                    */

  int           i, j;
  double        r1, phi, x0[2]={0e0, 0e0}, x1[2], x2[2], DA, alpha[2], beta[2];
  std::ofstream os;

  os.open("dynap.dat", std::ios::out);
  os << "# Dynamical Aperture:" << std::endl;
  os << "#    x      y" << std::endl;
  os << "#   [mm]   [mm]" << std::endl;
  os << "#" << std::endl;

  danot_(1);

  get_Map(); K = MapNorm(Map, g, A1, A0, Map_res, 1); get_ab(alpha, beta, 0);

  for (i = 0; i < 2; i++) {
    x_min[i] = 0.0; x_max[i] = 0.0;
  }

  DA = 0.0; r1 = r;
  for (i = 0; i < n_pts; i++) {
    phi = i*pi/(n_pts-1);
    if (i == 0)
      phi = 1e-3;
    else if (i == n_pts-1)
      phi -= 1e-3;
    get_r_stable(r1, phi, delta, n, eps, map);
    x2[X_] = r1*cos(phi); x2[Y_] = r1*sin(phi);
    for (j = 0; j < 2; j++) {
      x_min[j] = min(x2[j], x_min[j]); x_max[j] = max(x2[j], x_max[j]);
    }
    if (i == 0) {
      x0[X_] = x2[X_]; x0[Y_] = x2[Y_];
    } else
      DA += x1[X_]*x2[Y_] - x2[X_]*x1[Y_];
    os << std::fixed << std::setprecision(2)
	 << std::setw(8) << 1e3*x2[X_] << std::setw(8) << 1e3*x2[Y_]
       << std::endl;
    x1[X_] = x2[X_]; x1[Y_] = x2[Y_];
  }
  DA += x2[X_]*x0[Y_] - x0[X_]*x2[Y_];
  // factor of 2 from mid-plane symmetry
  DA = fabs(DA)/sqrt(beta[X_]*beta[Y_]);

  std::cout << std::endl;
  std::cout << std::fixed << std::setprecision(1)
       << "DA^ = " << std::setw(6) << 1e6*DA
       << ", x^ = " << std::setw(5) << 1e3*x_min[X_] << " - "
       << std::setw(4) << 1e3*x_max[X_] << " mm"
       << ", y^ = " << std::setw(5) << 1e3*x_min[Y_] << " - "
       << std::setw(4) << 1e3*x_max[Y_] << " mm" << std::endl;

  return DA;
} 


void get_rad(char *cav_name, double I_beam)
{
  int    j, h_rf, jj[ss_dim];
  double V_rf, f_rf, U0, phi0, delta_rf, alpha_c;

  get_cav(get_Fnum(cav_name), 1, h_rf, V_rf, f_rf);
  U0 = 1e9*E0*dE; phi0 = fabs(asin(U0/V_rf));

  for (j = 0; j < ss_dim; j++)
    jj[j] = 0;
  jj[delta_] = 1;
  alpha_c = Map[ct_][jj]/elem[n_elem-1].S;

  delta_rf = sqrt(-V_rf*cos(pi-phi0)*(2.0-(pi-2.0*(pi-phi0))*tan(pi-phi0))
	     /(alpha_c*pi*h_rf*1e9*E0));

  std::cout << std::endl;
  std::cout << std::fixed << std::setprecision(3)
       << "alphac           = " << alpha_c << std::endl;
  std::cout << std::fixed << std::setprecision(1)
       << "U0               = " << 1e-3*U0 << " keV" << std::endl;
  std::cout << std::fixed << std::setprecision(2)
       << "phi0             = 180 - " << phi0*180.0/pi << " deg" << std::endl;
  std::cout << std::endl;
  std::cout << std::fixed << std::setprecision(2)
       << "RF bucket height = " << 1e2*delta_rf <<"%" << std::endl;
}


void get_emittance(void)
{
  int    i;
  double U0, C;
  double alpha_c, alpha_z, beta_z, gamma_z, sigma_s, sigma_delta;

//  danot_(2);

  C = elem[n_elem-1].S;

  // compute one-turn map for damped system
  rad_on = true; emittance_on = true; cavity_on = true; totpath_on = false;

  get_COD(10, 1e-10, 0.0, true); K = MapNorm(Map, g, A1, A0, Map_res, 1);

  U0 = 1e9*E0*dE;

  // compute imaginary part of eigenvalues
  gettura_(nu_, rad_);

  // compute diffusion coefficients for the invariants (obtained from A1)
  // A1 = [Re(v_i) Im(v_i), , ], i = 1, 2, 3
  A1 += fixed_point; A1.propagate(1, n_elem); A1 -= A1.cst();

  for (i = 0; i < 3; i++) {
    // compute emittances
    eps_[i] = -D_[i]/(4.0*rad_[i]);
    // compute partition numbers
    part_numb_[i] = 2.0*(1.0+fixed_point[delta_])*rad_[i]/dE;
    // compute damping times
    tau_[i] = -C/(clight*rad_[i]);
    if (nu_[i] < 0.0) nu_[i] += 1.0;
  }

  // undamped system
  rad_on = false; emittance_on = false;

  get_COD(10, 1e-10, 0.0, true); K = MapNorm(Map, g, A1, A0, Map_res, 1);

  // Note, [x, p_x, y, p_y, ct, delta] for 6-dim dynamics

  alpha_c = Map[ct_][delta_]/C;

  // longitudinal alpha and beta
  alpha_z = -A1[ct_][ct_]*A1[delta_][ct_] - A1[ct_][delta_]*A1[delta_][delta_];
  beta_z = sqr(A1[ct_][ct_]) + sqr(A1[ct_][delta_]);
  gamma_z = (1.0+sqr(alpha_z))/beta_z;

  // bunch size
  sigma_s = sqrt(beta_z*eps_[Z_].cst());
  sigma_delta = sqrt(gamma_z*eps_[Z_].cst());

  std::cout << std::endl;
  std::cout << std::scientific << std::setprecision(3)
       << "momentum comp.:       alpha_c     = "
       << std::setw(9) << alpha_c << std::endl;
  std::cout << std::fixed << std::setprecision(1)
       << "dE per turn [keV]:    U0          = "
       << std::setw(9) << 1e-3*U0 << std::endl;
  std::cout << std::fixed << std::setprecision(3)
       << "part. numbers:        J_x         = "
	    << std::setw(9) << part_numb_[X_]
       << ",     J_y = " << std::setw(9) << part_numb_[Y_]
       << ",    J_z = " << std::setw(9) << part_numb_[Z_] << std::endl;
  std::cout << std::scientific << std::setprecision(3)
       << "emittance [m.rad]:    eps_x       = "
	    << std::setw(9) << eps_[X_].cst()
       << ",   eps_y = " << std::setw(9) << eps_[Y_].cst()
       << ",  eps_z = " << std::setw(9) << eps_[Z_].cst() << std::endl;
  std::cout << std::scientific << std::setprecision(3)
       << "bunch length [mm]:    sigma_s     = "
       << std::setw(9) << 1e3*sigma_s << std::endl;
  std::cout << std::scientific << std::setprecision(3)
       << "momentum spread:      sigma_delta = "
       << std::setw(9) << sigma_delta << std::endl;
  std::cout << std::endl;
  std::cout << std::fixed << std::setprecision(1)
       << "damping times [msec]:  tau_x      = " << std::setw(9) << 1e3*tau_[X_]
       << ",   tau_y = "  << std::setw(9)<< 1e3*tau_[Y_]
       << ",  tau_z = " << std::setw(9) << 1e3*tau_[Z_] << std::endl;
  std::cout << std::fixed << std::setprecision(5)
       << "fractional tunes:      nu_x       = " << std::setw(9) << nu_[X_]
       << ",    nu_y = " << std::setw(9) << nu_[Y_]
       << ",   nu_z = " << std::setw(9) << nu_[Z_] << std::endl;
  std::cout << std::fixed << std::setprecision(5)
       << "                       1-nu_x     = " << std::setw(9) << 1.0-nu_[X_]
       << ",  1-nu_y = " << std::setw(9) << 1.0-nu_[Y_]
       << ", 1-nu_z = " << std::setw(9) << 1.0-nu_[Z_] << std::endl;
}
