void h_zero(const int no, const int n_iter,
	    const double nu_x, const double nu_y, const bool prt_iter)
{
  bool          first = true;
  int           i, j, m, n;
  double        **A, *b, *w, **U, **V, *dbn, *bn, *bn_max, **A_inv;
  double        s_max, L;
  tps           h, h_re, h_im, K_re, K_im;
  ofstream      h_out, K_out;

  const bool    prt = true, lin_chrom = false, tune_shift = false;

  const int     m_max      = 100; // max no of constraints
  const double  s_cut      = 1e-8;
  const double  scl_ksi1   = 1e3, scl_ksi2      = 1e0, scl_ksi3  = 1e0;
  const double  scl_ksi4   = 1e0, scl_dnu_delta = 1e0;
  const double  scl_geom   = 1e0;
  const double  scl_chrom1 = 1e0, scl_chrom2    = 1e0;
  const double  scl_alphac2 = 0e-2;

  b = dvector(1, m_max); w = dvector(1, n_prm); dbn = dvector(1, n_prm);
  bn = dvector(1, n_prm); bn_max = dvector(1, n_prm);
  A = dmatrix(1, m_max, 1, n_prm); U = dmatrix(1, m_max, 1, n_prm);
  V = dmatrix(1, n_prm, 1, n_prm);
  A_inv = dmatrix(1, n_prm, 1, m_max);

  for (i = 1; i <= n_prm; i++) {
    L = get_L(prms[i-1], 1);
    if (bns[i-1] == Sext)
      bn_max[i] = b3L_max; 
    else if (bns[i-1] == Oct)
      bn_max[i] = b4L_max; 
    else {
      cout << "h_zero: undefined limit" << endl;
      exit(0);
    }
    // Note, Jacobian is in multipole strengths
    if (L != 0.0) bn_max[i] /= L;
    bn[i] = get_bn(prms[i-1], 1, bns[i-1]);
  }

  if (first) {
    // store initial sextupole strengths
    first = false;
    cout << endl;
    cout << "initial b3s:" << endl;
    cout << endl;
    for (i = 1; i <= n_prm; i++) {
      for (j = 1; j <= get_n_Kids(prms[i-1]); j++)
	b3s[i-1] = get_bn(prms[i-1], j, bns[i-1]);
      cout << scientific << setprecision(3)
	   << setw(11) << b3s[i-1] << setw(2) << bns[i-1] << endl;
    }
  } else if (true) {
    // initialize sextupoles
    for (i = 1; i <= n_prm; i++)
      for (j = 1; j <= get_n_Kids(prms[i-1]); j++)
	set_bn(prms[i-1], j, bns[i-1], b3s[i-1]);
  }

  n = 0;
  do {
    n++;
    for (i = 1; i <= n_prm; i++) {
      for (j = 1; j <= get_n_Kids(prms[i-1]); j++)
	set_bn_par(prms[i-1], j, bns[i-1], 7);
      
      danot_(no-1); get_Map_N();
      danot_(no);
      h = get_h(); h = h*Id; CtoR(h, h_re, h_im);
      if (!tune_shift)
	K_re = h_re;
      else {
	K = MapNorm(Map, g, A1, A0, Map_res, 1); K = K*Id; CtoR(K, K_re, K_im);
      }

      // Note, sine terms cancel for mirror symmetry
      m = 0;
      if (no >= 4) {
	// linear chromaticity
	A[++m][i] = scl_ksi1*h_ijklm_p(h_re, 1, 1, 0, 0, 1, 7);
	A[++m][i] = scl_ksi1*h_ijklm_p(h_re, 0, 0, 1, 1, 1, 7);

	// 1st order chromatic terms
	A[++m][i] = scl_chrom1*h_ijklm_p(h_re, 2, 0, 0, 0, 1, 7);
	A[++m][i] = scl_chrom1*h_ijklm_p(h_re, 0, 0, 2, 0, 1, 7);
	
	// 2nd order chromatic terms
        A[++m][i] = scl_chrom2*h_ijklm_p(h_re, 1, 0, 0, 0, 2, 7);

	// 1st order geometric terms
	A[++m][i] = scl_geom*h_ijklm_p_scl(h_re, 2, 1, 0, 0, 0, 7, nu_x, nu_y);
	A[++m][i] = scl_geom*h_ijklm_p_scl(h_re, 3, 0, 0, 0, 0, 7, nu_x, nu_y);
	A[++m][i] = scl_geom*h_ijklm_p_scl(h_re, 1, 0, 1, 1, 0, 7, nu_x, nu_y);
	A[++m][i] = scl_geom*h_ijklm_p_scl(h_re, 1, 0, 0, 2, 0, 7, nu_x, nu_y);
	A[++m][i] = scl_geom*h_ijklm_p_scl(h_re, 1, 0, 2, 0, 0, 7, nu_x, nu_y);
      }

      if (no >= 5) {
	// 2nd order geometric terms
	A[++m][i] = scl_geom*h_ijklm_p_scl(h_re, 3, 1, 0, 0, 0, 7, nu_x, nu_y);
	A[++m][i] = scl_geom*h_ijklm_p_scl(h_re, 2, 0, 1, 1, 0, 7, nu_x, nu_y);
	A[++m][i] = scl_geom*h_ijklm_p_scl(h_re, 4, 0, 0, 0, 0, 7, nu_x, nu_y);
	A[++m][i] = scl_geom*h_ijklm_p_scl(h_re, 2, 0, 0, 2, 0, 7, nu_x, nu_y);
	A[++m][i] = scl_geom*h_ijklm_p_scl(h_re, 2, 0, 2, 0, 0, 7, nu_x, nu_y);
	A[++m][i] = scl_geom*h_ijklm_p_scl(h_re, 0, 0, 4, 0, 0, 7, nu_x, nu_y);
	A[++m][i] = scl_geom*h_ijklm_p_scl(h_re, 1, 1, 2, 0, 0, 7, nu_x, nu_y);
	A[++m][i] = scl_geom*h_ijklm_p_scl(h_re, 0, 0, 3, 1, 0, 7, nu_x, nu_y);
      }
	
      if (no >= 6) {
	// 3rd order geometric terms
	A[++m][i] = scl_geom*h_ijklm_p_scl(h_re, 5, 0, 0, 0, 0, 7, nu_x, nu_y);
	A[++m][i] = scl_geom*h_ijklm_p_scl(h_re, 4, 1, 0, 0, 0, 7, nu_x, nu_y);
	A[++m][i] = scl_geom*h_ijklm_p_scl(h_re, 3, 2, 0, 0, 0, 7, nu_x, nu_y);
	A[++m][i] = scl_geom*h_ijklm_p_scl(h_re, 3, 0, 2, 0, 0, 7, nu_x, nu_y);
	A[++m][i] = scl_geom*h_ijklm_p_scl(h_re, 2, 1, 2, 0, 0, 7, nu_x, nu_y);
	A[++m][i] = scl_geom*h_ijklm_p_scl(h_re, 1, 0, 4, 0, 0, 7, nu_x, nu_y);
	A[++m][i] = scl_geom*h_ijklm_p_scl(h_re, 3, 0, 1, 1, 0, 7, nu_x, nu_y);
	A[++m][i] = scl_geom*h_ijklm_p_scl(h_re, 2, 1, 1, 1, 0, 7, nu_x, nu_y);
	A[++m][i] = scl_geom*h_ijklm_p_scl(h_re, 1, 0, 3, 1, 0, 7, nu_x, nu_y);
	A[++m][i] = scl_geom*h_ijklm_p_scl(h_re, 3, 0, 0, 2, 0, 7, nu_x, nu_y);
	A[++m][i] = scl_geom*h_ijklm_p_scl(h_re, 2, 1, 0, 2, 0, 7, nu_x, nu_y);
	A[++m][i] = scl_geom*h_ijklm_p_scl(h_re, 1, 0, 2, 2, 0, 7, nu_x, nu_y);
	A[++m][i] = scl_geom*h_ijklm_p_scl(h_re, 1, 0, 0, 4, 0, 7, nu_x, nu_y);
	A[++m][i] = scl_geom*h_ijklm_p_scl(h_re, 1, 0, 1, 3, 0, 7, nu_x, nu_y);
      }

      // amplitude dependent tune shifts
      if (no >= 5) {
	// linear
	A[++m][i] = scl_dnu/n_cell
	            *h_ijklm_p_scl(K_re, 2, 2, 0, 0, 0, 7, nu_x, nu_y);
	A[++m][i] = scl_dnu/n_cell
                    *h_ijklm_p_scl(K_re, 1, 1, 1, 1, 0, 7, nu_x, nu_y);
	A[++m][i] = scl_dnu/n_cell
	            *h_ijklm_p_scl(K_re, 0, 0, 2, 2, 0, 7, nu_x, nu_y);
      }

      if ((no >= 6) && (scl_dnu_delta != 0.0)) {
	A[++m][i] = scl_dnu_delta/n_cell
	            *h_ijklm_p_scl(K_re, 2, 2, 0, 0, 1, 7, nu_x, nu_y);
	A[++m][i] = scl_dnu_delta/n_cell
	            *h_ijklm_p_scl(K_re, 1, 1, 1, 1, 1, 7, nu_x, nu_y);
	A[++m][i] = scl_dnu_delta/n_cell
	            *h_ijklm_p_scl(K_re, 0, 0, 2, 2, 1, 7, nu_x, nu_y);
      }
	
      if (no >= 7) {
	// quadratic
	A[++m][i] = scl_dnu/n_cell
                    *h_ijklm_p_scl(K_re, 3, 3, 0, 0, 0, 7, nu_x, nu_y);
	A[++m][i] = scl_dnu/n_cell
	            *h_ijklm_p_scl(K_re, 1, 1, 2, 2, 0, 7, nu_x, nu_y);
	A[++m][i] = scl_dnu/n_cell
	            *h_ijklm_p_scl(K_re, 2, 2, 1, 1, 0, 7, nu_x, nu_y);
	A[++m][i] = scl_dnu/n_cell
	            *h_ijklm_p_scl(K_re, 0, 0, 3, 3, 0, 7, nu_x, nu_y);

	if (scl_dnu_delta != 0.0) {
	  A[++m][i] = scl_dnu_delta/n_cell
	              *h_ijklm_p_scl(K_re, 2, 2, 0, 0, 2, 7, nu_x, nu_y);
	  A[++m][i] = scl_dnu_delta/n_cell
	              *h_ijklm_p_scl(K_re, 1, 1, 1, 1, 2, 7, nu_x, nu_y);
	  A[++m][i] = scl_dnu_delta/n_cell
	              *h_ijklm_p_scl(K_re, 0, 0, 2, 2, 2, 7, nu_x, nu_y);
	}
      }

      if ((no >= 8) && (scl_dnu_delta != 0.0)) {
	A[++m][i] = scl_dnu_delta/n_cell
                    *h_ijklm_p_scl(K_re, 2, 2, 0, 0, 3, 7, nu_x, nu_y);
	A[++m][i] = scl_dnu_delta/n_cell
	            *h_ijklm_p_scl(K_re, 1, 1, 1, 1, 3, 7, nu_x, nu_y);
	A[++m][i] = scl_dnu_delta/n_cell
	            *h_ijklm_p_scl(K_re, 0, 0, 2, 2, 3, 7, nu_x, nu_y);

	A[++m][i] = scl_dnu_delta/n_cell
	            *h_ijklm_p_scl(K_re, 3, 3, 0, 0, 1, 7, nu_x, nu_y);
	A[++m][i] = scl_dnu_delta/n_cell
	            *h_ijklm_p_scl(K_re, 1, 1, 2, 2, 1, 7, nu_x, nu_y);
	A[++m][i] = scl_dnu_delta/n_cell
	            *h_ijklm_p_scl(K_re, 2, 2, 1, 1, 1, 7, nu_x, nu_y);
	A[++m][i] = scl_dnu_delta/n_cell
	            *h_ijklm_p_scl(K_re, 0, 0, 3, 3, 1, 7, nu_x, nu_y);
      }

      if (no >= 9) {
	// cubic
	A[++m][i] = scl_dnu/n_cell
	            *h_ijklm_p_scl(K_re, 4, 4, 0, 0, 0, 7, nu_x, nu_y);
	A[++m][i] = scl_dnu/n_cell
	            *h_ijklm_p_scl(K_re, 3, 3, 1, 1, 0, 7, nu_x, nu_y);
	A[++m][i] = scl_dnu/n_cell
	            *h_ijklm_p_scl(K_re, 2, 2, 2, 2, 0, 7, nu_x, nu_y);
	A[++m][i] = scl_dnu/n_cell
	            *h_ijklm_p_scl(K_re, 1, 1, 3, 3, 0, 7, nu_x, nu_y);
	A[++m][i] = scl_dnu/n_cell
	            *h_ijklm_p_scl(K_re, 0, 0, 4, 4, 0, 7, nu_x, nu_y);

	if (scl_dnu_delta != 0.0) {
	  A[++m][i] = scl_dnu_delta/n_cell
	              *h_ijklm_p_scl(K_re, 3, 3, 0, 0, 2, 7, nu_x, nu_y);
	  A[++m][i] = scl_dnu_delta/n_cell
	              *h_ijklm_p_scl(K_re, 1, 1, 2, 2, 2, 7, nu_x, nu_y);
	  A[++m][i] = scl_dnu_delta/n_cell
	              *h_ijklm_p_scl(K_re, 2, 2, 1, 1, 2, 7, nu_x, nu_y);
	  A[++m][i] = scl_dnu_delta/n_cell
	              *h_ijklm_p_scl(K_re, 0, 0, 3, 3, 2, 7, nu_x, nu_y);
	}
      }
      
      if ((no >= 5) && (scl_ksi2 != 0.0)) {
	// 2nd order chromaticity
	A[++m][i] = scl_ksi2/n_cell*h_ijklm_p(K_re, 1, 1, 0, 0, 2, 7);
	A[++m][i] = scl_ksi2/n_cell*h_ijklm_p(K_re, 0, 0, 1, 1, 2, 7);
      }

      if ((no >= 6) && (scl_ksi3 != 0.0)) {
	// 3rd order chromaticity
	A[++m][i] = scl_ksi3/n_cell*h_ijklm_p(K_re, 1, 1, 0, 0, 3, 7);
	A[++m][i] = scl_ksi3/n_cell*h_ijklm_p(K_re, 0, 0, 1, 1, 3, 7);
      }
	
      if ((no >= 7) && (scl_ksi4 != 0.0)) {
	// 4th order chromaticity
	A[++m][i] = scl_ksi4/n_cell*h_ijklm_p(K_re, 1, 1, 0, 0, 4, 7);
	A[++m][i] = scl_ksi4/n_cell*h_ijklm_p(K_re, 0, 0, 1, 1, 4, 7);
      }

      // 2nd order momentum compaction
      A[++m][i] = scl_alphac2*h_ijklm_p(Map[ct_], 0, 0, 0, 0, 2, 7)
	          /elem[n_elem-1].s;

      if (m > m_max) {
	cout << "h_zero: m_max exceeded " << m << "(" << m_max << ")" << endl;
	exit(0);
      }

      for (j = 1; j <= get_n_Kids(prms[i-1]); j++)
	clr_bn_par(prms[i-1], j, bns[i-1]);
    }
    
    m = 0;
    if (no >= 4) {
      // linear chromaticity
      b[++m] = -scl_ksi1*(ksi1[X_]*pi*4.0*Jx*delta
			     +h_ijklm(h_re, 1, 1, 0, 0, 1));
      b[++m] = -scl_ksi1*(ksi1[Y_]*pi*4.0*Jy*delta
			     +h_ijklm(h_re, 0, 0, 1, 1, 1));
      
      // 1st order chromatic terms
      b[++m] = -scl_chrom1*h_ijklm(h_re, 2, 0, 0, 0, 1);
      b[++m] = -scl_chrom1*h_ijklm(h_re, 0, 0, 2, 0, 1);
      
      // 2nd order chromatic terms
      b[++m] = -scl_chrom2*h_ijklm(h_re, 1, 0, 0, 0, 2);

      // 1st order geometric terms
      b[++m] = -scl_geom*h_ijklm_scl(h_re, 2, 1, 0, 0, 0, nu_x, nu_y);
      b[++m] = -scl_geom*h_ijklm_scl(h_re, 3, 0, 0, 0, 0, nu_x, nu_y);
      b[++m] = -scl_geom*h_ijklm_scl(h_re, 1, 0, 1, 1, 0, nu_x, nu_y);
      b[++m] = -scl_geom*h_ijklm_scl(h_re, 1, 0, 0, 2, 0, nu_x, nu_y);
      b[++m] = -scl_geom*h_ijklm_scl(h_re, 1, 0, 2, 0, 0, nu_x, nu_y);
    }
    
    if (no >= 5) {
      // 2nd order geometric terms
      b[++m] = -scl_geom*h_ijklm_scl(h_re, 3, 1, 0, 0, 0, nu_x, nu_y);
      b[++m] = -scl_geom*h_ijklm_scl(h_re, 2, 0, 1, 1, 0, nu_x, nu_y);
      b[++m] = -scl_geom*h_ijklm_scl(h_re, 4, 0, 0, 0, 0, nu_x, nu_y);
      b[++m] = -scl_geom*h_ijklm_scl(h_re, 2, 0, 0, 2, 0, nu_x, nu_y);
      b[++m] = -scl_geom*h_ijklm_scl(h_re, 2, 0, 2, 0, 0, nu_x, nu_y);
      b[++m] = -scl_geom*h_ijklm_scl(h_re, 0, 0, 4, 0, 0, nu_x, nu_y);
      b[++m] = -scl_geom*h_ijklm_scl(h_re, 1, 1, 2, 0, 0, nu_x, nu_y);
      b[++m] = -scl_geom*h_ijklm_scl(h_re, 0, 0, 3, 1, 0, nu_x, nu_y);
    }

    if (no >= 6) {
      // 3rd order geometric terms
      b[++m] = -scl_geom*h_ijklm_scl(h_re, 5, 0, 0, 0, 0, nu_x, nu_y);
      b[++m] = -scl_geom*h_ijklm_scl(h_re, 4, 1, 0, 0, 0, nu_x, nu_y);
      b[++m] = -scl_geom*h_ijklm_scl(h_re, 3, 2, 0, 0, 0, nu_x, nu_y);
      b[++m] = -scl_geom*h_ijklm_scl(h_re, 3, 0, 2, 0, 0, nu_x, nu_y);
      b[++m] = -scl_geom*h_ijklm_scl(h_re, 2, 1, 2, 0, 0, nu_x, nu_y);
      b[++m] = -scl_geom*h_ijklm_scl(h_re, 1, 0, 4, 0, 0, nu_x, nu_y);
      b[++m] = -scl_geom*h_ijklm_scl(h_re, 3, 0, 1, 1, 0, nu_x, nu_y);
      b[++m] = -scl_geom*h_ijklm_scl(h_re, 2, 1, 1, 1, 0, nu_x, nu_y);
      b[++m] = -scl_geom*h_ijklm_scl(h_re, 1, 0, 3, 1, 0, nu_x, nu_y);
      b[++m] = -scl_geom*h_ijklm_scl(h_re, 3, 0, 0, 2, 0, nu_x, nu_y);
      b[++m] = -scl_geom*h_ijklm_scl(h_re, 2, 1, 0, 2, 0, nu_x, nu_y);
      b[++m] = -scl_geom*h_ijklm_scl(h_re, 1, 0, 2, 2, 0, nu_x, nu_y);
      b[++m] = -scl_geom*h_ijklm_scl(h_re, 1, 0, 0, 4, 0, nu_x, nu_y);
      b[++m] = -scl_geom*h_ijklm_scl(h_re, 1, 0, 1, 3, 0, nu_x, nu_y);
    }

    // amplitude dependent tune shifts
    if (no >= 5) {
      // linear
      b[++m] = -scl_dnu/n_cell*h_ijklm_scl(K_re, 2, 2, 0, 0, 0, nu_x, nu_y);
      b[++m] = -scl_dnu/n_cell*h_ijklm_scl(K_re, 1, 1, 1, 1, 0, nu_x, nu_y);
      b[++m] = -scl_dnu/n_cell*h_ijklm_scl(K_re, 0, 0, 2, 2, 0, nu_x, nu_y);
    }
      
    if ((no >= 6) && (scl_dnu_delta != 0.0)) {
      b[++m] = -scl_dnu_delta/n_cell
	       *h_ijklm_scl(K_re, 2, 2, 0, 0, 1, nu_x, nu_y);
      b[++m] = -scl_dnu_delta/n_cell
	       *h_ijklm_scl(K_re, 1, 1, 1, 1, 1, nu_x, nu_y);
      b[++m] = -scl_dnu_delta/n_cell
	       *h_ijklm_scl(K_re, 0, 0, 2, 2, 1, nu_x, nu_y);
    }
    
    if (no >= 7) {
      // quadratic
      b[++m] = -scl_dnu/n_cell*h_ijklm_scl(K_re, 3, 3, 0, 0, 0, nu_x, nu_y);
      b[++m] = -scl_dnu/n_cell*h_ijklm_scl(K_re, 1, 1, 2, 2, 0, nu_x, nu_y);
      b[++m] = -scl_dnu/n_cell*h_ijklm_scl(K_re, 2, 2, 1, 1, 0, nu_x, nu_y);
      b[++m] = -scl_dnu/n_cell*h_ijklm_scl(K_re, 0, 0, 3, 3, 0, nu_x, nu_y);
      
      if (scl_dnu_delta != 0.0) {
	b[++m] = -scl_dnu_delta/n_cell
	         *h_ijklm_scl(K_re, 2, 2, 0, 0, 2, nu_x, nu_y);
	b[++m] = -scl_dnu_delta/n_cell
	         *h_ijklm_scl(K_re, 1, 1, 1, 1, 2, nu_x, nu_y);
	b[++m] = -scl_dnu_delta/n_cell
	         *h_ijklm_scl(K_re, 0, 0, 2, 2, 2, nu_x, nu_y);
      }
    }
      
    if ((no >= 8) && (scl_dnu_delta != 0.0)) {
	b[++m] = -scl_dnu_delta/n_cell
	         *h_ijklm_scl(K_re, 2, 2, 0, 0, 3, nu_x, nu_y);
	b[++m] = -scl_dnu_delta/n_cell
	         *h_ijklm_scl(K_re, 1, 1, 1, 1, 3, nu_x, nu_y);
	b[++m] = -scl_dnu_delta/n_cell
	         *h_ijklm_scl(K_re, 0, 0, 2, 2, 3, nu_x, nu_y);

	b[++m] = -scl_dnu_delta/n_cell
	         *h_ijklm_scl(K_re, 3, 3, 0, 0, 1, nu_x, nu_y);
	b[++m] = -scl_dnu_delta/n_cell
	         *h_ijklm_scl(K_re, 1, 1, 2, 2, 1, nu_x, nu_y);
	b[++m] = -scl_dnu_delta/n_cell
	         *h_ijklm_scl(K_re, 2, 2, 1, 1, 1, nu_x, nu_y);
	b[++m] = -scl_dnu_delta/n_cell
	         *h_ijklm_scl(K_re, 0, 0, 3, 3, 1, nu_x, nu_y);
    }

    if (no >= 9) {
      // cubic
      b[++m] = -scl_dnu/n_cell*h_ijklm_scl(K_re, 4, 4, 0, 0, 0, nu_x, nu_y);
      b[++m] = -scl_dnu/n_cell*h_ijklm_scl(K_re, 3, 3, 1, 1, 0, nu_x, nu_y);
      b[++m] = -scl_dnu/n_cell*h_ijklm_scl(K_re, 2, 2, 2, 2, 0, nu_x, nu_y);
      b[++m] = -scl_dnu/n_cell*h_ijklm_scl(K_re, 1, 1, 3, 3, 0, nu_x, nu_y);
      b[++m] = -scl_dnu/n_cell*h_ijklm_scl(K_re, 0, 0, 4, 4, 0, nu_x, nu_y);
      
      if (scl_dnu_delta != 0.0) {
	b[++m] = -scl_dnu_delta/n_cell
	         *h_ijklm_scl(K_re, 3, 3, 0, 0, 2, nu_x, nu_y);
	b[++m] = -scl_dnu_delta/n_cell
	         *h_ijklm_scl(K_re, 1, 1, 2, 2, 2, nu_x, nu_y);
	b[++m] = -scl_dnu_delta/n_cell
	         *h_ijklm_scl(K_re, 2, 2, 1, 1, 2, nu_x, nu_y);
	b[++m] = -scl_dnu_delta/n_cell
	         *h_ijklm_scl(K_re, 0, 0, 3, 3, 2, nu_x, nu_y);
      }
    }

    if ((no >= 5) && (scl_ksi2 != 0.0)) {
      // 2nd order chromaticity
      b[++m] = -scl_ksi2/n_cell*h_ijklm(K_re, 1, 1, 0, 0, 2);
      b[++m] = -scl_ksi2/n_cell*h_ijklm(K_re, 0, 0, 1, 1, 2);
    }
    
    if ((no >= 6) && (scl_ksi3 != 0.0)) {
      // 3rd order chromaticity
      b[++m] = -scl_ksi3/n_cell*h_ijklm(K_re, 1, 1, 0, 0, 3);
      b[++m] = -scl_ksi3/n_cell*h_ijklm(K_re, 0, 0, 1, 1, 3);
    }
    
    if ((no >= 7) && (scl_ksi4 != 0.0)) {
      // 4th order chromaticity
      b[++m] = -scl_ksi4/n_cell*h_ijklm(K_re, 1, 1, 0, 0, 4);
      b[++m] = -scl_ksi4/n_cell*h_ijklm(K_re, 0, 0, 1, 1, 4);
    }
    
    // 2nd order momentum compaction
    b[++m] = -scl_alphac2*h_ijklm(Map[ct_], 0, 0, 0, 0, 2)/elem[n_elem-1].s;
	
    if (lin_chrom) m = 2; // only correct linear chromaticity

    if (prt) {
      cout << endl;
      cout << n << " Ax = b:" << endl;
      cout << endl;
      for (i = 1; i <= m; i++) {
	cout  << setw(3) << i;
	for (j = 1; j <= n_prm; j++)
	  cout << scientific << setprecision(2) << setw(10) << A[i][j];
	cout << scientific << setprecision(2) << setw(10) << b[i] << endl;
      }
    }
    
    if (false) {
      for (i = 1; i <= m; i++)
	for (j = 1; j <= n_prm; j++)
	  U[i][j] = A[i][j];

      dsvdcmp(U, m, n_prm, w, V);

      s_max = -1e30;
      for (i = 1; i <= n_prm; i++)
	s_max = max(w[i], s_max);
      
      if (prt) {
	cout << endl;
	cout << "singular values:" << endl;
	cout << endl;
	for (i = 1; i <= n_prm; i++) {
	  cout << scientific << setprecision(3) << setw(11) << w[i];
	  if (w[i]/s_max < s_cut) {
	    w[i] = 0.0;
	    cout << " (zeroed)";
	  }
	  if (i % 8 == 0) cout << endl;
	}
	if (n_prm % 8 != 0) cout << endl;
      }

      dsvbksb(U, w, V, m, n_prm, b, dbn);
    } else
      SVD_lim(m, n_prm, A, b, bn_max, 1e-10, bn, dbn);

    if (prt) {
      m = 0;
      cout << endl;
      if (no >= 4) {
	cout << "linear chromaticity:" << endl;
	cout << scientific << setprecision(3)
	     << setw(11) << b[++m] << setw(11) << b[++m] << endl;
	cout << "1st order chromatic terms:" << endl;
	cout << scientific << setprecision(3)
	     << setw(11) << b[++m] << setw(11) << b[++m] << endl;
	cout << "2nd order chromatic terms:" << endl;
	cout << scientific << setprecision(3)
	     << setw(11) << b[++m] << endl;
	cout << "1st order geometric terms:" << endl;
	for (i = 1; i <= 5; i++)
	  cout << scientific << setprecision(3) << setw(11) << b[++m];
	cout << endl;
      }
      if (no >= 5) {
	cout << "2nd order geometric terms:" << endl;
	for (i = 1; i <= 8; i++)
	  cout << scientific << setprecision(3) << setw(11) << b[++m];
	cout << endl;
      }
      if (no >= 6) {
	cout << "3rd order geometric terms:" << endl;
	for (i = 1; i <= 14; i++) {
	  cout << scientific << setprecision(3) << setw(11) << b[++m];
	  if (i == 8) cout << endl;
	}
	cout << endl;
      }
      if (no >= 5) {
	cout << endl;
	cout << "amplitude dependent tune shifts" << endl;
	cout << "linear:" << endl;
	for (i = 1; i <= 3; i++)
	  cout << scientific << setprecision(3) << setw(11) << b[++m];
	cout << endl;
      }
      if ((no >= 6) && (scl_dnu_delta != 0.0)) {
	for (i = 1; i <= 3; i++)
	  cout << scientific << setprecision(3) << setw(11) << b[++m];
	cout << endl;
      }
      if (no >= 7) {
	cout << "quadratic:" << endl;
	for (i = 1; i <= 4; i++)
	  cout << scientific << setprecision(3) << setw(11) << b[++m];
	cout << endl;

	if (scl_dnu_delta != 0.0) {
	  for (i = 1; i <= 3; i++)
	    cout << scientific << setprecision(3) << setw(11) << b[++m];
	  cout << endl;
	}
      }
      if ((no >= 8) && (scl_dnu_delta != 0.0)) {
	  for (i = 1; i <= 7; i++)
	    cout << scientific << setprecision(3) << setw(11) << b[++m];
	  cout << endl;
      }
      if (no >= 9) {
	cout << "cubic:" << endl;
	for (i = 1; i <= 5; i++)
	  cout << scientific << setprecision(3) << setw(11) << b[++m];
	cout << endl;

	if (scl_dnu_delta != 0.0) {
	  for (i = 1; i <= 4; i++)
	    cout << scientific << setprecision(3) << setw(11) << b[++m];
	  cout << endl;
	}
      }
      if ((no >= 5) && (scl_ksi2 != 0.0)) {
	cout << "2nd order chromaticity:" << endl;
	for (i = 1; i <= 2; i++)
	  cout << scientific << setprecision(3) << setw(11) << b[++m];
	cout << endl;
      }
      if ((no >= 6) && (scl_ksi3 != 0.0)) {
	cout << "3rd order chromaticity:" << endl;
	for (i = 1; i <= 2; i++)
	  cout << scientific << setprecision(3) << setw(11) << b[++m];
	cout << endl;
      }
      if ((no >= 7) && (scl_ksi4 != 0.0)) {
	cout << "4th order chromaticity:" << endl;
	for (i = 1; i <= 2; i++)
	  cout << scientific << setprecision(3) << setw(11) << b[++m];
	cout << endl;
      }
      cout << "2nd order momentum compaction:" << endl;
      cout << scientific << setprecision(3) << setw(11) << b[++m] << endl;
   }
    
    for (i = 1; i <= n_prm; i++)
      for (j = 1; j <= get_n_Kids(prms[i-1]); j++)
	set_dbn(prms[i-1], j, bns[i-1], dbn[i]);
    bn[i] = get_bn(prms[i-1], 1, bns[i-1]);

    if (prt_iter || (n == n_iter)) {
      if (prt_iter) {
	sext_out << endl;
	sext_out << "n = " << n << ":" << endl;
      }
      for (i = 1; i <= n_prm; i++)
	for (j = 1; j <= get_n_Kids(prms[i-1]); j++)
	  sext_out << fixed << setprecision(7) 
		   << setw(6) << get_Name(prms[i-1]) << "(" << j << ") = "
		   << setw(11) << get_bnL(prms[i-1], 1, bns[i-1])
		   << setw(2) << bns[i-1] << endl;
    }
  } while (n < n_iter);

  if (prt_iter) {
    danot_(no-1); get_Map_N();
    danot_(no);
    h = get_h(); h = h*Id; CtoR(h, h_re, h_im);
    K = MapNorm(Map, g, A1, A0, Map_res, 1); K = K*Id; CtoR(K, K_re, K_im);
    file_wr(h_out, "h.dat"); h_out << h_re; h_out.close();
    file_wr(K_out, "K.dat"); K_out << K_re; K_out.close();
  }

  free_dvector(b, 1, m_max); free_dvector(w, 1, n_prm);
  free_dvector(dbn, 1, n_prm); free_dvector(bn, 1, n_prm);
  free_dvector(bn_max, 1, n_prm);
  free_dmatrix(A, 1, m_max, 1, n_prm);
  free_dmatrix(U, 1, m_max, 1, n_prm); free_dmatrix(V, 1, n_prm, 1, n_prm);
  free_dmatrix(A_inv, 1, n_prm, 1, m_max);
}


