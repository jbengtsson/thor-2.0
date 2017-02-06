/* Author:  Johan Bengtsson  */


long int          n_elem;
elem_type<double> elem[max_elem];
elem_type<tps>    elem_tps[max_elem];
Family            Families[max_Family];

bool lost;

double cl_rad, q_fluct;
double d_coeff1, d_coeff2, k_coeff1, k_coeff2;
double E0, dE, beta0, gamma0;
tps    I2, I4, I5, dcurly_H, dI4;
tps    D_[3]; // diff. coeff. for the linear invarient

bool rad_on    = false, H_exact        = false, totpath_on   = false;
bool cavity_on = false, quad_fringe_on = false, emittance_on = false;
bool IBS_on    = false;


/* initialize symplectic integrator */
void ini_si(void)
{

  std::cout << std::endl;
  std::cout << "ini_si" << std::endl;

  // compute 4th order symplectic integrator coefficients
  d_coeff1 = 1.0/(2.0*(2.0-pow(2.0, 1.0/3.0))); d_coeff2 = 0.5 - d_coeff1;
  k_coeff1 = 2.0*d_coeff1; k_coeff2 = 1.0 - 2.0*k_coeff1;
  // compute radiaton coefficients
  cl_rad = C_gamma*pow(E0, 3)/(2.0*pi);
//  q_fluct = 3.0*cu*c_gamma*hbc*pow(E0, 5)/(4.0*pi*pow(xmc2, 3.0));
  q_fluct =
    3.0*C_u*C_gamma*1e-9*h_bar*clight/(4.0*M_PI*cube(1e-9*m_e))*pow(E0, 5.0);
}

/* 2-dim rotation */
template<typename T>
void rot(double c, double s, ss_vect<T> &x)
{
  T  x1, px1;

  x1 = x[x_];
  x[x_] = c*x1 + s*x[y_]; x[y_] = -s*x1 + c*x[y_];
  px1 = x[px_];
  x[px_] = c*px1 + s*x[py_]; x[py_] = -s*px1 + c*x[py_];
}

/* Euclidian transformation (global to local coordinates) */
template<typename T>
void gtol(const elem_type<T> &elem, ss_vect<T> &x)
{
  x[px_] += elem.c1; x[py_] += elem.s1;
  x[x_] -= elem.dx[X_]; x[y_] -= elem.dx[Y_];
  rot(elem.droll[X_], elem.droll[Y_], x);
  x[px_] -= elem.c0;
}

/* Euclidian transformation (local to global coordinates) */
template<typename T>
void ltog(const elem_type<T> &elem, ss_vect<T> &x)
{
  x[px_] -= elem.c0;
  rot(elem.droll[X_], -elem.droll[Y_], x);
  x[x_] += elem.dx[X_]; x[y_] += elem.dx[Y_];
  x[px_] += elem.c1; x[py_] += elem.s1;
}

#if 1

template<typename T>
inline T get_ps(ss_vect<T> &x)
{
  T  ps, ps2;

  if (!H_exact)
    ps = 1.0 + x[delta_];
  else {
    ps2 = sqr(1.0+x[delta_]) - sqr(x[px_]) - sqr(x[py_]);
    if (ps2 >= 0.0)
      ps = sqrt(ps2);
    else {
      std::cout << "*** Speed of light exceeded!\n" << std::endl;
      exit(1);
    }
  }
  return ps;
}

#else

template<typename T>
inline T get_ps(const ss_vect<T> &ps)
{
  T p_s, p_s2;

  if (!H_exact)
    // Ultra relatistic approximation [cT, delta] -> [ct, p_t].
    p_s2 = sqr(1e0+ps[delta_]) - sqr(ps[px_]) - sqr(ps[py_]);
  else
    // delta -> p_t.
    p_s2 =
      1e0 + 2e0*ps[delta_]/beta0 + sqr(ps[delta_]) - sqr(ps[px_])
      - sqr(ps[py_]);
  if (p_s2 >= 0e0)
    p_s = sqrt(p_s2);
  else {
    printf("get_p_s: *** Speed of light exceeded!\n");
    p_s = NAN;
  }
  return(p_s);
}

#endif

// partial template-class specialization
// primary version
template<typename T>
class is_tps { };


// partial specialization
template<>
class is_tps<double> {
 public:
  static inline void get_ps(const ss_vect<double> &x,
			    elem_type<double> &elem) { }

  static inline double get_curly_H(const ss_vect<tps> &x)
    {
      std::cout << "get_curly_H: operation not defined for double" << std::endl;
      exit(1);
      return 0e0;
    }

  static inline double get_dI4(const double h, const double b2, const double L,
			       const ss_vect<tps> &x)
  {
    std::cout << "get_dI4: operation not defined for double" << std::endl;
    exit(1);
    return 0e0;
  }

  static inline void emittance(const double B2, const double u,
			       const double ps0, const ss_vect<double> &xp) { }
};


// partial specialization
template<>
class is_tps<tps> {
 public:
  static inline void get_ps(const ss_vect<tps> &x, elem_type<tps> &elem)
  { elem.A1 = x; }

  static inline tps get_curly_H(const ss_vect<tps> &A)
  {
    int          j;
    tps          curly_H[2];
    ss_vect<tps> eta;

    eta.zero();
    for (j = 0; j < 4; j++) {
      // Include parameter dependance.
      eta[j] = A[j][delta_] + h_ijklm_p(A[j], 0, 0, 0, 0, 1, 7)*tps(0e0, 7);
    }
 
    get_twoJ(2, eta, A, curly_H);

    return curly_H[X_];
  }

  static inline tps get_dI4(const ss_vect<tps> &A)
  {
    return A[x_][delta_] + h_ijklm_p(A[x_], 0, 0, 0, 0, 1, 7)*tps(0e0, 7);
  }

  static inline void emittance(const tps &B2, const tps &H_dL, const tps &ps0,
			       const ss_vect<tps> &A) {
    int           j;
    double        dens0;
    tps           v;
    ss_vect<tps>  A_inv;

    // Input is: A = [Re(v_i) Im(v_i), , ], i = 1, 2, 3.
    if (B2 > 0.0) {
      v = q_fluct*pow(B2.cst(), 1.5)*pow(ps0, 4)*H_dL; dens0 = v.cst();
      // compute "A^-1"
      A_inv = Inv(A);
      /* Add contribution to diffusion coeffs. of the invariants.  Diffusion
	 only occurs in the momentum components of the invariants.
	 The invariants are computed from A:
           J_i = z^T G_i z, i = 1, 2, 3  (quadratic form => 2nd rank tensor)
	   G_iab = (A_2i,a)^-1 (A_2i,b)^-1 + (A_2i+1,a^-1 (A_2i+1,b)^-1  */
      for (j = 0; j < 3; j++)
	D_[j] += (sqr(h_ijklm(A_inv[j*2], 0, 0, 0, 0, 1)
		      +h_ijklm_p(A_inv[j*2], 0, 0, 0, 0, 1, 7)*tps(0.0, 7))
		  +sqr(h_ijklm(A_inv[j*2+1], 0, 0, 0, 0, 1)
		       +h_ijklm_p(A_inv[j*2+1], 0, 0, 0, 0, 1, 7)*tps(0.0, 7)))
	  *dens0;
    }
  }
};


template<typename T>
T B2_perp(const double h_ref, const T B[], const ss_vect<T> u)
{
/*            -   - 2        -
   Calculate |B x e| , where e is a unit vector in the direction of
   propagation                                                       */
    
  T  x_n, e[3], B2;

  x_n = 1.0/(sqrt(sqr(1.0+u[x_]*h_ref)+sqr(u[px_])+sqr(u[py_])));
  e[X_] = u[px_]*x_n; e[Y_] = u[py_]*x_n; e[Z_] = (1.0+u[x_]*h_ref)*x_n;
  B2 = sqr(B[Y_]*e[Z_]-B[Z_]*e[Y_]) + sqr(B[X_]*e[Y_]-B[Y_]*e[X_])
       + sqr(B[Z_]*e[X_]-B[X_]*e[Z_]);

  return B2;
}


template<typename T, typename U>
void radiate(ss_vect<T> &x, const U L, double h_ref, T B[])
{
  T          ps0, ps1, H_dL, B2;
  ss_vect<T> u;

  // large ring: x' and y' approx. const. wrt radiation:  p_x = x'*(1+delta)
  u = x; ps0 = get_ps(x); u[px_] /= ps0; u[py_] /= ps0;

  // Note, H = -p_s => integral(H(s)) = path length
  H_dL = (1.0+u[x_]*h_ref+(sqr(u[px_])+sqr(u[py_]))/2.0)*L;
  B2 = B2_perp(h_ref, B, u);

  if (rad_on) {
    x[delta_] -= cl_rad*sqr(ps0)*B2*H_dL;
    ps1 = get_ps(x); x[px_] = u[px_]*ps1; x[py_] = u[py_]*ps1;
  }

  if (emittance_on) is_tps<T>::emittance(B2, H_dL, ps0, u);
}

template<typename T, typename U>
void thin_kick(const int order, const T an[], const T bn[], const T L,
	       const double h_bend, const double h_ref, const int thick,
	       ss_vect<U> &x)
{

/* The kick is given by

              e L       L delta    L x              e L
     Dp_x = - --- B_y + ------- - ----- ,    Dp_y = --- B_x
              p_0         rho     rho^2             p_0

    where

                           ====
                           \ 
      (B_y + iB_x) = B rho  >   (ia_n  + b_n ) (x + iy)^n-1
                           /
                           ====

    where

       e      1
      --- = -----
      p_0   B rho                                            */

  int	      j;
  U           BxoBrho, ByoBrho, ByoBrho1;
  U           B[3];
  ss_vect<U>  x0;

  x0 = x;
  if ((h_bend != 0.0) || ((1 <= order) && (order <= mpole_max))) {
    // compute field with Horner's rule
    if (order > 0) {
      ByoBrho = bn[order-1]; BxoBrho = an[order-1];
    } else {
      ByoBrho = 0.0; BxoBrho = 0.0;
    }
    for (j = order-1; j >= 1; j--) {
      ByoBrho1 = x0[x_]*ByoBrho - x0[y_]*BxoBrho + bn[j-1];
      BxoBrho = x0[y_]*ByoBrho + x0[x_]*BxoBrho + an[j-1]; ByoBrho = ByoBrho1;
    }

    if (rad_on || emittance_on) {
      B[X_] = BxoBrho; B[Y_] = ByoBrho + h_bend; B[Z_] = 0.0;
      radiate(x, L, h_ref, B);
    }

    if (h_ref != 0.0) {
      x[px_] -= L*(ByoBrho+(h_bend-h_ref)/2.0+h_ref*h_bend*x0[x_]
                -h_ref*x0[delta_]);
      x[ct_] += L*h_ref*x0[x_];
    } else
      x[px_] -= L*(ByoBrho+h_bend);
    x[py_] += L*BxoBrho;
  }
}

#if 1

template<typename T, typename U>
void drift_pass(const T L, ss_vect<U> &x)
{
  U u;

  if (!H_exact) {
    u = L/(1.0+x[delta_]);
    x[ct_] += u*(sqr(x[px_])+sqr(x[py_]))/(2.0*(1.0+x[delta_]));
  } else {
    u = L/get_ps(x); x[ct_] += u*(1.0+x[delta_]) - L;
  }
  x[x_] += x[px_]*u; x[y_] += x[py_]*u;
  if (totpath_on) x[ct_] += L;
}

#else

template <class T>
void drift_pass(const T L, ss_vect<T> &ps)
{
  T u, p_s, delta1;

    // [ct, p_t].
  if (false) {
    // Ultra relatistic approximation [cT, delta] -> [ct, p_t].
    u = L/get_ps(ps);
    ps[x_] += ps[px_]*u;
    ps[y_] += ps[py_]*u;
    ps[ct_] += u*(1e0+ps[delta_]) - L;
  } else {
    // delta -> p_t.
    p_s = sqrt(1e0+2e0*ps[delta_]/beta0+sqr(ps[delta_]));
    delta1 = p_s - 1e0;

    p_s = sqrt(sqr(1e0+delta1)-sqr(ps[px_])-sqr(ps[py_]));
    u = L/p_s;

    u = L/get_ps(ps);
    ps[x_] += ps[px_]*u;
    ps[y_] += ps[py_]*u;
    ps[ct_] += L*(1e0/beta0+ps[delta_])/p_s;
  }
  if (totpath_on) ps[ct_] += L;
}

#endif

/* phenomenological correction for magnet gap */
double get_psi(double h_bend, double phi, double gap)
{
  double  psi;

  const double k1 = 0.5, k2 = 0.0;

  if (phi == 0.0)
    psi = 0.0;
  else
    psi = k1*gap*h_bend*(1.0+sqr(sin(dtor(phi))))/cos(dtor(phi))
          *(1.0-k2*gap*h_bend*tan(dtor(phi)));

  return psi;
}


/* dipole fringe field using thin quadrupole fudge */
template<typename T>
void bend_HE_fringe(double h_bend, double phi, double gap, ss_vect<T> &x)
{
  x[px_] += h_bend*tan(dtor(phi))*x[x_];
  if (false)
    // warning: => diverging Taylor map (see SSC-141)
    x[py_] -= h_bend*tan(dtor(phi)-get_psi(h_bend, phi, gap))*x[y_]
              /(1.0+x[delta_]);
  else
    x[py_] -= h_bend*tan(dtor(phi)-get_psi(h_bend, phi, gap))*x[y_];
}


template<typename T>
void p_rot(double phi, ss_vect<T> &x)
{
  double      c, s;
  T           ps, p;
  ss_vect<T>  x1;

  c = cos(dtor(phi)); s = sin(dtor(phi));
  x1 = x; ps = get_ps(x); p = c*ps - s*x1[px_];
  x[x_] = x1[x_]*ps/p; x[px_] = s*ps + c*x1[px_];
  x[y_] += x1[x_]*x1[py_]*s/p;
  x[ct_] += (1.0+x1[delta_])*x1[x_]*s/p;
}


template<typename T>
void bend_fringe(double h_bend, ss_vect<T> &x)
{
  /* The Lie generator is (SSC-141)

              1                 y^2 p_x
     h = -/+ ---- ----------------------------------- ,
             2rho sqrt( (1+delta)^2 - p_x^2 - p_y^2 )

              -   -      -
     exp(:h:) x = x + :h:x + ...
  
  */

  double      coeff;
  T           u, ps, ps2, ps3;
  ss_vect<T>  x1;

  coeff = -h_bend/2.0;
  x1 = x; ps = get_ps(x); ps2 = sqr(ps); ps3 = ps2*ps;
  u = 1.0 + 4.0*coeff*x1[px_]*x1[y_]*x1[py_]/ps3;
  if (u >= 0.0) {
    x[y_] = 2.0*x1[y_]/(1.0+sqrt(u));
    x[x_] = x1[x_] - coeff*sqr(x1[y_])*(ps2+sqr(x1[px_]))/ps3;
    x[py_] = x1[py_] + 2.0*coeff*x1[px_]*x1[y_]/ps;
    x[ct_] = x1[ct_] - coeff*x1[px_]*sqr(x1[y_])*(1.0+x1[delta_])/ps3;
  } else {
    std::cout << "bend_fringe: *** Speed of light exceeded!" << std::endl;
    exit(0);
  }
}


/* quadrupole fringe field */
template<typename T>
void quad_fringe(T b2, ss_vect<T> &x)
{
  T  u, ps;

  u = b2/(12.0*(1.0+x[delta_])); ps = u/(1.0+x[delta_]);
  x[py_] /= 1.0 - 3.0*u*sqr(x[y_]); x[y_] -= u*cube(x[y_]);
  if (cavity_on) x[ct_] -= ps*cube(x[y_])*x[py_];
  x[px_] /= 1.0 + 3.0*u*sqr(x[x_]);
  if (cavity_on) x[ct_] += ps*cube(x[x_])*x[px_];
  x[x_] += u*cube(x[x_]); u = u*3.0; ps = ps*3.0;
  x[y_] = exp(-u*sqr(x[x_]))*x[y_]; x[py_] = exp(u*sqr(x[x_]))*x[py_];
  x[px_] += 2.0*u*x[x_]*x[y_]*x[py_];
  if (cavity_on) x[ct_] -= ps*sqr(x[x_])*x[y_]*x[py_];
  x[x_] = exp(u*sqr(x[y_]))*x[x_]; x[px_] = exp(-u*sqr(x[y_]))*x[px_];
  x[py_] -= 2.0*u*x[y_]*x[x_]*x[px_];
  if (cavity_on) x[ct_] += ps*sqr(x[y_])*x[x_]*x[px_];
}


/* multipole */
template<typename T>
void mpole_pass(const elem_type<T> &elem, ss_vect<T> &x)
{
  int     i;
  double  h_ref;
  T       L0, L1, L2, k1, k2;

  // Note, expanded
  gtol(elem, x);

  /* fringe fields */
  if (quad_fringe_on && (elem.mpole->bn[Quad-1] != 0.0))
    quad_fringe(elem.mpole->bn[Quad-1], x);

  if (!H_exact) {
    if (elem.mpole->h_bend != 0.0)
      bend_HE_fringe(elem.mpole->h_bend, elem.mpole->edge1,
		     elem.mpole->gap, x);
  } else {
    p_rot(elem.mpole->edge1, x); bend_fringe(elem.mpole->h_bend, x);
  }

  if (elem.L != 0.0) {
    // if (!H_exact ||
    // 	((elem.mpole->edge1 == 0.0) && (elem.mpole->edge2 == 0.0))) {
    if (!H_exact) {
      // polar coordinates
      h_ref = elem.mpole->h_bend; L0 = elem.L/elem.mpole->n_step;
    } else {
      // Cartesian coordinates
      h_ref = 0.0;
      if (elem.mpole->h_bend == 0.0)
	L0 = elem.L/elem.mpole->n_step;
      else
	// L0 = 1.0/elem.mpole->h_bend*(sin(dtor(elem.mpole->edge1))
	//      + sin(dtor(elem.mpole->edge2)))/elem.mpole->n_step;
	L0 =
	  2.0/elem.mpole->h_bend*sin(elem.L*elem.mpole->h_bend/2.0)
	  /elem.mpole->n_step;
    }
    switch (elem.mpole->method) {
    case Linear:
      std::cout << "*** matrix method not supported: " << elem.Name
		<< std::endl;
      exit(1);
      break;
    case Second:
      std::cout << "*** second order symplectic integrator not supported"
		<< std::endl;
      exit(1);
      break;
    case Fourth:
      L1 = d_coeff1*L0; L2 = d_coeff2*L0; k1 = k_coeff1*L0; k2 = k_coeff2*L0;
      dcurly_H = 0e0; dI4 = 0e0;
      for (i = 0; i < elem.mpole->n_step; i++) {
	if (emittance_on && (!cavity_on) && (elem.mpole->h_bend != 0e0)) {
	  dcurly_H += is_tps<tps>::get_curly_H(x);
	  dI4 += is_tps<tps>::get_dI4(x);
	}

	drift_pass(L1, x);
	thin_kick(elem.mpole->order, elem.mpole->an, elem.mpole->bn, k1,
		  elem.mpole->h_bend, h_ref, true, x);
	drift_pass(L2, x);
	thin_kick(elem.mpole->order, elem.mpole->an, elem.mpole->bn, k2,
		  elem.mpole->h_bend, h_ref, true, x);

	if (emittance_on && (!cavity_on) && (elem.mpole->h_bend != 0e0)) {
	  dcurly_H += 4e0*is_tps<tps>::get_curly_H(x);
	  dI4 += 4e0*is_tps<tps>::get_dI4(x);
	}

	drift_pass(L2, x);
	thin_kick(elem.mpole->order, elem.mpole->an, elem.mpole->bn, k1,
		  elem.mpole->h_bend, h_ref, true, x);
	drift_pass(L1, x);

	if (emittance_on && (!cavity_on) && (elem.mpole->h_bend != 0e0)) {
	  dcurly_H += is_tps<tps>::get_curly_H(x);
	  dI4 += is_tps<tps>::get_dI4(x);
	}
      }

      if (emittance_on && (!cavity_on) && (elem.mpole->h_bend != 0e0)) {
	dcurly_H /= 6e0*elem.mpole->n_step;
	dI4 *=
	  elem.mpole->h_bend
	  *(sqr(elem.mpole->h_bend)+2e0*elem.mpole->bn[Quad-1])
	  /(6e0*elem.mpole->n_step);
	I2 += elem.L*sqr(elem.mpole->h_bend); I4 += elem.L*dI4;
	I5 += elem.L*fabs(cube(elem.mpole->h_bend))*dcurly_H;
      }
      break;
    default:
      std::cout << "mpole_pass: undefined element %d\n" << elem.kind
		<< std::endl;
    }
  } else {
    k1 = 1.0;
    thin_kick(elem.mpole->order, elem.mpole->an, elem.mpole->bn,
	      k1, 0.0, 0.0, false, x);
  }

  /* fringe fields */
  if (!H_exact) {
    if (elem.mpole->h_bend != 0.0)
      bend_HE_fringe(elem.mpole->h_bend, elem.mpole->edge2,
		     elem.mpole->gap, x);
  } else {
    bend_fringe(-elem.mpole->h_bend, x); p_rot(elem.mpole->edge2, x);
  }
  
  if (quad_fringe_on && (elem.mpole->bn[Quad-1] != 0.0))
    quad_fringe(-elem.mpole->bn[Quad-1], x);

  // Note, expanded
  ltog(elem, x);
}


template<typename T>
void get_Axy(const wiggler_type *W, const T z, const ss_vect<T> &x,
	     T &AxoBrho, T &AyoBrho, T AxoBrhop[], T AyoBrhop[])
{
  int     i;
  double  ky, kz_n;
  T       cx, cz, sx, sz, chy, shy;
  
  AxoBrho = 0.0; AyoBrho = 0.0;

  for (i = 0; i < 3; ++i) {
    AxoBrhop[i] = 0.0; AyoBrhop[i] = 0.0;
  }

  for (i = 0; i < W->n_harm; i ++) {
    kz_n = W->harm[i]*2.0*M_PI/W->lambda; ky = sqrt(sqr(W->kxV[i])+sqr(kz_n));

    cx = cos(W->kxV[i]*x[x_]); sx = sin(W->kxV[i]*x[x_]);
    chy = cosh(ky*x[y_]); shy = sinh(ky*x[y_]); sz = sin(kz_n*z);

    AxoBrho += W->BoBrhoV[i]/kz_n*cx*chy*sz;
    AyoBrho += W->BoBrhoV[i]*W->kxV[i]/(ky*kz_n)*sx*shy*sz;

    // derivatives with respect to x
    AxoBrhop[X_] -= W->BoBrhoV[i]*W->kxV[i]/kz_n*sx*chy*sz;
    AyoBrhop[X_] += W->BoBrhoV[i]*sqr(W->kxV[i])/(ky*kz_n)*cx*shy*sz;

    // derivatives with respect to y
    AxoBrhop[Y_] += W->BoBrhoV[i]*ky/kz_n*cx*shy*sz;
    AyoBrhop[Y_] += W->BoBrhoV[i]*W->kxV[i]/kz_n*sx*chy*sz;

    if (rad_on) {
      cz = cos(kz_n*z);
      // derivatives with respect to z
      AxoBrhop[Z_] += W->BoBrhoV[i]*cx*chy*cz;
      AyoBrhop[Z_] += W->BoBrhoV[i]*W->kxV[i]/ky*sx*shy*cz;
    }
  }
}

template<typename T>
void wiggler_pass(const elem_type<T> &elem, ss_vect<T> &x)
{
  // first order symplectic integrator for wiggler using expanded Hamiltonian

  int         i;
  T           h, z, a11, a12, a21, a22, c11, c12, c21, c22, det, d1, d2, x2;
  T           B[3], AxoBrho, AyoBrho, AxoBrhop[3], AyoBrhop[3];
  T           hops0, ps0;
  
  h = elem.L/elem.wiggler->n_step; z = 0.0;
  for (i = 1; i <= elem.wiggler->n_step; ++i) {
    get_Axy(elem.wiggler, z, x, AxoBrho, AyoBrho, AxoBrhop, AyoBrhop);
    ps0 = 1.0 + x[delta_]; hops0 = h/ps0;
    a11 = hops0*AxoBrhop[X_]; a12 = hops0*AyoBrhop[X_];
    a21 = hops0*AxoBrhop[Y_]; a22 = hops0*AyoBrhop[Y_];
    det = 1.0 - a11 - a22 + a11*a22 - a12*a21;
    d1 = hops0*AxoBrho*AxoBrhop[X_]; d2 = hops0*AxoBrho*AxoBrhop[Y_];
    c11 = (1.0-a22)/det; c12 = a12/det; c21 = a21/det; c22 = (1.0-a11)/det;
    x2 = c11*(x[px_]-d1) + c12*(x[py_]-d2);
    x[py_] = c21*(x[px_]-d1) + c22*(x[py_]-d2); x[px_] = x2;
    x[x_] += hops0*(x[px_]-AxoBrho); x[y_] += hops0*x[py_];
    d1 = (x[px_]-AxoBrho)/ps0; d2 = (x[py_]-AyoBrho)/ps0;
    x[ct_] += h*(sqr(d1)+sqr(d2))/2.0;

    if (totpath_on) x[ct_] += h;

    if (rad_on || emittance_on) {
      B[X_] = -AyoBrhop[Z_]; B[Y_] = AxoBrhop[Z_];
      B[Z_] = AyoBrhop[X_] - AxoBrhop[Y_];
      radiate(x, h, 0.0, B);
    }

    z += h;
  }
}


template<typename T>
void get_Axy_Wu(const wiggler_type *W, const T z, const ss_vect<T> &x,
		T &AoBrho, T &dp, const bool hor) 
{
  int     i;
  double  ky, kz_n;
  T       cx, sx, sz, chy, shy;

  AoBrho = 0.0; dp = 0.0;

  for (i = 0; i < W->n_harm; i++) {
    kz_n = W->harm[i]*2.0*M_PI/W->lambda; ky = sqrt(sqr(W->kxV[i])+sqr(kz_n));

    cx  = cos(W->kxV[i]*x[x_]); sx = sin(W->kxV[i]*x[x_]);
    chy = cosh(ky*x[y_]); shy = sinh(ky*x[y_]); sz = sin(kz_n*z);

    if (hor) {
      // A_x/Brho
      AoBrho += W->BoBrhoV[i]/kz_n*cx*chy*sz;
      // dp_y
      if (W->kxV[i] == 0.0)
	dp += W->BoBrhoV[i]/kz_n*ky*x[x_]*shy*sz;
      else
	dp += W->BoBrhoV[i]/(W->kxV[i]*kz_n)*ky*sx*shy*sz;
    } else {
      // A_y/Brho
      AoBrho += W->BoBrhoV[i]*W->kxV[i]/(ky*kz_n)*sx*shy*sz;
      // dp_x
      dp += W->BoBrhoV[i]/kz_n*sqr(W->kxV[i]/ky)*cx*chy*sz;
    }
  }
}


template<typename T>
void wiggler_pass_Wu(const elem_type<T> &elem, ss_vect<T> &x)
{
  /* Second order symplectic integrator for insertion devices based on:

       Y.K. Wu, E. Forest, D.S. Robin "Explicit Symplectic Integrator for
       s-dependent Static Magnetic Field"                                    */

  int     i;
  T       h, z, hd, AxoBrho, AyoBrho, dpy, dpx;
  
  h = elem.L/elem.wiggler->n_step; z = 0.0;

  for (i = 1; i <= elem.wiggler->n_step; i++) {
    hd = h/(1.0+x[delta_]);

    // 1: half step in z
    z += 0.5*h;

    // 2: half drift in y
    get_Axy_Wu(elem.wiggler, z, x, AyoBrho, dpx, false);

    x[px_] -= dpx; x[py_] -= AyoBrho;
    x[y_] += 0.5*hd*x[py_];
    x[ct_] += sqr(0.5)*hd*sqr(x[py_])/(1.0+x[delta_]);
   
    get_Axy_Wu(elem.wiggler, z, x, AyoBrho, dpx, false);

    x[px_] += dpx; x[py_] += AyoBrho;

    // 3: full drift in x
    get_Axy_Wu(elem.wiggler, z, x, AxoBrho, dpy, true);

    x[px_] -= AxoBrho; x[py_] -= dpy; x[x_] += hd*x[px_];
    x[ct_] += 0.5*hd*sqr(x[px_])/(1.0+x[delta_]);

    if (totpath_on) x[ct_] += h;
   
    get_Axy_Wu(elem.wiggler, z, x, AxoBrho, dpy, true);

    x[px_] += AxoBrho; x[py_] += dpy;

    // 4: a half drift in y
    get_Axy_Wu(elem.wiggler, z, x, AyoBrho, dpx, false);

    x[px_] -= dpx; x[py_] -= AyoBrho;
    x[y_] += 0.5*hd*x[py_];
    x[ct_] += sqr(0.5)*hd*sqr(x[py_])/(1.0+x[delta_]);
   
    get_Axy_Wu(elem.wiggler, z, x, AyoBrho, dpx, false);

    x[px_] += dpx; x[py_] += AyoBrho;

    // 5: half step in z
    z += 0.5*h;
  }
}


template<typename T>
void cavity_pass(const elem_type<T> &elem, ss_vect<T> &x)
{
  T  delta;

  if (cavity_on && (elem.cavity->V_rf != 0.0)) {
    delta = -elem.cavity->V_rf/(E0*1e9)*sin(2.0*pi*elem.cavity->f_rf
	    /clight*x[ct_]);
    x[delta_] += delta;

    if (rad_on) dE -= is_double<T>::cst(delta);

    if (totpath_on) x[ct_] -= elem.cavity->h_rf/elem.cavity->f_rf*clight;
  }
}


template<typename T>
void kick_map_pass(const elem_type<T> &elem, ss_vect<T> &x)
{
  T       h;
  T       tx1, tz1;      /* thetax and thetaz retrieved from
			    interpolation routine */
  T       tx2, tz2;      /* thetax and thetaz retrieved from
			      interpolation routine */
  double  alpha0 = 0.0;  // 1/Brho
  double  alpha02= 0.0;  // alpha square
  int     n = 0, i = 0;
  bool    outoftable = false;
     
  n = elem.kick_map->n_step; alpha0 = clight/E0*1E-9*elem.kick_map->scl;
  alpha02= sqr(alpha0);
  
//  GtoL(X, elem->dS, elem->dT, 0.0, 0.0, 0.0);

  // (n+1) drifts, n kicks
  h = elem.L/(n+1);

  drift_pass(h, x);
  for (i = 1; i <= n; i++) {
    if (elem.kick_map->order == 1) {
      // first order kick
      if (elem.kick_map->method == 2)
        SplineInterpolation2(x[x_], x[y_], tx1, tz1, elem, outoftable);
      else
        LinearInterpolation2(x[x_], x[y_], tx1, tz1, elem, outoftable, 1);
      if (outoftable) {
	x[x_] = 1e30;
        return;
      }
      x[px_] += alpha0*tx1/n; x[py_] += alpha0*tz1/n;
    } else { 
      // second order kick
      if (elem.kick_map->method == 2)
        SplineInterpolation2(x[x_], x[y_], tx2, tz2, elem, outoftable);
      else
        LinearInterpolation2(x[x_], x[y_], tx2, tz2, elem, outoftable, 2);
      if (outoftable) {
	x[x_] = 1e30;
        return;
      }
      x[px_] += alpha02*tx2/(n*(1.0+x[delta_]));
      x[py_] += alpha02*tz2/(n*(1.0+x[delta_]));
    }  

    drift_pass(h, x);
  }
  
//  LtoG(X, elem->dS, elem->dT, 0.0, 0.0, 0.0);
}


/*double Kfunc(double **Cmat, int i, int j)
{
  double  lambdamat[3][3], Linv[3][3], lambda, arg;

  const int  Kdelta[3][3] = { {1, 0, 0}, {0, 1, 0}, {0, 0, 1} };

  lambdamat = Cmat - lambda ((1,0,0),(0,1,0),(0,0,1));
  Linv = Inv(Lambdamat);
  arg = (sqrt(lambda)/Det(Lambdamat))*(Kdelta[i][j]* Trace[Linv]-3*Linv[i][j]);

  return(2*sqr(pi)*qromb(arg, 0.0, 100.0)); //100 = infinity...
}*/


/*void add_IBS(const double Nb, const double L, const double **A1,
	     const double eps[], const double gamma_rel, double Deps[])
{*/
  /* A1 is passed, compute the invariants and emittances,
     The invariant matrices G for the uncoupled case are:

              |gamma alpha|
     G      = |           |
      x,y,z   |alpha beta |

  */

/*  const int  n_dim = 3;

  int     i, j;
  double  **Ainv, **Jmat, **boost, **G1, **G2, **G3;
  double  **Cmat1, **Cmat2, **Cmat3, **Cmat, **lambdamat, **Kmat;

  const double  tol = 1.E-9;

  const double  beta_rel = 1.0, P0 = E0*1e9*q_e/clight;

  const double  constA
    = Nb*sqr(r_e)*clight
      /(16.0*cube(pi*beta_rel)*eps[X_]*eps[Y_]*eps[Z_]*pow(gamma_rel, 4));

  Ainv = dmatrix(1, 2*n_dim, 1, 2*n_dim);
  Jmat = dmatrix(1, 2*n_dim, 1, 2*n_dim);
  boost = dmatrix(1, 2*n_dim, 1, 2*n_dim);
  G1 = dmatrix(1, 2*n_dim, 1, 2*n_dim);
  G2 = dmatrix(1, 2*n_dim, 1, 2*n_dim);
  G3 = dmatrix(1, 2*n_dim, 1, 2*n_dim);
  Cmat1 = dmatrix(1, 2*n_dim, 1, 2*n_dim);
  Cmat2 = dmatrix(1, 2*n_dim, 1, 2*n_dim);
  Cmat3 = dmatrix(1, 2*n_dim, 1, 2*n_dim);
  Cmat = dmatrix(1, 2*n_dim, 1, 2*n_dim);
  lambdamat = dmatrix(1, n_dim, 1, n_dim);
  Kmat = dmatrix(1, n_dim, 1, n_dim);

  Jmat = {
    {  0.0, 1.0,  0.0, 0.0, 0.0, 0.0 },
    { -1.0, 0.0,  0.0, 0.0, 0.0, 0.0 },
    {  0.0, 0.0,  0.0, 1.0, 0.0, 0.0 },
    {  0.0, 0.0, -1.0, 0.0, 0.0, 0.0 },
    {  0.0, 0.0,  0.0, 0.0, 0.0, 1.0 },
    {  0.0, 0.0,  0.0, 0.0,-1.0, 0.0 } 
  };

  boost = {
    { 1.0, 0.0, 0.0,           0.0,    0.0,    0.0          }, 
    { 0.0, 1.0, 0.0,           0.0,    0.0,    0.0          }, 
    { 0.0, 0.0, 1.0/gamma_rel, 0.0,    0.0,    0.0          }, 
    { 0.0, 0.0, 0.0,           1.0/P0, 0.0,    0.0          }, 
    { 0.0, 0.0, 0.0,           0.0,    1.0/P0, 0.0          }, 
    { 0.0, 0.0, 0.0,           0.0,    0.0,    gamma_rel/P0 }
  };

  // Construct the invariant matrices
  for (i = 0; i < 2*n_dim; i++)
    for(j = 0, j < 2*n_dim, j++ ) {
      G1[i][j] = Ainv[0][i]*Ainv[0][j] + Ainv[1][i]*Ainv[1][j];
      G2[i][j] = Ainv[2][i]*Ainv[2][j] + Ainv[3][i]*Ainv[3][j];
      G3[i][j] = Ainv[4][i]*Ainv[4][j] + Ainv[5][i]*Ainv[5][j];
    };
    
  // Transform the Invariant matrices into C.O.M. frame and into real
  // (not betatron) coordinates
//  M1 = (1/eps[X_])*Tp(boost)*G1*boost;
//  M2 = (1/eps[Y_])*Tp(boost)*G2*boost;
//  M3 = (1/eps[Z_])*Tp(boost)*G3*boost;

  // Build up C matrices out of momentum components of M1,2,3
  for(i = 1; i <= n_dim; i++)
    for(j = 1; j <= n_dim; j++) {
      Cmat1[i][j] = M1[2+i][2+j]; Cmat2[i][j] = M2[2+i][2+j];
      Cmat3[i][j] = M3[2+i][2+j];
    };
      
  Cmat = Cmat1+Cmat2+Cmat3;
  
  for(i = 1; i <= n_dim; i++)
    for(j = 1; j <= n_dim; j++)
      Kmat[i][j] = Kfunc(Cmat, i, j);
  
  // constA * Kmat = d Sigma/dt.  e.g. d <delta^2>/dt = constA Kmat[2,2].
  // Using Cmat, we can find the Deps.

  Deps[X_] += (L/Clight)*constA*Tr[Cmat1*Kmat];
  Deps[Y_] += (L/Clight)*constA*Tr[Cmat2*Kmat];
  Deps[Z_] += (L/Clight)*constA*Tr[Cmat3*Kmat];

  free_dmatrix(Ainv, 1, 2*n_dim, 1, 2*n_dim);
  free_dmatrix(Jmat, 1, 2*n_dim, 1, 2*n_dim);
  free_dmatrix(boost, 1, 2*n_dim, 1, 2*n_dim);
  free_dmatrix(G1, 1, 2*n_dim, 1, 2*n_dim);
  free_dmatrix(G2, 1, 2*n_dim, 1, 2*n_dim);
  free_dmatrix(G3, 1, 2*n_dim, 1, 2*n_dim);
  free_dmatrix(Cmat1, 1, 2*n_dim, 1, 2*n_dim);
  free_dmatrix(Cmat2, 1, 2*n_dim, 1, 2*n_dim);
  free_dmatrix(Cmat3, 1, 2*n_dim, 1, 2*n_dim);
  free_dmatrix(Cmat, 1, 2*n_dim, 1, 2*n_dim);
  free_dmatrix(lambdamat, 1, n_dim, 1, n_dim);
  free_dmatrix(Kmat, 1, n_dim, 1, n_dim);
}*/


/* symplectic integrator */
template<typename T>
bool si(const long int i0, const long int i1, ss_vect<T> &x,
	elem_type<T> elem[])
{
  long int  i;

  if (Check_Ampl(x)) return false;

  if (rad_on) dE = 0e0;

  if (emittance_on) {
    I2 = 0e0; I4 = 0e0; I5 = 0e0;

    for (i = 0; i < 3; i++)
      D_[i] = 0e0;
  }

  for (i = i0; i <= i1; i++) {
//    if (IBS_on) add_IBS();

    switch (elem[i-1].kind) {
    case Marker:
      break;
    case Drift:
      drift_pass(elem[i-1].L, x);
      break;
    case Mpole:
      mpole_pass(elem[i-1], x);
      break;
    case Wiggler:
      switch (elem[i-1].wiggler->method) {
      case 1:
	wiggler_pass(elem[i-1], x);
	break;
      case 2:
	wiggler_pass_Wu(elem[i-1], x);
	break;
      default:
	std::cout << "si: undefined method " << elem[i-1].wiggler->method
		  << std::endl;
      }
      break;
    case Cavity:
      cavity_pass(elem[i-1], x);
      break;
    case Kick_map:
      kick_map_pass(elem[i-1], x);
      break;
    default:
      std::cout << "si_templ: undefined element no "
	   << std::setw(3) << i << " type = " << elem[i-1].kind << std::endl;
      exit(0);
    }

    if (Check_Ampl(x)) {
//      cout << "si: particle lost at " << i << "(" << i1 << ")" << endl;
      return false;
    }

    // is_tps<T>::get_ps(x, elem[i-1]);
  }

  return true;
}
