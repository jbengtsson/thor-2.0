// E0 contain kinetic energy [GeV].
extern double E0, dE, beta0, gamma0;
extern tps    I2, I4, I5;
extern tps    D_[];

const double pi = M_PI;

const double clight  = 2.99792458e8;   // speed of light in vacuum
const double q_e     = 1.60217646e-19; // electron charge
const double m_e     = 0.51099906e6;   // electron rest mass [eV]
const double mu_0    = 4.0*M_PI*1e-7;  // permittivity of free space
const double eps_0   = 1.0/(sqr(clight)*mu_0); // permeability of free space
const double h_bar   = 6.58211899e-16; // reduced Planck constant [eV]
const double r_e     = q_e/(4.0*M_PI*eps_0*m_e); // classical electron radius
const double C_u     = 55.0/(24.0*sqrt(3.0));
const double C_gamma = 4.0*M_PI*r_e/(3.0*cube(1e-9*m_e));


// Element type codes.
enum elemkind { Marker = -1, Drift = 0, Mpole = 1, Cavity = 2,
		Thinkick = 3, Wiggler = 4, Undef = 5, Kick_map = 6 };

// Integration methods.
enum intmeth  { Linear = 0, Second = 2, Fourth = 4 };

// Multipole coefficients.
enum mpoles { Dip = 1, Quad = 2, Sext = 3, Oct = 4, Dec = 5, Dodec = 6 };

// Max multipole order.
/* const int mpole_max = Quad; */
/* const int mpole_max = Dodec; */
/* const int mpole_max = 12; */
const int mpole_max = 21;

template<typename T> class mpole_type {
 public:
  int    method, n_step;
  double dx_sys[2], dx_rms[2], dx_rnd[2];
  double droll_par, droll_sys, droll_rms, droll_rnd;
  double an_par[mpole_max], an_sys[mpole_max];
  double an_rms[mpole_max], an_rnd[mpole_max];
  double bn_par[mpole_max], bn_sys[mpole_max];
  double bn_rms[mpole_max], bn_rnd[mpole_max];
  T      an[mpole_max];
  T      bn[mpole_max];
  int    order, n_design;
  double edge1, edge2;
  double gap;
  double h_bend;
};

const int  n_harm_max = 10;

typedef struct {
  int    method, n_step;
  double lambda;
  int    n_harm;              // no of harmonics
  int    harm[n_harm_max];    // harmonic number
  double BoBrhoV[n_harm_max]; // B/Brho vertical
  double BoBrhoH[n_harm_max]; // B/Brho horizontal 
  double kxV[n_harm_max];     // kx 
  double kxH[n_harm_max];     // kx 
  double phi[n_harm_max];     // phi 
} wiggler_type;

typedef struct {
  double V_rf, f_rf;
  int    h_rf;
} cavity_type;

const int IDXMAX = 200, IDZMAX = 100;

typedef struct {
  int    method,                       /* interpolation method
					  (linear = 1, spline = 2) */
         n_step,                       // number of steps
         order;                        // first or second order kick map
  char   file_name[100];               // filename for insertion description
  double scl;                          // scale factor
  int    nx, nz;                       // hor/ver no of points
  double tabx[IDXMAX], tabz[IDZMAX], // hor/ver spacing
         thetax[IDZMAX][IDXMAX],
         thetax1[IDZMAX][IDXMAX],    // first order => 1
         thetaz[IDZMAX][IDXMAX],
         thetaz1[IDZMAX][IDXMAX],
         **tx, **tz,
         **f2x, **f2z,
         **tx1, **tz1,
         **f2x1, **f2z1,               // a voir
         *tab1, *tab2;                 // tab of x/z meshes from Radia
} kick_map_type;

const int name_length = 15;

/* define element structure */
template<typename T> class elem_type {
 public:
  char       Name[name_length];
  T          L;
  double     S;
  int        Fnum, Knum;
  double     dx[2];
  double     droll[2];  // cos(dr), sin(dr)
  double     droll_par;
  double     c0, c1, s1;
  int        kind;
  double     max_ampl[2][2];
  ss_vect<T> A1;
  double     Alpha[2], Beta[2], Nu[2], Eta[2], Etap[2];
  union {
    // nothing for drift
    mpole_type<T> *mpole;
    wiggler_type  *wiggler;
    cavity_type   *cavity;
    kick_map_type *kick_map;
  };
};


const int  max_Kids = 700;  // max no of kids

struct Family {
  char Name[name_length];
  int  n_Kids;
  int  Kids[max_Kids];
};


const int  max_elem   = 6000, // max no of elements
           max_Family = 1000; // max no of families
// RHIC
//const int  max_elem   = 7700, // max no of elements
//           max_Family = 2200; // max no of families


extern long int          n_elem;

extern elem_type<double> elem[];
extern elem_type<tps>    elem_tps[];
extern Family            Families[];

extern bool lost;

extern bool rad_on, H_exact, totpath_on, cavity_on, quad_fringe_on;
extern bool emittance_on, IBS_on;

void ini_si(void);
