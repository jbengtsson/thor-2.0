/* Author:	 Johan Bengtsson

   Definitions:  Polymorphic number class.              */

#include <sstream>
using namespace std;

#define sqr(x) ((x)*(x))
#define cube(x) ((x)*(x)*(x))

#define fract(x) ((x)-(int)(x))
#define nint(x) ((x) < 0 ? ((long)(x-0.5)) : ((long)(x+0.5))) 

#define min(a, b) ((a) < (b) ? (a) : (b))
#define max(a, b) ((a) > (b) ? (a) : (b))

#define sgn(n) ((n > 0) ? 1 : ((n < 0) ? -1 : 0)) 

#define dtor(x)  ((x)*pi/180.0)


const int max_str = 150;


template<typename T> class ss_vect;

// Polymorphic class for floating point and TPSA

class tps {
 public:
  tps(void);
  tps(const double);
  tps(const double, const int);
  tps(const tps &);
  ~tps(void);

  // initialize TPSA library
   friend void TPSAEps(const double);
  // trace level for TPSALib and LieLib
  friend void idprset(const int);

  const double cst(void) const;
  double operator[](const int) const;
  double operator[](const int []) const;
  void pook(const int [], const double);

  void exprt(double [], int [], int [], char []) const;
  void imprt(const int, double [], const int [], const int []);

  tps& operator=(const double);
  tps& operator+=(const double);
  tps& operator-=(const double);
  tps& operator*=(const double);
  tps& operator/=(const double);

  tps& operator=(const tps &);
  tps& operator+=(const tps &);
  tps& operator-=(const tps &);
  tps& operator*=(const tps &);
  tps& operator/=(const tps &);

  friend istream& operator>>(istream &, tps &);
  template<typename CharT, class Traits>
  friend basic_istream<CharT, Traits>&
    operator>>(basic_istream<CharT, Traits> &, ss_vect<tps> &);
  friend ostream& operator<<(ostream &, const tps &);

  friend double abs(const tps &);
  friend double abs2(const tps &);
  friend tps sqrt(const tps &);
//  friend tps sqr(const tps &);
  friend tps pow(const tps &, const int);
  friend tps exp(const tps &);
  friend tps log(const tps &);
  friend tps sin(const tps &);
  friend tps cos(const tps &);
  friend tps tan(const tps &);
  friend tps asin(const tps &);
  friend tps atan(const tps &);
  friend tps sinh(const tps &);
  friend tps cosh(const tps &);

  friend tps Der(const tps &, const int);
  friend tps pseudo_der(const tps &, const int);
  friend ss_vect<tps> FExpo(const tps &, const ss_vect<tps> &,
                            const int, const int, const int);
  friend tps LieExp(const tps &, const tps &);
  friend tps PB(const tps &, const tps &);
  friend tps Take(const tps &, const int);

  // R(nd2, nv) = P(nd2, nd2)*Q(nd2, nv)
  friend ss_vect<tps> operator*(const ss_vect<tps> &, const ss_vect<tps> &);
  // R(nv, nv) = P(nv, nv)*Q(nv, nv)
  friend ss_vect<tps> CCT(const ss_vect<tps> &, const ss_vect<tps> &);
  friend ss_vect<tps> MTREE(const ss_vect<tps> &);
  friend ss_vect<double> PPUSH(const ss_vect<tps> &, ss_vect<double> &);
  friend tps operator*(const tps &, const ss_vect<tps> &);

  friend ss_vect<tps> LieExp(const tps &, const ss_vect<tps> &);
  // Q(nv, nv) = P(nd2, nd2)^-1
  friend ss_vect<tps> Inv(const ss_vect<tps> &);
  // Q(nv, nv) = P(nv, nv)^-1
  friend ss_vect<tps> Inv_Ext(const ss_vect<tps> &);
  friend void PInv(const int, const tps [], const int, tps [], const int []);
  friend ss_vect<tps> PInv(const ss_vect<tps> &, const int []);
  friend void GoFix(const ss_vect<tps> &, ss_vect<tps> &,
		    ss_vect<tps> &, const int);
  friend tps MapNorm(const ss_vect<tps> &, tps &, ss_vect<tps> &,
		     ss_vect<tps> &, ss_vect<tps> &, const int);
  friend ss_vect<tps> MapNormF(const ss_vect<tps> &, ss_vect<tps> &,
			       ss_vect<tps> &, ss_vect<tps> &,
			       ss_vect<tps> &, const int, const int);
  friend ss_vect<tps> dHdJ(const tps &);
  friend void CtoR(const tps &, tps &, tps &);
  friend tps RtoC(const tps &, const tps &);
  friend tps LieFact_DF(const ss_vect<tps> &, ss_vect<tps> &);
  friend ss_vect<tps> FlowFact(const ss_vect<tps> &);
  friend tps Intd(const ss_vect<tps> &, const double);
 private:
  int     intptr; // index used by Fortran implementation
  double  r;      // floating-point calc. if intptr = 0
};


template<typename T> class elem_type;

// Abstract class for general beam dynamics
class channel {
  double  L, href, phi, phi1, phi2;
 public:
  virtual ~channel(void);
  virtual bool propagate(const long int, const long int) = 0;
};


// Class for single particle phase space dynamics

const int ps_dim = 6;         // phase space dimension
const int ss_dim = ps_dim+1; // state phase space dimension

// spatial components
enum spatial_index { X_ = 0, Y_ = 1, Z_ = 2 };

// phase space components
// (Note, e.g. spin components should be added here)
enum ps_index { x_ = 0, px_ = 1, y_ = 2, py_ = 3, delta_ = 4, ct_ = 5 };


template<typename T> class ss_vect {
 public:
  typedef T value_type;

  ss_vect(void);
// Let's the compiler synthetize the copy constructor
//  ss_vect(const T &a) { }
//  ss_vect(const ss_vect<T> &a) { }
  template<typename U>
    ss_vect(const ss_vect<U> &);

  ss_vect<double> cst(void) const;
  T& operator[](const int i) { return ss[i]; }
  const T& operator[](const int i) const { return ss[i]; }

  ss_vect<T>& operator*=(const double);
  ss_vect<T>& operator*=(const tps &);

  ss_vect<T>& operator=(const ss_vect<T> &);
  ss_vect<T>& operator+=(const ss_vect<T> &);
  ss_vect<T>& operator-=(const ss_vect<T> &);

  friend ss_vect<tps> operator+(const ss_vect<tps> &, const ss_vect<tps> &);

//  friend ss_vect<double> operator*(const ss_vect<tps> &,
//				   const ss_vect<double> &);
  // R(nd2, nv) = P(nd2, nd2)*Q(nd2, nv)
  friend ss_vect<tps> operator*(const ss_vect<tps> &, const ss_vect<tps> &);
  // R(nv, nv) = P(nv, nv)*Q(nv, nv)
  friend ss_vect<tps> CCT(const ss_vect<tps> &, const ss_vect<tps> &);
  friend ss_vect<tps> MTREE(const ss_vect<tps> &);
  friend ss_vect<double> PPUSH(const ss_vect<tps> &, ss_vect<double> &);
  friend tps operator*(const tps &, const ss_vect<tps> &);

  template<typename CharT, class Traits>
    friend basic_istream<CharT, Traits>&
    operator>>(basic_istream<CharT, Traits> &, ss_vect<T> &);

  template<typename CharT, class Traits>
    friend basic_ostream<CharT, Traits>&
    operator<<(basic_ostream<CharT, Traits> &, const ss_vect<T> &);

  void zero(void);
  void identity(void);

  friend ss_vect<tps> FExpo(const tps &, const ss_vect<tps> &,
                            const int, const int, const int);
  friend ss_vect<tps> LieExp(const tps &, const ss_vect<tps> &);
  // Q(nv, nv) = P(nd2, nd2)^-1
  friend ss_vect<tps> Inv(const ss_vect<tps> &);
  // Q(nv, nv) = P(nv, nv)^-1
  friend ss_vect<tps> Inv_Ext(const ss_vect<tps> &);
  friend ss_vect<tps> PInv(const ss_vect<tps> &, const int []);
  friend void GoFix(const ss_vect<tps> &, ss_vect<tps> &,
		    ss_vect<tps> &, const int);
  friend tps MapNorm(const ss_vect<tps> &, tps &, ss_vect<tps> &,
		      ss_vect<tps> &, ss_vect<tps> &, const int);
  friend ss_vect<tps> MapNormF(const ss_vect<tps> &, ss_vect<tps> &,
			       ss_vect<tps> &, ss_vect<tps> &,
			       ss_vect<tps> &, const int, const int);
  friend void dHdJ(const tps &, ss_vect<tps> &);
  friend void CtoR(const tps &, tps &, tps &);
  friend tps RtoC(const tps &, const tps &);
  friend tps LieFact_DF(const ss_vect<tps> &, ss_vect<tps> &);
  friend tps LieFact(const ss_vect &);
  friend ss_vect<tps> FlowFact(const ss_vect<tps> &);
  friend tps Intd(const ss_vect<tps> &, const double);

  bool propagate(const long int, const long int);
 private:
  // (Note, e.g. spin components should be added here)
  T  ss[ss_dim];
};


typedef ss_vect<double>  Vector; 
typedef Vector           Matrix[ss_dim];

/* pre-declare for:
     inline ss_vect<double> get_FS(const ss_vect<double> ps)    */
template<>
ss_vect<double> ss_vect<tps>::cst(void) const;

// partial template-class specialization
// primary version: is_double<>
template<typename T>
class is_double { };

// partial specialization
template<>
class is_double<double> {
 public:
  static inline double cst(const double x) { return x; }
};

// partial specialization
template<>
class is_double<tps> {
 public:
  static inline double cst(const tps &x) { return x.cst(); }
};

//template<typename T>
//inline T sqr(const T &a)
//{ return a*a; }


inline tps operator+(const tps &a, const tps &b)
{ return tps(a) += b; }

inline tps operator+(const tps &a, const double b)
{ return tps(a) += b; }

inline tps operator+(const double a, const tps &b)
{ return tps(a) += b; }

inline tps operator-(const tps &a, const tps &b)
{ return tps(a) -= b; }

inline tps operator-(const tps &a, const double b)
{ return tps(a) -= b; }

inline tps operator-(const double a, const tps &b)
{ return tps(a) -= b; }

inline tps operator*(const tps &a, const tps &b)
{ return tps(a) *= b; }

inline tps operator*(const tps &a, const double b)
{ return tps(a) *= b; }

inline tps operator*(const double a, const tps &b)
{ return tps(a) *= b; }

inline tps operator/(const tps &a, const tps &b)
{ return tps(a) /= b; }

inline tps operator/(const tps &a, const double b)
{ return tps(a) /= b; }

inline tps operator/(const double a, const tps &b)
{ return tps(a) /= b; }


inline tps operator+(const tps &x)
{ return tps(x); }

inline tps operator-(const tps &x)
{ return tps(x) *= -1.0; }


inline bool operator>(const tps &a, const tps &b)
{ return a.cst() > b.cst(); }

inline bool operator>(const tps &a, const double b)
{ return a.cst() > b; }

inline bool operator>(const double a, const tps &b)
{ return a > b.cst(); }


inline bool operator<(const tps &a, const tps &b)
{ return a.cst() < b.cst(); }

inline bool operator<(const tps &a, const double b)
{ return a.cst() < b; }

inline bool operator<(const double a, const tps &b)
{ return a < b.cst(); }


inline bool operator>=(const tps &a, const tps &b)
{ return a.cst() >= b.cst(); }

inline bool operator>=(const tps &a, const double b)
{ return a.cst() >= b; }

inline bool operator>=(const double a, const tps &b)
{ return a >= b.cst(); }


inline bool operator<=(const tps &a, const tps &b)
{ return a.cst() <= b.cst(); }

inline bool operator<=(const tps &a, const double b)
{ return a.cst() <= b; }

inline bool operator<=(const double a, const tps &b)
{ return a <= b.cst(); }


inline bool operator==(const tps &a, const tps &b)
{ return a.cst() == b.cst(); }

inline bool operator==(const tps &a, const double b)
{ return a.cst() == b; }

inline bool operator==(const double a, const tps &b)
{ return a == b.cst(); }


inline bool operator!=(const tps &a, const tps &b)
{ return a.cst() != b.cst(); }

inline bool operator!=(const tps &a, const double b)
{ return a.cst() != b; }

inline bool operator!=(const double a, const tps &b)
{ return a != b.cst(); }


template<typename T>
inline ss_vect<T>::ss_vect(void)
{
  int  i;

  for (i = 0; i < ss_dim; i++)
//    ss[i] = T();
    /* use operator=(const double r) instead of operator=(const tps &x)
       to avoid allocation of temporary variable */
    ss[i] = 0.0;
}

template<typename T>
template<typename U>
inline ss_vect<T>::ss_vect(const ss_vect<U> &a)
{
  int              i;

  for (i = 0; i < ss_dim; i++)
    ss[i] = a[i];
}
