/* Author:        Johan Bengtsson

    Definitions:  Interface to Fortran library for Truncated Power
		  Series Algebra.                                       */

//#include "mpi.h"

extern int  no_tps, ndpt_tps;

bool    ini_tps = false, header = false, res_basis = false, stable = false,
        debug_tpsa  = false;

const int  n_max    = 100;     // max iterations for LieExp

const int   name_len_for = 11; // name length in FORTRAN library is 10 + NUL
const char  tpsa_name[] = "tpsa       "; 

int  bufsize;  // Note, max no of monomials is (no+nv)!/(nv!*no!)


long int fact(long int n)
{
  if (n > 0)
    return n*fact(n-1);
  else if (n == 0)
    return 1;
  else {
    cout << "fact: neg. argument: " << n << endl;
    exit(1);
  }
}


long int nok(long int n, long int k)
{
  long int  j;
  double    u;

  u = 1.0;
  for (j = 0; j < k; j++)
    u *= (double)(n-j)/(double)(k-j);
  return (long int)(u+0.5);
}


void daexp(const int &intptr, double rbuf[], int ibuf1[], int ibuf2[],
	   char name[name_len_for])
{

  daexp_(intptr, rbuf, ibuf1, ibuf2, name); name[10] = '\0';
}



// Interface to Fortran TPSA library

extern "C" {
  // for g77 compability
//  void f_init(void);
//  void MAIN_() { cout << "call to MAIN_" << endl; }

  // for g95 compability
//  void g95_runtime_start(int argc, const char *argv[]); 
//  void g95_runtime_stop();
}
 

void TPSA_Ini(void)
{

  cout << endl;
  cout << scientific << setprecision(0)
       << "TPSA_Ini: no = " << no_tps << ", nv = " << ss_dim
       << ", eps = " << eps_tps << endl;

  // initialize g77 I/O
  // cout << "initilizing g77 i/o" << endl;
  // f_init();

  // initialize g95 I/O
  // cout << "initilizing g95 i/o" << endl;
  // g95_runtime_start(0, NULL); 

  // Initialize TPSA-lib.
  daini_(no_tps, ss_dim, 0);

  daeps_(eps_tps);

  // Initialize Lie-lib.
  lieinit_(no_tps, ss_dim, nd_tps, ndpt_tps, iref_tps, 0);

  bufsize = nok(no_tps+ss_dim, ss_dim);

  ini_tps = true;
}

void TPSAEps(const double eps)
{ daeps_(eps); eps_tps = eps; }

tps::tps(void) {

  if (!ini_tps) TPSA_Ini();
  intptr = 0;
  daall_(intptr, 1, tpsa_name, no_tps, ss_dim); dacon_(intptr, 0.0);
  if (debug_tpsa)
    cout << "tps(void):                        "
	 << ", intptr = " << intptr << endl;
}

tps::tps(const double r)
{

  if (!ini_tps) TPSA_Ini();
  intptr = 0; daall_(intptr, 1, tpsa_name, no_tps, ss_dim); dacon_(intptr, r);
  if (debug_tpsa)
    cout << "tps(const double r):              "
	 << ", intptr = " << intptr << endl;
}

tps::tps(const double r, const int i)
{

  if (!ini_tps) TPSA_Ini();
  intptr = 0; daall_(intptr, 1, tpsa_name, no_tps, ss_dim);
  if (i == 0)
    dacon_(intptr, r);
  else
    davar_(intptr, r, i);
  if (debug_tpsa)
    cout << "tps(const double r, const int i): "
	 << ", intptr = " << intptr << endl;
}

tps::tps(const tps &x) {

  if (!ini_tps) TPSA_Ini();
  intptr = 0;
  daall_(intptr, 1, tpsa_name, no_tps, ss_dim); dacop_(x.intptr, intptr);
  if (debug_tpsa)
    cout << "tps(const tps &x):                "
	 << ", intptr = " << intptr << endl;
}

tps::~tps(void) {

  if (debug_tpsa)
    cout << "~tps(void):                       "
	 << ", intptr = " << intptr << endl;

  dadal_(intptr, 1);
}


double tps::operator[](const int k) const
{
  int     i, jj[ss_dim];
  double  r;

  for (i = 0; i < ss_dim; i++)
    jj[i] = 0;
  jj[k] = 1;
  dapek_(intptr, jj, r);
  return(r);
}

double tps::operator[](const int jj[]) const
{
  double  r;

  dapek_(intptr, jj, r);
  return(r);
}

void tps::pook(const int jj[], const double r)
{ dapok_(intptr, jj, r); }

void tps::exprt(double rbuf[], int ibuf1[], int ibuf2[], char *name) const
{ daexp_(intptr, rbuf, ibuf1, ibuf2, name); }

void tps::imprt(const int n, double rbuf[],
		const int ibuf1[], const int ibuf2[])
{ rbuf[0] = n; daimp_(rbuf, ibuf1, ibuf2, intptr); }

tps& tps::operator=(const double r) { dacon_(intptr, r); return *this; }

tps& tps::operator+=(const double x)
{ dacad_(intptr, x, intptr); return *this; }

tps& tps::operator-=(const double x)
{ dacad_(intptr, -x, intptr); return *this; }

tps& tps::operator*=(const double x)
{ dacmu_(intptr, x, intptr); return *this; }

tps& tps::operator/=(const double x)
{ dacmu_(intptr, 1.0/x, intptr); return *this; }


tps& tps::operator=(const tps &x) { dacop_(x.intptr, intptr); return *this; }

tps& tps::operator+=(const tps &x)
{ daadd_(intptr, x.intptr, intptr); return *this; }

tps& tps::operator-=(const tps &x)
{ dasub_(intptr, x.intptr, intptr); return *this; }

tps& tps::operator*=(const tps &x)
{ damul_(intptr, x.intptr, intptr); return *this; }

tps& tps::operator/=(const tps &x)
{ dadiv_(intptr, x.intptr, intptr); return *this; }


tps sqrt(const tps &a)
{
  tps  b;

  dafun_("SQRT", a.intptr, b.intptr);
  return b;
}

tps pow(const tps &a, const int n)
{
  if (n < 0)
    return tps(pow(a, n+1)) /= a;
  else if (n == 0)
    return tps(1.0);
  else if (n == 1)
    return tps(a);
  else if (n > 1)
    return tps(pow(a, n-1)) *= a;
  else {
    cout << "pow: should never get here " << n << endl;
    exit(1);
  }
}

tps exp(const tps &a)
{
  tps  b;

  dafun_("EXP ", a.intptr, b.intptr);
  return b;
}

tps log(const tps &a)
{
  tps  b;

  dafun_("LOG ", a.intptr, b.intptr);
  return b;
}

#if false

tps sin(const tps &a)
{
  tps  b;

  dafun_("SIN ", a.intptr, b.intptr);
  return b;
}

tps cos(const tps &a)
{
  tps  b;

  dafun_("COS ", a.intptr, b.intptr);
  return b;
}

tps tan(const tps &a)
{
  tps  b;

  dafun_("TAN ", a.intptr, b.intptr);
  return b;
}

#else

tps sin_tps(const tps &a)
{
  int       n;
  long int  k;
  tps       r, b;

  b = 0.0; k = 1; r = a;
  for (n = 1; n <= no_tps; n += 2) {
    b += r/(double)k; k *= (n+1)*(n+2); r *= -sqr(a);
  }

  return b;
}

tps cos_tps(const tps &a)
{
  int       n;
  long int  k;
  tps       r, b;

  b = 0.0; k = 1; r = 1.0;
  for (n = 0; n <= no_tps; n += 2) {
    b += r/(double)k; k *= (n+1)*(n+2); r *= -sqr(a);
  }

  return b;
}

tps sin(const tps &a)
{
  // sin(a+b) = sin(a)*cos(b) + cos(a)*sin(b)
  double  cst;
  tps     b;

  cst = a.cst(); b = a - cst;

  return sin(cst)*cos_tps(b)+cos(cst)*sin_tps(b);
}

tps cos(const tps &a)
{
  // cos(a+b) = cos(a)*cos(b) - sin(a)*sin(b)
  double  cst;
  tps     b;

  cst = a.cst(); b = a - cst;

  return cos(cst)*cos_tps(b)-sin(cst)*sin_tps(b);
}

tps tan(const tps &a) { return sin(a)/cos(a); }

#endif

tps asin(const tps &a)
{
  tps  b;

  dafun_("ASIN", a.intptr, b.intptr);
  return b;
}

#if false

tps atan(const tps &a)
{
  tps  b;

  dafun_("ATAN", a.intptr, b.intptr);
  return b;
}

#else

tps atan(const tps &a)
{
  // arctan(a+b) to 8th order
  double  cst;
  tps     b, c;

  if (no_tps <= 8) {
    cst = a.cst(); b = a - cst;

    c = b/(1+sqr(cst)) - (cst*sqr(b))/pow(1.0+sqr(cst), 2) +
      ((-1+3.0*sqr(cst))*pow(b, 3))/(3.0*pow(1.0+sqr(cst), 3)) +
      ((cst-pow(cst, 3))*pow(b, 4))/pow(1+sqr(cst), 4) +
      ((1-10*sqr(cst)+5.0*pow(cst, 4))*pow(b, 5))/(5.0*pow(1+sqr(cst), 5))+
      ((-3*cst+10.0*pow(cst, 3)-3.0*pow(cst, 5))*pow(b, 6))/
      (3.0*pow(1+sqr(cst), 6)) +
      ((-1+21.0*sqr(cst)-35.0*pow(cst, 4)+7.0*pow(cst, 6))*pow(b, 7))/
      (7.0*pow(1+sqr(cst), 7)) +
      ((cst-7.0*pow(cst, 3)+7.0*pow(cst, 5)-pow(cst, 7))*pow(b, 8))/
      pow(1+sqr(cst), 8)
      + ((1.0-36.0*sqr(cst) + 126.0*pow(cst, 4) - 
        84.0*pow(cst, 6) + 9.0*pow(cst, 8))*pow(b, 9))/
      (9.0*pow(1.0+sqr(cst), 9))
      + atan(cst);
   } else {
    cout << "atan: only defined to 8th order (" << no_tps
	 << ")" << endl;
    exit(1);
  }

  return c;
}

#endif

/*int sgn(const double a)
{
  if (a > 0.0)
    return 1;
  else if (a == 0.0)
    return 0;
  else
    return -1;
}*/

tps atan2(const tps &b, const tps &a) {
  tps  c;

  if (a.cst() > 0.0)
    c = atan(b/a);
  else if (a.cst() == 0.0)
    if (b.cst() != 0.0)
      c = sgn(b.cst())*pi/2.0;
    else {
      cout << "atan2: 0/0 undefined" << endl;
      exit(1);
    }
  else
    if (b.cst() >= 0.0)
      c = atan(b/a) + pi;
    else
      c = atan(b/a) - pi;
  return c;
}

#if false

tps sinh(const tps &a)
{
  tps  b;

  dafun_("SINH", a.intptr, b.intptr);

  return b;
}

tps cosh(const tps &a)
{
  tps  b;

  dafun_("COSH", a.intptr, b.intptr);

  return b;
}

#else

tps sinh_tps(const tps &a)
{
  int       n;
  long int  k;
  tps       r, b;

  b = 0.0; k = 1; r = a;
  for (n = 1; n <= no_tps; n += 2) {
    b += r/(double)k; k *= (n+1)*(n+2); r *= sqr(a);
  }

  return b;
}

tps cosh_tps(const tps &a)
{
  int       n;
  long int  k;
  tps       r, b;

  b = 0.0; k = 1; r = 1.0;
  for (n = 0; n <= no_tps; n += 2) {
    b += r/(double)k; k *= (n+1)*(n+2); r *= sqr(a);
  }

  return b;
}

tps sinh(const tps &a)
{
  // sinh(a+b) = sinh(a)*cosh(b) + cosh(a)*sinh(b)
  double  cst;
  tps     b;

  cst = a.cst(); b = a - cst;

  return sinh(cst)*cosh_tps(b)+cosh(cst)*sinh_tps(b);
}

tps cosh(const tps &a)
{
  // cos(a+b) = cos(a)*cos(b) - sin(a)*sin(b)
  double  cst;
  tps     b;

  cst = a.cst(); b = a - cst;

  return cosh(cst)*cosh_tps(b)+sinh(cst)*sinh_tps(b);
}

#endif

const double tps::cst(void) const
{
  int     i, jj[ss_dim];
  double  r;

  for (i = 0; i < ss_dim; i++)
    jj[i] = 0;
  dapek_(intptr, jj, r);
  return r;
}

double abs(const tps &a)
{
  double  r;

  daabs_(a.intptr, r);
  return r;
}

double abs2(const tps &a)
{
  double  r;

  daabs2_(a.intptr, r);
  return r;
}


void idprset(const int level)
{
  idprset_(level);
}


tps Der(const tps &a, const int k)
{
  tps  b;

  dader_(k, a.intptr, b.intptr);
  return b;
}

tps pseudo_der(const tps &a, const int k)
{
  tps  b;

  datra_(k, a.intptr, b.intptr);
  return b;
}

tps LieExp(const tps &H, const tps &x)
{
  tps  y;

  exp1d_(H.intptr, x.intptr, y.intptr, eps_tps, n_max);
  return y;
}

ss_vect<tps> LieExp(const tps &H, const ss_vect<tps> &x)
{
  int           i, xintptrs[ps_dim], mapintptrs[ps_dim];
  ss_vect<tps>  map;

  for (i = 0; i < ps_dim; i++) {
    xintptrs[i] = x[i].intptr; mapintptrs[i] = map[i].intptr;
  }
  expnd2_(H.intptr, xintptrs, mapintptrs, eps_tps, n_max);
  return map;
}

ss_vect<tps> FExpo(const tps &H, const ss_vect<tps> &x,
                   const int k0, const int k1, const int k)
{
  int           i, xintptrs[ps_dim], mapintptrs[ps_dim];
  ss_vect<tps>  map;

  for (i = 0; i < ps_dim; i++) {
    xintptrs[i] = x[i].intptr; mapintptrs[i] = map[i].intptr;
  }
  fexpo_(H.intptr, xintptrs, mapintptrs, k0, k1, 1.0, k);
  for (i = ps_dim; i < ss_dim; i++)
    map[i] = tps(0.0, i+1);
  return map;
}

ss_vect<tps> MTREE(const ss_vect<tps> &x)
{
  int           i, xintptrs[ss_dim], yintptrs[ss_dim];
  ss_vect<tps>  y;

  for (i = 0; i < ss_dim; i++) {
    xintptrs[i] = x[i].intptr; yintptrs[i] = y[i].intptr;
  }
  etmtree_(xintptrs, yintptrs);
  return y;
}

ss_vect<double> PPUSH(const ss_vect<tps> &x, ss_vect<double> &y)
{
  int              i, xintptrs[ss_dim];
  ss_vect<double>  z;

  for (i = 0; i < ss_dim; i++)
    xintptrs[i] = x[i].intptr;
  etppush2_(xintptrs, y, z);
  return z;
}

tps operator*(const tps &x, const ss_vect<tps> &y)
{
  int  i, xintptrs[ps_dim], y1intptrs[ss_dim], zintptrs[ps_dim];
  tps  z;
  tps  y1[ss_dim];

  xintptrs[0] = x.intptr; zintptrs[0] = z.intptr;
  for (i = 0; i < ss_dim; i++) {
    y1[i] = (i < ps_dim)? y[i] : tps(0.0, i+1);
    y1intptrs[i] = y1[i].intptr;
  }
  dacct_(xintptrs, 1, y1intptrs, ss_dim, zintptrs, 1);
  return z;
}

ss_vect<tps> operator*(const ss_vect<tps> &x, const ss_vect<tps> &y)
{
  int           i, xintptrs[ps_dim], yintptrs[ps_dim], zintptrs[ps_dim];
  ss_vect<tps>  z;

  for (i = 0; i < ps_dim; i++) {
    xintptrs[i] = x[i].intptr; yintptrs[i] = y[i].intptr;
    zintptrs[i] = z[i].intptr;
  }
  etcct_(xintptrs, yintptrs, zintptrs);
  return z;
}

ss_vect<tps> CCT(const ss_vect<tps> &x, const ss_vect<tps> &y)
{
  int           i, xintptrs[ps_dim], yintptrs[ps_dim], zintptrs[ps_dim];
  ss_vect<tps>  z;

  for (i = 0; i < ps_dim; i++) {
    xintptrs[i] = x[i].intptr; yintptrs[i] = y[i].intptr;
    zintptrs[i] = z[i].intptr;
  }
  dacct_(xintptrs, ps_dim, yintptrs, ps_dim, zintptrs, ps_dim);
  return z;
}

ss_vect<tps> Inv(const ss_vect<tps> &x)
{
  int           i, xintptrs[ps_dim], yintptrs[ps_dim];
  ss_vect<tps>  y;

  for (i = 0; i < ps_dim; i++) {
    xintptrs[i] = x[i].intptr; yintptrs[i] = y[i].intptr;
  }
  etinv_(xintptrs, yintptrs);
  return y;
}

ss_vect<tps> Inv_Ext(const ss_vect<tps> &x)
{
  int           i, xintptrs[ps_dim], yintptrs[ps_dim];
  ss_vect<tps>  y;

  for (i = 0; i < ps_dim; i++) {
    xintptrs[i] = x[i].intptr; yintptrs[i] = y[i].intptr;
  }
  dainv_(xintptrs, ps_dim, yintptrs, ps_dim);
  return y;
}

ss_vect<tps> PInv(const ss_vect<tps> &x, const int jj[])
{
  int           i, xintptrs[ps_dim], yintptrs[ps_dim];
  ss_vect<tps>  y;

  for (i = 0; i < ps_dim; i++) {
    xintptrs[i] = x[i].intptr; yintptrs[i] = y[i].intptr;
  }
  etpin_(xintptrs, yintptrs, jj);
  return y;
}

void PInv(const int m, const tps x[], const int n, tps y[], const int jj[])
{
  int  i, xintptrs[m], yintptrs[n];

  for (i = 0; i < m; i++)
    xintptrs[i] = x[i].intptr;
  for (i = 0; i < n; i++)
    yintptrs[i] = y[i].intptr;
  dapin_(xintptrs, m, yintptrs, n, jj);
}

void GoFix(const ss_vect<tps> &xy, ss_vect<tps> &a1, ss_vect<tps> &a1inv,
	   const int nord)
{
  int  i, xyintptrs[ps_dim], a1intptrs[ps_dim], a1invintptrs[ps_dim];

  for (i = 0; i < ps_dim; i++) {
    xyintptrs[i] = xy[i].intptr; a1intptrs[i] = a1[i].intptr;
    a1invintptrs[i] = a1inv[i].intptr;
  }
  gofix_(xyintptrs, a1intptrs, a1invintptrs, nord);
}

tps MapNorm(const ss_vect<tps> &x, tps &g, ss_vect<tps> &a2,
	    ss_vect<tps> &a1, ss_vect<tps> &xy, const int nord)
{
  int  i, xintptrs[ps_dim], a2intptrs[ps_dim], a1intptrs[ps_dim];
  int  xyintptrs[ps_dim];
  tps  K;

  for (i = 0; i < ps_dim; i++) {
    xintptrs[i] = x[i].intptr; a2intptrs[i] = a2[i].intptr;
    a1intptrs[i] = a1[i].intptr; xyintptrs[i] = xy[i].intptr;
  }
  stable = mapnorm_(xintptrs, g.intptr, a2intptrs, a1intptrs, xyintptrs,
		    K.intptr, nord);
  return K;
}

ss_vect<tps> MapNormF(const ss_vect<tps> &x, ss_vect<tps> &g, ss_vect<tps> &a2,
		      ss_vect<tps> &a1, ss_vect<tps> &xy,
		      const int nord, const int kpmax)
{
  int   i, xintptrs[ps_dim], gintptrs[ps_dim], a2intptrs[ps_dim];
  int   a1intptrs[ps_dim], xyintptrs[ps_dim], Kintptrs[ps_dim];
  ss_vect<tps>  K;

  for (i = 0; i < ps_dim; i++) {
    xintptrs[i] = x[i].intptr; gintptrs[i] = g[i].intptr;
    a2intptrs[i] = a2[i].intptr; a1intptrs[i] = a1[i].intptr;
    xyintptrs[i] = xy[i].intptr;  Kintptrs[i] = K[i].intptr;
  }
  stable = mapnormf_(xintptrs, gintptrs, a2intptrs, a1intptrs, xyintptrs,
		     Kintptrs, nord, kpmax);
  return K;
}

ss_vect<tps> dHdJ(const tps &H)
{
  int           i, nuintptrs[ps_dim];
  ss_vect<tps>  nu;

  for (i = 0; i < ps_dim; i++)
    nuintptrs[i] = nu[i].intptr;
  dhdj_(H.intptr, nuintptrs);
  return nu;
}

void CtoR(const tps &a, tps &a_re, tps &a_im)
{
  ctor_(a.intptr, a_re.intptr, a_im.intptr);
}

tps RtoC(const tps &a_re, const tps &a_im)
{
  tps  a;

  rtoc_(a_re.intptr, a_im.intptr, a.intptr);
  return a;
}

tps LieFact_DF(const ss_vect<tps> &xy, ss_vect<tps> &x)
{
  /* Dragt-Finn factorization:

       M = M_lin exp(:h_3:) exp(:h_4:) ... 

  */
  int  i, xyintptrs[ps_dim], xintptrs[ps_dim];
  tps  H;

  for (i = 0; i < ps_dim; i++) {
    xyintptrs[i] = xy[i].intptr; xintptrs[i] = x[i].intptr;
  }
  liefact_(xyintptrs, xintptrs, H.intptr);
  return H;
}

tps LieFact(const ss_vect<tps> &xy)
{
  /* Single exponent Dragt-Finn factorization:

       M = exp(:h_2:) exp(:h_3:) exp(:h_4:) ... 

  */
  return Intd(FlowFact(xy), -1.0); 
}

ss_vect<tps> FlowFact(const ss_vect<tps> &xy)
{
  int           i, xyintptrs[ps_dim], Vintptrs[ps_dim];
  ss_vect<tps>  V;

  for (i = 0; i < ps_dim; i++) {
    xyintptrs[i] = xy[i].intptr; Vintptrs[i] = V[i].intptr;
  }
  flofacg_(xyintptrs, Vintptrs, eps_tps);
  return V;
}

tps Intd(const ss_vect<tps> &V, double eps)
{
  int  i, Vintptrs[ps_dim];
  tps  H;

  for (i = 0; i < ps_dim; i++)
    Vintptrs[i] = V[i].intptr;
  intd_(Vintptrs, H.intptr, eps);
  return H;
}

tps PB(const tps &a, const tps &b)
{
  tps  c;

  etpoi_(a.intptr, b.intptr, c.intptr);
  return c;
}

tps Take(const tps &H, const int n)
{
  tps  Hn;

  take_(H.intptr, n, Hn.intptr);
  return Hn;
}

istream& operator>>(istream &is, tps &a)
{
  int     i, n, no1, nv1;
  int     ibuf1[bufsize], ibuf2[bufsize], jj[ss_dim];
  double  rbuf[bufsize];
  char	  line[max_str];

  const bool  debug = true;

  if (debug) {
    darea77_(a.intptr, 8);
    return is;
  }

  is.getline(line, max_str); is.getline(line, max_str);
  sscanf(line, "%*17c %d %*6c %d", &no1, &nv1);
  if ((no1 <= no_tps) && (nv1 <= ss_dim)) {
    for (i = 0; i < 4; i++)
      is.getline(line, max_str);
    ibuf1[0] = no_tps; ibuf2[0] = ss_dim; n = 0;
    sscanf(line, "%d %le", &no1, &rbuf[n+1]);
    while (no1 >= 0) {
      n++; 
      for (i = 0; i < ss_dim; i++) {
	is.getline(line, max_str);
	sscanf(line, "%d", &jj[i]); 
      }
      hash_(no_tps, ss_dim, jj, ibuf1[n], ibuf2[n]);
      is.getline(line, max_str);
      sscanf(line, "%d %le", &no1, &rbuf[n+1]);
    } 
    is.getline(line, max_str);
    rbuf[0] = -no1;
    daimp_(rbuf, ibuf1, ibuf2, a.intptr);
  } else
    cout << "*** illegal no (" << no_tps << ") or nv ("
	 << ss_dim << ")" << endl;
  return is;
}


ostream& operator<<(ostream &os, const tps &a)
{
  char           name[name_len_for];
  int            i, j, ord, n, no;
  int            ibuf1[bufsize], ibuf2[bufsize], jj[ss_dim];
  double         rbuf[bufsize];
  ostringstream  s;

  const bool  debug = false;

  if (debug) {
    dapri77_(a.intptr, 7);
    return os;
  }

  daexp(a.intptr, rbuf, ibuf1, ibuf2, name);
  s << endl;
  
  name[10] = '\0'; i = 0;
  while ((i <= 10) && (name[i] != ' ')) {
    s << name[i]; i++;
  }
  n = (int) rbuf[0];
  s << ", NO = " << no_tps
    << ", NV = " << ss_dim << ", INA = " << a.intptr << endl;

  for (i = 1; i <= 66; i++)
    s << "-"; 
  s << endl;

  if (header) {
    s << endl;
    if (!res_basis) {
      s << "                                                        n"
	<< endl;
      s << "      ====     i  i   i  i  i   i  i     i             ===="
	<< endl;
      s << "      \\         1  2   3  4  5   6  7     n            \\   "
	<< endl;
      s << "  P =  |   a  x  p   y  p  d  ct  p ... p  ,    |I| =  |   i"
	<< endl;
      s << "      /     I     x      y         1     n             /     k"
	<< endl;
      s << "      ====                                             ===="
	<< endl;
      s << "       I                                               k=1"
	<< endl;
    } else {
      s << "                                                          n"
	<< endl;
      s << "      ====      i   i   i   i  i   i  i     i            ===="
	<< endl;
      s << "      \\        + 1 - 2 + 3 - 4  5   6  7     n           \\   "
	<< endl;
      s << "  P =  |   a  h   h   h   h   d  ct  p ... p  ,    |I| =  |   i"
	<< endl;
      s << "      /     I  x   x   y   y          1     n            /     k"
	<< endl;
      s << "      ====                                               ===="
	<< endl;
      s << "       I                                                 k=1"
	<< endl;
    }
  }
  
  if (n != 0) {
    s << endl;
    s << "   |I|         a              ";
    for (i = 1; i <= ss_dim; i++)
      s << "  i";
    s << endl;
    s << "                I              ";
    for (i = 1; i <= ss_dim; i++)
      s << setw(3) << i;
    s << endl;
    s << endl;
  } else
    s << "   ALL COMPONENTS ZERO " << endl;
  for (no = 0; no <= no_tps; no++) {
    for (i = 1; i <= n; i++) {
      dehash_(no_tps, ss_dim, ibuf1[i-1], ibuf2[i-1], jj);
      ord = 0;
      for (j = 0; j < ss_dim; j++) 
	ord += jj[j];
      if (ord == no)
	if (fabs(rbuf[i]) >= eps_tps) {
	  s << setw(5) << ord << scientific << setw(24)
	    << setprecision(16) << rbuf[i] << " ";
	  for (j = 0; j < ss_dim; j++)
	    s << setw(3) << jj[j];
	  s << endl;
	}
    }
  }
  if (n == 0) n = 1;
  s << setw(5) << -n
    << scientific << setw(24) << setprecision(16) << 0.0 << " ";
  for (j = 0; j < ss_dim; j++)
    s << setw(3) << 0;
  s << endl;

  return os << s.str();
}
