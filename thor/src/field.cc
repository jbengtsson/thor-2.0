/* Author:	 Johan Bengtsson

   Definitions:  Polymorphic number class.              */


// partial specialization
template<>
ss_vect<double> ss_vect<tps>::cst(void) const
{
  int              i;
  ss_vect<double>  x;

  for (i = 0; i < ss_dim; i++)
    x[i] = ss[i].cst();
  return x;
}


template<typename T>
ss_vect<T>& ss_vect<T>::operator=(const ss_vect<T> &x)
{
  int  i;

  for (i = 0; i < ss_dim; i++)
    ss[i] = x[i];
  return *this;
}

template<typename T>
ss_vect<T>& ss_vect<T>::operator+=(const ss_vect<T> &a)
{
  int  i;

  for (i = 0; i < ss_dim; i++)
    ss[i] += a[i];
  return *this;
}

template<typename T>
ss_vect<T>& ss_vect<T>::operator-=(const ss_vect<T> &a)
{
  int  i;

  for (i = 0; i < ss_dim; i++)
    ss[i] -= a[i];
  return *this;
}

template<typename T>
ss_vect<T>& ss_vect<T>::operator*=(const double a)
{
  int  i;

  for (i = 0; i < ss_dim; i++)
    ss[i] *= a;
  return *this;
}

template<>
ss_vect<tps>& ss_vect<tps>::operator*=(const tps &a)
{
  int  i;

  for (i = 0; i < ss_dim; i++)
    ss[i] *= a;
  return *this;
}


template<typename T>
ss_vect<T> operator+(const ss_vect<T> &x) { return ss_vect<T>(x); }

template<typename T>
ss_vect<T> operator-(const ss_vect<T> &x) { return ss_vect<T>(x) *= -1; }

// instantiate
template ss_vect<double> operator-(const ss_vect<double> &);
template ss_vect<tps> operator-(const ss_vect<tps> &);

ss_vect<double> operator+(const ss_vect<double> &a, const ss_vect<double> &b)
{ return ss_vect<double>(a) += b; }

ss_vect<tps> operator+(const ss_vect<tps> &a, const ss_vect<double> &b)
{ return ss_vect<tps>(a) += b; }

ss_vect<tps> operator+(const ss_vect<double> &a, const ss_vect<tps> &b)
{ return ss_vect<tps>(a) += b; }

ss_vect<tps> operator+(const ss_vect<tps> &a, const ss_vect<tps> &b)
{ return ss_vect<tps>(a) += b; }

ss_vect<double> operator-(const ss_vect<double> &a, const ss_vect<double> &b)
{ return ss_vect<double>(a) -= b; }

ss_vect<tps> operator-(const ss_vect<tps> &a, const ss_vect<double> &b)
{ return ss_vect<tps>(a) -= b; }

ss_vect<tps> operator-(const ss_vect<double> &a, const ss_vect<tps> &b)
{ return ss_vect<tps>(a) -= b; }

ss_vect<tps> operator-(const ss_vect<tps> &a, const ss_vect<tps> &b)
{ return ss_vect<tps>(a) -= b; }

ss_vect<double> operator*(const ss_vect<double> &a, const double b)
{ return ss_vect<double>(a) *= b; }

ss_vect<double> operator*(const double a, const ss_vect<double> &b)
{ return ss_vect<double>(b) *= a; }

ss_vect<tps> operator*(const ss_vect<tps> &a, const double b)
{ return ss_vect<tps>(a) *= b; }

ss_vect<tps> operator*(const double a, const ss_vect<tps> &b)
{ return ss_vect<tps>(b) *= a; }


template<typename T>
void ss_vect<T>::zero(void)
{
  int         i;

  for (i = 0; i < ss_dim; i++)
    ss[i] = 0.0;
}

// partial specialization
template<>
void ss_vect<tps>::identity(void)
{
  int           i;

  for (i = 0; i < ss_dim; i++)
   ss[i] = tps(0.0, i+1);
}


template<typename CharT, class Traits>
std::basic_istream<CharT, Traits>&
operator>>(std::basic_istream<CharT, Traits> &is, ss_vect<tps> &a)
{
  int                                 i;
  std::basic_istringstream<CharT, Traits>  s;

  for (i = 0; i < ss_dim; i++)
    darea77_(a[i].intptr, 7);
  return is;
}

// instantiate
template std::basic_istream<char, std::char_traits<char> >&
operator>>(std::basic_istream<char, std::char_traits<char> > &,
	   ss_vect<tps> &);

template<typename CharT, class Traits>
std::basic_ostream<CharT, Traits>&
operator<<(std::basic_ostream<CharT, Traits> &os, const ss_vect<double> &a)
{
  int                                 i;
  std::basic_ostringstream<CharT, Traits>  s;

  s.flags(os.flags()); s.imbue(os.getloc());
  for (i = 0; i < 6; i++)
    s << std::setprecision(os.precision()) << std::setw(os.width()) << a[i];
  return os << s.str();
}

// instantiate
template std::basic_ostream<char, std::char_traits<char> >&
operator<<(std::basic_ostream<char, std::char_traits<char> > &,
	   const ss_vect<double> &);

template<typename CharT, class Traits>
std::basic_ostream<CharT, Traits>&
operator<<(std::basic_ostream<CharT, Traits> &os, const ss_vect<tps> &a)
{
  int                                 i;
  std::basic_ostringstream<CharT, Traits>  s;

  s.flags(os.flags()); s.imbue(os.getloc());
  for (i = 0; i < ss_dim; i++)
    s << std::setprecision(os.precision()) << std::setw(os.width()) << a[i];
  return os << s.str();
}

// instantiate
template std::basic_ostream<char, std::char_traits<char> >&
operator<<(std::basic_ostream<char, std::char_traits<char> > &,
	   const ss_vect<tps> &);


template<typename T>
bool si(const long int, const long int, ss_vect<T> &, elem_type<T> []);

extern elem_type<double>  elem[];
extern elem_type<tps>     elem_tps[];

// partial specialization
template<>
bool ss_vect<double>::propagate(const long int i0, const long int i1)
{ return si(i0, i1, *this, elem); }

// partial specialization
template<>
bool ss_vect<tps>::propagate(const long int i0, const long int i1)
{ return si(i0, i1, *this, elem_tps); }
