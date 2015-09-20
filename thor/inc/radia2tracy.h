/* Author:  Johan Bengtsson  */

void Read_IDfile(char *fic_radia, double &L, int &pnx, int &pnz,
                 double tabx[IDXMAX],  double tabz[IDZMAX],
                 double thetax[IDZMAX][IDXMAX], double thetaz[IDZMAX][IDXMAX]);

template<typename T>
void LinearInterpolation2(T &X, T &Z, T &TX, T &TZ, const elem_type<T> &elem,
			  bool &out, int order);

template<typename T>
void SplineInterpolation2(T &X, T &Z, T &thetax, T &thetaz,
			  const elem_type<T> &elem, bool &out);

void Matrices4Spline(kick_map_type *kick_map);

template<typename T>
void spline(const double x[], const T y[], int const n,
	    double const yp1, const double ypn, T y2[]);

template<typename T, typename U>
void splint(const double xa[], const U ya[], const U y2a[],
	    const int n, const T &x, T &y);

template<typename T>
void splin2(const double x1a[], const double x2a[],
	    double **ya, double **y2a, const int m, const int n,
	    const T &x1, const T &x2, T &y);

void splie2(double x1a[], double x2a[], double **ya,
	    int m, int n, double **y2a);
