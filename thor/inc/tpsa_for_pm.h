extern bool  stable, debug_tpsa, header;

extern int  bufsize;  // Note, max no of monomials is (no+nv)!/(nv!*no!)


long int fact(long int n);

long int nok(long int n, long int k);

double getmat(const ss_vect<tps> &map, const int i, const int j);

void putmat(ss_vect<tps> &map, const int i, const int j, const double r);

void getlinmat(const int nv, const ss_vect<tps> &map, Matrix &mat);

void putlinmat(const int nv, const Matrix &mat, ss_vect<tps> &map);

void idprset(const int level);
