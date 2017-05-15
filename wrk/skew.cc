double get_an(const int Fnum, const int Knum, const int n)
{
  return elem[get_loc(Fnum, Knum)-1].mpole->an[n-1];
}


void set_dan(const int Fnum, const int Knum, const int n, const double dan)
{
  int  k;

  k = get_loc(Fnum, Knum) - 1;
  elem[k].mpole->an[n-1] += dan; elem_tps[k].mpole->an[n-1] += dan;
  if (n > elem[k].mpole->order) {
    elem[k].mpole->order = n; elem_tps[k].mpole->order = n;
  }
}


void set_dan(const int Fnum, const int n, const double dan)
{
  int  j;

  for (j = 1; j <= get_n_Kids(Fnum); j++)
    set_dan(Fnum, j, n, dan);
}


void set_an_par(const int Fnum, const int Knum, const int n, const int j)
{
  // set parameter dependence
  int     k;
  double  an;

  k = get_loc(Fnum, Knum) - 1;
  an = elem_tps[k].mpole->an[n-1].cst();
  elem_tps[k].mpole->an[n-1] = tps(an, j);
  if (n > elem_tps[k].mpole->order) elem_tps[k].mpole->order = n;
}


void set_an_par(const int Fnum, const int n, const int j)
{
  // set parameter dependence
  int  k;

  for (k = 1; k <= get_n_Kids(Fnum); k++)
    set_an_par(Fnum, k, n, j);
}


void clr_an_par(const int Fnum, const int Knum, const int n)
{
  // clear parameter dependence
  int     k;
  double  an;

  k = get_loc(Fnum, Knum) - 1;
  an = elem_tps[k].mpole->an[n-1].cst(); elem_tps[k].mpole->an[n-1] = an;
  // clear order
}


void clr_an_par(const int Fnum, const int n)
{
  // set parameter dependence
  int  k;

  for (k = 1; k <= get_n_Kids(Fnum); k++)
    clr_an_par(Fnum, k, n);
}


