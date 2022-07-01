#include "param_type.h"


void param_type::add_Fam(const std::string name, const int n,
			 const double bnL_min, const double bnL_max,
			 const double bnL_scl)
{
  const bool prt = false;

  this->name.push_back(name);
  this->Fnum.push_back(get_Fnum(name.c_str()));
  this->n.push_back(n);
  this->L.push_back(0e0);
  this->bnL_min.push_back(bnL_min);
  this->bnL_max.push_back(bnL_max);
  this->bnL_scl.push_back(bnL_scl);
  this->bnL.push_back(0e0);  
  this->n_prm = this->Fnum.size();

  if (prt) printf("add_prm: %2d\n", this->n_prm);
}


void param_type::create_Fam(const std::string Fname, const int n,
			    const double bnL_min, const double bnL_max,
			    const double bnL_scl,
			    const std::vector<int> &locs)
{
  const bool prt = false;

  this->name.push_back(Fname);
  this->Fnum.push_back(-1);
  this->n.push_back(n);
  this->L.push_back(0e0);
  this->bnL_min.push_back(bnL_min);
  this->bnL_max.push_back(bnL_max);
  this->bnL_scl.push_back(bnL_scl);
  this->bnL.push_back(0e0);
  this->locs.push_back(locs);
  this->n_prm = this->Fnum.size();

  if (prt) printf("add_prm: %2d\n", this->n_prm);
}


void param_type::ini_prm(void)
{
  int    k, Fnum;
  double bnL_ext;

  bool prt = false;

  if (prt) printf("\n");
  for (k = 0; k < n_prm; k++) {
    Fnum = (this->Fnum[k] > 0)? this->Fnum[k] : elem[locs[k][0]-1].Fnum;
    L[k] = get_L(Fnum, 1);
    bnL_ext = get_bnL(Fnum, 1, n[k]);
    if ((bnL_min[k] <= bnL_ext) && (bnL_ext <= bnL_max[k])) {
      bnL[k] = bnL_internal(bnL_ext, bnL_min[k], bnL_max[k]);
      if (prt) {
	printf("ini_prm: k = %1d Fnum = %1d\n", k, Fnum);
	print(k);
      }
    } else {
      printf("\nini_prm:\n outside range: %-8s %10.3e %10.3e [%10.3e, %10.3e]"
	     "\n",
	     name[k].c_str(), bnL_ext, bnL_ext, bnL_min[k], bnL_max[k]);
      exit(1);
    }
  }
}

void param_type::set_prm(void)
{
  int    k, Fnum;
  double bnL_ext;

  const bool prt = false;

  const int  n_prt = 6;

  if (prt) printf("\nset_prm:\n  ");
  for (k = 0; k <= n_prm; k++) {
    Fnum = (this->Fnum[k] > 0)? this->Fnum[k] : elem[locs[k][0]-1].Fnum;
    bnL_ext = bnL_bounded(bnL[k], bnL_min[k], bnL_max[k]);
    set_bnL(Fnum, n[k], bnL_ext, 0e0);
    if (prt) {
      printf(" %12.5e", bnL_ext);
      if (k % n_prt == 0) printf("\n  ");
    }
  }
  if (prt && (n_prm % n_prt != 0)) printf("\n");
}


void param_type::print(const int k) const
{
  double bnL_ext;

  // Bounded.
  bnL_ext = bnL_bounded(bnL[k], bnL_min[k], bnL_max[k])/L[k];
  printf(" %-8s %10.3e [%9.3e, %9.3e] %9.3e\n",
	 name[k].c_str(), bnL_ext, bnL_min[k], bnL_max[k], L[k]);
}


void param_type::print(void) const
{
  int k;

  printf("\n");
  for (k = 0; k < n_prm; k++)
    print(k);
}


double bnL_internal(const double bnL_bounded,
		   const double bnL_min, const double bnL_max)
{
#if 0
  return asin((2e0*(bnL_bounded-bnL_min))/(bnL_max-bnL_min)-1e0);
#else
  return bnL_bounded;
#endif
}


double bnL_bounded(const double bnL_internal,
		  const double bnL_min, const double bnL_max)
{
#if 0
  return bnL_min + (sin(bnL_internal)+1e0)*(bnL_max-bnL_min)/2e0;
#else
  return bnL_internal;
#endif
}
