#include "param_type.h"


void param_type::add_prm(const std::string name, const int n,
			 const double bn_min, const double bn_max,
			 const double bn_scl)
{
  const bool prt = false;

  this->name.push_back(name);
  this->Fnum.push_back(get_Fnum(name.c_str()));
  this->n.push_back(n);
  this->bn_min.push_back(bn_min);
  this->bn_max.push_back(bn_max);
  this->bn_scl.push_back(bn_scl);
  this->bn.push_back(0e0);  
  this->n_prm = this->Fnum.size();

  if (prt) printf("add_prm: %2d\n", this->n_prm);
}


void param_type::ini_prm(void)
{
  int    k;
  double bn_ext;

  bool prt = false;

  for (k = 0; k < n_prm; k++) {
    bn_ext = get_bn(Fnum[k], 1, n[k]);
    // Bounded.
    if ((bn_min[k] <= bn_ext) && (bn_ext <= bn_max[k])) {
      bn[k] = bn_internal(bn_ext, bn_min[k], bn_max[k]);
      if (prt) {
	printf("ini_prm:");
	print(k);
      }
  } else {
      printf("\nini_prm:\n outside range: %-8s %10.3e %10.3e [%10.3e, %10.3e]"
	     "\n",
	     name[k].c_str(), bn_ext, bn_ext, bn_min[k], bn_max[k]);
      exit(1);
    }
  }
}

void param_type::set_prm(void)
{
  int    i;
  double bn_ext;

  const bool prt = false;

  const int  n_prt = 6;

  if (prt) printf("\nset_prm:\n  ");
  for (i = 0; i <= n_prm; i++) {
    // Bounded.
    bn_ext = bn_bounded(bn[i], bn_min[i], bn_max[i]);
    set_bn(Fnum[i], n[i], bn_ext, 0e0);
    if (prt) {
      printf(" %12.5e", bn_ext);
      if (i % n_prt == 0) printf("\n  ");
    }
  }
  if (prt && (n_prm % n_prt != 0)) printf("\n");
}


void param_type::print(const int k) const
{
  double bn_ext;

  // Bounded.
  bn_ext = bn_bounded(bn[k], bn_min[k], bn_max[k]);
  printf(" %-8s %10.3e [%9.3e, %9.3e]\n",
	 name[k].c_str(), bn_ext, bn_min[k], bn_max[k]);
}


void param_type::print(void) const
{
  int k;

  printf("\n");
  for (k = 0; k < n_prm; k++)
    print(k);
}


double bn_internal(const double bn_bounded,
		   const double bn_min, const double bn_max)
{
  return asin((2e0*(bn_bounded-bn_min))/(bn_max-bn_min)-1e0);
}


double bn_bounded(const double bn_internal,
		  const double bn_min, const double bn_max)
{
  return bn_min + (sin(bn_internal)+1e0)*(bn_max-bn_min)/2e0;
}
