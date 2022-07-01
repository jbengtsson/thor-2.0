#include "thor_lib.h"

#include <vector>


struct param_type {
private:

public:
  std::vector<std::string>
    name;
  int
    n_prm;
  std::vector<double>
    L,
    bnL_min,
    bnL_max,
    bnL_scl,
    bnL;
  std::vector<int>
    Fnum,
    n;
  std::vector< std::vector<int> >
    locs;

  void add_Fam(const std::string Fname, const int n,
	       const double bnL_min, const double bnL_max,
	       const double bnL_scl);
  void create_Fam(const std::string Fname, const int n,
		  const double bnL_min, const double bnL_max,
		  const double bnL_scl, const std::vector<int> &locs);
  void ini_prm(void);
  void set_prm(void);
  void print(const int k) const;
  void print(void) const;
};

double bnL_internal(const double bnL_bounded, const double bnL_min,
		   const double bnL_max);
double bnL_bounded(const double bnL_internal, const double bnL_min,
		  const double bnL_max);

