#include <vector>


struct param_type {
private:

public:
  std::vector<std::string>
    name;
  int
    n_prm;
  std::vector<double>
    bn_min,
    bn_max,
    bn_scl,
    bn;
  std::vector<int>
    Fnum,
    n;

  void add_prm(const std::string Fname, const int n,
	       const double bn_min, const double bn_max,
	       const double bn_scl);
  void ini_prm(void);
  void set_prm(void);
  void print(const int k) const;
  void print(void) const;
};

double bn_internal(const double bn_bounded, const double bn_min,
		   const double bn_max);
double bn_bounded(const double bn_internal, const double bn_min,
		  const double bn_max);

