thor-2.0

Author: Johan Bengtsson

Self-Consistent Symplectic Integrator for charged particle beam dynamics,
based on TPSA (Truncated Power Series Algebra), aka PTC (Polymorphic Tracking
Code), originated 1994; by implementing a transparent polymorphic number object
with reference counting for FP/TPSA in C++.

The symplectic integrator for realistic modeling of magnetic lattices for
ring-based synchrotrons was initially implemented in Pascal, by the author,
with care taken for the software architecture and resulting records/modules
(-> "objects") to reflect the structure of the mathematical objects describing
the underlying beam dynamics model.


Requirements:

   GNU C/C++ and FORTRAN-95 compilers: gcc and gfortran.
   GNU Scientific Library GSL.
   GNU autoconf/automake environment and libtool.
   "Numerical Recipes in C": http://www.nr.com.

To install:

   mkdir git_repos
   cd git_repos
   git clone git@github.com:jbengtsson/thor-2.0.git
   cd shor-2.0
   ./make_thor.sh
