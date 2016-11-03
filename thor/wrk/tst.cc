#define NO 2

#include "thor_lib.h"
#include "tpsa_for.h"

int  no_tps   = NO,
     ndpt_tps = 5;


int main(int argc, char *argv[])
{

  danot_(no_tps-1);

  rad_on    = false; H_exact        = false; totpath_on   = false;
  cavity_on = false; quad_fringe_on = false; emittance_on = false;
  IBS_on    = false;

  rd_mfile(argv[1], elem); rd_mfile(argv[1], elem_tps);

  // initialize the symplectic integrator after energy has been defined
  ini_si();

  // disable from TPSALib- and LieLib log messages
  idprset(-1);
}
