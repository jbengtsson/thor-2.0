#define ORDER 10

#include "thor_lib.h"

#include "field.cc"
#include "tpsa_for_pm.cc"

#include "radia2tracy.cc"
#include "si.cc"
#include "rd_mfile.cc"
#include "prt_mfile.cc"

#include "tools.cc"

const int  nv_tps   = ss_dim; // no of variables
const int  nd_tps   = 3;      // degrees of freedom
const int  iref_tps = 0;      // resonance file: fort.7

double     eps_tps  = 1e-25;


// instantiate templates

template class ss_vect<double>; 

template class ss_vect<tps>;

template void rd_mfile(const char [], elem_type<double> []);

template void rd_mfile(const char [], elem_type<tps> []);
